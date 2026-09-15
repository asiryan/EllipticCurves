"""One feedback loop for exact points, certificates and rational search charts.

Local characters select a provably independent subset, not the entire useful
point pool. Inconclusive observations remain available as projection centres
and as candidates for exact division of relations. Chart coverage survives
basis changes. None of these bounded searches is a completeness theorem.
"""
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
import time

from bootstrap import save
from models import prepare_models, record_coverage, missing_intervals
from point_models import divide_basis, prepare_covers, refine_observed
from local_models import expand_models
from torsion_certificate import certify, torsion_points, translate_pool


def point_height(point):
    """Cheap exact size priority; not a canonical height or independence test."""
    x,y=map(Q,point)
    return max(abs(x.numerator).bit_length(),x.denominator.bit_length(),
               (2*abs(y.numerator).bit_length()+2)//3)


def reservoir(points, basis, count=32, offset=0):
    """Preserve small and rotating other observations, without assuming dependence."""
    certified=set(map(tuple,basis));pending=[p for p in points if tuple(p) not in certified]
    pending.sort(key=lambda p:(point_height(p),tuple(map(Q,p))))
    short=pending[:count//2];other=pending[count//2:]
    if other:
        i=offset%len(other);other=other[i:]+other[:i]
    return short+other[:count-len(short)]


def search(source, output, seconds=300, workers=4, anchors=2048, target=32,
           batch_size=8, job_seconds=2, import_run=None):
    # The archived implementations remain available for explicit comparisons.
    # They are never executed or combined by this solver.
    from seeded import ROOT, POLICY, clean, read, generate, diverse_vectors, order_models
    from seeded import boxes, batch_search, job_key, transported_models
    if import_run:raise ValueError('Unified search resumes its own checkpoint only')
    source,output=Path(source).resolve(),Path(output).resolve()
    if output==ROOT or not output.is_relative_to(ROOT):raise ValueError('Use a dedicated workspace directory')
    output.mkdir(parents=True,exist_ok=True)
    names=['seeded.py','unified.py','models.py','bootstrap.py','certificate.py',
           'torsion_certificate.py','point_models.py','local_models.py']
    config={'input_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
        'code_sha256':{n:hashlib.sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in names},
        'policy':POLICY,'anchors':anchors,'target':target,'workers':workers,
        'batch_size':batch_size,'job_seconds':job_seconds,'anchor_mode':'unified',
        'search_policy':{'shared_exact_observations':True,'relation_division':2,
            'retained_charts':True,'pending_reservoir':32,'adaptive_timeout':True}}
    config['code_sha256'].update({n:hashlib.sha256((ROOT/'tools/RankHunt'/n).read_bytes()).hexdigest()
        for n in ['blind_search.py','point_search.py','bounded_anchor_pool.py',
                  'anchor_diversity.py','point_arithmetic.py']})
    path=output/'checkpoint.json'
    if path.exists():
        state=read(path)
        if state['config']!=config:raise ValueError('Resume requires unchanged input, code and settings')
        found=state['found']
    else:
        found=clean(read(source))
        state={'config':config,'source':str(source),'generation':0,'wall_seconds':0,
            'models':[],'active':[],'coverage':{},'model_cache':{},'events':[],
            'preparations':[],'attempted_models':0,'finished_models':0,
            'timed_out_batches':0,'finished_slices':0,'refinements':[],
            'attempts':{},'reservoir_offset':0,'basis_prepared':None,'expansions':{},'sweep':0,
            'origin':'unified seeded search'}
    if not found['points']:raise ValueError('Known seed points are required')
    save(output/'basis.json',{'ainvs':found['ainvs'],'points':state.get('certified_basis',found['points'])})
    proof=certify(output/'basis.json')
    if not proof['all_selected_independent'] or not proof['points']:raise ValueError('No independent seed basis')
    state.setdefault('initial_lower_bound',proof['LowerBound'])
    known=set(map(tuple,found['points']));cache=state['model_cache'];coverage=state['coverage']
    charts={m['key']:m for m in state['models']}
    active=[charts[k] for k in state['active'] if k in charts]
    started=time.perf_counter();deadline=started+seconds;before=state['wall_seconds'];status='running'
    last_save=started;last_progress=started
    def left():return deadline-time.perf_counter()
    def signature():return hashlib.sha256(json.dumps(proof['points']).encode()).hexdigest()
    def checkpoint():
        nonlocal last_save
        state.update(found=found,models=list(charts.values()),active=[m['key'] for m in active],
            status=status,certified_basis=proof['points'],lower_bound=proof['LowerBound'],
            wall_seconds=before+time.perf_counter()-started)
        save(path,state)
        save(output/'points.json',{'ainvs':found['ainvs'],'points':proof['points'],
                                 'rank_lower_bound':proof['LowerBound']})
        last_save=time.perf_counter()
    def incorporate(observations, recertify=True):
        nonlocal proof
        new=[];old_bound=proof['LowerBound'];old_points=set(map(tuple,proof['points']))
        for observation in observations:
            point=clean({'ainvs':found['ainvs'],'points':[observation['point']]})['points'][0]
            if tuple(point) in known:continue
            observation={**observation,'point':point};known.add(tuple(point))
            found['points'].append(point);new.append(observation)
        if new:
            with (output/'discoveries.jsonl').open('a',encoding='utf-8') as stream:
                for observation in new:stream.write(json.dumps(observation)+'\n')
            if recertify:
                pending=reservoir(found['points'],proof['points'],32,state['reservoir_offset'])
                state['reservoir_offset']+=24
                candidates=proof['points']+[o['point'] for o in new]+pending
                candidates=list(dict.fromkeys(map(tuple,candidates)))
                save(output/'basis.json',{'ainvs':found['ainvs'],
                    'points':candidates})
                candidate=certify(output/'basis.json')
                if candidate['all_selected_independent'] and candidate['LowerBound']>=old_bound:
                    proof=candidate
                else:state['inconclusive_updates']=state.get('inconclusive_updates',0)+1
        if proof['LowerBound']>old_bound:
            selected={tuple(p) for p in proof['points']}-old_points
            event={'run':output.name,'lower_bound':proof['LowerBound'],
                'seconds':before+time.perf_counter()-started,'attempted_models':state['attempted_models'],
                'selected_new_points':[o for o in new if tuple(o['point']) in selected]}
            state['events'].append(event);state['generation']+=1
            print(json.dumps({k:v for k,v in event.items() if k!='selected_new_points'}),flush=True)
        return proof['LowerBound']>old_bound, bool(new)
    def register(models, producer, refresh=True):
        nonlocal active
        fresh=[];ordered=[]
        for model in models:
            if model['key'] not in charts:
                model=dict(model,producer=producer,preparation=len(state['preparations']),
                           uses_new_direction=None)
                charts[model['key']]=model;fresh.append(model)
            if refresh or model in fresh:ordered.append(charts[model['key']])
        if ordered:
            # Retain old charts and their coverage. Rotate the tail when there
            # are too many for one sweep, instead of deleting it on rank growth.
            active=ordered+[m for m in active if m['key'] not in {n['key'] for n in ordered}]
        return len(fresh)
    def refine():
        if left()<.08:return False
        pending=reservoir(found['points'],proof['points'],32,state['reservoir_offset'])
        if not pending:return False
        state['reservoir_offset']+=16
        extra,stats=refine_observed({'ainvs':found['ainvs'],'points':proof['points']+pending,
            'basis_size':len(proof['points'])},min(.3,max(.03,left()*.12)),cache,8)
        state['refinements'].append(stats)
        witnesses={tuple(clean({'ainvs':found['ainvs'],'points':[s['point']]})['points'][0]):s
                   for s in stats['steps']}
        grew,_=incorporate([{'point':p,'method':'relation_division','witness':witnesses[tuple(p)]}
                           for p in extra['points']])
        return grew
    def prepare_basis():
        preparing=time.perf_counter();sig=signature();basis={'ainvs':found['ainvs'],'points':proof['points']}
        count=max(anchors,len(proof['points']))
        pool={**basis,'approximate_heights':[point_height(p) for p in basis['points']],
              'vectors':[[int(i==j) for j in range(len(basis['points']))] for i in range(len(basis['points']))]}
        with tempfile.TemporaryDirectory(prefix='preparation-',dir=output) as temp:
            folder=Path(temp);save(folder/'basis.json',basis)
            try:
                pool=generate(folder/'basis.json',folder/'pool',count,count,
                    min(1,max(.03,left()*.2)),selector=diverse_vectors)
            except subprocess.TimeoutExpired:
                state['pool_timeouts']=state.get('pool_timeouts',0)+1
            pool=translate_pool(pool,2*count)
            models=prepare_models(pool,min(.8,max(.03,left()*.18)),cache)
            if not models and left()>.08:
                models=prepare_models(pool,min(.3,left()*.2),cache,minimal=False)
            register(order_models(models),'pointed',refresh=state.get('pool_basis')!=sig)
            if 'seed_covers' not in state and left()>.1:
                state['seed_covers']=True
                if torsion_points(tuple(found['ainvs'])):
                    cover=prepare_covers(basis,min(.15,left()*.08),cache,4)
                    register(cover,'cover')
        state['basis_prepared']=sig if len(models)>=len(pool['points']) else None
        state['pool_basis']=sig
        state['preparations'].append({'generation':state['generation'],'kind':'basis',
            'seconds':time.perf_counter()-preparing,'basis_lower_bound':proof['LowerBound'],
            'basis_points':basis['points'],'pool_points':pool['points'],'pool_vectors':pool.get('vectors'),
            'models':len(models)})
    def enrich():
        """New geometry is derived from the shared state only after a stall."""
        preparing=time.perf_counter();sig=signature();key=sig+':'+str(len(found['points']))
        attempt=state['expansions'].setdefault(key,{'calls':0,'sweep':-1})
        if attempt['sweep']==state['sweep'] or left()<.12:return False
        attempt.update(calls=attempt['calls']+1,sweep=state['sweep']);added=0
        pending=reservoir(found['points'],proof['points'],8,state['reservoir_offset'])
        if pending and left()>.1:
            state['reservoir_offset']+=4
            pool={'ainvs':found['ainvs'],'points':pending,
                  'approximate_heights':[point_height(p) for p in pending]}
            added+=register(order_models(prepare_models(pool,min(.2,left()*.08),cache)),'observed',refresh=False)
        if left()>.12 and torsion_points(tuple(found['ainvs'])):
            basis={'ainvs':found['ainvs'],'points':proof['points']}
            added+=register(prepare_covers(basis,min(.15,left()*.07),cache,4),'cover',refresh=False)
            if proof['LowerBound']<=3 and charts and min(m['coefficient_bits'] for m in charts.values())>80:
                with tempfile.TemporaryDirectory(prefix='isogeny-',dir=output) as temp:
                    models,records=transported_models(basis,Path(temp),min(8,anchors),
                        min(.35,max(.04,left()*.12)),cache)
                    added+=register(models,'dual_isogeny',refresh=False)
                    state.setdefault('isogeny_preparations',[]).append(records)
        if left()>.12 and charts:
            completed={json.loads(k)[0] for k,v in cache.get('local_models',{}).items() if v['complete']}
            usage=state.setdefault('local_root_attempts',{})
            seeds=sorted((m for m in charts.values() if not m.get('neighbour_origin') and m['key'] not in completed),
                key=lambda m:(usage.get(m['key'],0),m['coefficient_bits'],m['key']))
            if seeds:
                k=seeds[0]['key'];tries=usage.get(k,0);usage[k]=tries+1
                added+=register(expand_models(seeds,min(.35*2**min(tries,2),left()*.2),cache,
                    count=6,roots=1,beam=8,depth=32),'local_neighbour',refresh=False)
        state['preparations'].append({'generation':state['generation'],'kind':'enrichment',
            'seconds':time.perf_counter()-preparing,'models_added':added})
        return bool(added)

    search_boxes=[]
    for area in sorted({n*d for n,d in boxes()}):
        group=sorted((b for b in boxes() if b[0]*b[1]==area),key=lambda b:-b[1])
        search_boxes+=group[:1]+[b for b in group[1:] if b[1]==1]+[b for b in group[1:] if b[1]!=1]
    print(json.dumps({'run':output.name,'initial_lower_bound':proof['LowerBound'],'budget':seconds,
                      'algorithm':'unified'}),flush=True)
    checkpoint()
    try:
        if proof['LowerBound']<target and 'seed_division' not in state:
            refined,stats=divide_basis({'ainvs':found['ainvs'],'points':proof['points']},min(.2,max(.01,left()*.1)))
            state['seed_division']=stats
            if stats['steps']:
                incorporate([{'point':p,'method':'seed_division'} for p in refined['points']],recertify=False)
                save(output/'basis.json',refined);candidate=certify(output/'basis.json')
                if candidate['all_selected_independent'] and candidate['LowerBound']>=proof['LowerBound']:
                    proof=candidate
                else:state['seed_division']['certificate_inconclusive']=True
        with ThreadPoolExecutor(max_workers=workers) as executor:
            while left()>.05 and proof['LowerBound']<target:
                state['sweep']+=1
                if state['basis_prepared']!=signature():prepare_basis()
                if not active:
                    enrich()
                    if not active:status='model_preparation_incomplete';break
                if 'initial_features' not in state:
                    bits=sorted(m['coefficient_bits'] for m in active)
                    state['initial_features']={'quartic_bits_min':bits[0],
                        'quartic_bits_p10':bits[len(bits)//10],'model_count':len(active)}
                grew=False;worked=False;enriched=False
                # Bounded working set, with retained charts rejoining on later sweeps.
                width=max(anchors,len(proof['points']))+4
                sweep=active[:width];active=active[width:]+sweep if len(active)>width else active
                for bi,(n,d) in enumerate(search_boxes):
                    pending=[m for m in sweep if missing_intervals(coverage.get(m['key'],[]),n,d)]
                    for offset in range(0,len(pending),32):
                        if left()<=.05:break
                        group=pending[offset:offset+32]
                        chunks=[group[i:i+batch_size] for i in range(0,len(group),batch_size)]
                        futures=[]
                        for chunk in chunks:
                            retry=max(state['attempts'].get(job_key(m,n,d),0) for m in chunk)
                            waves=(len(chunks)+workers-1)//workers
                            limit=min(max(.02,left()/waves),job_seconds*len(chunk),.12*len(chunk)*2**min(retry,6))
                            futures.append(executor.submit(batch_search,found,chunk,n,d,limit,coverage))
                        observed=[];worked=True
                        for chunk,future in zip(chunks,futures):
                            result=future.result();observed+=result['observations']
                            state['attempted_models']+=len(chunk);state['finished_models']+=len(result['complete'])
                            state['timed_out_batches']+=int(result['timed_out']);state['finished_slices']+=len(result['slices'])
                            for m in chunk:
                                key=job_key(m,n,d)
                                if key not in result['complete']:state['attempts'][key]=state['attempts'].get(key,0)+1
                            for key,nn,lo,hi in result['slices']:
                                coverage[key]=record_coverage(coverage.get(key,[]),nn,lo,hi)
                        grew,new=incorporate(observed)
                        if new and not grew and left()>.08:grew=refine()
                        if grew:checkpoint();break
                        if time.perf_counter()-last_save>5:checkpoint()
                        if time.perf_counter()-last_progress>30:
                            print(json.dumps({'run':output.name,'lower_bound':proof['LowerBound'],
                                'seconds':round(before+time.perf_counter()-started,2),
                                'observed_points':len(found['points']),'retained_charts':len(charts)}),flush=True)
                            last_progress=time.perf_counter()
                    if grew or left()<=.05:break
                    # Try alternative exact charts before entering vastly larger boxes.
                    if n*d>=2**28 and (bi+1==len(search_boxes) or search_boxes[bi+1][0]*search_boxes[bi+1][1]>n*d):
                        if refine():grew=True;break
                        if enrich():enriched=True;break
                if not worked and not grew and not enriched:
                    if left()<=.05:break
                    if state['basis_prepared'] is None:continue
                    if any(missing_intervals(coverage.get(m['key'],[]),n,d)
                           for m in active for n,d in search_boxes):continue
                    if not enrich():status='frontier_exhausted';break
            if proof['LowerBound']>=target:status='target_reached'
            elif status=='running':status='budget_completed'
    except KeyboardInterrupt:status='interrupted'
    except Exception as error:
        status='error';state['error']=str(error);raise
    finally:checkpoint()
    print(json.dumps({'run':output.name,'status':status,'lower_bound':proof['LowerBound'],
                      'seconds':state['wall_seconds']}),flush=True)
    return {'ainvs':found['ainvs'],'points':proof['points'],'rank_lower_bound':proof['LowerBound']}
