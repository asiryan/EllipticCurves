"""Batch the preserved point-search formulas; keep checkpoints, not job files."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import time

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'RankHunt'))
from blind_search import (GP, ROOT, POLICY, boxes, clean,
                         diverse_vectors, order_models)
from bootstrap import save, checked_script, gp_process
from models import prepare_models, parity_vectors, parity_order, record_coverage, search_script
from torsion_certificate import certify
from point_models import divide_basis,prepare_covers,isogenous_sources


# Isolated module namespace: preserve the archived engine and other importers.
_pool_spec=importlib.util.spec_from_file_location('equation_anchor_pool',ROOT/'tools/RankHunt/bounded_anchor_pool.py')
_pool=importlib.util.module_from_spec(_pool_spec);_pool_spec.loader.exec_module(_pool)


def normalize_lattice_output(stdout):
    # GP writes tiny approximate heights as "1.23 E-95", which is not JSON.
    # Only the numeric lattice record is touched; exact point records are not.
    return '\n'.join(re.sub(r'(?<=[\d.])\s+[Ee]\s*([+-]?)\s*(\d+)',r'e\1\2',line)
                     if line.startswith('LATTICE ') else line for line in stdout.splitlines())


def _pool_gp(script,prefix,timeout):
    p=gp_process(script,timeout)
    if p.returncode or 'POOL_END' not in p.stdout or any('***' in l and 'Warning:' not in l for l in p.stderr.splitlines()):
        raise RuntimeError('Bounded anchor preparation failed: '+p.stderr)
    return normalize_lattice_output(p.stdout)


_pool.call_gp=_pool_gp
generate=_pool.generate


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def job_key(model, n, d):
    return f'{model["key"]}:{n}:{d}'


def transported_models(data, folder, count, timeout, cache):
    """Search neighbours using the same certified input span, then map to E."""
    from torsion_certificate import translate_pool
    deadline=time.perf_counter()+timeout;result=[];records=[]
    sources=isogenous_sources(data,min(.12,timeout),cache)
    for i,source in enumerate(sources):
        if deadline-time.perf_counter()<.05:break
        basis,division=divide_basis(source,min(.08,(deadline-time.perf_counter())*.2))
        path=folder/f'isogeny-{i}.json';save(path,basis)
        try:
            pool=generate(path,folder/f'isogeny-pool-{i}',count,count,
                          max(.05,deadline-time.perf_counter()),selector=diverse_vectors)
        except subprocess.TimeoutExpired:
            records.append({'alpha':source['transport_alpha'],'status':'pool_budget_exhausted'})
            continue
        pool=translate_pool(pool,count)
        models=order_models(prepare_models(pool,max(.05,(deadline-time.perf_counter())*.65),cache))
        covers=[]
        if deadline-time.perf_counter()>.03:
            covers=prepare_covers(basis,min(.1,deadline-time.perf_counter()),cache,4)
        local=[]
        for k in range(max(len(models),len(covers))):
            if k<len(models):local.append(models[k])
            if k<len(covers):local.append(covers[k])
        for model in local[:count]:
            model=dict(model,search_ainvs=source['ainvs'],transport_alpha=source['transport_alpha'],
                       transport_change=source['transport_change'],
                       pool_index=None)
            model['key']=hashlib.sha256(json.dumps([data['ainvs'],model['transport_alpha'],model['key']]).encode()).hexdigest()
            result.append(model)
        records.append({'alpha':source['transport_alpha'],'source':source,'divisions':division,
                        'models':len(local[:count])})
    return result,records


def batch_search(data, models, n, d, timeout, coverage=None):
    script = search_script(data, models, n, d, coverage)
    checked_script(script)
    started = time.perf_counter()
    timed_out = False
    try:
        p = gp_process(script,timeout)
        stdout, stderr, code = p.stdout, p.stderr, p.returncode
    except subprocess.TimeoutExpired as error:
        timed_out = True
        stdout, stderr, code = error.stdout or '', error.stderr or '', None
        if isinstance(stdout, bytes): stdout = stdout.decode(errors='replace')
        if isinstance(stderr, bytes): stderr = stderr.decode(errors='replace')
    if code not in (0, None) or any('***' in s and 'Warning:' not in s for s in stderr.splitlines()):
        raise RuntimeError('PARI point search failed: ' + stderr[-4000:])
    complete, observations, slices = [], [], []
    current = None
    for line in stdout.splitlines():
        if line.startswith('ANCHOR_BEGIN '):
            current = int(line.split()[1]) - 1
        elif line.startswith('ANCHOR_DONE '):
            complete.append(job_key(models[int(line.split()[1]) - 1], n, d))
        elif line.startswith('SLICE_DONE '):
            i,nn,lo,hi = json.loads(line[11:])
            if not (1<=i<=len(models) and nn==n and 1<=lo<=hi<=d):
                raise ValueError('Malformed completed search slice')
            slices.append((models[i-1]['key'],nn,lo,hi))
        elif line.startswith('POINT '):
            match = re.fullmatch(r'POINT \[(-?\d+(?:/\d+)?), (-?\d+(?:/\d+)?)\]', line)
            if match is None or current is None: raise ValueError('Malformed point output')
            point = clean({'ainvs':data['ainvs'], 'points':[match.groups()]})['points'][0]
            model=models[current]
            observations.append({'point':point, 'anchor':model['anchor'],
                                 'n':n, 'd':d, 'model':model['key'],
                                 'method':('dual_isogeny' if 'transport_alpha' in model else model.get('kind','pointed')),
                                 **({'transport_alpha':model['transport_alpha'],'anchor_ainvs':model['search_ainvs']}
                                    if 'transport_alpha' in model else {})})
    if not timed_out and ('SEARCH_END' not in stdout or len(complete) != len(models)):
        raise RuntimeError('Incomplete successful PARI batch')
    return {'complete':complete, 'observations':observations, 'slices':slices, 'timed_out':timed_out,
            'seconds':time.perf_counter() - started}


def independent_result(data, expected):
    path = Path(__file__).with_name('certificate.py')
    if hashlib.sha256(path.read_bytes()).hexdigest() != '9e0d0d2562fc53705e92a2eaa9a3f6e7c923f1cd3fd82b14df68b60268f4ad54':
        raise ValueError('Independent verifier changed')
    spec = importlib.util.spec_from_file_location('independent_certificate', path)
    verifier = importlib.util.module_from_spec(spec); spec.loader.exec_module(verifier)
    from torsion_certificate import torsion_points,select_basis,verify_certificate,singleton_certificate
    torsion=torsion_points(tuple(data['ainvs']))
    if torsion:
        result=select_basis(data,torsion,max_prime=2000)
        if result is None and expected==1:
            certificate=singleton_certificate(data)
        elif result is None or result['points']!=data['points']:
            raise ValueError('Independent verification did not retain the entire selected basis')
        else:certificate=result['certificate']
        claim=verify_certificate(data,certificate)
    else:
        certificate = verifier.build_certificate(data, max_prime=2000)
        claim = verifier.verify_certificate(data, certificate)
    if claim['rank_lower_bound'] != expected or not claim['all_points_independent_modulo_torsion']:
        raise ValueError('Independent verification did not confirm the selected basis')
    return {**data, 'rank_lower_bound':expected, 'certificate':certificate, 'verification':claim}


def search(source, output, seconds=300, workers=4, anchors=2048, target=32,
           batch_size=8, job_seconds=2, import_run=None, anchor_mode='adaptive'):
    if anchor_mode not in ('adaptive','fixed','frozen','parity','geometric'): raise ValueError('Unknown anchor policy')
    source, output = Path(source).resolve(), Path(output).resolve()
    if output == ROOT or not output.is_relative_to(ROOT): raise ValueError('Use a dedicated workspace directory')
    output.mkdir(parents=True, exist_ok=True)
    config = {'input_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
              'code_sha256':{name:hashlib.sha256(((Path(__file__).parent/name) if name in ('seeded.py','models.py','bootstrap.py','certificate.py','torsion_certificate.py','point_models.py') else (ROOT/'tools/RankHunt'/name)).read_bytes()).hexdigest()
                  for name in ['seeded.py','models.py','bootstrap.py','certificate.py','torsion_certificate.py','point_models.py','blind_search.py','point_search.py','bounded_anchor_pool.py',
                               'anchor_diversity.py','point_arithmetic.py']},
              'policy':POLICY, 'anchors':anchors, 'target':target, 'batch_size':batch_size,
              'job_seconds':job_seconds, 'workers':workers,'anchor_mode':anchor_mode,
              'search_policy':{'cached_inverse_maps':True,'denominator_slices':256,
                               'parity_balanced':anchor_mode=='parity',
                               'profile_budget_fraction':0.25,'profile_budget_max_seconds':2,
                               'small_point_division':[2,3,5] if anchor_mode=='geometric' else [],'division_depth':3,
                               'isogeny_cover_divisor_beam':256,'isogeny_cover_limit':4 if anchor_mode=='geometric' else 0,
                               'isogenous_neighbours':anchor_mode=='geometric',
                               'isogenous_trigger_quartic_bits':80,'isogenous_models_limit':8,
                               'isogenous_max_initial_lower_bound':3,'isogenous_model_slice_seconds':.12,
                               'box_order':'balanced, integer, remaining denominators' if anchor_mode=='geometric' else 'original'}}
    state_path = output/'checkpoint.json'
    if state_path.exists():
        state = read(state_path)
        if state['config'] != config: raise ValueError('Resume requires unchanged input, code and settings')
        found = state['found']; done = set(state['done']); models = state['models']
    else:
        found = clean(read(source)); done = set(); models = None
        state = {'config':config, 'source':str(source), 'generation':0, 'wall_seconds':0,
                 'attempted_models':0, 'finished_models':0, 'timed_out_batches':0,
                 'events':[], 'origin':'seeded search', 'imported_seconds':0}
        if import_run:
            legacy = Path(import_run).resolve(); old_config = read(legacy/'config.json')
            if old_config['input_sha256'] != config['input_sha256'] or old_config['anchors'] != anchors:
                raise ValueError('Import input or anchor policy differs')
            if old_config['policy'] != POLICY: raise ValueError('Import policy differs')
            old = read(legacy/'state.json'); found = clean(old['found']); done = set(old['done'])
            state.update(generation=old['generation'], imported_seconds=old['wall_seconds'],
                         origin='continue previous point search', imported_run=str(legacy))
            model_path = legacy/f'generation-{old["generation"]:03d}'/'models.json'
            # Old profiling files do not contain inverse maps; rebuild them.
            models = None
        if not found['points']: raise ValueError('Known seed points are required')
    save(output/'basis.json', {'ainvs':found['ainvs'],'points':state.get('certified_basis',found['points'])})
    proof = certify(output/'basis.json')
    if not proof['all_selected_independent'] or not proof['points']: raise ValueError('No independent seed basis')
    known = set(map(tuple, found['points']))
    if 'initial_lower_bound' not in state: state['initial_lower_bound'] = proof['LowerBound']
    initial_basis=clean(read(source))
    initial_points=set(map(tuple,initial_basis['points']))
    coverage = state.setdefault('coverage',{})
    model_cache = state.setdefault('model_cache',{})
    state.setdefault('finished_slices',0)
    state.setdefault('preparations',[])
    started = time.perf_counter(); deadline = started + seconds; before = state['wall_seconds']
    last_save = started; last_progress = started; status = 'running'
    def checkpoint():
        nonlocal last_save
        state.update(found=found, done=sorted(done), models=models, status=status,
                     certified_basis=proof['points'],
                     lower_bound=proof['LowerBound'], wall_seconds=before+time.perf_counter()-started)
        save(output/'checkpoint.json', state)
        save(output/'points.json', {'ainvs':found['ainvs'], 'points':proof['points'],
                                   'rank_lower_bound':proof['LowerBound']})
        last_save = time.perf_counter()
    checkpoint()
    print(json.dumps({'run':output.name, 'initial_lower_bound':proof['LowerBound'], 'budget':seconds}), flush=True)
    try:
        if anchor_mode=='geometric' and proof['LowerBound']<target and 'seed_division' not in state:
            refined,division=divide_basis({'ainvs':found['ainvs'],'points':proof['points']},
                                         min(.2,max(.01,(deadline-time.perf_counter())*.1)))
            state['seed_division']=division
            if division['steps']:
                save(output/'basis.json',refined);candidate=certify(output/'basis.json')
                if candidate['all_selected_independent'] and candidate['LowerBound']>=proof['LowerBound']:
                    proof=candidate;models=None
                    for p in proof['points']:
                        if tuple(p) not in known:found['points'].append(p);known.add(tuple(p))
                else:state['seed_division']['certificate_inconclusive']=True
            checkpoint()
        with ThreadPoolExecutor(max_workers=workers) as executor:
            while time.perf_counter() < deadline and proof['LowerBound'] < target:
                if models is None:
                    preparing=time.perf_counter()
                    save(output/'basis.json', {'ainvs':found['ainvs'], 'points':proof['points']})
                    with tempfile.TemporaryDirectory(prefix='preparation-', dir=output) as temp:
                        folder = Path(temp)
                        pool_basis=initial_basis if anchor_mode in ('fixed','frozen') else {'ainvs':found['ainvs'],'points':proof['points']}
                        if anchor_mode=='fixed':
                            pool={'ainvs':initial_basis['ainvs'],'points':initial_basis['points'],
                                  'approximate_heights':[0]*len(initial_basis['points'])}
                        else:
                            pool_source=folder/'basis.json';save(pool_source,pool_basis)
                            pool = generate(pool_source, folder/'pool', anchors, anchors,
                                            min(90,max(1,deadline-time.perf_counter())),
                                            selector=parity_vectors if anchor_mode=='parity' else diverse_vectors)
                            from torsion_certificate import translate_pool
                            pool=translate_pool(pool,anchors)
                        # A difficult later anchor must not consume the entire
                        # search budget after useful earlier models are ready.
                        profile_budget=min(2,max(.05,(deadline-time.perf_counter())*.25))
                        models = prepare_models(pool,profile_budget,model_cache)
                        if not models and deadline-time.perf_counter()>.1:
                            # Global minimal-model factorization can stall even
                            # though reduction without factorization is cheap.
                            fallback_budget=min(1,(deadline-time.perf_counter())*.25)
                            models=prepare_models(pool,fallback_budget,model_cache,minimal=False)
                        models = parity_order(models) if anchor_mode=='parity' else order_models(models)
                        if anchor_mode=='geometric' and anchors>=4 and deadline-time.perf_counter()>.1:
                            from torsion_certificate import torsion_points
                            if torsion_points(tuple(found['ainvs'])):
                                covers=prepare_covers({'ainvs':found['ainvs'],'points':proof['points']},
                                    min(.3,(deadline-time.perf_counter())*.1),model_cache,min(4,anchors//4))
                                neighbours=[]
                                if (proof['LowerBound']<=3 and models and
                                        min(m['coefficient_bits'] for m in models)>80 and deadline-time.perf_counter()>.2):
                                    neighbours,records=transported_models(pool_basis,folder,min(8,anchors//2),
                                        min(.4,(deadline-time.perf_counter())*.12),model_cache)
                                    state.setdefault('isogeny_preparations',[]).append(records)
                                # Distinct projection geometries share every new certified basis.
                                interleaved=[]
                                for i in range(max(len(models),len(covers),len(neighbours))):
                                    if i<len(models):interleaved.append(models[i])
                                    if i<len(neighbours):interleaved.append(neighbours[i])
                                    if i<len(covers):interleaved.append(covers[i])
                                models=interleaved[:anchors]
                    # Coefficients refer to an exactly checked independent basis.
                    # Mark only provable use of directions outside the initial span.
                    retained=initial_points<=set(map(tuple,pool_basis['points']))
                    extra=[i for i,p in enumerate(pool_basis['points']) if tuple(p) not in initial_points]
                    for model in models:
                        model['preparation']=len(state['preparations'])
                        model['uses_new_direction']=(None if model.get('kind')=='isogeny_cover' or 'transport_alpha' in model else
                            False if anchor_mode in ('fixed','frozen') else
                            any(pool['vectors'][model['pool_index']][i] for i in extra) if retained else None)
                    state['preparations'].append({'generation':state['generation'],
                        'basis_lower_bound':proof['LowerBound'],'seconds':time.perf_counter()-preparing,
                        'pool_points':pool['points'],'pool_vectors':pool.get('vectors'),
                        'torsion_translations':pool.get('torsion_translations'),
                        'basis_points':pool_basis['points'],
                        'models':[{'key':m['key'],'pool_index':m['pool_index'],
                                   'kind':m.get('kind','pointed'),
                                   'transport_alpha':m.get('transport_alpha'),
                                   'coefficient_bits':m['coefficient_bits'],
                                   'uses_new_direction':m['uses_new_direction']} for m in models]})
                    if not models:
                        models=None
                        status='model_preparation_incomplete'
                        break
                    checkpoint()
                if 'initial_features' not in state:
                    bits=sorted(m['coefficient_bits'] for m in models)
                    state['initial_features']={'quartic_bits_min':bits[0],
                        'quartic_bits_p10':bits[len(bits)//10],
                        'anchor_height_min':min(m['anchor_height'] for m in models),
                        'parity_classes':len({m['parity'] for m in models if 'parity' in m}),
                        'model_count':len(models)}
                improved = False
                search_boxes=list(boxes())
                if anchor_mode=='geometric':
                    areas=sorted({n*d for n,d in search_boxes});ordered=[]
                    for area in areas:
                        group=sorted((b for b in search_boxes if b[0]*b[1]==area),key=lambda b:-b[1])
                        ordered+=group[:1]+[b for b in group[1:] if b[1]==1]+[b for b in group[1:] if b[1]!=1]
                    search_boxes=ordered
                for n,d in search_boxes:
                    pending = [m for m in models if job_key(m,n,d) not in done]
                    for offset in range(0, len(pending), 32):
                        remaining = deadline-time.perf_counter()
                        if remaining <= .05: break
                        group = pending[offset:offset+32]
                        # An expensive first projection must not hide all the
                        # other geometries behind it in the same GP process.
                        transporting=any('transport_alpha' in m for m in models)
                        chunk_size=1 if transporting else batch_size
                        chunks = [group[i:i+chunk_size] for i in range(0,len(group),chunk_size)]
                        limit = min(job_seconds*batch_size, max(.05,remaining/((len(chunks)+workers-1)//workers)))
                        if transporting:limit=min(limit,.12)
                        futures = [executor.submit(batch_search,found,chunk,n,d,limit,coverage) for chunk in chunks]
                        observations = []; old_bound = proof['LowerBound']; new = []
                        for chunk,future in zip(chunks,futures):
                            result = future.result()
                            by_key={m['key']:m for m in chunk}
                            for observation in result['observations']:
                                model=by_key[observation['model']]
                                observation.update(preparation=model['preparation'],
                                                   uses_new_direction=model['uses_new_direction'])
                            observations.extend(result['observations'])
                            done.update(result['complete'])
                            state['attempted_models'] += len(chunk)
                            state['finished_models'] += len(result['complete'])
                            state['timed_out_batches'] += int(result['timed_out'])
                            for key,nn,lo,hi in result['slices']:
                                coverage[key]=record_coverage(coverage.get(key,[]),nn,lo,hi)
                            state['finished_slices'] += len(result['slices'])
                        for observation in observations:
                            point = tuple(observation['point'])
                            if point not in known:
                                known.add(point); found['points'].append(list(point)); new.append(observation)
                        if new:
                            old_points=set(map(tuple,proof['points']))
                            save(output/'basis.json', {'ainvs':found['ainvs'],
                                'points':proof['points']+[o['point'] for o in new]})
                            candidate = certify(output/'basis.json')
                            if candidate['LowerBound'] >= old_bound and candidate['all_selected_independent']:
                                proof = candidate
                            else:
                                # An inconclusive certificate cannot erase an
                                # earlier exact proof or imply Q-dependence.
                                state['inconclusive_updates']=state.get('inconclusive_updates',0)+1
                            # Only successful discoveries need permanent provenance.
                            with (output/'discoveries.jsonl').open('a',encoding='utf-8') as stream:
                                for observation in new: stream.write(json.dumps(observation)+'\n')
                        if proof['LowerBound'] > old_bound:
                            event = {'run':output.name, 'lower_bound':proof['LowerBound'],
                                     'seconds':before+time.perf_counter()-started,
                                     'attempted_models':state['attempted_models'],
                                     'selected_new_points':[o for o in new if tuple(o['point']) not in old_points
                                         and o['point'] in proof['points']]}
                            state['events'].append(event)
                            print(json.dumps({k:v for k,v in event.items() if k!='selected_new_points'}),flush=True)
                            state['generation'] += 1
                            if anchor_mode not in ('fixed','frozen'): models = None
                            improved = True; checkpoint(); break
                        if time.perf_counter()-last_save >= 5: checkpoint()
                        if time.perf_counter()-last_progress >= 30:
                            print(json.dumps({'run':output.name, 'lower_bound':proof['LowerBound'],
                                'seconds':round(before+time.perf_counter()-started,2),
                                'finished_models':state['finished_models']}),flush=True)
                            last_progress=time.perf_counter()
                    if improved or time.perf_counter() >= deadline-.05: break
                if not improved:
                    status = 'budget_completed' if time.perf_counter() >= deadline-.05 else 'pass_completed'
                    break
            if proof['LowerBound'] >= target: status='target_reached'
            elif status == 'running': status='budget_completed'
    except KeyboardInterrupt:
        status='interrupted'
    except Exception as error:
        status='error'; state['error']=str(error); raise
    finally:
        checkpoint()
    print(json.dumps({'run':output.name,'status':status,'lower_bound':proof['LowerBound'],
                      'seconds':state['wall_seconds']}),flush=True)
    return {'ainvs':found['ainvs'],'points':proof['points'],'rank_lower_bound':proof['LowerBound']}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description='Expand supplied points; an existing checkpoint resumes completed slices.')
    parser.add_argument('--input',required=True)
    parser.add_argument('--output',required=True)
    parser.add_argument('--seconds',type=float,default=10)
    parser.add_argument('--workers',type=int,default=2)
    parser.add_argument('--anchors',type=int,default=64)
    parser.add_argument('--target',type=int,default=32)
    parser.add_argument('--anchor-mode',choices=('adaptive','fixed','frozen','parity','geometric'),default='adaptive',
                        help='fixed: supplied anchors; frozen: their initial generated pool; adaptive/parity: rebuild after growth')
    args=parser.parse_args()
    if not (0<args.seconds<=7200 and 1<=args.workers<=24 and 1<=args.anchors<=4096 and 1<=args.target<=100):
        parser.error('Invalid bounded search settings')
    result=search(args.input,args.output,args.seconds,args.workers,args.anchors,args.target,
                  batch_size=4,anchor_mode=args.anchor_mode)
    save(Path(args.output)/'result.json',independent_result(result,result['rank_lower_bound']))
