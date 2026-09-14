"""Bounded search for new high-rank curves using the preserved point engine.

Select generated parameters by finite-field scores, then allocate increasing
point-search budgets by certified lower bounds. No score is a rank certificate.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path
import random
import shutil
import subprocess
import sys
import threading
import time

from blind_search import clean, certify
from point_search import ROOT, save as atomic_save


POLICY={'version':2,'candidate_order':'interleave full score, tail score, fixed shuffle',
        'promotion_order':'certified lower bound, tail score, full score, id',
        'random_seed':20260914,'full_rank_computation':False,
        'novelty':'exclude duplicate j and the known record j; global novelty not asserted'}
RECORD_A=['1','1','1',
    '-1284727764113567728281797636015784768866707681415849262157224232063',
    '560368321454261339256859338901915312332769858684945406858043869199456710681989058863306170127006181']


def read(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path,data):
    for attempt in range(6):
        try:return atomic_save(path,data)
        except PermissionError:
            if attempt==5:raise
            time.sleep(.01*2**attempt)


def j_invariant(data):
    a1,a2,a3,a4,a6=map(Q,data['ainvs'])
    b2=a1*a1+4*a2;b4=a1*a3+2*a4;b6=a3*a3+4*a6
    b8=a1*a1*a6+4*a2*a6-a1*a3*a4+a2*a3*a3-a4*a4
    disc=-b2*b2*b8-8*b4**3-27*b6*b6+9*b2*b4*b6
    if not disc:raise ValueError('Singular candidate')
    return (b2*b2-24*b4)**3/disc


def candidate_order(rows,count):
    by_score=sorted(rows,key=lambda r:(-r['score'],r['id']))
    by_tail=sorted(rows,key=lambda r:(-selection_tail(r),r['id']))
    shuffled=sorted(rows,key=lambda r:r['id']);random.Random(POLICY['random_seed']).shuffle(shuffled)
    result=[];seen=set()
    for group in zip(by_score,by_tail,shuffled):
        for row in group:
            if row['id'] not in seen:seen.add(row['id']);result.append(row)
    return result[:count]


def selection_tail(row):
    return row.get('confirmation_score',row['tail_score'])


def arithmetic_score(data,start,end,timeout=30):
    """Finite-field point counts only; no rank or L-function computation."""
    from blind_search import GP
    from point_search import vector
    if not 3<=start<end<=262139:raise ValueError('Invalid scoring interval')
    script=('default(realprecision,38);\nE=ellinit('+vector(data['ainvs'])+');\n'
        f's=0;forprime(p={start+1},{end},if(E.disc%p,s+=log((p+1-ellap(E,p))/p)));print(s);quit;\n')
    process=subprocess.run([str(GP),'-fq'],input=script,capture_output=True,text=True,timeout=timeout)
    if process.returncode or '***' in process.stderr:raise RuntimeError('Finite-field scoring failed: '+process.stderr)
    result=float(process.stdout.strip())
    if not math.isfinite(result):raise ValueError('Nonfinite score')
    return result


def promoted(rows,count):
    eligible=[r for r in rows if r.get('status') in {'budget_completed','target_reached'}]
    return sorted(eligible,key=lambda r:(-r['lower_bound'],-r['tail_score'],-r['score'],r['id']))[:count]


def adaptive_promoted(rows,count):
    """Favor certified growth while reserving slots for slower-starting curves."""
    eligible=[r for r in rows if r.get('status') not in {'error','pending'}]
    if not eligible:return []
    exploitation=sorted(eligible,key=lambda r:(-r['lower_bound'],-r.get('recent_gain_per_second',0),
        -selection_tail(r),r['id']))
    explore=0 if count==1 else max(1,count//3)
    chosen=exploitation[:max(1,count-explore)];seen={r['id'] for r in chosen}
    tail=sorted(eligible,key=lambda r:(-selection_tail(r),-r['score'],r['id']))
    cost=sorted(eligible,key=lambda r:(r.get('quartic_bits_p10',float('inf')),
        -r.get('gain_per_second',0),-selection_tail(r),r['id']))
    for pair in zip(tail,cost):
        for row in pair:
            if len(chosen)>=count:return chosen
            if row['id'] not in seen:chosen.append(row);seen.add(row['id'])
    return chosen


def prepare_multiscale(output,samples,heights,workers,seed):
    """One pooled input, with separate bounded samples at several heights."""
    output=Path(output).resolve()
    if output==ROOT or not output.is_relative_to(ROOT):raise ValueError('Use a dedicated workspace directory')
    output.mkdir(parents=True,exist_ok=True)
    dll=ROOT/'tools/RankHunt/bin/Release/net8.0/RankHunt.dll'
    config={'samples':samples,'heights':heights,'workers':workers,'seed':seed,'dll_sha256':digest(dll)}
    checkpoint=output/'generation.json'
    if checkpoint.exists():
        old=read(checkpoint)
        if old['config']!=config:raise ValueError('Generation settings changed')
        if old['status']=='completed':return
    save(checkpoint,{'config':config,'status':'running'})
    def one(index,height):
        count=samples//len(heights)+int(index<samples%len(heights))
        folder=output/f'h{height}'
        command=['dotnet',str(dll),'sample','--selection','bands','--samples',str(count),
            '--height',str(height),'--keep','256','--refine-keep','256','--final-keep','48',
            '--prime-bound','65521','--workers',str(workers),'--seed',str(seed+index),'--output',str(folder)]
        process=subprocess.run(command,capture_output=True,text=True,timeout=1200)
        if process.returncode:raise RuntimeError(process.stdout[-2000:]+'\n'+process.stderr)
        meta=read(folder/'complete.json')
        print(json.dumps({'generated_height':height,'primitive_draws':meta['PrimitiveDraws'],
                          'finalists':meta['exported'],'seconds':round(meta['seconds'],2)}),flush=True)
        return height,read(folder/'candidates.json'),meta
    results=[]
    with ThreadPoolExecutor(max_workers=2) as executor:
        futures=[executor.submit(one,i,h) for i,h in enumerate(heights)]
        for future in as_completed(futures):results.append(future.result())
    seen={str(j_invariant({'ainvs':RECORD_A}))};rows=[]
    for height,candidates,meta in sorted(results):
        for row in candidates:
            if row['j_invariant'] in seen:continue
            seen.add(row['j_invariant'])
            rows.append({**row,'file':f'h{height}/'+row['file'],'sampling_height':height})
    save(output/'candidates.json',rows)
    save(checkpoint,{'config':config,'status':'completed','primitive_draws':sum(r[2]['PrimitiveDraws'] for r in results),
                    'sampling_with_replacement':True,'distinct_finalists':len(rows)})
    print(f'Pooled {len(rows)} different candidate j-invariants.',flush=True)


def run_fast(args):
    from fast_search import search,independent_result
    check_engine()
    pool=Path(args.candidates).resolve();out=Path(args.output).resolve()
    if out==ROOT or not out.is_relative_to(ROOT):raise ValueError('Use a dedicated workspace output directory')
    record_j=j_invariant({'ainvs':RECORD_A});seen={record_j};source=[]
    for index,row in enumerate(read(pool/'candidates.json')):
        path=(pool/row['file']).resolve()
        if not path.is_relative_to(pool):raise ValueError('Candidate path escapes pool')
        data=clean(read(path));j=j_invariant(data)
        if j in seen:continue
        if not all(math.isfinite(float(row[k])) for k in ['score','tail_score']):raise ValueError('Nonfinite score')
        seen.add(j)
        source.append({**row,'id':f'c{index:03d}','input':str(path),'input_sha256':digest(path),
                       'j_invariant':str(j),'seed_lower_bound':None,'lower_bound':0,
                       'search_seconds':0,'status':'pending'})
    if args.confirmation_bound:
        print(f'Checking {len(source)} candidates on primes (65521,{args.confirmation_bound}]',flush=True)
        with ThreadPoolExecutor(max_workers=args.parallel_curves) as executor:
            futures={executor.submit(arithmetic_score,read(r['input']),65521,args.confirmation_bound):r for r in source}
            for future in as_completed(futures):
                futures[future]['confirmation_score']=future.result()
    selected=candidate_order(source,args.counts[0]);out.mkdir(parents=True,exist_ok=True)
    config={'inputs':{r['id']:r['input_sha256'] for r in selected},
        'counts':args.counts,'seconds':args.seconds,'target':args.target,'anchors':args.anchors,
        'workers':args.point_workers,'parallel':args.parallel_curves,'job_seconds':args.job_seconds,
        'batch_size':args.batch_size,
        'confirmation_bound':args.confirmation_bound,
        'controller_sha256':digest(Path(__file__)),'point_engine_sha256':digest(ROOT/'tools/RankHunt/fast_search.py')}
    checkpoint=out/'campaign.json'
    if checkpoint.exists():
        state=read(checkpoint)
        if state['config']!=config:raise ValueError('Campaign settings changed')
        if state['status']=='completed':return state
        rows=state['rows'];stage=state['stage'];active=state['active']
    else:rows=selected;stage=0;active=[r['id'] for r in rows]
    def persist(status='running'):
        save(checkpoint,{'config':config,'rows':rows,'stage':stage,'active':active,'status':status})
    persist()
    def work(row,budget):
        folder=out/row['id'];old=row['lower_bound'];used=row['search_seconds']
        if (folder/'checkpoint.json').exists():
            previous=read(folder/'checkpoint.json');used=previous['wall_seconds'];old=previous['lower_bound']
        if used<budget and old<args.target:
            search(row['input'],folder,max(1,budget-used),args.point_workers,args.anchors,args.target,
                   batch_size=args.batch_size,job_seconds=args.job_seconds)
        result=read(folder/'checkpoint.json');elapsed=result['wall_seconds']
        initial=result['initial_lower_bound'];bound=result['lower_bound']
        return {**row,'seed_lower_bound':initial,'lower_bound':bound,'search_seconds':elapsed,
            'status':result['status'],'gain_per_second':(bound-initial)/max(elapsed,.001),
            'recent_gain_per_second':(bound-max(initial,old))/max(elapsed-used,.001),
            **result.get('initial_features',{})}
    while stage<len(args.seconds) and active:
        print(json.dumps({'stage':stage,'curves':len(active),'cumulative_seconds':args.seconds[stage]}),flush=True)
        with ThreadPoolExecutor(max_workers=args.parallel_curves) as executor:
            futures=[executor.submit(work,r,args.seconds[stage]) for r in rows if r['id'] in active]
            for future in as_completed(futures):
                result=future.result();rows=[result if r['id']==result['id'] else r for r in rows]
                print(json.dumps({'finished_candidate':result['id'],'parameter':f'{result["u"]}/{result["v"]}',
                    'lower_bound':result['lower_bound'],'seconds':round(result['search_seconds'],2),
                    'gain_per_second':round(result['gain_per_second'],4),
                    'quartic_bits_p10':result.get('quartic_bits_p10')}),flush=True)
                persist()
        stage+=1
        if max(r['lower_bound'] for r in rows)>=args.target:break
        active=[r['id'] for r in adaptive_promoted(rows,args.counts[stage])] if stage<len(args.counts) else []
        persist()
    best=max(rows,key=lambda r:(r['lower_bound'],r['tail_score'],r['id']))
    result=independent_result(read(out/best['id']/'points.json'),best['lower_bound'])
    result.update(parameter={'u':best['u'],'v':best['v']},global_novelty_not_established=True)
    save(out/'best.json',result);persist('completed')
    print(json.dumps({'status':'completed','best_parameter':result['parameter'],
                      'independently_verified_lower_bound':best['lower_bound']}),flush=True)
    return result


def check_engine():
    manifest=read(ROOT/'results/point-engine-v3-20260914.manifest.json')
    for name,value in manifest['code_sha256'].items():
        if digest(ROOT/'tools/RankHunt'/name)!=value:raise ValueError('Preserved point engine changed: '+name)
    return manifest['code_sha256']


def worker_command(row,out,seconds,args):
    folder=out/'searches'/row['id']
    command=[sys.executable,str(ROOT/'tools/RankHunt/blind_search.py'),
        '--input',str(out/'inputs'/(row['id']+'.json')),'--output',str(folder),
        '--seconds',str(seconds),'--anchors',str(args.anchors),'--workers',str(args.point_workers),
        '--job-seconds',str(args.job_seconds),'--target',str(args.target)]
    if (folder/'state.json').exists():command.append('--resume')
    return command


def kill_worker(process):
    if process.poll() is not None:return
    if sys.platform=='win32':
        subprocess.run(['taskkill','/PID',str(process.pid),'/T','/F'],capture_output=True,timeout=15)
    else:process.kill()
    try:process.wait(timeout=10)
    except subprocess.TimeoutExpired:process.kill();process.wait(timeout=10)


def search_one(row,out,cumulative_seconds,args,stage,stop_event=None):
    folder=out/'searches'/row['id'];summary_path=folder/'summary.json'
    previous=read(summary_path) if summary_path.exists() else {}
    used=previous.get('wall_seconds',0)
    if used>=cumulative_seconds or previous.get('lower_bound',0)>=args.target:return {**row,**{
        'lower_bound':previous.get('lower_bound',row['lower_bound']),
        'search_seconds':used,'status':previous.get('status',row.get('status'))}}
    seconds=max(1,cumulative_seconds-used)
    command=worker_command(row,out,seconds,args)
    prefix=out/'processes'/f"{row['id']}-stage-{stage:02d}"
    save(prefix.with_suffix('.command.json'),{'argv':command,'cumulative_target_seconds':cumulative_seconds})
    start=time.perf_counter();last_bound=row['lower_bound'];interrupted=False
    with prefix.with_suffix('.stdout.txt').open('w') as stdout,prefix.with_suffix('.stderr.txt').open('w') as stderr:
        process=subprocess.Popen(command,stdout=stdout,stderr=stderr,
            creationflags=subprocess.CREATE_NO_WINDOW if sys.platform=='win32' else 0)
        save(prefix.with_suffix('.command.json'),{'argv':command,'pid':process.pid,
            'cumulative_target_seconds':cumulative_seconds})
        # A normal Windows read handle can block os.replace on summary.json.
        # Follow the append-only process log; read checkpoint JSON after exit.
        log_reader=prefix.with_suffix('.stdout.txt').open('r');pending=''
        try:
            while process.poll() is None:
                if (stop_event is not None and stop_event.is_set()) or time.perf_counter()-start>seconds+90:
                    interrupted=True;kill_worker(process);break
                pending+=log_reader.read()
                complete=pending.split('\n');pending=complete.pop()
                for line in complete:
                    if not line.startswith('{'):continue
                    try:event=json.loads(line)
                    except json.JSONDecodeError:continue
                    if event.get('lower_bound',0)>last_bound:
                        last_bound=event['lower_bound']
                        print(json.dumps({'candidate':row['id'],'lower_bound':last_bound,
                            'search_seconds':round(event['seconds'],3)}),flush=True)
                time.sleep(.25)
        except BaseException:
            kill_worker(process);raise
        finally:log_reader.close()
    result=read(summary_path) if summary_path.exists() else {}
    if interrupted and result:
        result['status']='interrupted';result['controller_timeout']=True;save(summary_path,result)
    status=result.get('status','error')
    if process.returncode and status not in {'error','interrupted'}:status='error'
    return {**row,'lower_bound':result.get('lower_bound',row['lower_bound']),
        'search_seconds':result.get('wall_seconds',used),'status':status,
        'last_process_exit_code':process.returncode,'last_process_wall_seconds':time.perf_counter()-start}


def run(args):
    engine_hashes=check_engine()
    pool=Path(args.candidates).resolve();out=Path(args.output).resolve()
    if out==ROOT or not out.is_relative_to(ROOT):raise ValueError('Use a dedicated output directory in this workspace')
    source=read(pool/'candidates.json');seen=set();ids=set();validated=[]
    record_j=j_invariant({'ainvs':RECORD_A})
    for index,original in enumerate(source):
        if not all(math.isfinite(float(original[k])) for k in ['score','tail_score']):raise ValueError('Nonfinite score')
        file=(pool/original['file']).resolve()
        if not file.is_relative_to(pool):raise ValueError('Candidate path escapes pool')
        data=clean(read(file));j=j_invariant(data)
        if 'j_invariant' in original and j!=Q(original['j_invariant']):raise ValueError('Incorrect candidate j-invariant')
        if j==record_j or j in seen:continue
        seen.add(j)
        row={**original,'id':f'c{index:03d}','source_file':str(file),'input_sha256':digest(file),
             'j_invariant':str(j),'prior_point_search_data_used':False}
        if row['id'] in ids:raise ValueError('Duplicate candidate id')
        ids.add(row['id']);validated.append(row)
    selected=candidate_order(validated,args.counts[0])
    if not selected:raise ValueError('No eligible new candidate curves')
    code_files=['record_hunt.py','CandidateSampler.cs','Hunt.cs','Program.cs','Structured302.cs','Data/icarm302-sections.json','audit_blind_search.py']
    code_hashes={name:digest(ROOT/'tools/RankHunt'/name) for name in code_files}
    config={'policy':POLICY,'engine_sha256':engine_hashes,'controller_sha256':code_hashes,
        'candidate_list_sha256':digest(pool/'candidates.json'),
        'selected_inputs':{r['id']:r['input_sha256'] for r in selected},
        'stages':list(zip(args.counts,args.seconds)),'point_workers':args.point_workers,
        'parallel_curves':args.parallel_curves,'anchors':args.anchors,'job_seconds':args.job_seconds,
        'target':args.target,'previously_studied':args.previously_studied,
        'dll_sha256':digest(ROOT/'tools/RankHunt/bin/Release/net8.0/RankHunt.dll')}
    # Normalize tuples to the same JSON representation used by resume.
    config=json.loads(json.dumps(config))
    if args.resume:
        if read(out/'config.json')!=config:raise ValueError('Resume requires unchanged inputs, code and policy')
        finished=read(out/'summary.json')
        if finished['status'] in {'campaign_completed','target_found'} and (out/'best-verified.json').exists():
            print('Campaign already completed and independently verified.',flush=True);return finished
        state=read(out/'state.json');rows=state['rows'];stage=state['stage'];active=state['active']
        elapsed_before=state['controller_wall_seconds']
    else:
        out.mkdir(parents=True,exist_ok=False)
        for name in ['inputs','searches','processes','code','audits']:(out/name).mkdir()
        save(out/'config.json',config);save(out/'candidate-source.json',source)
        for name in code_files:
            dest=out/'code'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/'tools/RankHunt'/name,dest)
        rows=[]
        for row in selected:
            dest=out/'inputs'/(row['id']+'.json');shutil.copy2(row['source_file'],dest)
            proof=certify(dest);save(out/'inputs'/(row['id']+'.certificate.json'),proof)
            if not proof['all_selected_independent'] or not proof['points']:raise ValueError('Unusable seed certificate')
            rows.append({**row,'seed_lower_bound':proof['LowerBound'],'lower_bound':proof['LowerBound'],
                         'status':'pending','search_seconds':0})
        stage=0;active=[r['id'] for r in rows];elapsed_before=0
    start=time.perf_counter();status='running'
    def checkpoint():
        ordered=sorted(rows,key=lambda r:(-r['lower_bound'],-r['tail_score'],r['id']))
        state={'rows':rows,'stage':stage,'active':active,'controller_wall_seconds':elapsed_before+time.perf_counter()-start}
        save(out/'state.json',state);save(out/'leaderboard.json',ordered)
        save(out/'summary.json',{'status':status,'target':args.target,'best_lower_bound':ordered[0]['lower_bound'],
            'best_candidate':ordered[0]['id'],'candidate_count':len(rows),'stage':stage,
            'controller_wall_seconds':state['controller_wall_seconds'],
            'aggregate_point_search_seconds':sum(r['search_seconds'] for r in rows),
            'record_j_excluded':True,'global_novelty_not_established':True})
    checkpoint()
    try:
        while stage<len(args.seconds) and active:
            print(f'Stage {stage}: {len(active)} curves, cumulative {args.seconds[stage]:g}s each',flush=True)
            tasks=[r for r in rows if r['id'] in active]
            stop_event=threading.Event();executor=ThreadPoolExecutor(max_workers=args.parallel_curves)
            try:
                futures=[executor.submit(search_one,r,out,args.seconds[stage],args,stage,stop_event) for r in tasks]
                for future in as_completed(futures):
                    result=future.result();rows=[result if r['id']==result['id'] else r for r in rows]
                    print(json.dumps({'finished_candidate':result['id'],'lower_bound':result['lower_bound'],
                        'seconds':round(result['search_seconds'],3),'status':result['status']}),flush=True)
                    checkpoint()
            except BaseException:
                stop_event.set();raise
            finally:executor.shutdown(wait=True,cancel_futures=True)
            stage+=1
            if max(r['lower_bound'] for r in rows)>=args.target:status='target_found';break
            active=[r['id'] for r in promoted(rows,args.counts[stage])] if stage<len(args.counts) else []
            checkpoint()
        if status=='running':status='campaign_completed'
    except KeyboardInterrupt:status='interrupted';raise
    except Exception as exc:status='error';save(out/'error.json',{'error':str(exc)});raise
    finally:checkpoint()
    # The best observed lower bound is not advertised as a result until the
    # independent point certificate and the complete search lineage pass.
    best=sorted(rows,key=lambda r:(-r['lower_bound'],-r['tail_score'],r['id']))[0]
    if (out/'searches'/best['id']/'summary.json').exists():
        from audit_blind_search import audit
        audit_dir=out/'audits'/f"{best['id']}-stage-{stage}"
        if not audit_dir.exists():audit(out/'searches'/best['id'],audit_dir,previously_studied=args.previously_studied)
        report=read(audit_dir/'audit.json')
        if report['certified_lower_bound']!=best['lower_bound']:raise ValueError('Final audit bound differs')
        save(out/'best-verified.json',{'candidate':best,'audit':str(audit_dir),
            'certified_lower_bound':best['lower_bound'],'target_reached':best['lower_bound']>=args.target,
            'different_j_from_reference':Q(best['j_invariant'])!=record_j,'global_novelty_not_established':True})
    print(f'DONE: {status}; best certified lower bound {best["lower_bound"]}',flush=True)
    return read(out/'summary.json')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--candidates',required=True);p.add_argument('--output',required=True)
    p.add_argument('--counts',type=int,nargs='+',default=[12,4,2])
    p.add_argument('--seconds',type=float,nargs='+',default=[15,60,240])
    p.add_argument('--point-workers',type=int,default=4);p.add_argument('--parallel-curves',type=int,default=2)
    p.add_argument('--anchors',type=int,default=2048);p.add_argument('--job-seconds',type=float,default=2)
    p.add_argument('--target',type=int,default=32);p.add_argument('--resume',action='store_true')
    p.add_argument('--previously-studied',action='store_true',help='Record prior work on these curves in the audit context')
    p.add_argument('--fast',action='store_true',help='Batch point searches with compact checkpoints and adaptive selection')
    p.add_argument('--confirmation-bound',type=int,default=262139,help='Later-prime check for --fast; 0 disables it')
    p.add_argument('--generate',action='store_true',help='Prepare a multiscale candidate pool before --fast')
    p.add_argument('--samples',type=int,default=10000000)
    p.add_argument('--heights',type=int,nargs='+',default=[3000,10000,30000,100000,1000000])
    p.add_argument('--generation-workers',type=int,default=4)
    p.add_argument('--generation-seed',type=int,default=20260916)
    p.add_argument('--batch-size',type=int,default=8,help='Pointed models per PARI process in --fast')
    args=p.parse_args()
    if not (len(args.counts)==len(args.seconds) and 1<=len(args.counts)<=8
        and all(1<=v<=1000 for v in args.counts) and all(1<=v<=86400 for v in args.seconds)
        and args.counts==sorted(args.counts,reverse=True) and args.seconds==sorted(set(args.seconds))
        and 1<=args.point_workers<=8 and 1<=args.parallel_curves<=4 and args.point_workers*args.parallel_curves<=32
        and 1<=args.anchors<=4096 and .05<=args.job_seconds<=60 and 1<=args.target<=100
        and 1<=args.batch_size<=32 and (args.confirmation_bound==0 or 65521<args.confirmation_bound<=262139)):p.error('Invalid bounded campaign settings')
    if args.generate:
        if not (args.fast and len(args.heights)==len(set(args.heights)) and 1<=len(args.heights)<=10
                and len(args.heights)<=args.samples<=100000000 and all(1<=h<=1000000000 for h in args.heights)
                and 1<=args.generation_workers<=8 and 0<=args.generation_seed<=2147483637):
            p.error('Invalid bounded generation settings')
        prepare_multiscale(args.candidates,args.samples,args.heights,args.generation_workers,args.generation_seed)
    (run_fast if args.fast else run)(args)
