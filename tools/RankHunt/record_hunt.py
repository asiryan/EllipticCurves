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
    by_tail=sorted(rows,key=lambda r:(-r['tail_score'],r['id']))
    shuffled=sorted(rows,key=lambda r:r['id']);random.Random(POLICY['random_seed']).shuffle(shuffled)
    result=[];seen=set()
    for group in zip(by_score,by_tail,shuffled):
        for row in group:
            if row['id'] not in seen:seen.add(row['id']);result.append(row)
    return result[:count]


def promoted(rows,count):
    eligible=[r for r in rows if r.get('status') in {'budget_completed','target_reached'}]
    return sorted(eligible,key=lambda r:(-r['lower_bound'],-r['tail_score'],-r['score'],r['id']))[:count]


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
    args=p.parse_args()
    if not (len(args.counts)==len(args.seconds) and 1<=len(args.counts)<=8
        and all(1<=v<=1000 for v in args.counts) and all(1<=v<=86400 for v in args.seconds)
        and args.counts==sorted(args.counts,reverse=True) and args.seconds==sorted(set(args.seconds))
        and 1<=args.point_workers<=8 and 1<=args.parallel_curves<=4 and args.point_workers*args.parallel_curves<=16
        and 1<=args.anchors<=4096 and .05<=args.job_seconds<=60 and 1<=args.target<=100):p.error('Invalid bounded campaign settings')
    run(args)
