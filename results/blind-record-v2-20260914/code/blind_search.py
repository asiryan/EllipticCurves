"""Bounded seeded point search with no reference points or recovery recipe.

This research heuristic accepts an integral curve and independently known seed
points. It proves lower bounds, not a complete Mordell-Weil group. Candidate
anchors, model priority and rectangular boxes depend only on observed points
and this fixed policy. No curve-specific coordinates or matrices are embedded.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import random
import re
import subprocess
import sys
import time

from point_arithmetic import on_curve, negate
from bounded_anchor_pool import generate
from point_search import ROOT, make_script, quartic_reduction_code, vector, save

GP=ROOT/'artifacts/native-validation/gp.exe'
DLL=ROOT/'tools/RankHunt/bin/Release/net8.0/RankHunt.dll'
POLICY={'version':2,'area_exponents':[20,24,28,32,36,40],
        'denominator_exponents':[0,4,8,12,16,20],
        'model_order':'interleave coefficient size, anchor height, fixed shuffle',
        'random_seed':20260914,'batch_size':32,
        'rebuild_on_certified_growth':True,'reference_access':False,
        'anchor_enumeration':'bounded best-first, at most 16*anchors vectors',
        'certificate_update':'previous independent basis plus newly observed points'}


def clean(data):
    a=list(map(Q,data['ainvs']))
    if len(a)!=5 or any(v.denominator!=1 for v in a):
        raise ValueError('This seeded search requires five integral Weierstrass coefficients')
    points=[]
    for p in data.get('points',[]):
        p=tuple(map(Q,p))
        if len(p)!=2 or not on_curve(a,p):raise ValueError('Invalid seed point')
        points.append(list(map(str,min(p,negate(a,p)))))
    points=sorted(set(map(tuple,points)),key=lambda p:(Q(p[0]),Q(p[1])))
    return {'ainvs':list(map(str,a)),'points':[list(p) for p in points]}


def file_guard(source, output):
    """Prevent this Python process from reading other workspace data files.

    This is an auditable data boundary, not an OS security sandbox. External
    programs receive generated arithmetic scripts / our own JSON only.
    """
    source=source.resolve();output=output.resolve();reads=set()
    def hook(event,args):
        if event!='open' or not isinstance(args[0],(str,bytes)):return
        path=Path(args[0]).resolve()
        if path.is_relative_to(ROOT) and path!=source and not path.is_relative_to(output):
            if path.suffix not in {'.py','.pyc','.dll','.deps.json','.runtimeconfig.json'}:
                raise PermissionError('Unexpected workspace data access: '+str(path))
        mode=args[1]
        if (isinstance(mode,str) and 'r' in mode) or mode is None:reads.add(str(path))
    sys.addaudithook(hook)
    return reads


def certify(path):
    p=subprocess.run(['dotnet',str(DLL),'basis','--input',str(path)],
        capture_output=True,text=True,timeout=60)
    if p.returncode:raise RuntimeError(p.stderr)
    return json.loads(p.stdout)


def gp_call(script, prefix, timeout):
    prefix.with_suffix('.gp').write_text(script)
    start=time.perf_counter();timed_out=False
    try:
        p=subprocess.run([str(GP),'-fq','-s','64M'],input=script,
            text=True,capture_output=True,timeout=timeout)
        stdout,stderr,code=p.stdout,p.stderr,p.returncode
    except subprocess.TimeoutExpired as exc:
        timed_out=True;code=None
        stdout=exc.stdout or '';stderr=exc.stderr or ''
        if isinstance(stdout,bytes):stdout=stdout.decode(errors='replace')
        if isinstance(stderr,bytes):stderr=stderr.decode(errors='replace')
    prefix.with_suffix('.stdout.txt').write_text(stdout)
    prefix.with_suffix('.stderr.txt').write_text(stderr)
    errors=any('***' in line and 'Warning:' not in line for line in stderr.splitlines())
    return stdout,{'seconds':time.perf_counter()-start,'timed_out':timed_out,
        'exit_code':code,'gp_errors':errors,'finished':code==0 and not errors and not timed_out}


def profile_models(pool, folder, timeout):
    a=pool['ainvs'];points=pool['points']
    script='default(parisizemax,536870912);\ndefault(realprecision,100);\nsetrand(20260914);\n'
    script+='E0=ellinit('+vector(a)+');E=ellminimalmodel(E0,&change);\n'
    script+='P=['+','.join(vector(p) for p in points)+'];P=vector(#P,i,ellchangepoint(P[i],change));\n'
    script+=('for(k=1,#P,x0=P[k][1];v0=2*P[k][2]+E.a1*x0+E.a3;'
        'D=x^4-2*(12*x0+E.b2)*x^2+32*v0*x+E.b2^2-8*E.b2*x0-48*x0^2-32*E.b4;'
        'den=denominator(content(D));F=den^2*D;'+quartic_reduction_code(True)+
        'print("MODEL ",[k,vector(5,j,Str(polcoef(C[1],j-1))),vector(3,j,Str(polcoef(C[2],j-1)))]));\n'
        'print("PROFILE_END");quit;\n')
    stdout,stats=gp_call(script,folder/'profile',timeout)
    models=[]
    for line in stdout.splitlines():
        if not line.startswith('MODEL '):continue
        index,p,q=json.loads(line[6:]);anchor=points[index-1]
        key=hashlib.sha256(json.dumps([anchor,p,q]).encode()).hexdigest()
        models.append({'key':key,'anchor':anchor,'pool_index':index,
            'coefficient_bits':max(abs(int(v)).bit_length() for v in p+q),
            'anchor_height':pool['approximate_heights'][index-1]})
    stats['models']=len(models)
    save(folder/'profile.json',stats)
    if stats['gp_errors']:raise RuntimeError('Model profiling failed; see stderr')
    return models


def order_models(models):
    rng=random.Random(POLICY['random_seed'])
    shuffled=list(models);rng.shuffle(shuffled)
    lists=[sorted(models,key=lambda m:(m['coefficient_bits'],m['anchor_height'],m['key'])),
           sorted(models,key=lambda m:(m['anchor_height'],m['coefficient_bits'],m['key'])),shuffled]
    result=[];seen=set()
    for group in zip(*lists):
        for m in group:
            if m['key'] not in seen:seen.add(m['key']);result.append(m)
    return result


def boxes():
    for area in POLICY['area_exponents']:
        for den in POLICY['denominator_exponents']:
            if 2*den<=area:yield 2**(area-den),2**den


def search_job(data, model, n, d, folder, timeout):
    script=make_script(data,'pointed',n,1,0,True,512,[model['anchor']],True,d)
    stdout,stats=gp_call(script,folder,timeout)
    points=re.findall(r'^POINT \[(-?\d+(?:/\d+)?), (-?\d+(?:/\d+)?)\]$',stdout,re.M)
    found=clean({'ainvs':data['ainvs'],'points':points})['points']
    stats.update({'model':model['key'],'numerator_bound':n,'denominator_bound':d,
        'points':found,'anchor':model['anchor'],'source_pool_index':model['pool_index'],
        'finished':stats['finished'] and 'SEARCH_END' in stdout})
    save(folder.with_suffix('.json'),stats)
    return stats


def run(args):
    source=Path(args.input).resolve();out=Path(args.output).resolve()
    if out==ROOT or not out.is_relative_to(ROOT):raise ValueError('Output must be a dedicated directory in this workspace')
    out.mkdir(parents=True,exist_ok=args.resume)
    module_names=['blind_search.py','point_search.py','bounded_anchor_pool.py','point_arithmetic.py']
    code_hashes={name:hashlib.sha256((Path(__file__).parent/name).read_bytes()).hexdigest() for name in module_names}
    input_hash=hashlib.sha256(source.read_bytes()).hexdigest()
    reads=file_guard(source,out)
    initial=clean(json.loads(source.read_text(encoding='utf-8-sig')))
    if not initial['points']:raise ValueError('Supply known seed points; equation-only bootstrap is outside this module')
    config={'policy':POLICY,'code_sha256':code_hashes,'input_sha256':input_hash,
        'anchors':args.anchors,'job_seconds':args.job_seconds,'workers':args.workers,'target':args.target}
    if args.resume:
        if json.loads((out/'config.json').read_text())!=config:raise ValueError('Resume requires unchanged input, code and policy')
        state=json.loads((out/'state.json').read_text())
        found=state['found'];done=set(state['done']);generation=state['generation'];total_jobs=state['total_jobs']
        elapsed_before=state['wall_seconds']
    else:
        save(out/'config.json',config);save(out/'input.json',initial)
        (out/'code').mkdir()
        for name in module_names:(out/'code'/name).write_bytes((Path(__file__).parent/name).read_bytes())
        found=initial;done=set();generation=0;total_jobs=0;elapsed_before=0
    save(out/'points.json',found);cert=certify(out/'points.json')
    if not cert['all_selected_independent'] or not cert['points']:
        raise ValueError('Current local certificate could not extract a nonempty independent seed basis')
    start=time.perf_counter();deadline=start+args.seconds
    summary={'reference_points_loaded':False,'reference_recipe_loaded':False,
        'seeded_search':True,'equation_only':False,'policy':POLICY,
        'initial_lower_bound':cert['LowerBound'],'lower_bound':cert['LowerBound'],
        'status':'running','events':[],'started_at':datetime.now().isoformat()}
    if args.resume:
        summary=json.loads((out/'summary.json').read_text());summary['status']='running'
    (out/'jobs').mkdir(exist_ok=True)
    def checkpoint():
        summary.update({'lower_bound':cert['LowerBound'],'point_count':len(found['points']),
            'generation':generation,'total_jobs':total_jobs,
            'wall_seconds':elapsed_before+time.perf_counter()-start})
        save(out/'points.json',found);save(out/'certificate.json',cert)
        save(out/'summary.json',summary)
        save(out/'state.json',{'found':found,'done':sorted(done),'generation':generation,
            'total_jobs':total_jobs,'wall_seconds':summary['wall_seconds']})
        save(out/'data-reads.json',sorted(reads))
    checkpoint()
    print(f'Output: {out}\nInitial lower bound {cert["LowerBound"]}; {args.seconds:g}s budget',flush=True)
    try:
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            while time.perf_counter()<deadline and cert['LowerBound']<args.target:
                folder=out/f'generation-{generation:03d}'
                if (folder/'models.json').exists():models=json.loads((folder/'models.json').read_text())
                else:
                    folder.mkdir(exist_ok=True)
                    basis={'ainvs':found['ainvs'],'points':cert['points']}
                    save(folder/'basis.json',basis)
                    if (folder/'pool/anchors.json').exists():pool=json.loads((folder/'pool/anchors.json').read_text())
                    else:pool=generate(folder/'basis.json',folder/'pool',args.anchors,args.anchors,min(60,max(1,deadline-time.perf_counter())))
                    models=order_models(profile_models(pool,folder,min(60,max(1,deadline-time.perf_counter()))))
                    save(folder/'models.json',models)
                if not models:raise RuntimeError('No usable models')
                print(f'Generation {generation}: {len(models)} models, bound {cert["LowerBound"]}',flush=True)
                improved=False
                for n,d in boxes():
                    candidates=[m for m in models if f'{m["key"]}:{n}:{d}' not in done]
                    for offset in range(0,len(candidates),POLICY['batch_size']):
                        remaining=deadline-time.perf_counter()
                        if remaining<=0:break
                        batch=candidates[offset:offset+POLICY['batch_size']]
                        # Bound pending work so a batch cannot overrun the total
                        # budget by many per-job limits.
                        limit=min(args.job_seconds,max(.05,remaining/(len(batch)/args.workers+1)))
                        futures=[]
                        for j,m in enumerate(batch):
                            prefix=out/'jobs'/f'{total_jobs+j:06d}'
                            futures.append((m,executor.submit(search_job,found,m,n,d,prefix,limit)))
                        old_points={tuple(p) for p in found['points']}
                        old_count=len(found['points']);old_bound=cert['LowerBound']
                        for m,future in futures:
                            result=future.result()
                            if result['gp_errors']:raise RuntimeError('GP error in point search; see job log')
                            found=clean({'ainvs':found['ainvs'],'points':found['points']+result['points']})
                            # Timeouts are recorded attempts, not exhausted boxes.
                            # A later invocation can use a fresh campaign/budget.
                            done.add(f'{m["key"]}:{n}:{d}')
                        total_jobs+=len(batch)
                        if len(found['points'])>old_count:
                            # Preserve the previous independent subgroup and
                            # avoid rechecking a growing archive of dependent
                            # observations at every batch. Nothing is discarded
                            # from the archive, and the final audit checks it all.
                            new_candidates=[p for p in found['points'] if tuple(p) not in old_points]
                            proof_input={'ainvs':found['ainvs'],'points':cert['points']+new_candidates}
                            save(out/'certificate-input.json',proof_input)
                            cert=certify(out/'certificate-input.json')
                            if cert['LowerBound']<old_bound:raise RuntimeError('Independent basis was lost')
                            event={'jobs':total_jobs,'n':n,'d':d,'new_distinct_points':len(found['points'])-old_count,
                                'old_lower_bound':old_bound,'lower_bound':cert['LowerBound'],
                                'seconds':elapsed_before+time.perf_counter()-start}
                            summary['events'].append(event)
                            print(json.dumps(event),flush=True)
                        checkpoint()
                        if cert['LowerBound']>old_bound:
                            if not cert['all_selected_independent']:raise RuntimeError('Unable to certify next search basis')
                            improved=True;generation+=1;checkpoint();break
                    if improved or time.perf_counter()>=deadline:break
                if not improved:
                    summary['status']='budget_completed' if time.perf_counter()>=deadline else 'policy_exhausted'
                    break
            if cert['LowerBound']>=args.target:summary['status']='target_reached'
            elif summary['status']=='running':summary['status']='budget_completed'
    except KeyboardInterrupt:
        summary['status']='interrupted'
    except Exception as exc:
        summary['status']='error';summary['error']=str(exc);raise
    finally:checkpoint()
    print(f'DONE: {summary["status"]}; bound {cert["LowerBound"]}; {summary["wall_seconds"]:.3f}s',flush=True)
    return summary


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',required=True)
    p.add_argument('--output',default=str(ROOT/'artifacts/blind-search'/datetime.now().strftime('%Y%m%d-%H%M%S')))
    p.add_argument('--seconds',type=float,default=600)
    p.add_argument('--job-seconds',type=float,default=2)
    p.add_argument('--anchors',type=int,default=2048)
    p.add_argument('--workers',type=int,default=2)
    p.add_argument('--target',type=int,default=31)
    p.add_argument('--resume',action='store_true')
    a=p.parse_args()
    if not (1<=a.seconds<=86400 and .05<=a.job_seconds<=60 and 1<=a.anchors<=4096 and 1<=a.workers<=8 and 1<=a.target<=100):p.error('Invalid bounded search settings')
    run(a)


if __name__=='__main__':main()
