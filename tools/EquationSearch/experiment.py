"""Frozen ICARM equation-only experiments: direct, one fixed anchor, adaptive anchors."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import html
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import time
import urllib.request

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'tools/RankHunt'))
from bootstrap import equation,save,discover,gp,prefix,vec,parse_points,select_basis,Q


def fetch(destination):
    url='https://elliptic-rank.icarm.cloud/curves'
    content=urllib.request.urlopen(url,timeout=30).read().decode()
    chosen={}
    for tr in re.findall(r'<tr\b.*?</tr>',content,re.S):
        attrs=dict(re.findall(r'data-([\w-]+)="([^"]*)"',tr))
        if 'rank' not in attrs: continue
        rank=int(attrs['rank'])
        if not 2<=rank<=15 or rank in chosen: continue
        match=re.search(r'\[\s*-?\d+(?:\s*,\s*-?\d+){4}\s*\]',html.unescape(tr))
        if match:
            a=list(map(str,json.loads(match.group())))
        else:
            # The table truncates large coefficients. Strip the public JSON
            # response immediately; no witness coordinates are saved or printed.
            detail=json.load(urllib.request.urlopen(url.replace('/curves','/curve/'+attrs['id']+'.json'),timeout=30))
            raw=detail.get('ainvs',detail.get('a_invariants'))
            if isinstance(raw,str): raw=json.loads(raw)
            if not isinstance(raw,list):raise ValueError('Missing coefficient field: '+str(list(detail)))
            a=list(map(str,raw))
        equation({'ainvs':a})
        chosen[rank]={'id':int(attrs['id']),'published_lower_bound':rank,'ainvs':a,
                      'url':url.replace('/curves','/curve/'+attrs['id'])}
    if set(chosen)!=set(range(2,16)): raise ValueError('Incomplete frozen sample')
    sample={'source':url,'html_sha256':hashlib.sha256(content.encode()).hexdigest(),
            'selection':'First row in conductor-sorted table for each published lower bound 2..15, before any search.',
            'published_point_coordinates_saved_or_passed_to_search':False,
            'truncated_rows_use_json_with_only_ainvs_retained':True,
            'curves':[chosen[r] for r in sorted(chosen)]}
    save(destination,sample);print(json.dumps(sample,indent=2))


def worker(args):
    from run import data_boundary
    source=Path(args.input).resolve();out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=True)
    reads=data_boundary(source,out);data=json.loads(source.read_text());started=time.perf_counter()
    if args.mode=='seed':
        data=equation(data);basis,state=discover(data,out,args.seconds,args.workers,2,seed_limit=1)
        result={**(basis or {'ainvs':data['ainvs'],'points':[],'rank_lower_bound':0}),
                'status':state['status'],'seconds':time.perf_counter()-started,'bootstrap':state}
    elif args.mode in ('unified','fixed','frozen','adaptive','geometric'):
        if set(data)!={'ainvs','points'} or len(data['points'])!=1: raise ValueError('Exactly one seed is required')
        from seeded import search,independent_result
        found=search(source,out,args.seconds,args.workers,args.anchors,args.target,batch_size=4,anchor_mode=args.mode)
        result=independent_result(found,found['rank_lower_bound'])
        result.update(seconds=time.perf_counter()-started,events=json.loads((out/'checkpoint.json').read_text())['events'])
    else:
        data=equation(data);from blind_search import boxes
        deadline=started+args.seconds;points={};events=[];result={**data,'points':[],'rank_lower_bound':0};jobs=[]
        for n,d in boxes():
            if time.perf_counter()>=deadline: break
            script=prefix(data['ainvs'])+'M=ellminimalmodel(E,&change);\n'
            script+=f'H=ellratpoints(M,[{n},{d}]);\n'
            script+='for(i=1,#H,W=ellchangepointinv(H[i],change);if(!ellisoncurve(E,W),error("Direct inverse"));print("POINT ",W));print("SEARCH_END");quit;\n'
            outtext,stats=gp(script,min(4,deadline-time.perf_counter()))
            jobs.append({'n':n,'d':d,**stats,'finished':'SEARCH_END' in outtext})
            for p in parse_points(outtext,list(map(Q,data['ainvs']))):points[p]=True
            if points:
                basis=select_basis({'ainvs':data['ainvs'],'points':list(points)})
                if basis and basis['rank_lower_bound']>result['rank_lower_bound']:
                    result=basis;events.append({'lower_bound':result['rank_lower_bound'],'seconds':time.perf_counter()-started})
            if result['rank_lower_bound']>=args.target:break
        result.update(seconds=time.perf_counter()-started,events=events,jobs=jobs)
    result.update(mode=args.mode,target=args.target,target_reached=result['rank_lower_bound']>=args.target)
    save(out/'result.json',result);save(out/'data-reads.json',sorted(reads))
    print(json.dumps({'mode':args.mode,'lower_bound':result['rank_lower_bound'],'seconds':result['seconds']}),flush=True)


def campaign(args):
    sample=json.loads(Path(args.sample).read_text());out=Path(args.output).resolve()
    if out==ROOT or not out.is_relative_to(ROOT):raise ValueError('Output must be a workspace subdirectory')
    if out.exists() and any(out.iterdir()):raise ValueError('Use a fresh experiment output')
    out.mkdir(parents=True,exist_ok=True)
    curves=[row for row in sample['curves'] if args.ranks is None or row['published_lower_bound'] in args.ranks]
    config={k:v for k,v in vars(args).items() if k not in ('command','output')}
    config['sample_sha256']=hashlib.sha256(Path(args.sample).read_bytes()).hexdigest()
    config['code_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(HERE.glob('*.py'))}
    save(out/'config.json',config)
    for row in curves:
        folder=out/str(row['id']);folder.mkdir()
        save(folder/'equation.json',{'ainvs':row['ainvs']})
    def call(row,mode,source,seconds):
        folder=out/str(row['id'])/mode
        cmd=[sys.executable,str(HERE/'experiment.py'),'worker','--input',str(source),'--output',str(folder),
             '--mode',mode,'--seconds',str(seconds),'--workers',str(args.workers),
             '--anchors',str(args.anchors),'--target',str(row['published_lower_bound'])]
        t=time.perf_counter()
        with (out/str(row['id'])/(mode+'.log')).open('w',encoding='utf-8') as stream:
            p=subprocess.run(cmd,stdout=stream,stderr=subprocess.STDOUT,timeout=seconds+120)
        if p.returncode:raise RuntimeError(f'Curve {row["id"]}, {mode}: see {mode}.log')
        result=json.loads((folder/'result.json').read_text());result['process_wall_seconds']=time.perf_counter()-t
        return result
    def curve(row):
        folder=out/str(row['id']);result={**row,'arms':{}}
        seed=call(row,'seed',folder/'equation.json',args.bootstrap_seconds)
        result['seed']=seed
        if seed['rank_lower_bound']==1:
            save(folder/'one-point.json',{'ainvs':row['ainvs'],'points':seed['points']})
        # Alternate order across curves; this is an exploratory shared-machine benchmark.
        modes=['direct','fixed','adaptive'] if row['id']%2 else ['adaptive','fixed','direct']
        for mode in modes:
            if mode!='direct' and seed['rank_lower_bound']!=1:continue
            found=call(row,mode,folder/('equation.json' if mode=='direct' else 'one-point.json'),args.seconds)
            result['arms'][mode]=found
            print(json.dumps({'id':row['id'],'published_bound':row['published_lower_bound'],'arm':mode,
                              'lower_bound':found['rank_lower_bound'],'seconds':round(found['seconds'],3),
                              'initial_points':0 if mode=='direct' else 1}),flush=True)
        return result
    started=time.perf_counter();results=[]
    with ThreadPoolExecutor(max_workers=args.parallel_curves) as executor:
        tasks={executor.submit(curve,row):row for row in curves}
        for future in as_completed(tasks):
            try:result=future.result()
            except Exception as error:result={**tasks[future],'error':str(error)}
            results.append(result)
            save(out/'results.json',{'config':config,'curves':sorted(results,key=lambda r:r['published_lower_bound']),
                                   'wall_seconds':time.perf_counter()-started})
    if any('error' in r for r in results):raise RuntimeError('Some experiments failed; inspect results.json')


def published_campaign(args):
    """Benchmark a frozen public sample, passing only one point to each worker."""
    sample_path=Path(args.sample).resolve();output=Path(args.output).resolve()
    if not output.is_relative_to(ROOT) or output.suffix!='.json':
        raise ValueError('Use a workspace JSON report')
    sample=json.loads(sample_path.read_text(encoding='utf-8-sig'))
    rows=sample['curves']
    retry_hash=None
    if args.retry_report:
        raw=Path(args.retry_report).read_bytes();previous=json.loads(raw)
        if previous['config']['sample_sha256']!=hashlib.sha256(sample_path.read_bytes()).hexdigest():
            raise ValueError('Retry requires the identical frozen sample')
        failed={t['id'] for t in previous['trials'] if not t['target_reached']}
        rows=[row for row in rows if row['id'] in failed]
        retry_hash=hashlib.sha256(raw).hexdigest()
    if args.ids is not None: rows=[row for row in rows if row['id'] in args.ids]
    def order(row):
        conductor=row.get('conductor')
        return (int(conductor) if conductor else float('inf'),row['id'])
    leaders=[min((r for r in rows if r['published_lower_bound']==rank),key=order)
             for rank in sorted({r['published_lower_bound'] for r in rows},reverse=True)]
    leader_ids={r['id'] for r in leaders}
    rows=leaders+([] if args.selection=='leaders' else
                  sorted((r for r in rows if r['id'] not in leader_ids),
                         key=lambda r:(-r['published_lower_bound'],order(r))))
    config={'sample_sha256':hashlib.sha256(sample_path.read_bytes()).hexdigest(),
            'code_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in HERE.glob('*.py')},
            'selection':args.selection,'ids':[r['id'] for r in rows],
            'seconds':args.seconds,'workers':args.workers,'anchors':args.anchors,
            'mode':args.mode,'parallel_trials':1,'external_watchdog_seconds':args.seconds+5,
            'initial_points':1,'published_point_index':0,'direct_search':False,'retry_report_sha256':retry_hash,
            'timing':'process_wall_seconds includes startup, model preparation, search and final verification; download excluded'}
    output.parent.mkdir(parents=True,exist_ok=True)
    if output.exists():
        report=json.loads(output.read_text())
        if report['config']!=config:raise ValueError('Resume requires the same sample, code and settings')
    else:
        report={'config':config,'source':sample['source'],'source_sha256':sample['source_sha256'],
                'scope':'Published single-seed benchmark; certifies lower bounds, not full ranks or general completeness.',
                'trials':[],'wall_seconds':0,'complete':False}
        save(output,report)
    done={t['id'] for t in report['trials']};before=report['wall_seconds'];started=time.perf_counter()
    for row in rows:
        if row['id'] in done:continue
        if len(row['points'])!=1:raise ValueError('Frozen sample must retain exactly one public point')
        temporary=tempfile.TemporaryDirectory(prefix='published-trial-',dir=output.parent)
        folder=Path(temporary.name).resolve()
        if folder==output.parent or not folder.is_relative_to(output.parent):
            raise ValueError('Temporary trial escaped workspace output')
        trial={'id':row['id'],'target':row['published_lower_bound'],'url':row['url'],
               'leader':row['id'] in leader_ids,'torsion':row.get('torsion'),
               'rank_lower_bound':None,'target_reached':False}
        try:
            source=folder/'one-point.json';out=folder/'search'
            save(source,{'ainvs':row['ainvs'],'points':row['points']})
            cmd=[sys.executable,str(HERE/'experiment.py'),'worker','--input',str(source),
                 '--output',str(out),'--mode',args.mode,'--seconds',str(args.seconds),
                 '--workers',str(args.workers),'--anchors',str(args.anchors),'--target',str(row['published_lower_bound'])]
            t=time.perf_counter();process=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
            try:
                stdout,stderr=process.communicate(timeout=args.seconds+5)
            except subprocess.TimeoutExpired:
                if sys.platform=='win32':
                    subprocess.run(['taskkill','/PID',str(process.pid),'/T','/F'],capture_output=True,timeout=5)
                else:process.kill()
                stdout,stderr=process.communicate(timeout=5)
                trial['status']='external_timeout'
            trial.update(process_wall_seconds=time.perf_counter()-t,returncode=process.returncode)
            state_path=out/'checkpoint.json'
            if state_path.exists():
                state=json.loads(state_path.read_text())
                trial.update(search_status=state['status'],checkpoint_lower_bound=state['lower_bound'],
                             engine_wall_seconds=state['wall_seconds'],
                             preparation_seconds=sum(p['seconds'] for p in state.get('preparations',[])),
                             attempted_models=state['attempted_models'],
                             events=[{k:e[k] for k in ('lower_bound','seconds','attempted_models')} for e in state['events']])
                trial['distinct_points_observed']=len(state['found']['points'])
                trial['engine_code_sha256']=state['config']['code_sha256']
                if args.mode=='unified':
                    trial['refinement_witnesses']=[s for r in state['refinements'] for s in r['steps']]
                    trial['refinement_seconds']=sum(r['seconds'] for r in state['refinements'])
                    trial['retained_charts']=len(state['models'])
                trial['isogeny_preparations']=state.get('isogeny_preparations',[])
                trial['selected_point_methods']=[o.get('method','pointed') for e in state['events']
                                                 for o in e.get('selected_new_points',[])]
            if process.returncode==0 and (out/'result.json').exists():
                result=json.loads((out/'result.json').read_text())
                trial.update(result=result,rank_lower_bound=result['rank_lower_bound'],
                             target_reached=result['target_reached'],status=state['status'])
                reads=json.loads((out/'data-reads.json').read_text())
                unexpected=[p for p in reads if Path(p).is_relative_to(ROOT) and Path(p)!=source
                            and not Path(p).is_relative_to(out) and Path(p).suffix not in
                            {'.py','.pyc','.dll','.deps.json','.runtimeconfig.json'}]
                if unexpected:raise ValueError('Unexpected benchmark data reads: '+str(unexpected))
                trial['unexpected_workspace_data_reads']=len(unexpected)
                # Keep witnesses and certificates once; concise events above suffice.
                trial['result'].pop('events',None)
            else:
                trial.setdefault('status','seed_not_certified' if 'No independent seed basis' in stderr else 'error')
                trial['error']=(stdout+'\n'+stderr)[-3000:]
        except Exception as error:
            trial.update(status='error',error=str(error))
        finally:
            if folder==output.parent or not folder.is_relative_to(output.parent):
                raise ValueError('Unsafe temporary trial cleanup')
            temporary.cleanup()
        report['trials'].append(trial);report['wall_seconds']=before+time.perf_counter()-started
        save(output,report)
        print(json.dumps({k:trial.get(k) for k in ('id','target','rank_lower_bound','status','process_wall_seconds')}),flush=True)
    report['complete']=True;report['wall_seconds']=before+time.perf_counter()-started
    save(output,report)
    print(json.dumps({'complete':True,'trials':len(report['trials']),
                      'targets_reached':sum(t['target_reached'] for t in report['trials']),
                      'wall_seconds':report['wall_seconds']}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    f=sub.add_parser('fetch');f.add_argument('--output',default=str(HERE/'icarm-equations.json'))
    w=sub.add_parser('worker');w.add_argument('--input',required=True);w.add_argument('--output',required=True)
    w.add_argument('--mode',choices=('seed','direct','unified','fixed','frozen','adaptive','geometric'),required=True)
    w.add_argument('--seconds',type=float,default=30);w.add_argument('--workers',type=int,default=6)
    w.add_argument('--anchors',type=int,default=512);w.add_argument('--target',type=int,required=True)
    c=sub.add_parser('run');c.add_argument('--sample',default=str(HERE/'icarm-equations.json'));c.add_argument('--output',required=True)
    c.add_argument('--seconds',type=float,default=30);c.add_argument('--bootstrap-seconds',type=float,default=10)
    c.add_argument('--workers',type=int,default=6);c.add_argument('--parallel-curves',type=int,default=4)
    c.add_argument('--anchors',type=int,default=512);c.add_argument('--ranks',type=int,nargs='+')
    s=sub.add_parser('published');s.add_argument('--sample',required=True);s.add_argument('--output',required=True)
    s.add_argument('--selection',choices=('all','leaders'),default='all')
    s.add_argument('--mode',choices=('unified','fixed','frozen','adaptive','geometric'),default='unified')
    s.add_argument('--seconds',type=float,default=3);s.add_argument('--workers',type=int,default=2)
    s.add_argument('--anchors',type=int,default=16)
    s.add_argument('--retry-report',help='Retest only failures from a report using this identical sample')
    s.add_argument('--ids',type=int,nargs='+',help='Optional explicit subset, recorded in report configuration')
    args=p.parse_args()
    if args.command=='fetch':fetch(args.output)
    elif args.command=='worker':worker(args)
    elif args.command=='published':
        if not (0<args.seconds<=60 and 1<=args.workers<=24 and 1<=args.anchors<=4096):
            p.error('Invalid bounded benchmark settings')
        published_campaign(args)
    else:campaign(args)
