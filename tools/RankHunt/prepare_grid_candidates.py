"""Reuse parameter scores from a grid; regenerate all seed points from formulas."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time

from point_search import ROOT,save
from record_hunt import j_invariant,RECORD_A


def prepare(source,output):
    start=time.perf_counter();source=Path(source).resolve();output=Path(output).resolve()
    if output==ROOT or not output.is_relative_to(ROOT):raise ValueError('Use a dedicated workspace output directory')
    rows=json.loads((source/'candidates.json').read_text());output.mkdir(parents=True,exist_ok=False)
    seen={j_invariant({'ainvs':RECORD_A})};result=[]
    for row in rows:
        if row.get('Control',False):continue
        u,v=int(row['U']),int(row['V'])
        name=f'curve_{u}_{v}.json';target=output/name
        p=subprocess.run(['dotnet',str(ROOT/'tools/RankHunt/bin/Release/net8.0/RankHunt.dll'),'export',
            '--family','icarm302-17','--u',str(u),'--v',str(v),'--output',str(target)],
            capture_output=True,text=True,timeout=60)
        if p.returncode:raise RuntimeError(p.stderr)
        data=json.loads(target.read_text());j=j_invariant(data)
        if j in seen:continue
        seen.add(j)
        result.append({'id':f'{u}_{v}','u':u,'v':v,'file':name,'score':row['Score'],
            'screening_score':row['ScreeningScore'],'tail_score':row['ValidationScore'],
            'seed_lower_bound':data['rank_lower_bound'],'j_invariant':str(j),
            'seed_origin':'fresh export of the 17 generic family sections',
            'old_point_files_loaded':False,'published_novelty_verified':False})
    save(output/'candidates.json',result)
    save(output/'complete.json',{'kind':'previous grid parameter scores with freshly generated seeds',
        'parameter_source':str(source),'source_sha256':hashlib.sha256((source/'candidates.json').read_bytes()).hexdigest(),
        'exported_curves':len(result),'old_point_files_loaded':False,'record_j_excluded':True,
        'seconds':time.perf_counter()-start})
    print(f'Prepared {len(result)} curves from parameter scores in {time.perf_counter()-start:.3f}s.',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input',required=True);p.add_argument('--output',required=True)
    args=p.parse_args();prepare(args.input,args.output)
