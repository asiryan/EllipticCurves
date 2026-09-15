"""Replay the saved, reference-guided search recipe without loading missing witnesses.

This is a previously tuned recipe for ICARM302, not an equation-only solver.
Missing points must be returned by an actual bounded search before they can
participate in any anchor or coordinate restoration.
"""
import argparse
from datetime import datetime
from fractions import Fraction as Q
import json
from pathlib import Path
from types import SimpleNamespace
import time

from point_arithmetic import combine
from point_search import ROOT, run, save
from record_recovery import cli
from save_record_reproduction import point_set, read
from search_translated_point import restore


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bundle',default=str(ROOT/'results/record-reproduction-20260914'))
    p.add_argument('--output',default=str(ROOT/'artifacts/record-replay'/datetime.now().strftime('%Y%m%d-%H%M%S')))
    p.add_argument('--timeout',type=float,default=180)
    args=p.parse_args()
    if not 0<args.timeout<=600:p.error('Invalid per-search timeout')
    bundle=Path(args.bundle);out=Path(args.output);out.mkdir(parents=True,exist_ok=False)
    current=read(bundle/'seed/input17.json')
    current_path=out/'initial17.json';save(current_path,current)
    summary={'reference_guided_precomputed_recipe':True,'blind_recovery':False,'operations':[]}
    start=time.perf_counter()
    for i,op in enumerate(read(bundle/'operations.json')):
        if op['kind']=='search':
            anchor=read(bundle/op['directory']/'anchor.json')
            basis={'ainvs':current['ainvs'],'points':anchor['known_basis']}
            if not point_set(basis)<=point_set(current):raise ValueError('Recipe requests an unseen anchor point')
            a=list(map(Q,current['ainvs']))
            point=combine(a,[tuple(map(Q,q)) for q in basis['points']],anchor['integer_coefficients'])
            if point!=tuple(map(Q,anchor['points'][0])):raise ValueError('Invalid anchor recipe')
            anchor_path=out/f'{i:02d}-anchor.json'
            save(anchor_path,{'ainvs':current['ainvs'],'points':[list(map(str,point))]})
            print(f"Search {i}: N={op['height']}, D={op['denominator_height']}",flush=True)
            prefix=out/f'{i:02d}-search'
            current=run(SimpleNamespace(input=str(current_path),output=str(prefix),
                mode='pointed',height=op['height'],denominator_height=op['denominator_height'],
                anchors=1,effort=0,minimal=True,quartic_minimal=True,
                anchor_input=str(anchor_path),stack_mb=512,timeout=args.timeout,
                gp=str(ROOT/'artifacts/native-validation/gp.exe')))
            current_path=prefix.with_suffix('.json')
            summary['operations'].append({'operation':i,'seconds':current['seconds'],
                'lower_bound':current['exact_certificate']['LowerBound'],'finished':current['finished']})
            if not current['finished']:
                summary['status']='search_incomplete';break
        elif op['kind']=='restore':
            recipe=read(bundle/op['file'])
            point=restore(current,recipe['observed_representative'],recipe['known_points'],recipe['known_coefficients'])
            current={'ainvs':current['ainvs'],'points':current['points']+[point],
                'source':str(current_path.resolve()),'restored_published_index':recipe['published_index']}
            current_path=out/f'{i:02d}-restored.json';save(current_path,current)
            summary['operations'].append({'operation':i,'restored_index':recipe['published_index']})
        else:raise ValueError('Unknown recipe operation')
        save(out/'summary.json',summary)
    else:
        summary['status']='recipe_completed'
    summary['certificate']=json.loads(cli('verify','--input',current_path))
    summary['result']=str(current_path.resolve());summary['seconds']=time.perf_counter()-start
    save(out/'summary.json',summary)
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':
    main()
