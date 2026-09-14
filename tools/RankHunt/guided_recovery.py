"""Reference-guided reproduction, not blind discovery of the published points.

Published witnesses guide model selection and rectangular search bounds.
Each actual point search receives only previously obtained points and an anchor
constructed from them. No witness coordinate is injected as a search result.
"""
import argparse
import contextlib
from datetime import datetime
from fractions import Fraction as Q
import json
from pathlib import Path
from types import SimpleNamespace
import time

from point_search import ROOT, run, save
from point_anchor_pool import generate
from point_search_audit import audit
from record_recovery import cli
from point_arithmetic import negate


def rounded_bound(n):
    p=10**(len(str(n))-1)
    return next(m*p for m in [1,2,5,10] if m*p>=n)


def matches(data, reference):
    if data['ainvs']!=reference['ainvs']:
        raise ValueError('Reference curve mismatch')
    a=list(map(Q,data['ainvs']))
    actual=[tuple(map(Q,p)) for p in data['points']]
    available={p:p for p in actual}
    available.update({negate(a,p):negate(a,p) for p in actual})
    matched=[];points=[]
    for i,p in enumerate(reference['points']):
        q=tuple(map(Q,p))
        if q in available:
            matched.append(i+1)
            points.append(list(map(str,available[q])))
    return matched,points


def add_catalog(catalog, report, pool, source):
    # Older reports retain successive height improvements; newer reports also
    # retain area improvements. Either gives valid bounds, never a guarantee
    # that this is the globally best model.
    for row in report['preimages']:
        anchor,index,_,h,_,_=row
        h=Q(h);dn=rounded_bound(h.denominator)
        nm=rounded_bound(max(abs(h.numerator),h.denominator))
        area=nm*dn
        old=catalog.get(str(index))
        if old is None or area<int(old['area']):
            catalog[str(index)]={'target_index':index,'anchor':pool['points'][anchor-1],
                'anchor_index':anchor,'anchor_pool':str(source),'numerator_bound':nm,
                'denominator_bound':dn,'area':area}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',required=True)
    p.add_argument('--reference',default=str(ROOT/'tests/Fixtures/icarm-302.json'))
    p.add_argument('--output',default=str(ROOT/'artifacts/guided-recovery'/datetime.now().strftime('%Y%m%d-%H%M%S')))
    p.add_argument('--rounds',type=int,default=16)
    p.add_argument('--anchors',type=int,default=2048)
    p.add_argument('--timeout',type=float,default=120)
    p.add_argument('--max-area',type=int,default=10**13)
    p.add_argument('--seed-catalog',action='append',default=[],help='Existing audit directory with audit.json')
    args=p.parse_args()
    if not (1<=args.rounds<=32 and 1<=args.anchors<=4096 and 0<args.timeout<=600 and 1<=args.max_area<=10**15):
        p.error('Invalid bounded recovery limits')
    directory=Path(args.output).resolve();directory.mkdir(parents=True,exist_ok=False)
    reference=json.loads(Path(args.reference).read_text())
    found=json.loads(Path(args.input).read_text())
    matched,points=matches(found,reference)
    catalog={}
    for source in args.seed_catalog:
        report=json.loads((Path(source)/'audit.json').read_text())
        pool=json.loads(Path(report['seed_input']).read_text())
        add_catalog(catalog,report,pool,report['seed_input'])
    summary={'reference_guided_model_selection_and_bounds':True,'blind_recovery':False,
        'search_input_includes_unrecovered_witnesses':False,'initial_input':str(Path(args.input).resolve()),
        'initial_matched_indices':matched,'matched_indices':matched,'runs':[],
        'status':'running','started_at':datetime.now().isoformat()}
    start=time.perf_counter()
    print(f'Output: {directory}\nStarting with {len(matched)} exactly matched points',flush=True)
    for number in range(args.rounds):
        stage=directory/f'round-{number:02d}';stage.mkdir()
        basis=stage/'basis.json'
        save(basis,{'ainvs':found['ainvs'],'points':points,'matched_indices':matched,
                    'source':'Previously obtained points, with reference sign and order'})
        cert=json.loads(cli('verify','--input',basis))
        if cert['LowerBound']!=len(points):
            raise RuntimeError('Matched basis independence not certified')
        save(stage/'basis.certificate.json',cert)
        # Fresh pools use only the current independent matched basis.
        pool=generate(basis,stage/'pool',args.anchors,args.anchors,60)
        with (stage/'audit-console.txt').open('w') as stream, contextlib.redirect_stdout(stream):
            report=audit(stage/'pool/anchors.json',args.reference,stage/'audit',120,100000,True)
        add_catalog(catalog,report,pool,stage/'pool/anchors.json')
        save(directory/'catalog.json',catalog)
        candidates=[r for k,r in catalog.items() if int(k) not in matched]
        if not candidates:
            summary['status']='no_model';break
        selected=min(candidates,key=lambda r:r['area'])
        if selected['area']>args.max_area:
            summary['status']='search_area_limit';summary['next_model']=selected;break
        anchor=stage/'anchor.json'
        save(anchor,{'ainvs':found['ainvs'],'points':[selected['anchor']],
                    'selection':selected,'reference_guided_selection':True})
        input_path=stage/'search-input.json';save(input_path,found)
        print(f"ROUND {number}: target #{selected['target_index']}, N={selected['numerator_bound']}, D={selected['denominator_bound']}",flush=True)
        summary['active_round']=number;save(directory/'summary.json',summary)
        with (stage/'search-console.txt').open('w') as stream, contextlib.redirect_stdout(stream):
            result=run(SimpleNamespace(input=str(input_path),output=str(stage/'search'),
                mode='pointed',height=selected['numerator_bound'],denominator_height=selected['denominator_bound'],
                anchors=1,effort=0,minimal=True,quartic_minimal=True,anchor_input=str(anchor),stack_mb=512,
                timeout=args.timeout,gp=str(ROOT/'artifacts/native-validation/gp.exe')))
        previous=set(matched)
        found=result;matched,points=matches(found,reference)
        row={'round':number,'selected_model':selected,'seconds':result['seconds'],
            'finished':result['finished'],'timed_out':result['timed_out'],
            'new_distinct_points':result['new_distinct_up_to_sign'],
            'lower_bound':result['exact_certificate']['LowerBound'],
            'new_matched_indices':sorted(set(matched)-previous),'matched_count':len(matched)}
        summary['runs'].append(row);summary['matched_indices']=matched
        summary['best_lower_bound']=result['exact_certificate']['LowerBound']
        summary['best_points_file']=str(stage/'search.json')
        summary['seconds']=time.perf_counter()-start
        save(directory/'matched-points.json',{'ainvs':found['ainvs'],'points':points,'matched_indices':matched,
             'reference_guided_reproduction':True,'all_points_already_present_in_search_input_or_output_up_to_sign':True})
        save(directory/'summary.json',summary)
        print(f"RESULT: bound {row['lower_bound']}, {len(matched)}/31 exact matches, new indices {row['new_matched_indices']}, {row['seconds']}s",flush=True)
        if len(matched)==len(reference['points']):
            summary['status']='all_reference_points_reproduced';break
        if selected['target_index'] not in matched:
            summary['status']='search_timeout' if result['timed_out'] else 'predicted_point_not_recovered'
            break
    else:
        summary['status']='round_limit'
    summary.pop('active_round',None);summary['seconds']=time.perf_counter()-start
    save(directory/'summary.json',summary)
    print(f"DONE: {summary['status']}; {len(matched)}/31 matches; {summary['seconds']:.3f}s",flush=True)


if __name__=='__main__':
    main()
