"""Search a reference-guided coset representative, then undo its known translate."""
import argparse
from fractions import Fraction as Q
import json
from pathlib import Path
from types import SimpleNamespace

from guided_recovery import matches
from point_arithmetic import add, combine, negate, on_curve
from point_search import ROOT, run, save
from record_recovery import cli


def restore(found, representative, known, coefficients):
    """Require actual observed points before any exact group operation."""
    a=list(map(Q,found['ainvs']))
    actual={tuple(map(Q,p)) for p in found['points']}
    actual|={negate(a,p) for p in list(actual)}
    r=tuple(map(Q,representative))
    basis=[tuple(map(Q,p)) for p in known]
    if r not in actual or any(p not in actual for p in basis):
        raise ValueError('A translated restoration requires previously observed points')
    result=add(a,r,negate(a,combine(a,basis,coefficients)))
    if result is None or not on_curve(a,result):
        raise ValueError('Invalid restored point')
    return list(map(str,result))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',required=True)
    p.add_argument('--diagnostic',required=True,help='Directory with catalog.json and witnesses.json')
    p.add_argument('--target',type=int,required=True,help='Published witness index')
    p.add_argument('--output',required=True)
    p.add_argument('--timeout',type=float,default=120)
    p.add_argument('--max-area',type=int,default=10**13)
    args=p.parse_args()
    if not 0<args.timeout<=600 or not 0<args.max_area<=10**15:p.error('Invalid limits')
    diag=Path(args.diagnostic);directory=Path(args.output)
    witness=json.loads((diag/'witnesses.json').read_text())
    catalog=json.loads((diag/'catalog.json').read_text())
    selected=min((v for v in catalog.values() if v['translate']['published_index']==args.target),key=lambda v:v['area'])
    if selected['area']>args.max_area:raise ValueError('Search area exceeds limit')
    directory.mkdir(parents=True,exist_ok=False)
    save(directory/'anchor.json',{'ainvs':witness['ainvs'],'points':[selected['anchor']],
                                 'selection':selected,'reference_guided_selection':True})
    result=run(SimpleNamespace(input=args.input,output=str(directory/'search'),mode='pointed',
        height=selected['numerator_bound'],denominator_height=selected['denominator_bound'],
        anchors=1,effort=0,minimal=True,quartic_minimal=True,
        anchor_input=str(directory/'anchor.json'),stack_mb=512,timeout=args.timeout,
        gp=str(ROOT/'artifacts/native-validation/gp.exe')))
    known=json.loads(Path(witness['known_basis']).read_text())
    rep=witness['points'][selected['target_index']-1]
    restored=restore(result,rep,known['points'],selected['translate']['known_coefficients'])
    reference=json.loads(Path(witness['reference']).read_text())
    if list(map(Q,restored))!=list(map(Q,reference['points'][args.target-1])):
        raise RuntimeError('Restored point differs from the published witness')
    derived={'ainvs':result['ainvs'],'points':result['points']+[restored],
             'source':str((directory/'search.json').resolve()),
             'restoration':{'published_index':args.target,'observed_representative':rep,
                'known_basis':str(Path(witness['known_basis']).resolve()),
                'known_coefficients':selected['translate']['known_coefficients'],
                'restored_point':restored},'reference_guided_reproduction':True}
    save(directory/'derived.json',derived)
    cert=json.loads(cli('verify','--input',directory/'derived.json'))
    save(directory/'derived.certificate.json',cert)
    indices,points=matches(derived,reference)
    save(directory/'matched.json',{'ainvs':derived['ainvs'],'points':points,'matched_indices':indices})
    print(json.dumps({'restored_published_index':args.target,'lower_bound':cert['LowerBound'],
                      'matched_count':len(indices),'derived_file':str(directory/'derived.json')},indent=2),flush=True)


if __name__=='__main__':
    main()
