"""Audit and preserve the completed reference-guided campaign, without searches."""
import argparse
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import shutil

from guided_recovery import matches
from point_arithmetic import combine, negate, on_curve
from point_search import ROOT, make_script, save
from record_recovery import prepare, cli
from search_translated_point import restore


def read(path):
    return json.loads(Path(path).read_text())


def point_set(data):
    a=list(map(Q,data['ainvs']))
    result=set()
    for p in data['points']:
        p=tuple(map(Q,p))
        if not on_curve(a,p):raise ValueError('Off-curve point in saved data')
        result.add(min(p,negate(a,p)))
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',default=str(ROOT/'artifacts/record-point-audit/translated-target31/derived.json'))
    parser.add_argument('--output',default=str(ROOT/'results/record-reproduction-20260914'))
    args=parser.parse_args()
    output=Path(args.output);output.mkdir(parents=True,exist_ok=True)
    seed=output/'seed';seed.mkdir(exist_ok=True)
    seed_path,_=prepare(seed)
    starting=read(seed_path)
    shutil.copy2(ROOT/'tools/RankHunt/Data/icarm302-record-basis.json',output/'basis-change.json')
    reference=read(ROOT/'tests/Fixtures/icarm-302.json')
    operations=[];searches=[]
    # The driver saves snapshots of the previous result as search-input.json.
    # Resolve those identical snapshots to the original GP output and script.
    originals={}
    for root in [ROOT/'artifacts/record-point-audit',ROOT/'artifacts/guided-recovery']:
        for gp in root.rglob('*.gp'):
            candidate=gp.with_suffix('.json')
            if candidate.is_file():
                d=read(candidate)
                if isinstance(d,dict) and d.get('mode')=='pointed':
                    originals[json.dumps(d,sort_keys=True)]=candidate

    def audit(path):
        path=Path(path);data=read(path)
        if data['ainvs']!=starting['ainvs']:raise ValueError('Curve changed')
        if point_set(data)==point_set(starting):
            return data
        if data.get('mode')=='pointed' and not path.with_suffix('.gp').exists():
            path=originals[json.dumps(data,sort_keys=True)]
        parent_path=Path(data['source']);parent=audit(parent_path)
        if 'restoration' in data:
            r=data['restoration'];known=read(r['known_basis'])
            actual=restore(parent,r['observed_representative'],known['points'],r['known_coefficients'])
            if actual!=r['restored_point']:raise ValueError('Restoration changed')
            expected={*point_set(parent),*point_set({'ainvs':data['ainvs'],'points':[actual]})}
            if point_set(data)!=expected:raise ValueError('Unexplained derived point')
            name=f'{len(operations):02d}-restore-{r["published_index"]}'
            save(output/(name+'.json'),{'ainvs':data['ainvs'],'known_points':known['points'],**r,
                'exact_reconstruction_verified':True})
            operations.append({'kind':'restore','file':name+'.json','published_index':r['published_index']})
            return data
        if data.get('mode')!='pointed' or not data['finished'] or data['timed_out']:
            raise ValueError('Unexpected or unfinished search in successful lineage')
        a=list(map(Q,data['ainvs']));anchor=read(data['anchor_source'])
        selected=anchor.get('selection')
        if selected:
            pool=read(selected['anchor_pool']);j=selected['anchor_index']-1
            if anchor['points']!=[pool['points'][j]]:raise ValueError('Anchor pool mismatch')
            basis_path=Path(pool['source'])
            basis=read(basis_path) if basis_path.is_file() else starting
            coefficients=pool['vectors'][j]
        else:
            basis=starting;coefficients=anchor['vector']
        if not point_set(basis)<=point_set(parent):raise ValueError('Anchor basis contains unseen points')
        reconstructed=combine(a,[tuple(map(Q,p)) for p in basis['points']],coefficients)
        if reconstructed!=tuple(map(Q,anchor['points'][0])):raise ValueError('Anchor not a known-point combination')
        script=path.with_suffix('.gp').read_text()
        expected_script=make_script(parent,'pointed',data['height'],data['anchors'],0,
            data['minimal_model_preprocessing'],data['pari_stack_limit_mb'],
            anchor['points'],data['quartic_minimal_preprocessing'],data['denominator_height'])
        if script!=expected_script:raise ValueError('Saved GP script differs from regenerated known-point-only script')
        forbidden=r'\b(?:ellrank|ellrankinit|ell2cover|ellanalyticrank|ellgenerators)\s*\('
        if re.search(forbidden,script):raise ValueError('Forbidden rank computation')
        stdout=path.with_suffix('.stdout.txt').read_text()
        emitted=re.findall(r'^POINT \[(-?\d+(?:/\d+)?), (-?\d+(?:/\d+)?)\]$',stdout,re.M)
        if point_set(data)!=point_set(parent)|point_set({'ainvs':data['ainvs'],'points':emitted}):
            raise ValueError('Saved search contains points absent from input and stdout')
        before,_=matches(parent,reference);after,_=matches(data,reference)
        number=len(searches)+1;name=f'{len(operations):02d}-search-{number:02d}'
        folder=output/name;folder.mkdir(exist_ok=True)
        for suffix in ['.gp','.stdout.txt','.stderr.txt','.json']:
            shutil.copy2(path.with_suffix(suffix),folder/('search'+suffix))
        save(folder/'input.json',parent)
        save(folder/'anchor.json',{'ainvs':data['ainvs'],'points':anchor['points'],
            'known_basis':basis['points'],'integer_coefficients':coefficients,
            'exact_reconstruction_verified':True})
        row={'kind':'search','directory':name,'original_result':str(path.resolve()),
            'height':data['height'],'denominator_height':data['denominator_height'],
            'seconds':data['seconds'],'lower_bound':data['exact_certificate']['LowerBound'],
            'new_direct_matches':sorted(set(after)-set(before)),
            'new_distinct_points':data['new_distinct_up_to_sign'],'gp_output_points_checked':len(emitted),
            'script_regenerated_from_known_points_only':True,'anchor_combination_verified':True}
        searches.append(row);operations.append(row)
        return data

    final=audit(Path(args.input))
    indices,points=matches(final,reference)
    if indices!=list(range(1,32)) or points!=reference['points']:
        raise ValueError('Not all 31 published coordinates match exactly in order')
    result={'ainvs':final['ainvs'],'points':points,'published_indices':indices,
        'reference_guided_reproduction':True,'blind_discovery':False,
        'all_coordinates_match_published_order_and_sign':True}
    save(output/'points31.json',result)
    certificate=json.loads(cli('verify','--input',output/'points31.json'))
    save(output/'csharp-certificate.json',certificate)
    if certificate['LowerBound']!=31 or certificate['hypotheses']:raise ValueError('C# certificate failed')
    verifier_path=ROOT/'artifacts/rank-package-audit/original/elliptic_rank_search/certificate.py'
    digest=hashlib.sha256(verifier_path.read_bytes()).hexdigest()
    if digest!='9e0d0d2562fc53705e92a2eaa9a3f6e7c923f1cd3fd82b14df68b60268f4ad54':
        raise ValueError('Independent verifier changed')
    spec=importlib.util.spec_from_file_location('audited_certificate',verifier_path)
    verifier=importlib.util.module_from_spec(spec);spec.loader.exec_module(verifier)
    proof=verifier.build_certificate(result,max_prime=2000)
    claim=verifier.verify_certificate(result,proof)
    if claim['rank_lower_bound']!=31 or not claim['all_points_independent_modulo_torsion'] or claim['conditional_assumptions']:
        raise ValueError('Independent certificate failed')
    save(output/'python-certificate.json',proof);save(output/'python-verification.json',claim)
    shutil.copy2(verifier_path,output/'independent_certificate_verifier.py')
    save(output/'operations.json',operations)
    unsuccessful=read(ROOT/'artifacts/guided-recovery/20260914-143421/round-10/search.json')
    blind_attempt=read(ROOT/'artifacts/record-recovery/20260914-142724/pointed.json')
    summary={'reference_guided':True,'blind_recovery':False,'initial_generic_points':17,
        'published_points_matched_exactly':31,'certified_lower_bound':31,'exact_rank_not_claimed':True,
        'successful_search_calls':len(searches),
        'successful_search_seconds':round(sum(r['seconds'] for r in searches),3),
        'guided_timeout_seconds':unsuccessful['seconds'],
        'earlier_128_anchor_unsuccessful_search_seconds':blind_attempt['seconds'],
        'times_exclude_model_selection_diagnostics_and_development':True,
        'all_search_scripts_regenerated_and_output_provenance_checked':True,
        'all_anchors_are_exact_integer_combinations_of_previously_observed_points':True,
        'reference_used_for_basis_alignment_model_selection_search_bounds_and_coset_choice':True,
        'unrecovered_reference_coordinates_injected_into_searches':False,
        'new_record_not_claimed':True,'csharp':certificate,'independent_python':claim,
        'independent_verifier_sha256':digest,'operations':'operations.json'}
    save(output/'summary.json',summary)
    hashes={str(p.relative_to(output)):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(output.rglob('*')) if p.is_file() and p.name!='sha256.json'}
    save(output/'sha256.json',hashes)
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':
    main()
