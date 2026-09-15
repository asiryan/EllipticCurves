"""Check a completed seeded search's exact lineage; never load a point reference."""
import argparse
from bisect import bisect_right
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import shutil
import time

from blind_search import clean
from point_arithmetic import combine
from point_search import ROOT, make_script, save


def read(path):return json.loads(Path(path).read_text())


def key_points(points):return {tuple(p) for p in points}


def integer_inverse(matrix):
    """Invert a unimodular matrix with exact arithmetic, rejecting other inputs."""
    n=len(matrix)
    if not n or any(len(row)!=n or any(Q(v).denominator!=1 for v in row) for row in matrix):
        raise ValueError('Invalid integral basis change')
    rows=[list(map(Q,row))+[Q(i==j) for j in range(n)] for i,row in enumerate(matrix)]
    for i in range(n):
        pivot=next((j for j in range(i,n) if rows[j][i]),None)
        if pivot is None:raise ValueError('Singular basis change')
        rows[i],rows[pivot]=rows[pivot],rows[i]
        scale=rows[i][i];rows[i]=[v/scale for v in rows[i]]
        for j in range(n):
            if j!=i and rows[j][i]:
                scale=rows[j][i];rows[j]=[x-scale*y for x,y in zip(rows[j],rows[i])]
    inverse=[row[n:] for row in rows]
    if any(v.denominator!=1 for row in inverse for v in row):
        raise ValueError('Basis change is not unimodular')
    return [[int(v) for v in row] for row in inverse]


def audit(run_dir, output, previously_studied=True):
    start=time.perf_counter();summary=read(run_dir/'summary.json')
    if summary['status']=='running':raise ValueError('Audit a completed or stopped campaign')
    config=read(run_dir/'config.json');initial=clean(read(run_dir/'input.json'))
    output.mkdir(parents=True,exist_ok=False)
    for name,digest in config['code_sha256'].items():
        if hashlib.sha256((run_dir/'code'/name).read_bytes()).hexdigest()!=digest:
            raise ValueError('Frozen source changed: '+name)
    external_data=[]
    for name in read(run_dir/'data-reads.json'):
        path=Path(name).resolve()
        if path.is_relative_to(ROOT) and not path.is_relative_to(run_dir.resolve()) and path.suffix not in {'.py','.pyc','.dll'}:
            external_data.append(path)
    if len(external_data)!=1 or hashlib.sha256(external_data[0].read_bytes()).hexdigest()!=config['input_sha256']:
        raise ValueError('Read log does not contain exactly the declared input data')
    if clean(read(external_data[0]))!=initial:raise ValueError('Saved seed differs from declared input')
    # Retain the exact bytes named by input_sha256, including ignored metadata.
    # input.json is the normalized arithmetic input and can have another hash.
    shutil.copy2(external_data[0],output/'declared-input.json')
    # The saved source can differ after future development. It remains the
    # authoritative version executed by this campaign.
    if hashlib.sha256((run_dir/'code/point_search.py').read_bytes()).hexdigest()!=hashlib.sha256((ROOT/'tools/RankHunt/point_search.py').read_bytes()).hexdigest():
        raise ValueError('Use the saved point-search version to audit its GP scripts')
    for name in ['input.json','config.json','summary.json','certificate.json','data-reads.json']:
        shutil.copy2(run_dir/name,output/name)
    known=key_points(initial['points']);a=list(map(Q,initial['ainvs']))
    starts=[0]+[e['jobs'] for e in summary['events'] if e['lower_bound']>e['old_lower_bound']]
    generations={};short_lattices={};verified_anchors=set();jobs_checked=0;printed_points=0;timeouts=0
    point_origin={p:'seed' for p in known}
    for path in sorted((run_dir/'jobs').glob('*.json')):
        number=int(path.stem)
        if number>=summary['total_jobs']:continue
        job=read(path);g=bisect_right(starts,number)-1
        if g not in generations:
            folder=run_dir/f'generation-{g:03d}'
            basis=clean(read(folder/'basis.json'));pool=read(folder/'pool/anchors.json')
            if not key_points(basis['points'])<=known:raise ValueError('Generation uses unseen basis points')
            # Preserve the exact original basis ordering used by vector entries.
            basis=read(folder/'basis.json')
            generations[g]=(basis,pool)
            lattice=folder/'pool/lattice.stdout.txt'
            if lattice.exists():
                lines=lattice.read_text().splitlines()
                transform,_=json.loads(next(l[8:] for l in lines if l.startswith('LATTICE ')))
                short=[tuple(map(Q,json.loads(l[6:]))) for l in lines if l.startswith('BASIS ')]
                original=[tuple(map(Q,p)) for p in basis['points']];n=len(original)
                if len(transform)!=n or len(short)!=n:raise ValueError('Wrong LLL basis size')
                inverse=integer_inverse(transform)
                for j,p in enumerate(short):
                    if combine(a,original,[transform[i][j] for i in range(n)])!=p:
                        raise ValueError('LLL point change is not exact')
                short_lattices[g]=(short,inverse)
        basis,pool=generations[g];index=job['source_pool_index']-1
        if pool['points'][index]!=job['anchor']:raise ValueError('Wrong anchor in job')
        anchor_key=(g,index)
        if anchor_key not in verified_anchors:
            coefficients=pool['vectors'][index]
            if len(coefficients)!=len(basis['points']) or any(Q(v).denominator!=1 for v in coefficients):
                raise ValueError('Invalid anchor coefficient vector')
            if g in short_lattices:
                short,inverse=short_lattices[g]
                coordinates=[sum(x*y for x,y in zip(row,coefficients)) for row in inverse]
                p=combine(a,short,coordinates)
            else:p=combine(a,[tuple(map(Q,p)) for p in basis['points']],coefficients)
            if p!=tuple(map(Q,job['anchor'])):raise ValueError('Anchor is not its claimed known-point combination')
            verified_anchors.add(anchor_key)
        script=path.with_suffix('.gp').read_text()
        expected=make_script(initial,'pointed',job['numerator_bound'],1,0,True,512,[job['anchor']],True,job['denominator_bound'])
        if script!=expected:raise ValueError('GP script differs from known-point-only regeneration')
        if re.search(r'\b(?:ellrank|ellrankinit|ell2cover|ellanalyticrank|ellgenerators|read|system)\s*\(',script):
            raise ValueError('Unexpected operation in a search script')
        stdout=path.with_suffix('.stdout.txt').read_text()
        points=re.findall(r'^POINT \[(-?\d+(?:/\d+)?), (-?\d+(?:/\d+)?)\]$',stdout,re.M)
        actual=clean({'ainvs':initial['ainvs'],'points':points})['points']
        if actual!=job['points']:raise ValueError('Job result differs from actual POINT output')
        for p in key_points(actual)-known:point_origin[p]=number
        known.update(map(tuple,actual));jobs_checked+=1;printed_points+=len(points)
        timeouts+=int(job['timed_out'])
    final=clean(read(run_dir/'points.json'))
    if known!=key_points(final['points']):raise ValueError('Final archive differs from seed plus observed outputs')
    selected=read(run_dir/'certificate.json')
    witness=clean({'ainvs':initial['ainvs'],'points':selected['points']})
    if not key_points(witness['points'])<=known:raise ValueError('Certificate contains unseen points')
    verifier_path=ROOT/'artifacts/rank-package-audit/original/elliptic_rank_search/certificate.py'
    if hashlib.sha256(verifier_path.read_bytes()).hexdigest()!='9e0d0d2562fc53705e92a2eaa9a3f6e7c923f1cd3fd82b14df68b60268f4ad54':raise ValueError('Verifier changed')
    spec=importlib.util.spec_from_file_location('independent_certificate',verifier_path)
    verifier=importlib.util.module_from_spec(spec);spec.loader.exec_module(verifier)
    proof=verifier.build_certificate(witness,max_prime=2000)
    claim=verifier.verify_certificate(witness,proof)
    if claim['rank_lower_bound']!=summary['lower_bound'] or not claim['all_points_independent_modulo_torsion']:
        raise ValueError('Independent proof did not match')
    save(output/'independent-points.json',witness)
    save(output/'python-certificate.json',proof);save(output/'python-verification.json',claim)
    shutil.copy2(verifier_path,output/'independent_certificate_verifier.py')
    save(output/'point-origins.json',[{'point':p,'first_observed_job':point_origin[tuple(p)]} for p in witness['points']])
    save(output/'all-found-points.json',final)
    shutil.copytree(run_dir/'code',output/'code')
    # Keep the full attempt corpus in artifacts; preserve direct witness jobs
    # and their exact anchor expressions in the compact results bundle.
    origins=sorted({point_origin[tuple(p)] for p in witness['points'] if point_origin[tuple(p)]!='seed'})
    (output/'witness-jobs').mkdir()
    for number in origins:
        prefix=run_dir/'jobs'/f'{number:06d}';job=read(prefix.with_suffix('.json'))
        g=bisect_right(starts,number)-1;basis,pool=generations[g];j=job['source_pool_index']-1
        for suffix in ['.gp','.stdout.txt','.stderr.txt','.json']:
            shutil.copy2(prefix.with_suffix(suffix),output/'witness-jobs'/(prefix.name+suffix))
        save(output/'witness-jobs'/(prefix.name+'.anchor-proof.json'),
            {'basis':basis,'coefficients':pool['vectors'][j],'anchor':job['anchor'],
             'basis_origins':[point_origin[tuple(clean({'ainvs':basis['ainvs'],'points':[p]})['points'][0])] for p in basis['points']]})
    report={'campaign':str(run_dir.resolve()),'status':summary['status'],
        'seed_count':len(initial['points']),'certified_lower_bound':summary['lower_bound'],
        'wall_seconds':summary['wall_seconds'],'jobs_checked':jobs_checked,
        'used_anchor_combinations_exactly_checked':len(verified_anchors),
        'lll_basis_changes_exactly_checked':len(short_lattices),
        'raw_point_outputs_checked':printed_points,'timed_out_attempts':timeouts,
        'archive_is_exact_seed_union_search_outputs':True,
        'all_used_anchors_come_from_previously_observed_points':True,
        'all_search_scripts_regenerated':True,'reference_points_or_recipe_loaded':False,
        'seeded_not_equation_only':True,'curve_previously_used_in_method_development':previously_studied,
        'no_full_rank_computation':True,'independent_verification':claim,
        'audit_seconds':time.perf_counter()-start}
    save(output/'audit.json',report)
    hashes={str(p.relative_to(output)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(output.rglob('*')) if p.is_file()}
    save(output/'sha256.json',hashes)
    print(json.dumps(report,indent=2),flush=True)
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',required=True);p.add_argument('--output',required=True)
    args=p.parse_args();audit(Path(args.run),Path(args.output))
