"""Preserve completed curve-search campaigns and an independently checked best curve."""
import argparse
from collections import Counter
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess

from point_search import ROOT,save,vector
from point_arithmetic import on_curve
from record_hunt import j_invariant,RECORD_A


def read(path):return json.loads(Path(path).read_text())
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def preserve(runs,output):
    output=Path(output).resolve()
    if output==ROOT or not output.is_relative_to(ROOT):raise ValueError('Use a dedicated workspace output directory')
    output.mkdir(parents=True,exist_ok=False);reports=[];best=None
    for index,name in enumerate(runs,1):
        run=Path(name).resolve();summary=read(run/'summary.json');verified=read(run/'best-verified.json')
        if summary['status'] not in {'campaign_completed','target_found'}:raise ValueError('Campaign is not complete')
        audit_dir=Path(verified['audit']);audit=read(audit_dir/'audit.json')
        if audit['certified_lower_bound']!=verified['certified_lower_bound']:raise ValueError('Conflicting result')
        for name,value in read(audit_dir/'sha256.json').items():
            if digest(audit_dir/name)!=value:raise ValueError('Audit artifact changed: '+name)
        config=read(run/'config.json')
        for name,value in config['controller_sha256'].items():
            if digest(run/'code'/name)!=value:raise ValueError('Controller snapshot changed')
        for name,value in config['selected_inputs'].items():
            if digest(run/'inputs'/(name+'.json'))!=value:raise ValueError('Candidate input changed')
        dest=output/f'campaign-{index}';dest.mkdir()
        for name in ['config.json','summary.json','leaderboard.json','best-verified.json','candidate-source.json']:
            shutil.copy2(run/name,dest/name)
        shutil.copytree(run/'inputs',dest/'inputs');shutil.copytree(run/'code',dest/'code')
        shutil.copytree(audit_dir,dest/'best-audit')
        rows=read(run/'leaderboard.json')
        report={'source':str(run),'tested_curves':len(rows),'best_lower_bound':verified['certified_lower_bound'],
            'distribution':dict(sorted(Counter(r['lower_bound'] for r in rows).items())),
            'controller_wall_seconds':summary['controller_wall_seconds'],
            'aggregate_point_search_seconds':summary['aggregate_point_search_seconds'],
            'best_candidate':verified['candidate'],'audit_seconds':audit['audit_seconds']}
        reports.append(report)
        if best is None or verified['certified_lower_bound']>best[0]['certified_lower_bound']:
            best=(verified,audit_dir,run)
    verified,audit_dir,run=best;source=read(audit_dir/'independent-points.json');bound=verified['certified_lower_bound']
    original=[tuple(map(Q,p)) for p in source['points']];a=list(map(Q,source['ainvs']))
    script='default(parisizemax,536870912);\nE0=ellinit('+vector(a)+');E=ellminimalmodel(E0,&m);\n'
    script+='P=['+','.join(vector(p) for p in original)+'];\n'
    script+='print("AINVS ",vector(5,j,Str(E[j])));print("CHANGE ",vector(4,j,Str(m[j])));\n'
    script+='for(i=1,#P,R=ellchangepoint(P[i],m);if(!ellisoncurve(E,R),error("Off curve"));print("POINTJSON ",vector(2,j,Str(R[j]))));print("MINIMAL_END");quit;\n'
    (output/'minimal-model.gp').write_text(script)
    p=subprocess.run([str(ROOT/'artifacts/native-validation/gp.exe'),'-fq','-s','64M'],
        input=script,text=True,capture_output=True,timeout=30)
    (output/'minimal-model.stdout.txt').write_text(p.stdout);(output/'minimal-model.stderr.txt').write_text(p.stderr)
    errors=any('***' in line and 'Warning:' not in line for line in p.stderr.splitlines())
    if p.returncode or 'MINIMAL_END' not in p.stdout or errors:raise ValueError('Minimal model conversion failed')
    lines=p.stdout.splitlines()
    b=list(map(Q,json.loads(next(l[6:] for l in lines if l.startswith('AINVS ')))))
    u,r,s,t=map(Q,json.loads(next(l[7:] for l in lines if l.startswith('CHANGE '))))
    points=[tuple(map(Q,json.loads(l[10:]))) for l in lines if l.startswith('POINTJSON ')]
    a1,a2,a3,a4,a6=a
    expected=[(a1+2*s)/u,(a2-s*a1+3*r-s*s)/u**2,(a3+r*a1+2*t)/u**3,
        (a4-s*a3+2*r*a2-(t+r*s)*a1+3*r*r-2*s*t)/u**4,
        (a6+r*a4+r*r*a2+r**3-t*a3-r*t*a1-t*t)/u**6]
    if b!=expected or len(points)!=len(original):raise ValueError('Wrong model change')
    for p,q in zip(original,points):
        if not on_curve(b,q) or (u*u*q[0]+r,u**3*q[1]+s*u*u*q[0]+t)!=p:
            raise ValueError('Exact point transport failed')
    data={'ainvs':list(map(str,b)),'points':[list(map(str,p)) for p in points],
        'rank_lower_bound':bound,'parameter':{k:verified['candidate'][k] for k in ['u','v']},
        'source_campaign':str(run),'model_change_from_family':list(map(str,[u,r,s,t])),
        'global_novelty_not_established':True}
    if j_invariant(data)==j_invariant({'ainvs':RECORD_A}):raise ValueError('Reference curve was not excluded')
    verifier=audit_dir/'independent_certificate_verifier.py'
    if digest(verifier)!='9e0d0d2562fc53705e92a2eaa9a3f6e7c923f1cd3fd82b14df68b60268f4ad54':raise ValueError('Unexpected verifier')
    spec=importlib.util.spec_from_file_location('saved_certificate',verifier)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    proof=module.build_certificate(data,max_prime=2000);claim=module.verify_certificate(data,proof)
    if claim['rank_lower_bound']!=bound or not claim['all_points_independent_modulo_torsion']:raise ValueError('Final independent proof failed')
    save(output/'best-curve.json',data);save(output/'best-certificate.json',proof);save(output/'best-verification.json',claim)
    shutil.copy2(verifier,output/'independent_certificate_verifier.py')
    summary={'campaigns':reports,'best_lower_bound':bound,'target':32,'target_reached':bound>=32,
        'independent_verification':claim,'different_j_from_reference':True,'global_novelty_not_established':True,
        'exact_model_change_checked':True,'no_full_rank_computation':True}
    save(output/'summary.json',summary)
    hashes={p.relative_to(output).as_posix():digest(p) for p in sorted(output.rglob('*')) if p.is_file()}
    save(output/'sha256.json',hashes)
    print(json.dumps({'best_lower_bound':bound,'target_reached':bound>=32,'ainvs':data['ainvs'],
        'independent_verification':claim,'output':str(output)},indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--runs',nargs='+',required=True);p.add_argument('--output',required=True)
    args=p.parse_args();preserve(args.runs,args.output)
