"""Diagnostic reference translates modulo known points, never search input.

Small representatives of a reference coset can be easier to recover than the
published representative itself. After a search actually finds R = Q + sum(cP),
the original Q can be recovered by subtracting the already known sum(cP).
"""
import argparse
from fractions import Fraction as Q
import json
import math
from pathlib import Path
import random
import subprocess

from guided_recovery import matches
from point_arithmetic import combine, add, on_curve
from point_search import ROOT, vector, save


def short_coset_vectors(gram, cross, center, trials):
    """Babai rounding and bounded greedy restarts; heights are only a heuristic."""
    n=len(gram);lower=[[0.0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1):
            v=gram[i][j]-sum(lower[i][k]*lower[j][k] for k in range(j))
            lower[i][j]=math.sqrt(v) if i==j else v/lower[j][j]
    first=[0]*n
    for i in reversed(range(n)):
        first[i]=round(center[i]-sum(lower[j][i]*(first[j]-center[j]) for j in range(i+1,n))/lower[i][i])
    rng=random.Random(20260914);answers=set()
    for trial in range(trials):
        z=first.copy()
        if trial:
            for _ in range(rng.randint(1,8)):
                z[rng.randrange(n)]+=rng.choice([-1,1])
        gradient=[cross[i]+sum(gram[i][j]*z[j] for j in range(n)) for i in range(n)]
        for _ in range(200):
            changes=[round(-gradient[i]/gram[i][i]) for i in range(n)]
            deltas=[2*changes[i]*gradient[i]+changes[i]**2*gram[i][i] for i in range(n)]
            i=min(range(n),key=lambda j:deltas[j])
            if deltas[i]>=-1e-8:break
            step=changes[i];z[i]+=step
            for j in range(n):gradient[j]+=step*gram[j][i]
        answers.add(tuple(z))
    return answers


def generate(source, reference, directory, candidates=2048, per_target=16):
    known=json.loads(Path(source).read_text())
    ref=json.loads(Path(reference).read_text())
    matched,_=matches(known,ref)
    missing=[i for i in range(1,len(ref['points'])+1) if i not in matched]
    directory.mkdir(parents=True,exist_ok=False)
    a=list(map(Q,known['ainvs']))
    basis=[tuple(map(Q,p)) for p in known['points']]
    targets=[tuple(map(Q,ref['points'][i-1])) for i in missing]
    script='default(realprecision,100);\ndefault(parisizemax,536870912);\n'
    script+='E=ellinit('+vector(a)+');P=['+','.join(vector(p) for p in basis)+'];\n'
    script+='W=['+','.join(vector(p) for p in targets)+'];\n'
    script+='G=ellheightmatrix(E,concat(P,W));n=#P;\n'
    script+=('B=matrix(n,n,i,j,G[i,j]);T=qflllgram(B);H=T~*B*T;'
             'print("BASIS ",[vector(n,i,Vec(T[i,])),vector(n,i,Vec(H[i,]))]);\n'
             'for(k=1,#W,c=T~*vector(n,i,G[i,n+k])~;'
             'print("COSET ",[k,Vec(c),Vec(matsolve(H,-c)),G[n+k,n+k]]));\n'
             'print("TRANSLATES_END");quit;\n')
    (directory/'diagnostic.gp').write_text(script)
    proc=subprocess.run([str(ROOT/'artifacts/native-validation/gp.exe'),'-fq','-s','64M'],
        input=script,text=True,capture_output=True,timeout=60)
    (directory/'stdout.txt').write_text(proc.stdout)
    (directory/'stderr.txt').write_text(proc.stderr)
    if proc.returncode or 'TRANSLATES_END' not in proc.stdout or any('***' in l and 'Warning:' not in l for l in proc.stderr.splitlines()):
        raise RuntimeError(proc.stderr)
    transform,gram=json.loads(next(l[6:] for l in proc.stdout.splitlines() if l.startswith('BASIS ')))
    rows=[]
    for line in proc.stdout.splitlines():
        if not line.startswith('COSET '):continue
        k,cross,center,base_height=json.loads(line[6:])
        for z in short_coset_vectors(gram,cross,center,candidates):
            v=[sum(transform[i][j]*z[j] for j in range(len(z))) for i in range(len(z))]
            height=base_height+2*sum(z[i]*cross[i] for i in range(len(z)))+sum(z[i]*gram[i][j]*z[j] for i in range(len(z)) for j in range(len(z)))
            rows.append([k,v+[1],height])
    chosen=[]; points=[]
    for k,index in enumerate(missing,1):
        options=sorted((r for r in rows if r[0]==k),key=lambda r:r[2])
        seen=set()
        for _,v,height in options:
            key=tuple(v)
            if key in seen:continue
            seen.add(key)
            p=add(a,targets[k-1],combine(a,basis,v[:-1]))
            if p is None or not on_curve(a,p):
                raise RuntimeError('Invalid translated reference')
            points.append(list(map(str,p)))
            chosen.append({'published_index':index,'known_coefficients':v[:-1],
                           'approximate_height':height})
            if len(seen)>=per_target:break
    result={'ainvs':known['ainvs'],'points':points,'translates':chosen,
        'known_basis':str(Path(source).resolve()),'reference':str(Path(reference).resolve()),
        'diagnostic_only':True,'published_targets_used':True,'not_search_results':True}
    save(directory/'witnesses.json',result)
    print(json.dumps({'missing':missing,'translated_witnesses':len(points),
        'best_heights':{i:min((r['approximate_height'] for r in chosen if r['published_index']==i),default=None) for i in missing}},indent=2),flush=True)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',required=True)
    p.add_argument('--reference',default=str(ROOT/'tests/Fixtures/icarm-302.json'))
    p.add_argument('--output',required=True)
    p.add_argument('--candidates',type=int,default=2048)
    p.add_argument('--per-target',type=int,default=16)
    args=p.parse_args()
    if not 1<=args.candidates<=4096 or not 1<=args.per_target<=64:p.error('Invalid limits')
    generate(args.input,args.reference,Path(args.output),args.candidates,args.per_target)
