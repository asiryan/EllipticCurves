"""Strictly bounded exploration of the height lattice of known points.

No qfminim: its storage limit does not bound its enumeration work. Here both
the number of discovered vectors and their coefficient sizes are capped.
Approximate heights affect priority only; every point identity is checked.
"""
from fractions import Fraction as Q
import heapq
import json
from pathlib import Path
import subprocess
import time

from point_arithmetic import combine, on_curve
from point_search import ROOT, vector, save


def bounded_vectors(gram, limit, count):
    n=len(gram);seen={};heap=[]
    def insert(v,score):
        if not any(v) or len(seen)>=limit:return
        if next(x for x in v if x)<0:v=tuple(-x for x in v)
        if v not in seen:
            seen[v]=score;heapq.heappush(heap,(score,v))
    for i in range(n):insert(tuple(int(i==j) for j in range(n)),gram[i][i])
    expanded=0
    while heap and len(seen)<limit:
        score,v=heapq.heappop(heap);expanded+=1
        gradient=[sum(gram[i][j]*v[j] for j in range(n) if v[j]) for i in range(n)]
        for i in range(n):
            for step in [-1,1]:
                w=list(v);w[i]+=step
                if abs(w[i])>3 or sum(map(abs,w))>12:continue
                insert(tuple(w),score+2*step*gradient[i]+gram[i][i])
                if len(seen)>=limit:break
            if len(seen)>=limit:break
    rows=sorted((score,v) for v,score in seen.items())[:count]
    return rows,{'vectors_considered':len(seen),'vectors_expanded':expanded,
                 'vector_limit':limit,'coefficient_bound':3,'l1_bound':12}


def call_gp(script,prefix,timeout):
    prefix.with_suffix('.gp').write_text(script)
    try:
        p=subprocess.run([str(ROOT/'artifacts/native-validation/gp.exe'),'-fq','-s','64M'],
            input=script,text=True,capture_output=True,timeout=timeout)
    except subprocess.TimeoutExpired as exc:
        for suffix,value in [('.stdout.txt',exc.stdout),('.stderr.txt',exc.stderr)]:
            value=value or '';value=value.decode(errors='replace') if isinstance(value,bytes) else value
            prefix.with_suffix(suffix).write_text(value)
        raise
    prefix.with_suffix('.stdout.txt').write_text(p.stdout);prefix.with_suffix('.stderr.txt').write_text(p.stderr)
    if p.returncode or 'POOL_END' not in p.stdout or any('***' in l and 'Warning:' not in l for l in p.stderr.splitlines()):
        raise RuntimeError('Bounded anchor preparation failed: '+p.stderr)
    return p.stdout


def generate(source,directory,count=2048,candidates=2048,timeout=60,selector=None):
    start=time.perf_counter();deadline=start+timeout
    data=json.loads(Path(source).read_text());a=list(map(Q,data['ainvs']))
    points=[tuple(map(Q,p)) for p in data['points']]
    if not points or any(not on_curve(a,p) for p in points):raise ValueError('Invalid known basis')
    directory.mkdir(parents=True,exist_ok=False)
    script='default(realprecision,100);\ndefault(parisizemax,536870912);\n'
    script+='E=ellinit('+vector(a)+');P=['+','.join(vector(p) for p in points)+'];\n'
    script+=('G=ellheightmatrix(E,P);T=qflllgram(G);if(abs(matdet(T))!=1,error("Bad LLL change"));H=T~*G*T;'
        'print("LATTICE ",[vector(#P,i,Vec(T[i,])),vector(#P,i,Vec(H[i,]))]);\n'
        'for(k=1,#P,S=[0];for(j=1,#P,if(T[j,k],S=elladd(E,S,ellmul(E,P[j],T[j,k]))));'
        'print("BASIS ",vector(2,j,Str(S[j]))));\nprint("POOL_END");quit;\n')
    stdout=call_gp(script,directory/'lattice',max(.1,deadline-time.perf_counter()))
    transform,gram=json.loads(next(l[8:] for l in stdout.splitlines() if l.startswith('LATTICE ')))
    short_basis=[tuple(map(Q,json.loads(l[6:]))) for l in stdout.splitlines() if l.startswith('BASIS ')]
    n=len(points)
    for j,p in enumerate(short_basis):
        if combine(a,points,[transform[i][j] for i in range(n)])!=p:raise ValueError('Exact LLL basis change failed')
    rows,stats=(selector or bounded_vectors)(gram,max(n,candidates*16),count)
    script='default(parisizemax,536870912);\nE=ellinit('+vector(a)+');\n'
    script+='P=['+','.join(vector(p) for p in short_basis)+'];\n'
    script+='V=['+','.join(vector(v) for _,v in rows)+'];\n'
    script+=('for(k=1,#V,S=[0];for(j=1,#P,if(V[k][j],S=elladd(E,S,ellmul(E,P[j],V[k][j]))));'
        'if(!ellisoncurve(E,S),error("Anchor off curve"));print("ANCHOR ",vector(2,j,Str(S[j]))));\n'
        'print("POOL_END");quit;\n')
    stdout=call_gp(script,directory/'anchors',max(.1,deadline-time.perf_counter()))
    anchors=[json.loads(l[7:]) for l in stdout.splitlines() if l.startswith('ANCHOR ')]
    if len(anchors)!=len(rows):raise ValueError('Incomplete anchor output')
    vectors=[]
    for p,(_,v) in zip(anchors,rows):
        if combine(a,short_basis,v)!=tuple(map(Q,p)):raise ValueError('Exact bounded anchor reconstruction failed')
        vectors.append([sum(transform[i][j]*v[j] for j in range(n)) for i in range(n)])
    report={'ainvs':data['ainvs'],'points':anchors,'vectors':vectors,
        'approximate_heights':[r[0] for r in rows],'source':str(Path(source).resolve()),
        'source_point_count':n,'selected_anchors':len(anchors),'new_independent_points':0,
        'published_witness_list_loaded':False,'selection_is_bounded_and_heuristic':True,
        'exact_lll_basis_checked':True,'enumeration':stats,'seconds':time.perf_counter()-start}
    save(directory/'anchors.json',report)
    print(f"Prepared {len(anchors)} anchors from {stats['vectors_considered']} bounded candidates in {report['seconds']:.3f}s",flush=True)
    return report
