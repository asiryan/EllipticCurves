"""Build search anchors as bounded short combinations of known points only."""
import argparse
from fractions import Fraction as Q
import json
from pathlib import Path
import subprocess
import time

from point_search import ROOT, vector, save
from point_arithmetic import combine, on_curve


def generate(source, directory, count=128, candidates=2048, timeout=60):
    data=json.loads(Path(source).read_text(encoding="utf-8-sig"))
    a=list(map(Q,data['ainvs']))
    points=[tuple(map(Q,p)) for p in data['points']]
    if not points or any(not on_curve(a,p) for p in points):
        raise ValueError('Invalid input points')
    directory.mkdir(parents=True,exist_ok=False)
    script='default(realprecision,80);\ndefault(parisizemax,536870912);\n'
    script+='E=ellinit('+vector(a)+');P=['+','.join(vector(p) for p in points)+'];\n'
    script+='G=ellheightmatrix(E,P);T=qflllgram(G);if(abs(matdet(T))!=1,error("nonunimodular LLL"));G2=T~*G*T;\n'
    script+=f'V=qfminim(G2,2*vecmax(vector(#P,i,G2[i,i])),{candidates},2);print("ENUMERATED ",V[1]);\n'
    script+=r'''{
for(k=1,matsize(V[3])[2],v=T*V[3][,k];S=[0];
  for(j=1,#P,if(v[j],S=elladd(E,S,ellmul(E,P[j],v[j]))));
  if(!ellisoncurve(E,S),error("anchor off curve"));
  print("ANCHOR ",[vector(2,j,Str(S[j])),Vec(v),v~*G*v]));
}
print("ANCHORS_END");quit;
'''
    (directory/'anchors.gp').write_text(script)
    start=time.perf_counter()
    p=subprocess.run([str(ROOT/'artifacts/native-validation/gp.exe'),'-fq','-s','64M'],
        input=script,text=True,capture_output=True,timeout=timeout)
    (directory/'stdout.txt').write_text(p.stdout)
    (directory/'stderr.txt').write_text(p.stderr)
    if p.returncode or 'ANCHORS_END' not in p.stdout or any('***' in l and 'Warning:' not in l for l in p.stderr.splitlines()):
        raise RuntimeError(p.stderr)
    rows=[json.loads(l[7:]) for l in p.stdout.splitlines() if l.startswith('ANCHOR ')]
    rows.sort(key=lambda r:r[2])
    selected=rows[:count]
    for coords, coeff, _ in selected:
        if combine(a,points,coeff)!=tuple(map(Q,coords)):
            raise ValueError('Independent exact anchor reconstruction failed')
    report={'ainvs':data['ainvs'],'points':[r[0] for r in selected],
        'vectors':[r[1] for r in selected],'approximate_heights':[r[2] for r in selected],
        'source':str(Path(source).resolve()),'source_point_count':len(points),
        'candidate_vectors_stored':len(rows),'selected_anchors':len(selected),
        'new_independent_points':0,'published_witness_list_loaded':False,
        'selection_is_bounded_and_heuristic':True,'seconds':time.perf_counter()-start}
    save(directory/'anchors.json',report)
    print(f"Prepared {len(selected)} verified anchors from {len(rows)} candidate combinations")
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',required=True)
    p.add_argument('--output',required=True)
    p.add_argument('--count',type=int,default=128)
    p.add_argument('--candidates',type=int,default=2048)
    p.add_argument('--timeout',type=float,default=60)
    a=p.parse_args()
    if not (1<=a.count<=a.candidates<=4096 and 0<a.timeout<=600):
        p.error('Invalid bounded anchor search limits')
    generate(a.input,Path(a.output),a.count,a.candidates,a.timeout)
