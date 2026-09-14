"""Diagnose known witnesses on pointed quartics; never supply them to a search."""
import argparse
import json
from pathlib import Path
import subprocess
import time

from point_search import ROOT, vector, save, quartic_reduction_code


def audit(source, witness, directory, timeout=120, height=100000, quartic_minimal=False):
    data = json.loads(Path(source).read_text(encoding="utf-8-sig"))
    target = json.loads(Path(witness).read_text(encoding="utf-8-sig"))
    if data["ainvs"] != target["ainvs"]:
        raise ValueError("Audit inputs must use exactly the same curve model")
    directory.mkdir(parents=True, exist_ok=False)
    script = "\n".join([
        'default(realprecision,100); default(parisizemax,536870912);',
        f'E=ellinit({vector(data["ainvs"])});',
        'P=['+','.join(vector(p) for p in data["points"])+'];',
        'W=['+','.join(vector(p) for p in target["points"])+'];',
        'check(b,msg)=if(!b,error(msg));',
        'for(i=1,#P,check(ellisoncurve(E,P[i]),"seed off curve"));',
        'for(i=1,#W,check(ellisoncurve(E,W[i]),"witness off curve"));',
        'roundtrips=0;besth=vector(#W,j,0);bestcost=vector(#W,j,0);',
        r'''{
for(k=1,#P,
  x0=P[k][1];v0=2*P[k][2]+E.a1*x0+E.a3;
  D=x^4-2*(12*x0+E.b2)*x^2+32*v0*x+E.b2^2-8*E.b2*x0-48*x0^2-32*E.b4;
  den=denominator(content(D));F=den^2*D; REDUCE_QUARTIC
  for(j=1,#W,for(si=0,1,
    Q=if(si,ellneg(E,W[j]),W[j]);
    if(Q[1]==x0,print("EXCEPTION ",[k,j,si,"anchor_abscissa"]);next);
    slope=(2*Q[2]+E.a1*Q[1]+E.a3-v0)/(Q[1]-x0);
    sq=8*Q[1]-(slope^2-E.b2-4*x0);
    check(sq^2==subst(D,x,slope),"inverse discriminant");
    dd=m[2][2,1]*slope-m[2][1,1];
    if(dd==0,
      z=(den*sq*m[2][2,1]^2-polcoef(m[3],2))/m[1];
      check(z^2+polcoef(C[2],2)*z==polcoef(C[1],4),"inverse infinity");
      backt=m[2][1,1]/m[2][2,1];
      backsq=(m[1]*z+polcoef(m[3],2))/m[2][2,1]^2/den;
      xx=(backt^2-E.b2-4*x0+backsq)/8;
      yy=(v0+backt*(xx-x0)-E.a1*xx-E.a3)/2;
      check([xx,yy]==Q,"exact infinity round trip");
      print("INFINITY_PREIMAGE ",[k,j,si,Str(z)]);next);
    h=(m[2][1,2]-m[2][2,2]*slope)/dd;
    dd=m[2][2,1]*h+m[2][2,2];
    z=(den*sq*dd^2-subst(m[3],x,h))/m[1];
    check(z^2+subst(C[2],x,h)*z==subst(C[1],x,h),"reduced inverse");
    backt=(m[2][1,1]*h+m[2][1,2])/dd;
    backsq=(m[1]*z+subst(m[3],x,h))/dd^2/den;
    xx=(backt^2-E.b2-4*x0+backsq)/8;
    yy=(v0+backt*(xx-x0)-E.a1*xx-E.a3)/2;
    check([xx,yy]==Q,"exact round trip");
    roundtrips++;hh=max(abs(numerator(h)),denominator(h));hc=hh*denominator(h);
    EMIT_PREIMAGE
  ));
);
}
print("ROUND_TRIPS ",roundtrips);print("AUDIT_END");quit;
''']).replace('REDUCE_QUARTIC',quartic_reduction_code(quartic_minimal))
    compact=len(data['points'])>128
    emit='print("PREIMAGE ",[k,j,si,Str(h),Str(z),Str(hh)]);'
    if compact:
        emit='if(besth[j]==0 || hh<besth[j] || hc<bestcost[j],besth[j]=if(besth[j],min(besth[j],hh),hh);bestcost[j]=if(bestcost[j],min(bestcost[j],hc),hc);'+emit+');'
    script=script.replace('EMIT_PREIMAGE',emit)
    (directory/"audit.gp").write_text(script, encoding="ascii")
    start=time.perf_counter()
    proc=subprocess.run([str(ROOT/'artifacts/native-validation/gp.exe'),'-fq','-s','64M'],
                        input=script, text=True, capture_output=True, timeout=timeout)
    (directory/'stdout.txt').write_text(proc.stdout)
    (directory/'stderr.txt').write_text(proc.stderr)
    if proc.returncode or 'AUDIT_END' not in proc.stdout or any('***' in line and 'Warning:' not in line for line in proc.stderr.splitlines()):
        raise RuntimeError(proc.stderr)
    rows=[json.loads(line[9:]) for line in proc.stdout.splitlines() if line.startswith('PREIMAGE ')]
    exceptions=[json.loads(line[10:]) for line in proc.stdout.splitlines() if line.startswith('EXCEPTION ')]
    infinities=[json.loads(line[18:]) for line in proc.stdout.splitlines() if line.startswith('INFINITY_PREIMAGE ')]
    best=[]
    areas=[]
    for index in range(1,len(target['points'])+1):
        candidates=[r for r in rows if r[1]==index]
        r=min(candidates,key=lambda r:int(r[5])) if candidates else None
        best.append({'witness_index':index,'best_anchor':r[0] if r else None,
                     'sign_flipped':bool(r[2]) if r else None,'minimum_height':r[5] if r else None,
                     'inside_search_box':int(r[5])<=height if r else False})
        from fractions import Fraction
        if candidates:
            r=min(candidates,key=lambda r:int(r[5])*Fraction(r[3]).denominator)
            h=Fraction(r[3])
            areas.append({'witness_index':index,'best_anchor':r[0],
                'numerator_bound':str(max(abs(h.numerator),h.denominator)),
                'denominator_bound':str(h.denominator),
                'search_area':str(max(abs(h.numerator),h.denominator)*h.denominator)})
    report={'diagnostic_only':True,'quartic_minimal_preprocessing':quartic_minimal,'seed_input':str(Path(source).resolve()),
            'witness_reference':str(Path(witness).resolve()),'seconds':time.perf_counter()-start,
            'seeds':len(data['points']),'witnesses':len(target['points']),
            'exact_round_trips':int(next(l[12:] for l in proc.stdout.splitlines() if l.startswith('ROUND_TRIPS '))),
            'all_preimages_saved':not compact,'search_height':height,
            'exact_infinity_round_trips':len(infinities),'infinity_preimages':infinities,
            'witnesses_inside_search_box':sum(r['inside_search_box'] for r in best),
            'minimum_required_heights':best,'minimum_search_areas':areas,'exceptions':exceptions,'preimages':rows}
    save(directory/'audit.json',report)
    print(json.dumps({k:v for k,v in report.items() if k not in ('exceptions','preimages')},indent=2))
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',required=True)
    p.add_argument('--witness',default=str(ROOT/'tests/Fixtures/icarm-302.json'))
    p.add_argument('--output',required=True)
    p.add_argument('--timeout',type=float,default=120)
    p.add_argument('--height',type=int,default=100000)
    p.add_argument('--quartic-minimal',action='store_true')
    a=p.parse_args()
    audit(a.input,a.witness,Path(a.output),a.timeout,a.height,a.quartic_minimal)
