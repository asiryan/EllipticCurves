"""Exact construction and bounded search of quadratic sections of ICARM #302.

Only arithmetic coefficient data is read from the audited package. PARI handles
polynomial arithmetic; every constructed bisection is verified by identities.
"""
import argparse
from fractions import Fraction
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools/RankStructureAudit"))
from audit import poly, read

LATTICE = ROOT / "tools/RankHunt/Data/icarm302-lattice.json"


def save(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2) + "\n", encoding="utf-8")


def run_gp(script, prefix, timeout=120):
    prefix.parent.mkdir(parents=True, exist_ok=True)
    prefix.with_suffix(".gp").write_text(script, encoding="ascii")
    start = time.perf_counter()
    proc = subprocess.run([str(ROOT / "artifacts/native-validation/gp.exe"), "-fq", "-s", "64M"],
                          input=script, capture_output=True, text=True, timeout=timeout)
    prefix.with_suffix(".stdout.txt").write_text(proc.stdout, encoding="utf-8")
    prefix.with_suffix(".stderr.txt").write_text(proc.stderr, encoding="utf-8")
    if proc.returncode or any("***" in line and "Warning:" not in line for line in proc.stderr.splitlines()):
        raise RuntimeError(proc.stderr)
    if "HUNT_END" not in proc.stdout:
        raise RuntimeError("PARI did not finish")
    return proc.stdout, time.perf_counter()-start


def prelude():
    data = read(ROOT / "tools/RankHunt/Data/icarm302-sections.json")
    points = ",".join(f'[{poly(p["x_coeffs"])},{poly(p["y_coeffs"])}]' for p in data["sections"])
    return "\n".join([
        "T='T; w='w; z='z;",
        'check(b,msg)=if(!b,error(msg));',
        "L=446667*T^2+471466*T+239031;",
        "p=318552*T^2+368554*T-72570; q=733413*T^2-45082*T-14960;",
        "D=5*(7174492962*T^4-7114589515*T^3-22069002960*T^2+3909144679*T-205134150);",
        "EE=882769396002*T^4+811447034567*T^3-1174040743*T^2-32493137198*T-2386325360;",
        "B=p*q*(L+p+q)-p*EE-q*D;",
        "A2=D+EE+L^2/4; A4=D*EE+L*B/2; A6=B^2/4;",
        f"seeds=[{points}];",
        "seeds=vector(17,i,[seeds[i][1],seeds[i][2]-(L*seeds[i][1]+B)/2]);",
        'coeffs(f)=vector(poldegree(f,T)+1,i,Str(polcoef(f,i-1,T)));',
        "padd(P,Q)={my(m,X,Y);if(#P==1,return(Q));if(#Q==1,return(P));if(P[1]==Q[1],if(P[2]+Q[2]==0,return([0]));m=(3*P[1]^2+2*A2*P[1]+A4)/(2*P[2]),m=(Q[2]-P[2])/(Q[1]-P[1]));X=m^2-A2-P[1]-Q[1];Y=-P[2]+m*(P[1]-X);return([X,Y]);};",
    ]) + "\n"


def construct(output, limit):
    lattice = read(LATTICE)
    vectors = lattice["bisection_vectors"][:limit]
    gram = lattice["gram"]
    for v in vectors:
        assert len(v) == 17 and all(x in (-1, 0, 1) for x in v)
        assert sum(v[i]*int(gram[i][j])*v[j] for i in range(17) for j in range(17)) == 10
    script = prelude() + "vectors=" + json.dumps(vectors) + ";\n" + r'''
makebis(v)={
  my(tau=[0],X,Y,dx,ll,NN,MM,aa,cc,qu,qv,delta,ff,hh=1,sf=1,unit,rr);
  for(j=1,17,if(v[j],tau=padd(tau,[seeds[j][1],v[j]*seeds[j][2]])));
  check(#tau==2,"zero trace");X=tau[1];Y=tau[2];
  dx=denominator(X);check(issquare(dx/pollead(dx),&ll),"trace denominator square");
  check(poldegree(ll,T)==3,"trace degree");
  NN=X*ll^2;MM=Y*ll^3;
  aa=lift(Mod(-MM,ll^2)/Mod(NN,ll^2));
  cc=(MM+aa*NN)/ll^2;
  qu=A2+(NN-aa^2)/ll^2;qv=A4+(2*aa*cc+NN*qu)/ll^2;
  check(ll^2*qu-NN==ll^2*A2-aa^2,"quadratic factor 2");
  check(ll^2*qv-NN*qu==ll^2*A4+2*aa*cc,"quadratic factor 1");
  check(-NN*qv==ll^2*A6-cc^2,"quadratic factor 0");
  delta=qu^2-4*qv;check(type(delta)=="t_POL","polynomial discriminant");
  ff=factor(delta);
  for(j=1,matsize(ff)[1],sf*=ff[j,1]^(ff[j,2]\2);if(ff[j,2]%2,hh*=ff[j,1]));
  unit=delta/(hh*sf^2);check(poldegree(unit,T)==0,"constant unit");unit=polcoef(unit,0,T);
  if(issquare(abs(unit),&rr),hh*=sign(unit);sf*=rr,hh*=numerator(unit)*denominator(unit);sf/=denominator(unit));
  check(delta==hh*sf^2,"squarefree discriminant");check(poldegree(hh,T)==2,"conic degree");
  return([coeffs(hh),coeffs(ll),coeffs(aa),coeffs(cc),coeffs(qu),coeffs(qv),coeffs(sf)]);
};
for(i=1,#vectors,bb=makebis(vectors[i]);print("BISECTION ",[i-1,bb]));
print("HUNT_END");quit;
'''
    stdout, seconds = run_gp(script, output / "construct", 120)
    rows = []
    for line in stdout.splitlines():
        if line.startswith("BISECTION "):
            index, values = json.loads(line[len("BISECTION "):])
            rows.append({"index": index, "trace_vector": vectors[index], **dict(zip(
                ["h", "l", "a", "c", "quadratic_u", "quadratic_v", "sqrt_factor"], values))})
    if len(rows) != len(vectors):
        raise RuntimeError("Missing bisections")
    save(output / "bisections.json", {"count": len(rows), "seconds": seconds, "bisections": rows})
    print(json.dumps({"constructed": len(rows), "seconds": seconds}, indent=2))


def parametrize(output):
    rows = read(output / "bisections.json")["bisections"]
    rows = list({tuple(r["h"]): r for r in reversed(rows)}.values())
    rows.sort(key=lambda r: r["index"])
    hs = "[" + ",".join("[" + ",".join(str(Fraction(c)) for c in r["h"]) + "]" for r in rows) + "]"
    script = "T='T;\nH=" + hs + ";\n" + r'''
check(b,msg)=if(!b,error(msg));
for(i=1,#H,h=H[i];Q=[2*h[3],h[2],0;h[2],2*h[1],0;0,0,-2];s=qfsolve(Q);if(type(s)!="t_COL",print("OBSTRUCTION ",[i-1,s]);next);M=qfparam(Q,s,3);F=M*[T^2,T,1]~;check(h[3]*F[1]^2+h[2]*F[1]*F[2]+h[1]*F[2]^2==F[3]^2,"parametrization identity");print("PARAM ",[i-1,vector(3,j,vector(3,k,Str(M[j,k]))),vector(3,j,Str(s[j]))]));
print("HUNT_END");quit;
'''
    stdout, seconds = run_gp(script, output / "parametrize", 120)
    params, obstructions = [], []
    for line in stdout.splitlines():
        if line.startswith("PARAM "):
            index, matrix, point = json.loads(line[6:])
            params.append({"index": rows[index]["index"], "matrix_descending": matrix, "conic_point": point})
        elif line.startswith("OBSTRUCTION "):
            obstructions.append(json.loads(line[12:]))
    save(output / "parametrizations.json", {"count": len(params), "seconds": seconds,
        "parametrizations": params, "obstructions": obstructions})
    print(json.dumps({"parametrized": len(params), "obstructions": len(obstructions), "seconds": seconds}, indent=2))


def intersections(output, limit):
    rows = read(output / "bisections.json")["bisections"]
    gram = [[int(c) for c in row] for row in read(LATTICE)["gram"]]
    ids = [p["index"] for p in read(output / "parametrizations.json")["parametrizations"]]
    pairs = []
    for pos, i in enumerate(ids):
        v = rows[i]["trace_vector"]
        for j in ids[pos+1:]:
            w = rows[j]["trace_vector"]
            dot = sum(v[a]*gram[a][b]*w[b] for a in range(17) for b in range(17))
            if abs(dot) == 7:
                pairs.append([i, j, 1 if dot > 0 else -1])
    pairs = pairs[:limit]
    script = "T='T; C=vector(" + str(len(rows)) + ");\n"
    for r in rows:
        script += f'C[{r["index"]+1}]=[' + ",".join(poly(r[k]) for k in
                    ["h", "l", "a", "c", "quadratic_u", "quadratic_v", "sqrt_factor"]) + "];\n"
    script += "pairs=" + json.dumps(pairs) + ";\n" + r'''
for(k=1,#pairs,pp=pairs[k];i=pp[1]+1;j=pp[2]+1;sg=pp[3];AA=C[i];BB=C[j];den=AA[3]*BB[2]-sg*BB[3]*AA[2];if(den==0,next);xx=(AA[4]*BB[2]-sg*BB[4]*AA[2])/den;R=gcd(numerator(xx^2+AA[5]*xx+AA[6]),numerator(xx^2+BB[5]*xx+BB[6]));ff=factor(R);for(a=1,matsize(ff)[1],if(poldegree(ff[a,1],T)==1,tt=-polcoef(ff[a,1],0,T)/polcoef(ff[a,1],1,T);if(subst(AA[2]*BB[2]*den,T,tt)==0,next);if(issquare(subst(AA[1],T,tt),&s1)&&issquare(subst(BB[1],T,tt),&s2),print("MEETING ",[pp,Str(tt),Str(s1),Str(s2)])))));
print("HUNT_END");quit;
'''
    stdout, seconds = run_gp(script, output / "intersections", 120)
    points = [json.loads(line[8:]) for line in stdout.splitlines() if line.startswith("MEETING ")]
    save(output / "intersections.json", {"pairs_checked": len(pairs), "seconds": seconds, "meetings": points})
    print(json.dumps({"pairs_checked": len(pairs), "rational_meetings": len(points), "seconds": seconds,
                      "first": points[:3]}, indent=2))


def advance(output, limit):
    rows = read(output / "bisections.json")["bisections"]
    params = {r["index"]: r for r in read(output / "parametrizations.json")["parametrizations"]}
    meetings = read(output / "intersections.json")["meetings"][:limit]
    script = "w='w; jobs=vector(" + str(len(meetings)) + ");\n"
    for k, (pair, t, _, _) in enumerate(meetings, 1):
        i, j, _ = pair
        matrix = "[" + ";".join(",".join(str(Fraction(c)) for c in r) for r in params[i]["matrix_descending"]) + "]"
        h = "[" + ",".join(str(Fraction(c)) for c in rows[j]["h"]) + "]"
        script += f"jobs[{k}]=[{i},{j},{Fraction(t)},{matrix},{h}];\n"
    script += r'''
check(b,msg)=if(!b,error(msg));
stepquartic(F,x0,y0)={my(b,c,d3,d4,dx,x1,y1);if(y0==0,return([]));b=subst(deriv(F),w,x0)/(2*y0);c=(subst(deriv(F,2),w,x0)/2-b^2)/(2*y0);d3=subst(deriv(F,3),w,x0)/6;d4=polcoef(F,4,w);if(c^2==d4,return([]));dx=(d3-2*b*c)/(c^2-d4);x1=x0+dx;y1=y0+b*dx+c*dx^2;check(y1^2==subst(F,w,x1),"quartic step");return([x1,y1]);};
for(k=1,#jobs,J=jobs[k];V=J[4]*[w^2,w,1]~;N=V[1];DD=V[2];h=J[5];F=h[1]*DD^2+h[2]*N*DD+h[3]*N^2;check(poldegree(F,w)==4 && poldisc(F)!=0,"smooth quartic");ff=factor(N-J[3]*DD);for(bi=1,matsize(ff)[1],if(poldegree(ff[bi,1],w)!=1,next);x0=-polcoef(ff[bi,1],0,w)/polcoef(ff[bi,1],1,w);check(issquare(subst(F,w,x0),&y0),"quartic seed square");P=stepquartic(F,x0,y0);if(#P==0,next);dn=subst(DD,w,P[1]);if(dn==0,next);tt=subst(N,w,P[1])/dn;print("ADVANCED ",[k-1,J[1],J[2],Str(tt),Str(x0),Str(y0),Str(P[1]),Str(P[2]),vector(5,a,Str(polcoef(F,a-1,w))) ])));
print("HUNT_END");quit;
'''
    # GP's deriv(f, variable) differs from repeated derivative notation.
    script = script.replace("deriv(F,2)", "deriv(deriv(F))").replace("deriv(F,3)", "deriv(deriv(deriv(F)))")
    stdout, seconds = run_gp(script, output / "advance", 120)
    candidates = []
    for line in stdout.splitlines():
        if not line.startswith("ADVANCED "): continue
        k, i, j, t, x0, y0, x1, y1, quartic = json.loads(line[9:])
        q = Fraction(t)
        candidates.append({"pair": [i, j], "parameter": [str(q.numerator), str(q.denominator)],
            "meeting_index": k, "quartic_coefficients": quartic, "quartic_seed": [x0, y0],
            "quartic_point": [x1, y1], "source": "osculating parabola from a rational bisection intersection"})
    unique = {tuple(c["parameter"]): c for c in candidates}
    save(output / "advanced.json", {"seconds": seconds, "raw_count": len(candidates), "candidates": list(unique.values())})
    print(json.dumps({"advanced": len(candidates), "unique_parameters": len(unique), "seconds": seconds,
                      "first_parameters": [c["parameter"] for c in list(unique.values())[:2]]}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["construct", "parametrize", "intersections", "advance"])
    parser.add_argument("--output", type=Path, default=ROOT / "artifacts/bisection-hunt")
    parser.add_argument("--limit", type=int, default=318)
    args = parser.parse_args()
    if not 1 <= args.limit <= 318:
        parser.error("limit must be between 1 and 318")
    if args.mode == "construct": construct(args.output, args.limit)
    elif args.mode == "parametrize": parametrize(args.output)
    elif args.mode == "intersections": intersections(args.output, args.limit)
    else: advance(args.output, args.limit)
