"""Find smaller parameters via the Jacobian of two compatible conics.

The quartic-to-cubic map and its inverse are checked exactly for every point.
PARI's auxiliary rank output is only a source of points; final ranks of the
target curves are certified by RankHunt bisections.
"""
import argparse
from fractions import Fraction
import json
from pathlib import Path
import subprocess

from bisection_hunt import ROOT, read, save, run_gp


def script_for(candidate, matrix, radius, effort, combinations, max_digits):
    coeff = "+".join(f"({Fraction(c)})*w^{i}" for i, c in enumerate(candidate["quartic_coefficients"]))
    x0, y0 = map(Fraction, candidate["quartic_seed"])
    mat = "[" + ";".join(",".join(str(Fraction(c)) for c in row) for row in matrix) + "]"
    return f"w='w;F={coeff};w0={x0};q0={y0};M={mat};radius={radius};effort={effort};box={int(combinations=='box')};maxdigits={max_digits};setrand(20260914);\n" + r'''
check(b,msg)=if(!b,error(msg));
f=subst(F,w,w+w0);check(polcoef(f,0)==q0^2,"quartic constant");
aa=polcoef(f,4);bb=polcoef(f,3);cc=polcoef(f,2);dd=polcoef(f,1);pp=dd/(2*q0);
E=ellinit([0,cc,0,bb*dd-4*aa*q0^2,bb^2*q0^2+aa*dd^2-4*aa*cc*q0^2]);
X0=pp^2-cc;Z0=-X0/(2*q0);Y0=dd*Z0-q0*bb;P0=[X0,Y0];
check(ellisoncurve(E,P0),"auxiliary seed");
proofY='proofY;proofX='proofX;proofZ=-proofX/(2*q0);proofU=(proofY-dd*proofZ+q0*bb)/(2*q0*(proofZ^2-aa));proofV=q0+pp*proofU+proofZ*proofU^2;
proofCubic=proofX^3+cc*proofX^2+(bb*dd-4*aa*q0^2)*proofX+bb^2*q0^2+aa*dd^2-4*aa*cc*q0^2;
check(Mod(numerator(proofV^2-subst(f,w,proofU)),proofY^2-proofCubic)==0,"symbolic inverse identity");
Emin=ellminimalmodel(E,&ch);Pmin=ellchangepoint(P0,ch);
print("AUXILIARY ",vector(5,i,Str(Emin[i])));print("SEED ",vector(2,i,Str(Pmin[i])));
R=ellrank(Emin,effort,[Pmin]);print("AUX_BOUNDS ",[R[1],R[2]]);print("GENERATORS ",#R[4]);
gens=R[4];if(#gens==0 && ellorder(Emin,Pmin)==0,gens=[Pmin]);
tors=elltors(Emin);TP=[[0]];for(i=1,#tors[3],old=TP;TP=concat(vector(tors[2][i],j,vector(#old,k,elladd(Emin,old[k],ellmul(Emin,tors[3][i],j-1))))));
halvings=0;for(i=1,#gens,for(step=1,8,foundhalf=0;for(j=1,#TP,target=elladd(Emin,gens[i],TP[j]);if(ellisdivisible(Emin,target,2,&halfpoint),check(ellmul(Emin,halfpoint,2)==target,"auxiliary exact half");gens[i]=halfpoint;halvings++;foundhalf=1;break));if(!foundhalf,break)));
print("HALVINGS ",halvings);
PTS=concat([Pmin],gens);for(i=1,#TP,if(#TP[i]==2,PTS=concat(PTS,[TP[i]])));print("AUX_POINTS ",vector(#PTS,i,vector(2,j,Str(PTS[i][j]))));
emit(P)={my(Q,X,Y,Z,u,v,ww,V,tt);if(#P==1,return());Q=ellchangepointinv(P,ch);X=Q[1];Y=Q[2];Z=-X/(2*q0);if(Z^2==aa,return());u=(Y-dd*Z+q0*bb)/(2*q0*(Z^2-aa));v=q0+pp*u+Z*u^2;ww=w0+u;check(v^2==subst(F,w,ww),"quartic inverse");V=M*[ww^2,ww,1]~;if(V[2]==0,return());tt=V[1]/V[2];if(max(#Str(abs(numerator(tt))),#Str(denominator(tt)))>maxdigits,return());print("PARAMETER ",[Str(tt),Str(ww),Str(v)]);};
if(box && #gens<=4,forvec(v=vector(#gens,i,[-radius,radius]),P=[0];for(i=1,#gens,P=elladd(Emin,P,ellmul(Emin,gens[i],v[i])));emit(P)),for(i=1,#gens,for(k=-radius,radius,if(k,emit(ellmul(Emin,gens[i],k)))));for(i=1,#gens,for(j=i+1,#gens,for(a=-radius,radius,if(a,for(b=-radius,radius,if(b,emit(elladd(Emin,ellmul(Emin,gens[i],a),ellmul(Emin,gens[j],b))))))))));
print("HUNT_END");quit;
'''


def main(args):
    data = read(args.input / "advanced.json")["candidates"]
    params = {p["index"]: p["matrix_descending"] for p in read(args.input / "parametrizations.json")["parametrizations"]}
    cases, found = [], {}
    for index in range(args.start, min(args.start+args.limit, len(data))):
        candidate = data[index]
        script = script_for(candidate, params[candidate["pair"][0]], args.radius, args.effort, args.combinations, args.max_parameter_digits)
        prefix = args.output / f"aux_{index:04}"
        try:
            stdout, seconds = run_gp(script, prefix, args.timeout)
        except subprocess.TimeoutExpired as exc:
            stdout = exc.stdout or b""
            stdout = stdout.decode(errors="replace") if isinstance(stdout,bytes) else stdout
            prefix.with_suffix(".stdout.txt").write_text(stdout, encoding="utf-8")
            cases.append({"index": index, "status": "timeout", "seconds": args.timeout})
            print(f"{index}: auxiliary rank search timed out", flush=True)
            continue
        count = 0
        for line in stdout.splitlines():
            if not line.startswith("PARAMETER "): continue
            t, w, v = json.loads(line[10:]); t = Fraction(t)
            key = (str(t.numerator),str(t.denominator))
            found[key] = {"pair": candidate["pair"], "parameter": list(key), "auxiliary_index": index,
                          "quartic_point": [w,v], "source": "small combinations of points on the auxiliary elliptic curve"}
            count += 1
        meta = {line.split(" ",1)[0]:json.loads(line.split(" ",1)[1]) for line in stdout.splitlines()
                if line.startswith(("AUXILIARY ","SEED ","AUX_BOUNDS ","GENERATORS ","HALVINGS ","AUX_POINTS "))}
        auxiliary = {"ainvs": meta["AUXILIARY"], "points": meta.pop("AUX_POINTS"),
                     "source": "auxiliary elliptic curve of the conic pair", "pair": candidate["pair"]}
        save(args.output / f"auxiliary_curve_{index:04}.json", auxiliary)
        validation = subprocess.run(["dotnet", str(ROOT / "tools/RankHunt/bin/Release/net8.0/RankHunt.dll"),
            "verify", "--input", str(args.output / f"auxiliary_curve_{index:04}.json")],capture_output=True,text=True,timeout=20)
        if validation.returncode: raise RuntimeError(validation.stderr)
        meta["exact_auxiliary_lower_bound"] = json.loads(validation.stdout)
        cases.append({"index": index, "status": "finished", "seconds": seconds, "raw_parameters": count, **meta})
        save(args.output / "advanced.json", {"candidates": list(found.values())})
        save(args.output / "auxiliary-summary.json", {"settings": {k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()}, "cases": cases, "unique_parameters": len(found)})
        print(f"{index}: {count} parameters, {seconds:.3f}s, auxiliary {meta.get('AUX_BOUNDS')}", flush=True)
    save(args.output / "advanced.json", {"candidates": list(found.values())})
    save(args.output / "auxiliary-summary.json", {"settings": {k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()}, "cases": cases, "unique_parameters": len(found)})
    # Keep the certifier input self-contained.
    save(args.output / "bisections.json", read(args.input / "bisections.json"))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", type=Path, default=ROOT / "artifacts/bisection-hunt")
    p.add_argument("--output", type=Path, default=ROOT / "artifacts/bisection-hunt/auxiliary")
    p.add_argument("--start", type=int, default=0)
    p.add_argument("--limit", type=int, default=4)
    p.add_argument("--radius", type=int, default=2)
    p.add_argument("--effort", type=int, default=1)
    p.add_argument("--combinations", choices=["pairwise","box"], default="pairwise")
    p.add_argument("--max-parameter-digits", type=int, default=200)
    p.add_argument("--timeout", type=float, default=30)
    a = p.parse_args()
    if not (a.start >= 0 and 1 <= a.limit <= 254 and 1 <= a.radius <= 8 and 0 <= a.effort <= 5 and 0 < a.timeout <= 120 and 1<=a.max_parameter_digits<=1000):
        p.error("Invalid finite search limits")
    main(a)
