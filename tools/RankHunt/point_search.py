"""Bounded PARI searches, retaining exact points and checking rank with C#.

pointed: project lines through known points to quartics, reduce, search, map back.
covers: compute actual locally soluble 2-covers using PARI ell2cover.
rank: PARI ellrank with supplied points, keeping its upper bound informational.
Every GP line is generated from validated rational numbers and bounded integers.
"""
import argparse
from fractions import Fraction
import json
from pathlib import Path
import re
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]


def save(path, data):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def vector(values):
    return "[" + ",".join(str(Fraction(v)) for v in values) + "]"


def quartic_reduction_code(minimal=False):
    if not minimal:
        return 'C=hyperellred(F,&m);'
    return ('C0=hyperellminimalmodel(F,&m0);C=hyperellred(C0,&m1);' +
            'md=m1[2][2,1]*x+m1[2][2,2];' +
            'mh=m0[1]*m1[3]+subst(m0[3],x,(m1[2][1,1]*x+m1[2][1,2])/md)*md^2;' +
            'm=[m0[1]*m1[1],m0[2]*m1[2],mh];')


def make_script(data, mode, height, anchors, effort, minimal=False, stack_mb=512,
                anchor_points=None, quartic_minimal=False, denominator_height=None):
    a = [Fraction(x) for x in data["ainvs"]]
    if len(a) != 5:
        raise ValueError("Expected five Weierstrass coefficients")
    points = [[Fraction(x), Fraction(y)] for x, y in (anchor_points if anchor_points is not None else data["points"])]
    search_denominator = height if denominator_height is None else denominator_height
    lines = [f"default(parisizemax,{stack_mb*1048576});", "default(realprecision,100);", "setrand(20260914);",
             f"E=ellinit({vector(a)});", "print(\"VERSION \",version());",
             "print(\"MODEL \",[E.a1,E.a2,E.a3,E.a4,E.a6]);",
             "P=[" + ",".join(vector(p) for p in points) + "];",
             'for(i=1,#P,if(!ellisoncurve(E,P[i]),error("Input point off curve")));']
    emit = 'if(!ellisoncurve(E,W),error("Mapped point off curve"));print("POINT ",W);'
    if minimal and mode != "direct":
        lines += ['Eoriginal=E;E=ellminimalmodel(Eoriginal,&modelchange);',
                  'P=vector(#P,i,ellchangepoint(P[i],modelchange));',
                  'print("MINIMAL_MODEL ",[E.a1,E.a2,E.a3,E.a4,E.a6]);']
        emit = ('if(!ellisoncurve(E,W),error("Mapped point off curve"));' +
                'W=ellchangepointinv(W,modelchange);' +
                'if(!ellisoncurve(Eoriginal,W),error("Minimal model inverse failed"));print("POINT ",W);')
    if mode == "pointed":
        # D(t) is the discriminant of the residual quadratic after intersecting
        # V=v0+t*(x-x0) with V^2=4*x^3+b2*x^2+2*b4*x+b6 and removing P.
        lines += [f"for(k=1,min(#P,{anchors})," +
            'print("ANCHOR_BEGIN ",k);' +
            'x0=P[k][1];v0=2*P[k][2]+E.a1*x0+E.a3;' +
            'D=x^4-2*(12*x0+E.b2)*x^2+32*v0*x+E.b2^2-8*E.b2*x0-48*x0^2-32*E.b4;' +
            'den=denominator(content(D));F=den^2*D;' +
            quartic_reduction_code(quartic_minimal) +
            'infq=polcoef(C[2],2);infp=polcoef(C[1],4);' +
            'if(m[2][2,1]!=0 && issquare(infq^2+4*infp,&infr),' +
            'inz=Set([(-infq+infr)/2,(-infq-infr)/2]);' +
            'for(infi=1,#inz,slope=m[2][1,1]/m[2][2,1];' +
            'square=(m[1]*inz[infi]+polcoef(m[3],2))/m[2][2,1]^2/den;' +
            'if(square^2!=subst(D,x,slope),error("Quartic infinity transformation failed"));' +
            'xx=(slope^2-E.b2-4*x0+square)/8;' +
            'yy=(v0+slope*(xx-x0)-E.a1*xx-E.a3)/2;W=[xx,yy];' + emit + '));' +
            f'H=hyperellratpoints(C,[{height},{search_denominator}]);' +
            'for(j=1,#H,h=H[j][1];z=H[j][2];dd=m[2][2,1]*h+m[2][2,2];' +
            'if(dd==0,next);slope=(m[2][1,1]*h+m[2][1,2])/dd;' +
            'square=(m[1]*z+subst(m[3],x,h))/dd^2/den;' +
            'if(square^2!=subst(D,x,slope),error("Quartic transformation failed"));' +
            'xx=(slope^2-E.b2-4*x0+square)/8;' +
            'yy=(v0+slope*(xx-x0)-E.a1*xx-E.a3)/2;W=[xx,yy];' + emit +
            ');print("ANCHOR_DONE ",k," ",#H));']
    elif mode == "covers":
        lines += ['print("COVERS_BEGIN");C=ell2cover(E);print("COVERS_READY ",#C);',
            f'for(k=1,min(#C,{anchors}),print("COVER_BEGIN ",k);R=C[k][1];M=C[k][2];' +
            f'H=hyperellratpoints(R,[{height},{height}]);' +
            'for(j=1,#H,if(H[j][2]==0,next);W=substvec(M,[x,y],H[j]);' + emit +
            ');print("COVER_DONE ",k," ",#H));']
    elif mode == "rank":
        lines += [f'print("RANK_BEGIN");R=ellrank(E,{effort},P);print("PARI_BOUNDS ",[R[1],R[2],R[3]]);',
                  'for(k=1,#R[4],W=R[4][k];' + emit + ');']
    elif mode == "halve":
        relations = data["relations"]
        if any(not isinstance(j, int) or isinstance(j, bool) or not 0 <= j < len(points)
               for relation in relations for j in relation):
            raise ValueError("Invalid relation indices")
        lines += ['L=[' + ','.join(vector([j+1 for j in relation]) for relation in relations) + '];',
                  'for(k=1,#L,R=[0];for(j=1,#L[k],R=elladd(E,R,P[L[k][j]]));' +
                  'if(R==[0],print("EXACT_RELATION ",k);next);' +
                  'if(ellisdivisible(E,R,2,&W),if(ellmul(E,W,2)!=R,error("Incorrect half"));' + emit + '));']
    elif mode == "direct":
        lines += ['E0=E;E=ellminimalmodel(E0,&change);',
                  f'H=ellratpoints(E,[{height},1]);',
                  'for(k=1,#H,W=ellchangepointinv(H[k],change);if(!ellisoncurve(E0,W),error("Direct map failed"));print("POINT ",W));']
    else:
        raise ValueError("Unknown mode")
    lines += ['print("SEARCH_END");quit;']
    return "\n".join(lines) + "\n"


def run(args):
    data = json.loads(Path(args.input).read_text(encoding="utf-8"))
    prefix = Path(args.output).resolve()
    prefix.parent.mkdir(parents=True, exist_ok=True)
    anchor_source = getattr(args, "anchor_input", None)
    anchor_points = None
    if anchor_source:
        if args.mode != "pointed":
            raise ValueError("Separate anchors apply only to pointed search")
        anchor_data = json.loads(Path(anchor_source).read_text(encoding="utf-8-sig"))
        if list(map(Fraction,anchor_data["ainvs"])) != list(map(Fraction,data["ainvs"])):
            raise ValueError("Anchor curve differs from the search curve")
        anchor_points = anchor_data["points"]
    quartic_minimal = getattr(args,"quartic_minimal",False)
    denominator_height = getattr(args,"denominator_height",None)
    if denominator_height is not None and (args.mode != "pointed" or not 1 <= denominator_height <= args.height):
        raise ValueError("A denominator bound requires pointed mode and 1 <= bound <= height")
    script = make_script(data, args.mode, args.height, args.anchors, args.effort, args.minimal, args.stack_mb,
                         anchor_points, quartic_minimal, denominator_height)
    prefix.with_suffix(".gp").write_text(script, encoding="utf-8")
    start = time.perf_counter()
    timed_out = False
    try:
        process = subprocess.run([str(Path(args.gp).resolve()), "-fq", "-s", "64M"], input=script,
                                 text=True, capture_output=True, timeout=args.timeout)
        stdout, stderr, exit_code = process.stdout, process.stderr, process.returncode
    except subprocess.TimeoutExpired as exc:
        timed_out = True
        stdout, stderr = exc.stdout or b"", exc.stderr or b""
        stdout = stdout.decode(errors="replace") if isinstance(stdout, bytes) else stdout
        stderr = stderr.decode(errors="replace") if isinstance(stderr, bytes) else stderr
        exit_code = None
    prefix.with_suffix(".stdout.txt").write_text(stdout, encoding="utf-8")
    prefix.with_suffix(".stderr.txt").write_text(stderr, encoding="utf-8")
    found = re.findall(r"^POINT \[(-?\d+(?:/\d+)?), (-?\d+(?:/\d+)?)\]$", stdout, re.M)
    a1, _, a3, _, _ = [Fraction(a) for a in data["ainvs"]]
    def normalize(p):
        x, y = map(Fraction, p)
        return str(x), str(min(y, -y-a1*x-a3))
    initial = {normalize(p) for p in data["points"]}
    merged = sorted(initial | {normalize(p) for p in found}, key=lambda p: (Fraction(p[0]), Fraction(p[1])))
    result = {"ainvs": [str(Fraction(a)) for a in data["ainvs"]], "points": merged,
              "source": str(Path(args.input).resolve()), "mode": args.mode,
              "minimal_model_preprocessing": args.minimal or args.mode == "direct",
              "pari_stack_limit_mb": args.stack_mb,
              "height": args.height, "anchors": args.anchors, "timeout_seconds": args.timeout,
              "seconds": round(time.perf_counter()-start, 3), "timed_out": timed_out,
              "process_exit_code": exit_code, "finished": "SEARCH_END" in stdout and exit_code == 0,
              "initial_distinct_up_to_sign": len(initial), "new_distinct_up_to_sign": len(merged)-len(initial),
              "new_points_are_not_automatically_independent": True}
    result["quartic_minimal_preprocessing"] = quartic_minimal
    result["denominator_height"] = denominator_height if denominator_height is not None else args.height
    if anchor_source:
        result["anchor_source"] = str(Path(anchor_source).resolve())
        result["available_anchor_count"] = len(anchor_points)
    bounds = re.search(r"^PARI_BOUNDS (\[.*\])$", stdout, re.M)
    if bounds:
        result["pari_reported_bounds"] = json.loads(bounds[1])
        result["pari_upper_bound_independently_certified"] = False
    # PARI may return 0 even after reporting an error and proceeding to quit.
    result["gp_errors_detected"] = any("***" in line and "Warning:" not in line for line in stderr.splitlines())
    result["finished"] = result["finished"] and not result["gp_errors_detected"]
    path = prefix.with_suffix(".json")
    save(path, result)
    verification = subprocess.run(["dotnet", str(ROOT / "tools/RankHunt/bin/Release/net8.0/RankHunt.dll"),
        "verify", "--input", str(path)], capture_output=True, text=True, timeout=60)
    if verification.returncode != 0:
        raise RuntimeError("C# exact certificate failed: " + verification.stderr)
    result["exact_certificate"] = json.loads(verification.stdout)
    save(path, result)
    print(json.dumps({k: v for k, v in result.items() if k not in ("ainvs", "points", "source")}, indent=2), flush=True)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True, help="Output prefix without extension")
    parser.add_argument("--gp", default=str(ROOT / "artifacts/native-validation/gp.exe"))
    parser.add_argument("--mode", choices=["pointed", "covers", "rank", "direct", "halve"], default="pointed")
    parser.add_argument("--height", type=int, default=10000)
    parser.add_argument("--anchors", type=int, default=32)
    parser.add_argument("--effort", type=int, default=1)
    parser.add_argument("--minimal", action="store_true", help="Minimize first and map all points back to the input model")
    parser.add_argument("--quartic-minimal", action="store_true", help="Minimize each pointed quartic before reduction")
    parser.add_argument("--anchor-input", help="Separate JSON of pointed-search anchors on the same curve")
    parser.add_argument("--denominator-height", type=int, help="Separate maximum denominator for pointed search")
    parser.add_argument("--stack-mb", type=int, default=512)
    parser.add_argument("--timeout", type=float, default=30)
    arguments = parser.parse_args()
    if not (1 <= arguments.height <= 10**10 and 1 <= arguments.anchors <= 1024 and
            0 <= arguments.effort <= 10 and 0 < arguments.timeout <= 3600 and 64 <= arguments.stack_mb <= 2048):
        parser.error("Invalid search limit")
    run(arguments)
