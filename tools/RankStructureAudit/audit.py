"""Independently check the archive's exact identities using PARI, without SymPy.

The archive is input data: no script from it is executed. Mathematical strings
are parsed as a restricted arithmetic AST before translation to GP. All outputs
are written outside the original archive tree.
"""
import ast
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]
ORIGINAL = ROOT / "artifacts/rank-structure-audit/original/elliptic_rank_structure"
STRUCTURE = ORIGINAL / "rank_structure"
OUTPUT = ROOT / "artifacts/rank-structure-audit/reproduced"


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def arithmetic(text, variables=None):
    """Translate integer arithmetic to GP, or evaluate it exactly with Fraction."""
    def visit(n):
        if isinstance(n, ast.Constant) and type(n.value) is int:
            return str(n.value) if variables is None else Fraction(n.value)
        if isinstance(n, ast.Name) and n.id in {"T", "U", "w"}:
            return n.id if variables is None else variables[n.id]
        if isinstance(n, ast.UnaryOp) and isinstance(n.op, (ast.UAdd, ast.USub)):
            v = visit(n.operand)
            if variables is None:
                return f"({'+' if isinstance(n.op, ast.UAdd) else '-'}{v})"
            return v if isinstance(n.op, ast.UAdd) else -v
        if isinstance(n, ast.BinOp) and isinstance(n.op, (ast.Add, ast.Sub, ast.Mult, ast.Div, ast.Pow)):
            a, b = visit(n.left), visit(n.right)
            if isinstance(n.op, ast.Pow):
                if not isinstance(n.right, ast.Constant) or type(n.right.value) is not int or not 0 <= n.right.value <= 100:
                    raise ValueError("Require a small nonnegative integer exponent")
            if variables is None:
                op = {ast.Add: "+", ast.Sub: "-", ast.Mult: "*", ast.Div: "/", ast.Pow: "^"}[type(n.op)]
                return f"({a}{op}{b})"
            if isinstance(n.op, ast.Add): return a + b
            if isinstance(n.op, ast.Sub): return a - b
            if isinstance(n.op, ast.Mult): return a * b
            if isinstance(n.op, ast.Div): return a / b
            return a ** int(b)
        raise ValueError(f"Non-arithmetic archive expression: {ast.dump(n)}")
    return visit(ast.parse(text, mode="eval").body)


def poly(coeffs):
    return "+".join(f"({Fraction(c)})*T^{i}" for i, c in enumerate(coeffs))


def evaluate(coeffs, t):
    result = Fraction(0)
    for c in reversed(coeffs):
        result = result * t + Fraction(c)
    return result


def section_points(basis, n, d):
    t = Fraction(n, d)
    return [[str(evaluate(p["x_coeffs"], t) * d**4), str(evaluate(p["y_coeffs"], t) * d**6)] for p in basis]


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    manifest = read(STRUCTURE / "PACKAGE_MANIFEST.json")
    for name, digest in manifest.items():
        path = (ORIGINAL / name).resolve()
        assert path.is_relative_to(ORIGINAL.resolve())
        assert hashlib.sha256(path.read_bytes()).hexdigest() == digest, name
    data = read(STRUCTURE / "general_sections43.json")
    clean = read(STRUCTURE / "general_sections43_unseeded.json")
    bc = read(STRUCTURE / "rank18_base_change.json")
    bi = read(STRUCTURE / "bisection_0.json")
    lines = [
        "default(parisizemax,536870912);",
        "U='U; T='T; w='w;",
        'check(b,msg)=if(!b,error(msg));',
        "L=446667*T^2+471466*T+239031;",
        "p=318552*T^2+368554*T-72570; q=733413*T^2-45082*T-14960;",
        "D=5*(7174492962*T^4-7114589515*T^3-22069002960*T^2+3909144679*T-205134150);",
        "EE=882769396002*T^4+811447034567*T^3-1174040743*T^2-32493137198*T-2386325360;",
        "B=p*q*(L+p+q)-p*EE-q*D;",
        "b2=L^2+4*(D+EE); b4=L*B+2*D*EE; b6=B^2; b8=(b2*b6-b4^2)/4;",
        "disc=-b2^2*b8-8*b4^3-27*b6^2+9*b2*b4*b6; c4=b2^2-24*b4;",
        'check(poldegree(disc,T)==24,"discriminant degree");',
        'check(poldegree(gcd(disc,deriv(disc,T)),T)==0,"squarefree discriminant");',
        'check(poldegree(gcd(c4,disc),T)==0,"c4 gcd");',
        "deg(f)=if(f==0,-100,poldegree(f,T));",
        "sectioncheck(P)=check(P[2]^2-(L*P[1]+B)*P[2]-P[1]*(P[1]+D)*(P[1]+EE)==0,\"section identity\");",
    ]
    for label, basis in [("P", data["basis"]), ("Q", clean["basis"])]:
        assert len(basis) == 17
        points = ",".join(f'[{arithmetic(p["x"])},{arithmetic(p["y"])}]' for p in basis)
        lines += [f"{label}=[{points}];", f"for(i=1,17,sectioncheck({label}[i]));"]
        for i, pt in enumerate(basis, 1):
            for j, axis in enumerate(["x", "y"], 1):
                lines.append(f'check({label}[{i}][{j}]==({poly(pt[axis+"_coeffs"])}),"coefficient data");')
                lines.append(f'check(deg({label}[{i}][{j}])<={4 if j == 1 else 6},"degree bound");')
    lines += [
        "pair(A,C)={my(dx=A[1]-C[1],dy=A[2]-C[2]);return(2-deg(gcd(dx,dy))-min(4-deg(dx),6-deg(dy)));};",
        "G=matrix(17,17,i,j,if(i==j,4,pair(P[i],P[j])));",
        "expected=" + "[" + ";".join(",".join(str(int(x)) for x in row) for row in data["gram"]) + "];",
        'check(G==expected,"height pairing"); check(matdet(G)==1092,"Gram determinant");',
        'mins=vector(17,k,matdet(matrix(k,k,i,j,G[i,j]))); check(vecmin(mins)>0,"positive definiteness");',
        'print("MINORS=",mins);',
    ]
    for name, value in {"phi": bc["T_of_w"], "uofw": bc["U_of_w"], "hh": bc["conic_h"],
                        "ll": bi["l"], "aa": bi["a"], "cc": bi["c"], "qu": bi["quadratic_u"],
                        "qv": bi["quadratic_v"], "tx": bi["tau_x"],
                        "XX": bc["extra_x_over_conic"], "yy": bc["extra_y_over_conic"],
                        "sf": "5*(" + bi["square_factor"] + ")/4024036"}.items():
        lines.append(f"{name}={arithmetic(value)};")
    lines += [
        'check(uofw^2==subst(hh,T,phi),"conic parametrization");',
        'check(poldisc(hh)!=0,"smooth conic");',
        'check(max(poldegree(numerator(phi),w),poldegree(denominator(phi),w))==2,"base change degree");',
        "A2=D+EE+L^2/4; A4=D*EE+L*B/2; A6=B^2/4; NN=tx*ll^2;",
        'check(ll^2*qu-NN==ll^2*A2-aa^2,"bisection coefficient 2");',
        'check(ll^2*qv-NN*qu==ll^2*A4+2*aa*cc,"bisection coefficient 1");',
        'check(-NN*qv==ll^2*A6-cc^2,"bisection coefficient 0");',
        'check(Mod(XX^2+qu*XX+qv,U^2-hh)==0,"extra x on quadratic");',
        'check(ll*(yy-(L*XX+B)/2)==aa*XX-cc,"extra y");',
        'check(XX==(-qu+sf*U)/2,"portable x formula");',
        'check(subst(hh,T,164518/924945)%19==3,"record misses conic");',
        'check(kronecker(3,19)==-1,"nonresidue");',
        'coeffs(f)=vector(poldegree(f,T)+1,i,Str(polcoef(f,i-1,T)));',
    ]
    for name in ["ll", "aa", "cc", "qu", "sf"]:
        lines.append(f'print("COEFF_{name}=",coeffs({name}));')
    lines += ['print("ALL_EXACT_IDENTITIES_PASSED");', "quit;"]
    gp_script = "\n".join(lines) + "\n"
    (OUTPUT / "identities.gp").write_text(gp_script, encoding="ascii")
    started = time.perf_counter()
    proc = subprocess.run([str(ROOT / "artifacts/native-validation/gp.exe"), "-q", "-f"],
                          input=gp_script, capture_output=True, text=True, timeout=90)
    seconds = time.perf_counter() - started
    (OUTPUT / "identities.stdout.txt").write_text(proc.stdout, encoding="utf-8")
    (OUTPUT / "identities.stderr.txt").write_text(proc.stderr, encoding="utf-8")
    assert proc.returncode == 0 and not any("***" in line and "Warning:" not in line for line in proc.stderr.splitlines()), proc.stderr
    assert "ALL_EXACT_IDENTITIES_PASSED" in proc.stdout
    values = dict(line.split("=", 1) for line in proc.stdout.splitlines() if "=" in line)
    for tag in ["specialization_-58_237", "specialization_164518_924945"]:
        d = read(STRUCTURE / (tag + ".json"))
        assert section_points(data["basis"], *d["parameter"]) == d["points"], tag
    # Check the asserted 17 + 14 decomposition, not merely independence of an
    # arbitrary list of 31 points on the same curve. The model change of scale
    # 3 is independently verified in RankPackageAudit.
    specialization = read(STRUCTURE / "specialization_164518_924945.json")
    record = read(ORIGINAL / "rank_search/record302.json")
    decomposition = read(STRUCTURE / "record_decomposition.json")
    a = list(map(Fraction, specialization["ainvs"]))
    b = list(map(Fraction, record["ainvs"]))
    scale = 3
    r = (scale**2 * (b[0]**2 + 4*b[1]) - (a[0]**2 + 4*a[1])) / 12
    s = (scale*b[0] - a[0]) / 2
    shift = (scale**3*b[2] - a[2] - r*a[0]) / 2
    mapped = [[str((Fraction(x)-r)/scale**2), str((Fraction(y)-s*(Fraction(x)-r)-shift)/scale**3)]
              for x, y in specialization["points"]]
    indices = decomposition["additional_published_point_indices_1_based"]
    assert len(set(indices)) == 14 and all(1 <= j <= 31 for j in indices)
    assert decomposition["ainvs"] == record["ainvs"]
    assert decomposition["points"] == mapped + [record["points"][j-1] for j in indices]
    d = read(STRUCTURE / "rank18_curve.json")
    t = arithmetic(bc["T_of_w"], {"w": Fraction(d["w"])})
    uu = arithmetic(bc["U_of_w"], {"w": Fraction(d["w"])})
    assert [t.numerator, t.denominator] == d["parameter"]
    points = section_points(data["basis"], t.numerator, t.denominator)
    points.append([str(arithmetic(bc[k], {"T": t, "U": uu}) * t.denominator**power)
                   for k, power in [("extra_x_over_conic", 4), ("extra_y_over_conic", 6)]])
    assert points == d["points"], "rank18 specialization"
    portable = {
        "source": "elliptic_rank_structure_20260914.zip, independently checked by tools/RankStructureAudit/audit.py",
        "source_sections_sha256": hashlib.sha256((STRUCTURE / "general_sections43.json").read_bytes()).hexdigest(),
        "sections": [{k: p[k] for k in ("x_coeffs", "y_coeffs")} for p in data["basis"]],
        "base_change": {name: json.loads(values["COEFF_" + gp]) for name, gp in
                        [("l", "ll"), ("a", "aa"), ("c", "cc"), ("quadratic_u", "qu"), ("sqrt_factor", "sf")]},
    }
    (OUTPUT / "icarm302-sections.json").write_text(json.dumps(portable, indent=2) + "\n", encoding="utf-8")
    clean_specialization = read(STRUCTURE / "specialization_-58_237.json")
    clean_specialization["points"] = section_points(clean["basis"], -58, 237)
    (OUTPUT / "unseeded-specialization.json").write_text(json.dumps(clean_specialization, indent=2) + "\n", encoding="utf-8")
    dll = ROOT / "tools/RankHunt/bin/Release/net8.0/RankHunt.dll"
    certs = {}
    for tag in ["specialization_-58_237", "specialization_164518_924945", "record_decomposition", "rank18_curve", "unseeded-specialization"]:
        path = (OUTPUT if tag == "unseeded-specialization" else STRUCTURE) / (tag + ".json")
        proc = subprocess.run(["dotnet", str(dll), "verify", "--input", str(path)], capture_output=True, text=True, timeout=30)
        assert proc.returncode == 0, proc.stderr
        certs[tag] = json.loads(proc.stdout)
        assert certs[tag]["LowerBound"] == (31 if tag == "record_decomposition" else 18 if tag == "rank18_curve" else 17)
    result = {"archive_manifest_files_verified": len(manifest), "exact_section_identities": 34,
              "discriminant_degree": 24, "squarefree_discriminant": True, "c4_coprime_to_discriminant": True,
              "gram_determinant": 1092, "positive_principal_minors": json.loads(values["MINORS"]),
              "base_change_identity_verified": True, "extra_section_identity_verified": True,
              "published_specialization_coordinates_match": True,
              "record_decomposition_17_plus_14_checked": True,
              "record_not_on_conic": {"prime": 19, "quadratic_nonresidue": 3},
              "independent_csharp_checks": certs, "symbolic_gp_seconds": seconds,
              "section_discovery_replayed": False, "saturation_computed": False}
    (OUTPUT / "audit.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
