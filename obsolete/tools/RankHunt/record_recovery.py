"""Search for 31 independent points on ICARM302 from 17 generic sections.

The published witness list is never loaded. This is a seeded point-search
experiment on a fixed, known curve, not an equation-only or parameter search.
Only point search and exact lower-bound certificates run; there is no full
rank computation or global Selmer-cover construction in this workflow.
"""
import argparse
from datetime import datetime
from fractions import Fraction as Q
import hashlib
import json
from math import isqrt
from pathlib import Path
import subprocess
import time
from types import SimpleNamespace

from point_search import ROOT, run, save
from point_arithmetic import combine, on_curve


RECORD = list(map(Q, [1, 1, 1,
    -1284727764113567728281797636015784768866707681415849262157224232063,
    560368321454261339256859338901915312332769858684945406858043869199456710681989058863306170127006181]))
DLL = ROOT / "tools/RankHunt/bin/Release/net8.0/RankHunt.dll"


def cli(*args):
    proc = subprocess.run(["dotnet", str(DLL), *map(str, args)],
                          capture_output=True, text=True, timeout=60)
    if proc.returncode:
        raise RuntimeError(proc.stderr or proc.stdout)
    return proc.stdout


def align_basis(points):
    change = json.loads((ROOT / "tools/RankHunt/Data/icarm302-record-basis.json").read_text())
    forward = [[int(c) for c in row] for row in change["forward"]]
    inverse = [[int(c) for c in row] for row in change["inverse"]]
    n = len(points)
    if len(forward) != n or len(inverse) != n or any(len(r) != n for r in forward+inverse):
        raise ValueError("Invalid basis matrix size")
    if any(sum(forward[i][k]*inverse[k][j] for k in range(n)) != int(i==j)
           for i in range(n) for j in range(n)):
        raise ValueError("The basis change is not an integer inverse")
    original = [tuple(map(Q,p)) for p in points]
    aligned = [combine(RECORD,original,row) for row in inverse]
    if any(p is None or not on_curve(RECORD,p) for p in aligned):
        raise ValueError("Invalid aligned point")
    if [combine(RECORD,aligned,row) for row in forward] != original:
        raise ValueError("Exact inverse basis change failed")
    return [list(map(str,p)) for p in aligned], change["published_indices"]


def prepare(directory, align=True):
    family_path = directory / "family17.json"
    cli("export", "--family", "icarm302-17", "--u", "164518", "--v", "924945",
        "--output", family_path)
    family = json.loads(family_path.read_text())
    a1, a2, a3, a4, a6 = map(Q, family["ainvs"])
    b1, b2, b3, _, _ = RECORD
    u = Q(3)
    r = (u*u*(b1*b1+4*b2)-(a1*a1+4*a2))/12
    s = (u*b1-a1)/2
    t = (u**3*b3-a3-r*a1)/2
    transformed = [(a1+2*s)/u, (a2-s*a1+3*r-s*s)/u**2,
        (a3+r*a1+2*t)/u**3,
        (a4-s*a3+2*r*a2-(t+r*s)*a1+3*r*r-2*s*t)/u**4,
        (a6+r*a4+r*r*a2+r**3-t*a3-r*t*a1-t*t)/u**6]
    if transformed != RECORD:
        raise ValueError("The family-to-record model change failed")
    points = []
    for x, y in family["points"]:
        x, y = Q(x), Q(y)
        xx, yy = (x-r)/u**2, (y-s*(x-r)-t)/u**3
        if yy*yy+b1*xx*yy+b3*yy != xx**3+b2*xx*xx+RECORD[3]*xx+RECORD[4]:
            raise ValueError("A mapped generic section is off the record curve")
        points.append([str(xx), str(yy)])
    if len(points) != 17:
        raise ValueError("Expected exactly 17 generic sections")
    save(directory / "generic17.json", {"ainvs":list(map(str,RECORD)),"points":points})
    indices = None
    if align:
        points, indices = align_basis(points)
    source = directory / "input17.json"
    save(source, {"ainvs": list(map(str, RECORD)), "points": points,
        "source": "17 generic sections at T=164518/924945, mapped to the published model",
        "model_change": {"u": str(u), "r": str(r), "s": str(s), "t": str(t)},
        "published_witness_list_loaded": False,
        "reference_used_for_basis_alignment": align,
        "aligned_published_indices": indices})
    cert = json.loads(cli("verify", "--input", source))
    if cert["LowerBound"] != 17:
        raise ValueError("Initial independent rank certificate is not 17")
    save(directory / "input17.certificate.json", cert)
    return source, cert


def check_bisections(source, directory):
    data = json.loads(source.read_text())
    unique = {tuple(b["h"]): b for b in data["bisections"]}
    t = Q(164518, 924945)
    hits = []
    for h, b in unique.items():
        value = sum(Q(c)*t**i for i, c in enumerate(h))
        if (value >= 0 and isqrt(value.numerator)**2 == value.numerator
                and isqrt(value.denominator)**2 == value.denominator):
            hits.append(b["index"])
    report = {"source": str(source), "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "parameter": str(t), "total_bisections": len(data["bisections"]),
        "distinct_square_conditions": len(unique), "rational_conditions_at_record": hits,
        "published_witness_list_loaded": False}
    save(directory / "bisection-check.json", report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default=str(ROOT / "artifacts/record-recovery" / datetime.now().strftime("%Y%m%d-%H%M%S")))
    parser.add_argument("--pointed-seconds", type=float, default=60)
    parser.add_argument("--height", type=int, default=100000)
    parser.add_argument("--anchors", type=int, default=64)
    parser.add_argument("--anchor-input", help="Optional separate anchor pool from known-point combinations")
    parser.add_argument("--quartic-minimal", action="store_true")
    parser.add_argument("--original-basis", action="store_true", help="Keep the original generic section basis")
    parser.add_argument("--bisections", default=str(ROOT / "artifacts/bisection-hunt/bisections.json"),
                        help="Optional existing pool to diagnose; missing file skips this check")
    args = parser.parse_args()
    if not (0 < args.pointed_seconds <= 3600 and 1 <= args.height <= 10**10
            and 1 <= args.anchors <= 1024):
        parser.error("Invalid point-search timeout, height or anchor limit")
    directory = Path(args.output).resolve()
    directory.mkdir(parents=True, exist_ok=False)
    start = time.perf_counter()
    source, initial_cert = prepare(directory, align=not args.original_basis)
    summary = {"started_at": datetime.now().isoformat(), "target_lower_bound": 31,
        "input": str(source), "input_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "published_witness_list_loaded": False, "parameter_was_supplied": True,
        "reference_used_for_basis_alignment": not args.original_basis,
        "workflow": "point_search_and_lower_bound_certificate_only",
        "initial_certificate": initial_cert, "best_lower_bound": 17,
        "phases": [], "status": "running"}
    bisections = Path(args.bisections).resolve()
    if bisections.exists():
        summary["bisection_check"] = check_bisections(bisections, directory)
    save(directory / "summary.json", summary)
    print(f"Output: {directory}\nInitial 17 generic points certified. Published witnesses are not input.", flush=True)
    prefix = directory / "pointed"
    print(f"BEGIN pointed: budget {args.pointed_seconds:g}s, current lower bound 17", flush=True)
    summary["active_phase"] = "pointed"
    save(directory / "summary.json", summary)
    result = run(SimpleNamespace(input=str(source), output=str(prefix), mode="pointed",
        height=args.height, anchors=args.anchors, minimal=True, stack_mb=512,
        anchor_input=args.anchor_input, quartic_minimal=args.quartic_minimal,
        effort=0, timeout=args.pointed_seconds,
        gp=str(ROOT / "artifacts/native-validation/gp.exe")))
    cert = result["exact_certificate"]
    summary["phases"].append({"mode": "pointed", "result": str(prefix.with_suffix('.json')),
        **{k: result[k] for k in ["seconds", "timed_out", "finished", "gp_errors_detected",
            "initial_distinct_up_to_sign", "new_distinct_up_to_sign", "exact_certificate"]}})
    summary["best_lower_bound"] = max(summary["best_lower_bound"], cert["LowerBound"])
    summary["best_points_file"] = str(prefix.with_suffix(".json"))
    summary.pop("active_phase", None)
    summary["status"] = "target_reached" if summary["best_lower_bound"] >= 31 else "budget_completed"
    summary["seconds"] = round(time.perf_counter()-start, 3)
    save(directory / "summary.json", summary)
    print(f"DONE: lower bound {summary['best_lower_bound']}; {summary['seconds']}s", flush=True)


if __name__ == "__main__":
    main()
