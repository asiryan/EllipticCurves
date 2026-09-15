"""Bounded specialization checks and a durable summary of this audit.

Run audit.py, build RankHunt, then run RankHunt enrich before this script.
No open-ended rank search is started here.
"""
import hashlib
import importlib.util
import json
from math import gcd
from pathlib import Path
import subprocess
import time

from audit import ROOT, OUTPUT, STRUCTURE, read


def save(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def main():
    dll = ROOT / "tools/RankHunt/bin/Release/net8.0/RankHunt.dll"
    samples = OUTPUT / "samples"
    samples.mkdir(exist_ok=True)
    rows = []
    start = time.perf_counter()
    for n, d in [(n, d) for d in range(1, 4) for n in range(-3, 4) if gcd(n, d) == 1]:
        path = samples / f"w_{n}_{d}.json"
        proc = subprocess.run(["dotnet", str(dll), "export", "--family", "icarm302-18",
                               "--u", str(n), "--v", str(d), "--output", str(path)],
                              capture_output=True, text=True, timeout=20)
        if proc.returncode:
            raise RuntimeError(proc.stderr)
        data = read(path)
        rows.append({"w": f"{n}/{d}", "lower_bound": data["rank_lower_bound"],
                     "points": len(data["points"]), "file": str(path.relative_to(ROOT))})
    total = time.perf_counter() - start
    # Use the already audited, unchanged verifier from the first package to
    # preserve standalone certificates. C# verification above is independent.
    verifier_path = ROOT / "artifacts/rank-package-audit/original/elliptic_rank_search/certificate.py"
    expected = "9e0d0d2562fc53705e92a2eaa9a3f6e7c923f1cd3fd82b14df68b60268f4ad54"
    if hashlib.sha256(verifier_path.read_bytes()).hexdigest() != expected:
        raise RuntimeError("Certificate verifier differs from the previously audited version")
    spec = importlib.util.spec_from_file_location("audited_certificate", verifier_path)
    verifier = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(verifier)
    destination = ROOT / "results/rank-structure-20260914"
    checks = {}
    for name, source in [
        ("rank18", OUTPUT / "our-rank18.json"),
        ("candidate-5037-2518", ROOT / "artifacts/rank-structure-audit/enriched-h10000/curve_-5037_2518.json"),
    ]:
        data = read(source)
        proof = verifier.build_certificate(data, max_prime=1009)
        checks[name] = verifier.verify_certificate(data, proof)
        save(destination / (name + ".json"), data)
        save(destination / (name + ".certificate.json"), proof)
    assert checks["rank18"]["rank_lower_bound"] == 18
    assert checks["candidate-5037-2518"]["rank_lower_bound"] == 17
    attempts = {}
    for name in ["top1", "rank18"]:
        path = ROOT / f"artifacts/rank-structure-audit/search/{name}.json"
        if path.exists():
            d = read(path)
            attempts[name] = {k: v for k, v in d.items() if k not in {"ainvs", "points", "source"}}
    result = {"specializations": rows, "sample_seconds_including_cli_startup": total,
              "independent_python_certificates": checks, "bounded_point_search_attempts": attempts}
    save(OUTPUT / "experiments.json", result)
    save(destination / "experiments.json", result)
    save(destination / "audit.json", read(OUTPUT / "audit.json"))
    save(destination / "enrichment.json", read(ROOT / "artifacts/rank-structure-audit/enriched-h10000/enrichment.json"))
    print(json.dumps({"samples": len(rows), "bounds": sorted(set(r["lower_bound"] for r in rows)),
                      "seconds": total, "certificate_checks": checks}, indent=2))


if __name__ == "__main__":
    main()
