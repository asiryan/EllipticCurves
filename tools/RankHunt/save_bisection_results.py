"""Preserve this campaign's results and independently check its certificates."""
import hashlib
import importlib.util
import json
from pathlib import Path

from bisection_hunt import ROOT, read, save


def main():
    campaign = ROOT / "artifacts/bisection-hunt"
    output = ROOT / "results/rank19-20260914"
    verifier_path = ROOT / "artifacts/rank-package-audit/original/elliptic_rank_search/certificate.py"
    if hashlib.sha256(verifier_path.read_bytes()).hexdigest() != "9e0d0d2562fc53705e92a2eaa9a3f6e7c923f1cd3fd82b14df68b60268f4ad54":
        raise RuntimeError("The independently audited Python verifier has changed")
    spec = importlib.util.spec_from_file_location("audited_certificate", verifier_path)
    verifier = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(verifier)
    checks = {}
    for name, path, lower in [
        ("rank19", campaign / "aux-small/certified/curve_0012.json", 19),
        ("auxiliary", campaign / "aux-box/auxiliary_curve_0251.json", 2),
    ]:
        data = read(path)
        proof = verifier.build_certificate(data, max_prime=2000)
        claim = verifier.verify_certificate(data, proof)
        if claim["rank_lower_bound"] < lower:
            raise RuntimeError("Independent certificate did not reach the claimed lower bound")
        save(output / (name + ".json"), data)
        save(output / (name + ".certificate.json"), proof)
        checks[name] = claim
    save(output / "proof-checks.json", checks)
    old = read(campaign / "certified/summary.json")
    fast = read(campaign / "certified-fast/summary.json")
    if old["count"] != fast["count"]:
        raise RuntimeError("Different benchmark input counts")
    for row in old["rows"]:
        if read(campaign / "certified" / row["file"]) != read(campaign / "certified-fast" / row["file"]):
            raise RuntimeError("Optimization changed the exact output: " + row["file"])
    phases = {name: read(campaign / (name + ".json"))["seconds"]
              for name in ["bisections", "parametrizations", "intersections", "advanced"]}
    summary = {"construction_seconds": phases, "initial_certification_seconds": old["seconds"],
        "optimized_certification_seconds": fast["seconds"], "speedup_in_this_comparison": old["seconds"]/fast["seconds"],
        "all_exact_outputs_match": True, "main_parameter_count": fast["count"],
        "main_best_certified_lower_bound": fast["best_lower_bound"], "main_results": fast["rows"],
        "auxiliary_search": read(campaign / "aux-small/auxiliary-summary.json"),
        "compact_parameter": read(output / "rank19.json")["parameter"],
        "no_upper_rank_bound_proved": True, "novelty_not_established": True}
    save(output / "search-summary.json", summary)
    pair = [283, 292]
    data = read(campaign / "bisections.json")
    parametrization = next(p for p in read(campaign / "parametrizations.json")["parametrizations"] if p["index"] == pair[0])
    original = read(campaign / "advanced.json")["candidates"][251]
    if original["pair"] != pair:
        raise RuntimeError("Pair provenance changed")
    save(output / "pair19.json", {"bisections": [b for b in data["bisections"] if b["index"] in pair],
        "parametrization": parametrization, "quartic_coefficients": original["quartic_coefficients"],
        "quartic_seed": original["quartic_seed"], "section_data_sha256": hashlib.sha256(
            (ROOT / "tools/RankHunt/Data/icarm302-sections.json").read_bytes()).hexdigest(),
        "generic_rank_lower_bound_over_auxiliary_curve": 19,
        "auxiliary_rational_rank_lower_bound": checks["auxiliary"]["rank_lower_bound"]})
    attempts = {}
    for name in ["small18-rank", "compact19-pointed"]:
        d = read(campaign / "further" / (name + ".json"))
        attempts[name] = {k: v for k,v in d.items() if k not in {"ainvs", "points"}}
    save(output / "further-search.json", attempts)
    print(json.dumps({"certificate_checks": checks, "main_parameters": fast["count"],
        "speedup": summary["speedup_in_this_comparison"], "results": str(output)}, indent=2))


if __name__ == "__main__":
    main()
