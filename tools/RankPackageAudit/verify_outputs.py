"""Audit the supplied archive and independently reproduced outputs.

Run from the repository root after the commands in docs/rank-package-audit-20260914.md.
The supplied certificate verifier is tested on valid and deliberately corrupted
inputs. The separate C# executable checks the mathematics with the project library.
"""
import copy
import csv
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import sys

ROOT = Path("artifacts/rank-package-audit")
ORIGINAL = ROOT / "original/elliptic_rank_search"
OUTPUT = ROOT / "reproduced"
sys.path.insert(0, str(ORIGINAL.resolve()))
import certificate


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def read_tsv(path):
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def compare_scores(original_name, reproduced_name):
    expected = {(r["u"], r["v"]): r for r in read_tsv(ORIGINAL / original_name)}
    actual = {(r["u"], r["v"]): r for r in read_tsv(OUTPUT / reproduced_name)}
    require(expected.keys() == actual.keys(), f"Candidate sets: {original_name}")
    error = max(abs(float(expected[key]["score_65521"]) - float(row["full_score"]))
                for key, row in actual.items())
    require(error < 1e-10, f"Score mismatch: {original_name}")
    return {"rows_including_control": len(actual), "max_score_error": error}


manifest = read_json(ORIGINAL / "MANIFEST.json")
for name, expected_hash in manifest.items():
    digest = hashlib.sha256((ORIGINAL / name).read_bytes()).hexdigest()
    require(digest == expected_hash, f"Modified archive input: {name}")

record = read_json(ORIGINAL / "record302.json")
cert = read_json(ORIGINAL / "record302_certificate.json")
result = {"archive_manifest_files_verified": len(manifest),
          "record": certificate.verify_certificate(record, cert)}
for label, data_path, cert_path in [
    ("t0_archived", ORIGINAL / "gp_driver_t0.json", ORIGINAL / "gp_driver_t0.certificate.json"),
    ("t0_fresh_pari", OUTPUT / "gp_t0.json", OUTPUT / "gp_t0.certificate.json"),
]:
    result[label] = certificate.verify_certificate(read_json(data_path), read_json(cert_path))

invalid = []
c = copy.deepcopy(cert)
c["rank_lower_bound"] += 1
invalid.append(("inflated_lower_bound", record, c))
c = copy.deepcopy(cert)
bits = c["independent_rows"][0]["bits"]
c["independent_rows"][0]["bits"] = ("1" if bits[0] == "0" else "0") + bits[1:]
invalid.append(("changed_character_bit", record, c))
c = copy.deepcopy(cert)
c["independent_rows"].append(copy.deepcopy(c["independent_rows"][0]))
invalid.append(("duplicate_matrix_row", record, c))
c = copy.deepcopy(cert)
c["torsion_witness"]["prime"] = 4
invalid.append(("composite_witness_prime", record, c))
c = copy.deepcopy(cert)
c["input_sha256"] = "0" * 64
invalid.append(("wrong_input_hash", record, c))
d = copy.deepcopy(record)
d["points"][0][1] = str(Fraction(d["points"][0][1]) + 1)
invalid.append(("point_off_curve", d, cert))
rejections = {}
for label, data, proof in invalid:
    try:
        certificate.verify_certificate(data, proof)
    except ValueError as exc:
        rejections[label] = str(exc)
    else:
        raise AssertionError(f"Invalid certificate accepted: {label}")
result["invalid_inputs_rejected"] = rejections

direct = read_json(OUTPUT / "csharp_t0_equation_only.json")
direct_cert = certificate.build_certificate(direct, max_prime=1009)
result["csharp_equation_only_points"] = certificate.verify_certificate(direct, direct_cert)
require(result["csharp_equation_only_points"]["rank_lower_bound"] >= 12,
        "Independent Python certificate of points found by C#")
(OUTPUT / "csharp_t0_equation_only.certificate.json").write_text(
    json.dumps(direct_cert, indent=2) + "\n", encoding="utf-8")

crt = read_tsv(OUTPUT / "crt_candidates.tsv")
require(crt == read_tsv(ORIGINAL / "crt_candidates.tsv"), "All CRT rows must agree")
constraints_checked = 0
for row in crt:
    u, v = int(row["u"]), int(row["v"])
    require(v > 0 and math.gcd(u, v) == 1 and max(abs(u), v) <= 1000000,
            "CRT primitive vector and height")
    for condition in row["reasons"].removeprefix("projective_CRT:").split(","):
        p, t = condition.split("=")
        p = int(p)
        require(v % p == 0 if t == "inf" else (u-int(t)*v) % p == 0,
                "CRT projective congruence")
        constraints_checked += 1
result["crt"] = {"identical_candidate_rows": len(crt),
                 "congruences_checked_independently": constraints_checked,
                 "record_was_generated": any((r["u"], r["v"]) == ("164518", "924945") for r in crt)}
result["deep_rescore"] = compare_scores("candidates_deep.tsv", "csharp_deep.tsv")
result["crt_rescore"] = compare_scores("crt_rescored.tsv", "csharp_crt_rescored.tsv")
(OUTPUT / "output_audit.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
print(json.dumps(result, indent=2))
