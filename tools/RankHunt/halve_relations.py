"""Try bounded exact halving of local Kummer relations, preserving every point.

This can improve a certificate without increasing the rational span of the input
points. It is not a computation of full saturation or an upper rank bound.
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys
from point_search import ROOT, save

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input", required=True)
parser.add_argument("--output-dir", required=True)
parser.add_argument("--rounds", type=int, default=4)
parser.add_argument("--timeout", type=float, default=30)
args = parser.parse_args()
if not 1 <= args.rounds <= 8:
    parser.error("Use 1 to 8 rounds")
directory = Path(args.output_dir).resolve()
directory.mkdir(parents=True, exist_ok=True)
current = Path(args.input).resolve()
history = []
for iteration in range(1, args.rounds+1):
    relations = directory / f"relations_{iteration}.json"
    subprocess.run(["dotnet", str(ROOT / "tools/RankHunt/bin/Release/net8.0/RankHunt.dll"),
                    "relations", "--input", str(current), "--output", str(relations)], check=True)
    if not json.loads(relations.read_text())["relations"]:
        break
    prefix = directory / f"halves_{iteration}"
    subprocess.run([sys.executable, str(ROOT / "tools/RankHunt/point_search.py"), "--input", str(relations),
                    "--output", str(prefix), "--mode", "halve", "--minimal", "--timeout", str(args.timeout)], check=True)
    current = prefix.with_suffix(".json")
    result = json.loads(current.read_text())
    history.append({"round": iteration, "lower_bound": result["exact_certificate"]["LowerBound"],
                    "new_points": result["new_distinct_up_to_sign"], "finished": result["finished"]})
    save(directory / "history.json", history)
    if not result["finished"] or not result["new_distinct_up_to_sign"]:
        break
print(f"Retained result: {current}")
