"""Resume a bounded batch of point searches; every result has a C# certificate."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import subprocess
import sys

from point_search import ROOT, save


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--mode", choices=["pointed", "covers", "rank", "direct"], default="pointed")
    parser.add_argument("--height", type=int, default=100000)
    parser.add_argument("--timeout", type=float, default=30)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--limit", type=int, default=48)
    parser.add_argument("--minimal", action="store_true")
    args = parser.parse_args()
    if not (1 <= args.workers <= 4 and 1 <= args.limit <= 10000):
        parser.error("Invalid worker/count limit")
    directory = Path(args.input_dir).resolve()
    output = Path(args.output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    ranking = directory / "candidates.json"
    if ranking.exists():
        paths = [directory / f"curve_{r['U']}_{r['V']}.json"
                 for r in json.loads(ranking.read_text()) if not r["Control"]][:args.limit]
    else:
        paths = sorted(directory.glob("curve_*.json"))[:args.limit]
    manifest = {"mode": args.mode, "height": args.height, "timeout": args.timeout,
                "inputs": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    if args.minimal:
        manifest["minimal"] = True
    manifest_path = output / "batch.json"
    if manifest_path.exists() and json.loads(manifest_path.read_text()) != manifest:
        raise ValueError("Changed inputs or search settings; choose a different output directory")
    save(manifest_path, manifest)

    def job(path):
        prefix = output / path.stem
        target = prefix.with_suffix(".json")
        if target.exists():
            result = json.loads(target.read_text())
            if "exact_certificate" in result:
                return path.stem, result, True
        command = [sys.executable, str(ROOT / "tools/RankHunt/point_search.py"),
            "--input", str(path), "--output", str(prefix), "--mode", args.mode,
            "--height", str(args.height), "--anchors", "64", "--timeout", str(args.timeout)]
        if args.minimal:
            command.append("--minimal")
        process = subprocess.run(command, text=True, capture_output=True)
        if process.returncode:
            return path.stem, {"error": process.stderr[-3000:]}, False
        return path.stem, json.loads(target.read_text()), False

    rows = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for future in as_completed([pool.submit(job, path) for path in paths]):
            name, result, cached = future.result()
            row = {"curve": name, "mode": args.mode, "lower_bound": result.get("exact_certificate", {}).get("LowerBound"),
                   "new_points": result.get("new_distinct_up_to_sign"), "seconds": result.get("seconds"),
                   "timed_out": result.get("timed_out"), "finished": result.get("finished"), "error": result.get("error")}
            rows.append(row)
            rows.sort(key=lambda r: (-(r["lower_bound"] or 0), r["curve"]))
            save(output / "leaderboard.json", rows)
            print(f"{len(rows)}/{len(paths)} {name}: lower={row['lower_bound']}, new={row['new_points']}, "
                  f"timeout={row['timed_out']}, finished={row['finished']}, cached={cached}", flush=True)


if __name__ == "__main__":
    main()
