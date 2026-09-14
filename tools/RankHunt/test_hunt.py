"""Small integration checks for exact maps, halving and checkpoint recovery.

Run after building RankHunt. Uses the local GP executable and temporary files
under artifacts; it never starts the large search campaign.
"""
import json
from fractions import Fraction
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from point_search import ROOT, make_script

DLL = ROOT / "tools/RankHunt/bin/Release/net8.0/RankHunt.dll"


class HuntTests(unittest.TestCase):
    def setUp(self):
        parent = ROOT / "artifacts/rank-hunt-tests"
        parent.mkdir(parents=True, exist_ok=True)
        self.temp = tempfile.TemporaryDirectory(dir=parent)
        self.directory = Path(self.temp.name)

    def tearDown(self):
        if not self.directory.resolve().is_relative_to((ROOT / "artifacts/rank-hunt-tests").resolve()):
            raise RuntimeError("Temporary directory escaped the test workspace")
        self.temp.cleanup()

    def cli(self, *args):
        process = subprocess.run(["dotnet", str(DLL), *map(str, args)], text=True, capture_output=True, timeout=40)
        self.assertEqual(process.returncode, 0, process.stderr)
        return process.stdout

    def source(self, a, points):
        path = self.directory / "input.json"
        path.write_text(json.dumps({"ainvs": list(map(str, a)),
                                   "points": [list(map(str, p)) for p in points]}))
        return path

    def search(self, source, mode, *extra):
        prefix = self.directory / mode
        process = subprocess.run([sys.executable, str(ROOT / "tools/RankHunt/point_search.py"),
            "--input", str(source), "--output", str(prefix), "--mode", mode, "--height", "100",
            "--timeout", "15", *extra], capture_output=True, text=True, timeout=30)
        self.assertEqual(process.returncode, 0, process.stderr)
        result = json.loads(prefix.with_suffix(".json").read_text())
        self.assertTrue(result["finished"], prefix.with_suffix(".stderr.txt").read_text())
        a1, a2, a3, a4, a6 = map(Fraction, result["ainvs"])
        for p in result["points"]:
            x, y = map(Fraction, p)
            self.assertEqual(y*y+a1*x*y+a3*y, x**3+a2*x*x+a4*x+a6)
        return result

    def test_pointed_map_with_rational_anchor(self):
        source = self.source([0, 0, 0, -25, 4], [["-109/25", "686/125"]])
        result = self.search(source, "pointed", "--minimal")
        self.assertGreater(result["new_distinct_up_to_sign"], 0)

    def test_actual_two_covers_map_back(self):
        source = self.source([0, 0, 0, -25, 4], [[0, 2]])
        result = self.search(source, "covers", "--minimal")
        self.assertGreaterEqual(result["exact_certificate"]["LowerBound"], 2)

    def test_even_point_gets_an_exact_half(self):
        source = self.source([0, 0, 1, -1, 0], [[1, 0]])  # twice (0,0)
        relations = self.directory / "relations.json"
        self.cli("relations", "--input", source, "--output", relations)
        self.assertEqual(json.loads(relations.read_text())["relations"], [[0]])
        result = self.search(relations, "halve", "--minimal")
        self.assertIn(["0", "-1"], result["points"])  # sign-normalized (0,0)
        self.assertEqual(result["exact_certificate"]["ImageDimension"], 1)

    def test_equation_strings_cannot_be_gp_commands(self):
        with self.assertRaises(ValueError):
            make_script({"ainvs": ["0", "0", "0", "1;quit", "1"], "points": []}, "pointed", 10, 1, 0)

    def test_grid_checkpoint_resumes_same_result(self):
        folder = self.directory / "grid"
        arguments = ("grid", "--height", "4", "--keep", "32", "--refine-keep", "16", "--final-keep", "4",
                     "--prime-bound", "16382", "--workers", "2", "--output", folder)
        self.cli(*arguments)
        before = json.loads((folder / "candidates.json").read_text())
        checkpoint = json.loads((folder / "grid-checkpoint.json").read_text())
        self.assertEqual(checkpoint["PrimitiveCount"], 23)
        # Replay a saved completed grid, not the expensive full campaign.
        self.cli(*arguments)
        self.assertEqual(before, json.loads((folder / "candidates.json").read_text()))
        self.assertTrue(json.loads((folder / "run.json").read_text())["resumed_grid"])

    def test_recovered_sections_match_independent_fixtures(self):
        # The doubled parameter also checks that homogeneous model scaling uses
        # the reduced fraction, not the raw CLI denominator.
        for family, n, d, fixture, lower in [
            ("icarm302-17", -116, 474, "icarm302-structured17.json", 17),
            ("icarm302-18", 0, 3, "icarm302-basechange18.json", 18),
        ]:
            output = self.directory / (family + ".json")
            self.cli("export", "--family", family, "--u", n, "--v", d, "--output", output)
            actual = json.loads(output.read_text())
            expected = json.loads((ROOT / "tests/Fixtures" / fixture).read_text())
            self.assertEqual(actual["ainvs"], expected["ainvs"])
            self.assertEqual(actual["points"], expected["points"])
            self.assertEqual(actual["rank_lower_bound"], lower)

    def test_generic_rank_is_not_assumed_at_specialization(self):
        output = self.directory / "t0.json"
        self.cli("export", "--family", "icarm302-17", "--u", 0, "--v", 1, "--output", output)
        actual = json.loads(output.read_text())
        self.assertEqual(len(actual["points"]), 17)
        # This is a finite-character lower bound, not an upper-bound assertion.
        self.assertEqual(actual["rank_lower_bound"], 12)

    def test_enrichment_preserves_existing_points(self):
        folder = self.directory / "original"
        folder.mkdir()
        path = folder / "curve_-58_237.json"
        self.cli("export", "--family", "icarm302", "--u", -58, "--v", 237, "--output", path)
        before = json.loads(path.read_text())
        (folder / "candidates.json").write_text(json.dumps([{"U": -58, "V": 237, "Control": False}]))
        out = self.directory / "enriched"
        self.cli("enrich", "--input-dir", folder, "--output", out)
        after = json.loads((out / path.name).read_text())
        self.assertEqual(before, json.loads(path.read_text()))
        self.assertTrue(set(map(tuple, before["points"])) <= set(map(tuple, after["points"])))
        self.assertEqual(before["rank_lower_bound"], 3)
        self.assertEqual(after["rank_lower_bound"], 17)


if __name__ == "__main__":
    unittest.main()
