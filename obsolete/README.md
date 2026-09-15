# Historical research archive

The maintained standalone search algorithm is now in
[elliptic-rank-search](../elliptic-rank-search/README.md).

This archive preserves the previous working files, including uncommitted
improvements to `unified.py` and `test_unified.py` present at extraction time.
No archived result was discarded or rewritten to remove its historical context.

## Contents

- `tools/EquationSearch/`: the previous combined-workspace implementation, experiments, and tests.
- `tools/RankHunt/`: earlier search engines, construction experiments, and their original library-dependent checks.
- `tools/RankPackageAudit/` and `tools/RankStructureAudit/`: historical audit harnesses.
- `results/`: prior run packages, proof files, summaries, and submission material.
- `docs/`: research reports describing those earlier experiments.
- `moved-files.sha256.json`: original paths, new paths, and SHA256 hashes for all 783 tracked files moved here.

Every tracked file was hashed before and after the move; the hashes match.
Generated or ignored files inside the moved directories were also retained.

The old documents and programs retain their original relative paths and runtime
assumptions. They are historical records, not the supported way to run the
extracted algorithm. Use the new directory for current instructions and tests.
The main library, its applications, and its tests remain at the repository root.
