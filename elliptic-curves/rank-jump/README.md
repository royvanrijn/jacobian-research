# Rank-jump understanding

Retrospective analysis, pinned to `7471645`; independent of the active curve search.

- [Findings, paired studies, mechanisms and missing implications](ANALYSIS.md)
- [Per-fibre quotient-block reports](FIBRE_REPORTS.md)
- [Comparison panel (CSV)](../../artifacts/generated-results/elliptic-curves/rank_jump_comparison_panel_v1.csv)
- [Frozen small experiment](EXPERIMENT.json)

Replay from the repository root, using Python 3.12 and its standard library:

```sh
python3 -m unittest discover -s elliptic-curves/rank-jump -p 'test_*.py'
python3 elliptic-curves/rank-jump/retrospective.py check
python3 elliptic-curves/rank-jump/blocks.py check
python3 elliptic-curves/rank-jump/cover_experiment.py check
python3 elliptic-curves/rank-jump/render.py --check
```

The portable `rank_jump_*inputs_v1.json` files retain the exact projections
needed for replay. `capture` and `build` refuse to overwrite outputs; they
are construction commands, not replay requirements. Capture additionally
needs the pinned git history and hash-matching historical transcripts or
bundles. No command runs a point search, parameter sweep or full descent.
`MATH_STATUS.json`, generated `STATUS.md`, search policies and live outputs
are unchanged.
