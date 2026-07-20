#!/usr/bin/env python3
"""Regression checks for ordinary-to-Laurent Newton support translation."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import sympy as sp

from jcsearch.newton import monomial_expansion, stage_d_families
from jcsearch.translated import TranslatedTwoDivisorChart, x, y

chart = TranslatedTwoDivisorChart()
for a in range(9):
    for b in range(9):
        assert sp.cancel(chart.pullback(monomial_expansion(a, b)) - x**a*y**b) == 0

families = stage_d_families()
assert set(families) == {
    "total_ratio_2_3", "weighted_3_2_ratio_2_3",
    "weighted_2_3_ratio_2_3", "corner_2_7_square_cube",
    "corner_7_2_square_cube",
}

print("PASS exact monomial translation for 0<=a,b<=8")

result_path = ROOT / "results" / "stage_d_newton_bands.json"
if result_path.exists():
    records = json.loads(result_path.read_text())
    exact = [r for r in records if any(run.get("prime") == 0 and run.get("empty")
                                       for run in r["runs"])]
    timed_out = [r for r in records if "timeout" in r["failure_class"]]
    assert len(exact) == 4
    assert [r["family"] for r in timed_out] == ["total_ratio_2_3"]
    assert all(r["identity_keller_control_in_kernels"] for r in records)
    assert all(r["target_bracket_exponent_reachable"] for r in records)
    print("PASS cached Stage D outcomes: four exact emptiness results and one timeout")
else:
    print("SKIP cached Stage D outcome audit (run stage_d_newton_bands.py first)")
