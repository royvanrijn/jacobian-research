from __future__ import annotations

import importlib.util
from pathlib import Path
import sys


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "cas"
    / "run_fermigier_rank20_mwrank_descent.py"
)
spec = importlib.util.spec_from_file_location("r20mwrank", SCRIPT)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
assert spec.loader is not None
spec.loader.exec_module(module)


def test_exact_rank20_from_external_lower_bound():
    lower, classification = module.classify_bounds(0, 20)
    assert lower == 20
    assert classification == "M0_exact_rank20"


def test_residual_interval():
    lower, classification = module.classify_bounds(0, 22)
    assert lower == 20
    assert classification == "M2_residual_rank_interval"


def test_mwrank_finds_new_lower_direction():
    lower, classification = module.classify_bounds(21, 22)
    assert lower == 21
    assert classification == "M3_rank_at_least21"


def test_upper_bound_below_known_rank_is_a_contradiction():
    try:
        module.classify_bounds(0, 19)
    except ArithmeticError:
        pass
    else:
        raise AssertionError("expected contradiction")


class FakeMwrankCurve:
    def selmer_rank(self):
        return 22

    def rank(self):
        return 0

    def rank_bound(self):
        return 22

    def certain(self):
        return False


def test_collect_result_keeps_raw_and_combined_intervals_separate():
    result = module.collect_result(FakeMwrankCurve(), two_torsion_rank=0, elapsed_seconds=1.5)
    assert result.selmer_rank == 22
    assert result.mwrank_lower == 0
    assert result.mwrank_upper == 22
    assert result.effective_lower == 20
    assert result.classification == "M2_residual_rank_interval"
