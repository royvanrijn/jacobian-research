from __future__ import annotations

import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "cas" / "run_fermigier_rank20_pari_descent.py"
spec = importlib.util.spec_from_file_location("r20pari", SCRIPT)
module = importlib.util.module_from_spec(spec)
import sys
sys.modules[spec.name] = module
assert spec.loader is not None
spec.loader.exec_module(module)


def test_parse_efforts_deduplicates_and_preserves_order():
    assert module.parse_efforts("0,1,2,1") == (0, 1, 2)


def test_exact_rank20():
    lower, classification = module.classify_bounds(18, 20)
    assert lower == 20
    assert classification == "P0_exact_rank20"


def test_unresolved_residual_interval():
    lower, classification = module.classify_bounds(18, 22)
    assert lower == 20
    assert classification == "P2_residual_rank_interval"


def test_new_rank_direction():
    lower, classification = module.classify_bounds(21, 22)
    assert lower == 21
    assert classification == "P3_rank_at_least21"


def test_contradiction_with_certified_rank():
    try:
        module.classify_bounds(18, 19)
    except ArithmeticError:
        pass
    else:
        raise AssertionError("expected contradiction")
