#!/usr/bin/env python3
"""Focused exact quotient replay on the rational two-section component."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
SCRIPT = CAS / "screen_mestre_transverse_component_independence.py"


def load() -> object:
    sys.path.insert(0, str(CAS))
    spec = importlib.util.spec_from_file_location("mestre_component_independence", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class MestreTransverseComponentIndependenceTest(unittest.TestCase):
    def test_small_exact_grid_has_no_spurious_gain(self) -> None:
        screen = load()
        result = screen.screen(
            root_height=2, parameter_height=2, prime_bound=101
        )
        self.assertEqual(result["candidate_pair_count"], 24)
        self.assertEqual(result["admissible_specialization_count"], 24)
        self.assertEqual(result["strict_quotient_rank_gain_count"], 0)
        self.assertEqual(result["best_observed"]["visible_rank"], 9)
        self.assertEqual(result["best_observed"]["augmented_rank"], 9)


if __name__ == "__main__":
    unittest.main()
