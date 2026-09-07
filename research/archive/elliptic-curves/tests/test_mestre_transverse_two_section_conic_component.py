#!/usr/bin/env python3
"""Focused checks for the second exact conic-rational component."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
SCRIPT = CAS / "verify_mestre_transverse_two_section_conic_component.py"


def load() -> object:
    sys.path.insert(0, str(CAS))
    spec = importlib.util.spec_from_file_location("mestre_conic_component", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class MestreTransverseTwoSectionConicComponentTest(unittest.TestCase):
    def test_exact_component_and_seed_intersection(self) -> None:
        verifier = load()
        result = verifier.replay()
        self.assertEqual(result["admissible_exact_sample_count"], 1922)
        self.assertTrue(result["all_recursive_residuals_vanish"])
        self.assertTrue(result["split_six_root_parameterization"]["all_root_product_coefficients_match"])
        self.assertEqual(
            result["seed_finite_intersection"]["common_affine_quartic_point"],
            ["4932/455", "-107740485691272/438652175"],
        )


if __name__ == "__main__":
    unittest.main()
