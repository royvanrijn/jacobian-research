#!/usr/bin/env python3
"""Focused checks for the rank-seven rational two-section component."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
PROBE = CAS / "probe_mestre_transverse_two_section.py"
VERIFY = CAS / "verify_mestre_transverse_two_section_component.py"
AUDIT = CAS / "audit_mestre_transverse_two_section_specialization.py"


def load(name: str, path: Path):
    sys.path.insert(0, str(CAS))
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


class MestreTransverseTwoSectionTest(unittest.TestCase):
    def test_rank_seven_hensel_branches(self) -> None:
        probe = load("mestre_transverse_probe", PROBE)
        result = probe.run((17, 19, 23), 4)
        self.assertEqual(result["exact_jacobian_rank"], 7)
        self.assertEqual(result["transverse_free_coordinate"], "c4")
        self.assertTrue(
            all(row["jacobian_rank"] == 7 for row in result["finite_field_tangent_checks"])
        )
        self.assertTrue(
            all(row["all_seven_residuals_zero_mod_prime_power"]
                for row in result["hensel_continuations"])
        )

    def test_exact_component_identity(self) -> None:
        verifier = load("mestre_transverse_verify", VERIFY)
        result = verifier.replay()
        self.assertEqual(result["admissible_exact_sample_count"], 301)
        self.assertTrue(result["all_recursive_residuals_vanish"])
        self.assertEqual(result["leading_invariant"], "D=16*(z-36)^2/9")
        self.assertTrue(
            result["split_six_root_parameterization"]["all_root_product_coefficients_match"]
        )
        self.assertEqual(
            result["seed_finite_intersection"]["common_affine_quartic_point"],
            ["33/5", "936/25"],
        )

    def test_finite_reduction_non_promotion_audit(self) -> None:
        audit = load("mestre_transverse_audit", AUDIT)
        result = audit.replay()
        self.assertEqual(result["visible_certificate"]["rank"], 9)
        self.assertEqual(result["augmented_certificate"]["rank"], 9)


if __name__ == "__main__":
    unittest.main()
