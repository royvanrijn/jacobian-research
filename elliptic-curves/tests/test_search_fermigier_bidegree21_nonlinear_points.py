#!/usr/bin/env python3
"""Lightweight checks for the bounded degree-32 rational-point sieve."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "elliptic-curves/cas/search_fermigier_bidegree21_nonlinear_points.py"
ARTIFACT = (
    ROOT
    / "artifacts/generated-results/"
    "elliptic_fermigier_bidegree21_p13_r20e1_nonlinear_points_h1024.json"
)
EXPECTED_SCRIPT_SHA256 = (
    "1dbe3cb7b95991d671f4df64109012afe36ac568556224298e8d8892941b3044"
)
EXPECTED_ARTIFACT_SHA256 = (
    "dd281569a1da8eb1c07a635faecb8b9f27269751c2639fb6b94f3a1bada46310"
)

SPEC = importlib.util.spec_from_file_location("fermigier_nonlinear_points", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class FermigierNonlinearPointSieveTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = json.loads(ARTIFACT.read_text())

    def test_pinned_files_and_result_digest(self) -> None:
        self.assertEqual(sha256_file(SCRIPT), EXPECTED_SCRIPT_SHA256)
        self.assertEqual(sha256_file(ARTIFACT), EXPECTED_ARTIFACT_SHA256)
        self.assertEqual(self.data["result_sha256"], MODULE.result_digest(self.data))

    def test_declared_box_is_complete_and_exact(self) -> None:
        search = self.data["search_region"]
        self.assertEqual(search["height"], 1024)
        self.assertEqual(search["sieve_primes"], [17, 19, 23, 29, 31])
        self.assertGreater(search["crt_modulus"], 2 * search["height"])
        self.assertEqual(search["C_D_pairs_scanned"], (2 * 1024 + 1) * 1024)
        self.assertEqual(search["primitive_exact_candidates"], 21819)
        self.assertEqual(search["exact_homogeneous_evaluations"], 21819)
        self.assertEqual(search["exact_hits"], 0)
        self.assertTrue(self.data["scope"]["affine_box_complete"])

    def test_factor_boundary_and_known_intersections(self) -> None:
        factor = self.data["nonlinear_factor"]
        self.assertEqual(factor["total_degree"], 32)
        self.assertEqual(factor["term_count"], 561)
        self.assertEqual(
            factor["primitive_coefficient_sha256"],
            "5c60ea4247ddc7eb99f1cc6726c592e569be3ee7c1de0b74892e0bad252d6eda",
        )
        boundary = self.data["projective_boundary"]
        self.assertEqual(boundary["factor_signature_over_QQ_after_C_equals_1"], [[32, 1]])
        self.assertEqual(boundary["rational_projective_point_count"], 0)
        intersections = self.data["known_linear_component_intersections"]
        self.assertEqual(intersections["distinct_rational_intersection_points"], 0)
        self.assertEqual(len(intersections["lines"]), 5)
        self.assertTrue(
            all(item["rational_affine_intersection_count"] == 0 for item in intersections["lines"])
        )

    def test_negative_outcome_is_not_overclaimed(self) -> None:
        outcome = self.data["outcome"]
        self.assertEqual(outcome["affine_degree32_hits_in_box"], 0)
        self.assertEqual(outcome["valid_genus_at_most_one_points_found"], 0)
        self.assertEqual(outcome["new_sections"], 0)
        self.assertEqual(outcome["new_specializations"], 0)
        self.assertFalse(outcome["target_met"])
        self.assertFalse(self.data["scope"]["all_rational_points_on_degree32_component_classified"])

    def test_exact_homogenized_evaluator(self) -> None:
        coefficients = {(32, 0): 2, (1, 1): -5, (0, 0): 3}
        C, K, D = 2, -3, 5
        expected = 2 * C**32 - 5 * C * K * D**30 + 3 * D**32
        self.assertEqual(MODULE.evaluate_homogeneous(coefficients, C, K, D), expected)


if __name__ == "__main__":
    unittest.main()
