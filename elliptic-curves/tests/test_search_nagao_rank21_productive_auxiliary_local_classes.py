#!/usr/bin/env python3
"""Regression tests for the auxiliary finite-class lift screen."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
if str(CAS) not in sys.path:
    sys.path.insert(0, str(CAS))

from search_nagao_rank21_productive_auxiliary_local_classes import (  # noqa: E402
    EXPECTED_INPUT_SHA256,
    FiniteCurve,
    LOCAL_PRIMES,
    TARGET_RESIDUES,
    load_input,
    parameter_residues,
)


SCRIPT = CAS / "search_nagao_rank21_productive_auxiliary_local_classes.py"
ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "elliptic_nagao_rank21_productive_auxiliary_local_classes.json"
)
INPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "elliptic_nagao_rank21_productive_auxiliary_orbit.json"
)


class ProductiveAuxiliaryLocalClassTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = json.loads(ARTIFACT.read_text(encoding="utf-8"))

    def test_input_and_finite_curve_orders(self) -> None:
        self.assertEqual(hashlib.sha256(INPUT.read_bytes()).hexdigest(), EXPECTED_INPUT_SHA256)
        _, auxiliary, basis = load_input(ROOT)
        self.assertEqual(len(basis), 8)
        curves = tuple(FiniteCurve(auxiliary.weierstrass_coefficients, prime) for prime in LOCAL_PRIMES)
        self.assertEqual(tuple(len(curve.points) for curve in curves), (20, 48, 90))
        expected_intersections = (
            {3, 10},
            {17, 20},
            {12, 31, 40, 43, 52, 71},
        )
        for curve, targets, expected in zip(
            curves, TARGET_RESIDUES, expected_intersections, strict=True
        ):
            residues = parameter_residues(curve, auxiliary)
            image = {value for value in residues if value is not None}
            self.assertEqual(targets & image, expected)

    def test_complete_state_and_exact_lift_counts(self) -> None:
        screen = self.data["finite_screen"]
        self.assertEqual(screen["complete_product_state_count"], 86400)
        self.assertEqual(screen["complete_product_expected_count"], 20 * 48 * 90)
        self.assertEqual(screen["selected_target_state_count"], 192)
        self.assertEqual(screen["state_counts_after_each_basis_coordinate"][-1], 86400)
        lifts = self.data["exact_lifts"]
        self.assertEqual(lifts["distinct_nonsingular_parameter_count"], 182)
        self.assertEqual(lifts["proxy_below_190_count"], 0)
        self.assertAlmostEqual(lifts["minimum_log_radical_upper_proxy"], 1239.138185968238)

    def test_every_stored_lift_has_forced_valuations(self) -> None:
        for record in self.data["exact_lifts"]["top_32"]:
            valuations = dict(record["small_prime_valuations"])
            self.assertGreaterEqual(valuations.get("13", valuations.get(13, 0)), 4)
            self.assertGreaterEqual(valuations.get("37", valuations.get(37, 0)), 2)
            self.assertGreaterEqual(valuations.get("83", valuations.get(83, 0)), 2)

    def test_artifact_hash_and_scope(self) -> None:
        self.assertFalse(self.data["target_hit"])
        self.assertEqual(self.data["conclusion"]["exact_conductor_call_count"], 0)
        self.assertEqual(
            self.data["reproduction"]["script_sha256"],
            hashlib.sha256(SCRIPT.read_bytes()).hexdigest(),
        )


if __name__ == "__main__":
    unittest.main()
