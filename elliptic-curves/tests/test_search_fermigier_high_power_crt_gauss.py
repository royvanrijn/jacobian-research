#!/usr/bin/env python3
"""Focused replay checks for the canonical-u high-power CRT lane."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "elliptic-curves/cas/search_fermigier_high_power_crt_gauss.py"
ARTIFACT = ROOT / "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_high_power_crt_gauss.json"
SPEC = importlib.util.spec_from_file_location("fermigier_high_power", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class HighPowerCrtGaussTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = json.loads(ARTIFACT.read_text())

    def test_exact_population_constants(self) -> None:
        self.assertEqual(self.data["population"]["raw_sign_quotiented_count"], 31_524)
        self.assertEqual(self.data["population"]["fresh_count"], 31_520)
        self.assertEqual(
            self.data["population"]["fresh_parameter_sha256"],
            MODULE.EXPECTED_FRESH_DIGEST,
        )
        self.assertEqual(self.data["population"]["broad_neighbourhood_intersection_count"], 0)

    def test_prior_seed_is_calibration_only(self) -> None:
        self.assertIn("673709/29965", self.data["prior_exclusion"]["intersection"])
        self.assertEqual(
            self.data["seed_calibration"]["point_search"],
            "not run unless the exact conductor gate completes below target",
        )

    def test_exact_local_constraints(self) -> None:
        for row in self.data["selected"]:
            homogeneous = int(row["radical"]["homogenized_discriminant_factor"])
            for prime, minimum in zip(MODULE.PRIMES, (2, 2, 2), strict=True):
                self.assertEqual(homogeneous % prime**minimum, 0)
            self.assertEqual(row["point_search"], "not run; exact completed conductor below target required")

    def test_gate_is_conductor_first(self) -> None:
        self.assertEqual(self.data["outcome"]["point_search_calls"], 0)
        self.assertEqual(self.data["outcome"]["rank_calls"], 0)
        self.assertEqual(self.data["outcome"]["strict_target_feasible_fresh_fibres"], 0)
        for row in self.data["selected"]:
            conductor = row["conductor"]
            if conductor["status"] == "completed":
                self.assertFalse(conductor["below_strict_log_conductor_target"])

    def test_result_digest(self) -> None:
        self.assertEqual(self.data["result_sha256"], MODULE.result_digest(self.data))


if __name__ == "__main__":
    unittest.main()
