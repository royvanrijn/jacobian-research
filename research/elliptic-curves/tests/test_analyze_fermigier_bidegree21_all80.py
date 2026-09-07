#!/usr/bin/env python3
"""Focused replay tests for the pinned all-80 bidegree-(2,1) audit."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "elliptic-curves"), str(ROOT / "elliptic-curves/cas")]
SCRIPT = ROOT / "elliptic-curves/cas/analyze_fermigier_bidegree21_all80.py"
ARTIFACT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves/elliptic_fermigier_bidegree21_all80.json"
)
EXPECTED_SCRIPT_SHA256 = (
    "d34caeb6b34bbad14c7c7cdf29436c455e25c957080c967d7f17ba0759d28fa2"
)
EXPECTED_ARTIFACT_SHA256 = (
    "2c3aa7a8fc57ad7160397506e8db47bb07ea8c988bab87c9e51b1529000301f5"
)

SPEC = importlib.util.spec_from_file_location("fermigier_bidegree21_all80", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class FermigierBidegree21All80Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = json.loads(ARTIFACT.read_text())

    def test_pinned_files_and_result_digest(self) -> None:
        self.assertEqual(sha256_file(SCRIPT), EXPECTED_SCRIPT_SHA256)
        self.assertEqual(sha256_file(ARTIFACT), EXPECTED_ARTIFACT_SHA256)
        stable = dict(self.data)
        digest = stable.pop("result_sha256")
        stable.pop("generated_at_utc")
        self.assertEqual(
            digest,
            hashlib.sha256(
                json.dumps(stable, sort_keys=True, separators=(",", ":")).encode()
            ).hexdigest(),
        )

    def test_population_and_exact_component_histogram(self) -> None:
        self.assertEqual(self.data["population"]["pair_count"], 80)
        classification = self.data["rational_component_classification"]
        self.assertEqual(classification["completed_pair_count"], 80)
        histogram = {item["type"]: item for item in classification["histogram"]}
        self.assertEqual(set(histogram), {
            "degree_drop_plus",
            "degree_drop_minus",
            "numerator_denominator_cancellation",
            "E22_anchor_pole",
            "rank20_anchor_pole",
        })
        for label in ("degree_drop_plus", "degree_drop_minus"):
            self.assertEqual(histogram[label]["count"], 80)
            self.assertEqual(histogram[label]["factor_signature"], [[8, 1]])
            self.assertEqual(histogram[label]["squareclass_genus"], 3)
        for label in (
            "numerator_denominator_cancellation",
            "E22_anchor_pole",
            "rank20_anchor_pole",
        ):
            self.assertEqual(histogram[label]["count"], 80)
            self.assertEqual(histogram[label]["factor_signature"], [[1, 4], [6, 1]])
            self.assertEqual(histogram[label]["squareclass_genus"], 2)

    def test_every_residual_has_a_degree_32_irreducibility_witness(self) -> None:
        residual = self.data["residual_irreducibility"]
        self.assertEqual(residual["completed_pair_count"], 80)
        self.assertEqual(residual["irreducible_pair_count"], 80)
        self.assertEqual(residual["witness_prime_histogram"], [
            {"count": 79, "prime": 101},
            {"count": 1, "prime": 103},
        ])
        self.assertEqual(residual["qq_fallback_pairs"], [])
        self.assertEqual(residual["unresolved_pairs"], [])
        for record in residual["records"]:
            witness = record["irreducibility_witness"]
            self.assertTrue(record["residual_irreducible_over_QQ"])
            self.assertTrue(witness["irreducible_degree_32_witness"])
            self.assertEqual(witness["residual_total_degree"], 32)
            self.assertEqual(witness["factor_signature"], [[32, 1]])

    def test_endpoint_exact_algebra_replay(self) -> None:
        pair = next(
            item
            for item in self.data["population"]["pairs"]
            if item["left"] == "P22" and item["right"] == "R20E8"
        )
        record = MODULE.exact_component_classification((
            pair["left"],
            Fraction(pair["left_x"]),
            pair["right"],
            Fraction(pair["right_x"]),
        ))
        self.assertEqual(record["valid_genus_at_most_one_components"], [])
        self.assertTrue(
            record["exact_identity_checks"]["cancellation_divisibility_by_cT_plus_1"]
        )

    def test_modular_witness_replay_in_both_characteristics(self) -> None:
        expected = {
            ("P13", "R20E1"): 101,
            ("P15", "R20E8"): 103,
        }
        pairs = {
            (item["left"], item["right"]): item
            for item in self.data["population"]["pairs"]
        }
        for key, prime in expected.items():
            item = pairs[key]
            witness = MODULE.modular_attempt(
                Fraction(item["left_x"]),
                Fraction(item["right_x"]),
                prime,
                timeout=60,
            )
            self.assertTrue(witness["irreducible_degree_32_witness"])
            self.assertEqual(witness["factor_signature"], [[32, 1]])

    def test_negative_outcome_is_explicit(self) -> None:
        self.assertTrue(self.data["scope"]["all_80_independent_pairs_classified"])
        self.assertIsNone(self.data["breakthrough"])
        self.assertEqual(self.data["outcome"]["valid_genus_at_most_one_components"], 0)
        self.assertFalse(self.data["outcome"]["target_met"])


if __name__ == "__main__":
    unittest.main()
