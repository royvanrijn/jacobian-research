#!/usr/bin/env python3
"""Pinned tests for the H=200000 genuine simultaneous-square search."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
PROGRAM = ROOT / "elliptic-curves"
sys.path.insert(0, str(CAS))
sys.path.insert(0, str(PROGRAM))
SCRIPT = CAS / "search_fermigier_exceptional_pair_simultaneous_h200000.py"
ARTIFACT = ROOT / (
    "artifacts/generated-results/elliptic-curves/"
    "elliptic_fermigier_exceptional_pair_simultaneous_h200000.json"
)
EXPECTED_SCRIPT_SHA256 = (
    "563583f301f75106c64b71de6f57c9963d19cf79826041a8137dcd52589049fd"
)
EXPECTED_ARTIFACT_SHA256 = (
    "66f23d9ab4a931ef9e3123c021f56a184505f26c1fad6b647301e7f0d1d52fbe"
)

SPEC = importlib.util.spec_from_file_location("fermigier_pair_h200000", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class FermigierPairSimultaneousH200000Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = json.loads(ARTIFACT.read_text())

    def test_pinned_files_and_stable_digest(self) -> None:
        self.assertEqual(sha256(SCRIPT), EXPECTED_SCRIPT_SHA256)
        self.assertEqual(sha256(ARTIFACT), EXPECTED_ARTIFACT_SHA256)
        self.assertEqual(
            self.data["result_sha256"],
            "6076c5c3d5ac3db63e341ec5ca7533e50a8e995561fffa3261c8978d7a8bbe24",
        )
        self.assertEqual(
            self.data["result_sha256"], MODULE.stable_result_digest(self.data)
        )
        self.assertEqual(self.data["source"]["script_sha256"], EXPECTED_SCRIPT_SHA256)

    def test_complete_box_and_separate_factor_semantics(self) -> None:
        outcome = self.data["outcome"]
        self.assertTrue(outcome["all_80_direction_searches_completed"])
        self.assertEqual(self.data["search_box"]["projective_height_bound"], 200_000)
        self.assertTrue(self.data["search_box"]["one_pass"])
        self.assertEqual(self.data["search_box"]["retries"], 0)
        self.assertEqual(outcome["direction_count"], 80)
        self.assertEqual(outcome["fiber_product_pair_count"], 3_160)
        self.assertFalse(outcome["product_square_surrogate_used"])
        self.assertTrue(all(row["search"]["status"] == "completed" for row in self.data["directions"]))
        self.assertTrue(all(row["both_factors_checked_separately"] for row in self.data["pair_results"]))

    def test_projective_local_bitsets_replay(self) -> None:
        self.assertEqual(self.data["local_sieve"]["primes"], [13, 23, 37, 41, 43])
        self.assertTrue(self.data["local_sieve"]["built_before_any_bounded_search"])
        transport, _ = MODULE.load_transport(ROOT)
        directions = MODULE.build_directions(transport)
        _, sieve = MODULE.build_local_sieve(directions)
        self.assertEqual(sieve["pair_count"], 3_160)
        self.assertEqual(
            sieve["pair_manifest_sha256"],
            self.data["local_sieve"]["pair_manifest_sha256"],
        )
        self.assertEqual(
            sieve["direction_bitsets_hex"],
            self.data["local_sieve"]["direction_bitsets_hex"],
        )

    def test_every_pair_has_only_the_two_anchor_squares(self) -> None:
        anchors = ["39508/39", "28917/10"]
        self.assertEqual(len(self.data["pair_results"]), 3_160)
        for row in self.data["pair_results"]:
            self.assertEqual(row["simultaneous_square_parameters"], anchors)
            self.assertEqual(row["beyond_anchor_parameters"], [])
        outcome = self.data["outcome"]
        self.assertEqual(outcome["signed_simultaneous_parameter_count_including_anchors"], 2)
        self.assertEqual(outcome["signed_third_parameter_count"], 0)
        self.assertEqual(outcome["new_third_parameter_count"], 0)
        self.assertEqual(self.data["third_parameter_certifications"], [])
        self.assertFalse(outcome["target_met"])

    def test_two_single_cover_points_are_exact_and_rejected(self) -> None:
        expected = {
            ("P14__R20E2", "75804/2587"),
            ("P17__R20E4", "130924/9385"),
        }
        rows = self.data["individual_beyond_anchor_incidences"]
        self.assertEqual(
            {(row["direction_id"], row["canonical_parameter_T"]) for row in rows},
            expected,
        )
        self.assertTrue(all(not row["prior_parameter"] for row in rows))
        self.assertTrue(
            all(
                row["rejection_reason"]
                == "no second affine cover is square at this parameter"
                for row in rows
            )
        )
        directions = {
            row["direction_id"]: tuple(Fraction(value) for value in row["polynomial_low_to_high"])
            for row in self.data["directions"]
        }
        for row in rows:
            parameter = Fraction(row["signed_parameter_T"])
            ordinate = Fraction(row["ordinate"])
            self.assertEqual(
                MODULE.evaluate_polynomial(directions[row["direction_id"]], parameter),
                ordinate * ordinate,
            )

    def test_prior_snapshot_is_pinned(self) -> None:
        snapshot = self.data["prior_parameter_snapshot"]
        self.assertEqual(snapshot["source_count"], 40)
        self.assertEqual(snapshot["unique_prior_parameter_count"], 1_808)
        self.assertEqual(
            snapshot["prior_parameter_sha256"],
            "1774a4c7957a66da428de5a0e2b3327574f32afe244216082734e48a0e6e23a8",
        )


if __name__ == "__main__":
    unittest.main()
