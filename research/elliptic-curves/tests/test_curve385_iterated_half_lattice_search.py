#!/usr/bin/env python3
"""Structural regression for the frozen curve-385 iterative recovery."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / "artifacts/generated-results/elliptic-curves"
BLIND = ART / "curve385_iterated_half_lattice_blind_v1.json"
VERIFIED = ART / "curve385_iterated_half_lattice_verification_v1.json"

BLIND_SHA256 = "356001898f738f607d984e081663a015825e11de0c606d35055af156eb2d7502"
VERIFIED_SHA256 = "b281556f5d08250f67b69b2c62a640ac17ba4d03325e4402e85c7d60882c3ae5"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def multiply(left, right):
    return [
        [
            sum(left[row][inner] * right[inner][column] for inner in range(len(right)))
            for column in range(len(right[0]))
        ]
        for row in range(len(left))
    ]


class Curve385IteratedHalfLatticeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.blind = json.loads(BLIND.read_text())
        cls.verified = json.loads(VERIFIED.read_text())

    def test_frozen_bytes(self) -> None:
        self.assertEqual(digest(BLIND), BLIND_SHA256)
        self.assertEqual(digest(VERIFIED), VERIFIED_SHA256)
        self.assertEqual(
            self.verified["phase_boundary"][
                "blind_artifact_sha256_before_public_fixture_import"
            ],
            BLIND_SHA256,
        )

    def test_blind_transition_and_limits(self) -> None:
        self.assertEqual(self.blind["status"], "STOPPED_AT_DECLARED_LIFT_LIMIT")
        self.assertFalse(self.blind["blindness_boundary"]["public_rank29_fixture_loaded"])
        initial = self.blind["initial_transition"]
        self.assertEqual((initial["rank_before"], initial["rank_after"]), (17, 20))
        self.assertEqual(initial["selected_new_direction_count"], 3)
        self.assertEqual(
            initial["discovered_group_saturation"]["status"],
            "PASS_BASIS_EQUALS_DISCOVERED_GROUP",
        )
        self.assertEqual(len(self.blind["iterations"]), 1)
        iteration = self.blind["iterations"][0]
        self.assertEqual((iteration["basis_rank_before"], iteration["basis_rank_after"]), (20, 29))
        self.assertEqual(iteration["ranking"]["all_lifts_including_zero_word"], 344)
        self.assertEqual(iteration["ranking"]["nonzero_quotient_word_lifts"], 301)
        self.assertEqual(iteration["searched_new_chart_count"], 301)
        self.assertEqual(iteration["bounded_complete_count"], 301)
        self.assertEqual(iteration["timeout_count"], 0)
        self.assertEqual(iteration["pari_failure_count"], 0)
        self.assertEqual(iteration["new_independent_direction_count"], 9)
        self.assertEqual(iteration["finite_index_saturation_event_count"], 0)
        self.assertEqual(
            iteration["discovered_group_saturation"]["status"],
            "PASS_BASIS_EQUALS_DISCOVERED_GROUP",
        )
        self.assertEqual(self.blind["stop"]["next_nonzero_lift_count"], 176_085)

    def test_public_subgroup_equality(self) -> None:
        self.assertEqual(
            self.verified["status"], "PASS_BLIND_M29_EQUALS_DISPLAYED_PUBLIC_M29"
        )
        transition = self.verified["transition"]
        self.assertEqual(transition["initial_deep43_blind_gain"], 3)
        self.assertEqual(transition["first_lift_round_blind_gain"], 9)
        self.assertEqual(transition["blind_quotient_dimension_recovered"], 12)
        self.assertEqual(transition["unrecovered_displayed_public_quotient_dimension"], 0)
        equality = self.verified["subgroup_equality"]
        left = equality["blind_basis_in_public_basis_rows"]
        right = equality["public_basis_in_blind_basis_rows"]
        identity = [[int(row == column) for column in range(29)] for row in range(29)]
        self.assertEqual(multiply(left, right), identity)
        self.assertEqual(multiply(right, left), identity)
        self.assertEqual(abs(equality["blind_to_public_determinant"]), 1)
        self.assertEqual(abs(equality["public_to_blind_determinant"]), 1)


if __name__ == "__main__":
    unittest.main()
