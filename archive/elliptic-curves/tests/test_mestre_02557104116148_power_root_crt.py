#!/usr/bin/env python3
"""Pinned checkpoint tests for the completed projective p-adic Mestre lane."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(CAS))

import search_mestre_02557104116148_direct_rational as direct  # noqa: E402
from search_mestre_root_tuple_scale import point_digest, point_on_short_curve  # noqa: E402
from search_mestre_root_tuple_scale_max200 import (  # noqa: E402
    mod3_independence_certificate,
)


Q = Fraction
ARTIFACT = (
    ROOT
    / "artifacts/generated-results/elliptic_mestre_02557104116148_power_root_crt.json"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class MestrePowerRootCheckpointTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = json.loads(ARTIFACT.read_text())

    def test_complete_projective_population(self) -> None:
        self.assertEqual(
            sha256(ARTIFACT),
            "f21b5be50a8ca46f4fc8771ada88d260d3e58bab99a739a8c810da1d39e48699",
        )
        self.assertEqual(
            self.data["result_sha256"],
            "979d7dccb228863b88d984528bb5f4c3cae02bb7ad7dea583a589b803339289e",
        )
        audit = self.data["lattice_population_audit"]
        self.assertEqual(audit["complete_projective_branch_combinations"], 14_336)
        self.assertEqual(audit["bounded_vectors_visited"], 1_146_880)
        self.assertEqual(audit["inside_completed_direct_box"], 1_424)
        self.assertEqual(audit["inside_completed_farey_annulus"], 284)
        self.assertEqual(audit["above_height_cap"], 7_916)
        self.assertEqual(audit["raw_retained_vector_instances"], 192_752)
        self.assertEqual(audit["unique_reduced_positive_parameters"], 48_188)
        self.assertEqual(
            audit["candidate_pair_lines_sha256"],
            "c37418dba25aa32603796cc158dc01ae1f427943cbc3b8d6414921c1b97f7f70",
        )
        profiles = self.data["complete_projective_p_adic_profiles"]
        self.assertEqual(
            [(row["prime"], row["complete_maximal_ball_count"]) for row in profiles],
            [(2, 2), (3, 2), (5, 8), (7, 4), (11, 4), (13, 14), (19, 2)],
        )
        p3 = next(row for row in profiles if row["prime"] == 3)
        self.assertEqual(
            p3["infinity"]["level_root_counts_after_u_divisible_by_p_filter"],
            [1, 3, 0, 0, 0, 0, 0, 0],
        )

    def test_exact_rank16_leader_certificate(self) -> None:
        record = next(
            row for row in self.data["selected_records"]
            if row["parameter"] == "32047/460"
        )
        conductor = record["conductor_phase"]
        self.assertEqual(
            conductor["conductor"],
            "8338729221410559102486702655689159810354695267392082947939199856143910",
        )
        self.assertEqual(conductor["root_number"], 1)
        self.assertTrue(conductor["below_strict_log_conductor_target_numerically"])
        stage = record["point_stages"]["H1000000"]
        self.assertEqual(stage["stable_numerical_rank"], 16)
        certificate = stage["finite_reduction_attempt"]
        self.assertEqual(certificate["certified_algebraic_rank_lower_bound"], 16)
        self.assertEqual(certificate["combined_exact_rank_over_F3"], 16)
        self.assertEqual(
            certificate["certificate_primes"],
            [71, 73, 101, 127, 131, 149, 257, 263, 277, 281, 347, 349, 353, 359, 367],
        )
        points = tuple(
            (Q(point["x"]), Q(point["y"])) for point in stage["numerical_subset"]
        )
        coefficients = direct.family_coefficients(Q(32_047, 460))
        self.assertEqual(point_digest(points), certificate["point_sha256"])
        self.assertTrue(all(point_on_short_curve(coefficients, point) for point in points))
        replay = mod3_independence_certificate(
            coefficients, points, prime_bound=certificate["certificate_prime_bound"]
        )
        self.assertEqual(json.loads(json.dumps(replay)), certificate)


if __name__ == "__main__":
    unittest.main()
