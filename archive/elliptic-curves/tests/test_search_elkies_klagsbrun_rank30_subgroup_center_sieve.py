#!/usr/bin/env python3

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from math import gcd
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
GENERATED = ROOT / "artifacts" / "generated-results"
sys.path.insert(0, str(CAS))

from search_elkies_klagsbrun_rank30 import exact_linear_combination  # noqa: E402
from search_elkies_klagsbrun_rank30_companion_center_sieve import (  # noqa: E402
    allowed_offset_mask,
    exact_nonoverlap_proof,
    load_companion_centers,
)
from search_elkies_klagsbrun_rank30_denominator_sieve import (  # noqa: E402
    build_offset_residue_masks,
    homogeneous_square_value,
)
from search_elkies_klagsbrun_rank30_subgroup_center_sieve import (  # noqa: E402
    exact_prior_center_separation,
    generate_subgroup_centers,
)


Q = Fraction
SCRIPT = CAS / "search_elkies_klagsbrun_rank30_subgroup_center_sieve.py"
ARTIFACT = (
    GENERATED / "elliptic_elkies_klagsbrun_rank30_subgroup_center_sieve.json"
)
EXPECTED_FULL_SHA256 = (
    "079e43b25d8e61df848b02f9ba336a98eb66f68a8859262c0660110ae5ef4e0c"
)
EXPECTED_SELECTED_SHA256 = (
    "584c037eda54b870bcbc7a50e42a01b0f08787aab6f4792b249cf6db6bb6441a"
)
EXPECTED_PRIMITIVE_COUNT = 118_154_598_712
EXPECTED_SURVIVOR_COUNT = 555_231
EXPECTED_SURVIVOR_SHA256 = (
    "467e408ad6fbf47d00bfe70375d9f6d6aa66e373ef91a0ef86bfb7e965032700"
)


class Rank30SubgroupCenterSieveTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.population, cls.selected, cls.selection = generate_subgroup_centers()

    def test_all_841_centers_and_blind_selection_replay_exactly(self) -> None:
        self.assertEqual(len(self.population), 841)
        self.assertEqual(len(self.selected), 128)
        self.assertEqual(
            self.selection["full_population_sha256"], EXPECTED_FULL_SHA256
        )
        self.assertEqual(
            self.selection["selected_sha256"], EXPECTED_SELECTED_SHA256
        )
        self.assertEqual(
            self.selection["prior_public_and_companion_x_overlap_count"], 0
        )
        self.assertEqual(
            [center.bit_height for center in self.selected],
            sorted(center.bit_height for center in self.selected),
        )
        for center in self.selected:
            self.assertEqual(exact_linear_combination(center.relation), center.point)

    def test_exact_nonoverlap_with_every_prior_region(self) -> None:
        direct = exact_nonoverlap_proof(
            self.selected, denominator_min=50_001, offset_radius=16_384
        )
        self.assertTrue(direct["all_exact_nonoverlap_checks_passed"])
        self.assertEqual(
            direct["prior_deep_x_pair_charts"]["exact_center_chart_pair_count"],
            128 * 406,
        )
        old, _ = load_companion_centers()
        prior = exact_prior_center_separation(
            self.selected,
            old,
            denominator_min=50_001,
            offset_radius=16_384,
        )
        self.assertTrue(prior["all_prior_center_boxes_disjoint"])
        self.assertTrue(prior["prior_public_boxes"]["passed"])
        self.assertTrue(prior["prior_companion_boxes"]["passed"])

    def test_modular_mask_matches_exact_values(self) -> None:
        center = self.selected[31]
        radius = 7
        prime = 19
        denominator = next(
            value
            for value in range(50_001, 50_100)
            if gcd(value, center.denominator_root) == 1
        )
        masks = build_offset_residue_masks(radius, (prime,))
        allowed = allowed_offset_mask(center, denominator, prime, masks)
        for offset in range(-radius, radius + 1):
            numerator = (
                center.numerator * denominator**2
                + offset * center.denominator_root**2
            )
            root = center.denominator_root * denominator
            value = homogeneous_square_value(numerator, root) % prime
            is_residue = any(square * square % prime == value for square in range(prime))
            self.assertEqual(bool(allowed & (1 << (offset + radius))), is_residue)

    def test_generated_artifact_is_complete_and_negative(self) -> None:
        data = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        self.assertEqual(
            data["reproduction"]["script_sha256"],
            hashlib.sha256(SCRIPT.read_bytes()).hexdigest(),
        )
        self.assertEqual(
            data["center_population"]["full_population_sha256"],
            EXPECTED_FULL_SHA256,
        )
        self.assertEqual(
            data["center_population"]["selected_sha256"],
            EXPECTED_SELECTED_SHA256,
        )
        result = data["search_result"]
        self.assertTrue(result["search_complete"])
        self.assertFalse(result["wall_cap_reached"])
        self.assertEqual(result["completed_center_count"], 128)
        self.assertEqual(
            result["processed_primitive_candidate_count"],
            EXPECTED_PRIMITIVE_COUNT,
        )
        self.assertEqual(
            result["declared_primitive_candidate_count"],
            EXPECTED_PRIMITIVE_COUNT,
        )
        self.assertEqual(
            result["modular_survivor_count_after_primitivity"],
            EXPECTED_SURVIVOR_COUNT,
        )
        self.assertEqual(
            result["modular_survivor_manifest_sha256"],
            EXPECTED_SURVIVOR_SHA256,
        )
        self.assertEqual(
            result["exact_nonsquare_count_after_sieve"],
            EXPECTED_SURVIVOR_COUNT,
        )
        self.assertEqual(result["exact_square_abscissa_count"], 0)
        self.assertFalse(result["rank30_target_hit"])


if __name__ == "__main__":
    unittest.main()
