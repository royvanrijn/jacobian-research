#!/usr/bin/env python3

from fractions import Fraction
from hashlib import sha256
import json
from math import gcd
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
GENERATED = ROOT / "artifacts" / "generated-results"
sys.path.insert(0, str(CAS))

from elkies_klagsbrun_rank29 import PUBLISHED_POINTS  # noqa: E402
from search_elkies_klagsbrun_rank30 import (  # noqa: E402
    exact_linear_combination,
    point_negate,
)
from search_elkies_klagsbrun_rank30_companion_center_sieve import (  # noqa: E402
    SOURCE_ARTIFACT_HASHES,
    allowed_offset_mask,
    exact_nonoverlap_proof,
    load_companion_centers,
    normalized_abscissa,
)
from search_elkies_klagsbrun_rank30_denominator_sieve import (  # noqa: E402
    build_offset_residue_masks,
    homogeneous_square_value,
    map_square_abscissa,
)


Q = Fraction
SCRIPT = CAS / "search_elkies_klagsbrun_rank30_companion_center_sieve.py"
ARTIFACT = (
    GENERATED / "elliptic_elkies_klagsbrun_rank30_companion_center_sieve.json"
)


class ElkiesKlagsbrunCompanionCenterSieveTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.centers, cls.inventory = load_companion_centers()

    def test_pinned_sources_and_all_32_centers_replay_exactly(self) -> None:
        self.assertEqual(len(self.inventory), len(SOURCE_ARTIFACT_HASHES))
        for record in self.inventory:
            expected = SOURCE_ARTIFACT_HASHES[record["path"]]
            self.assertEqual(record["sha256"], expected)
            self.assertEqual(
                sha256((ROOT / record["path"]).read_bytes()).hexdigest(), expected
            )
        self.assertEqual([item["point_record_count"] for item in self.inventory], [0, 0, 0, 57])
        self.assertEqual(len(self.centers), 32)
        public_x = {point[0] for point in PUBLISHED_POINTS}
        self.assertEqual(len({center.x for center in self.centers}), 32)
        for center in self.centers:
            self.assertNotIn(center.x, public_x)
            self.assertEqual(center.denominator_root**2, center.x.denominator)
            self.assertEqual(exact_linear_combination(center.relation), center.point)

    def test_exact_nonoverlap_gate_covers_every_center_and_x_pair(self) -> None:
        proof = exact_nonoverlap_proof(
            self.centers, denominator_min=3163, offset_radius=16_384
        )
        self.assertTrue(proof["all_exact_nonoverlap_checks_passed"])
        self.assertTrue(proof["public_center_box_separation"]["passed"])
        self.assertTrue(proof["companion_box_pairwise_separation"]["passed"])
        self.assertTrue(proof["prior_deep_x_offset_charts"]["passed"])
        pair_proof = proof["prior_deep_x_pair_charts"]
        self.assertTrue(pair_proof["passed"])
        self.assertEqual(pair_proof["exact_center_chart_pair_count"], 32 * 406)
        self.assertEqual(pair_proof["failure_count"], 0)
        self.assertGreater(
            Q(pair_proof["minimum_clearance_record"]["exact_clearance_ratio"]), 1
        )

    def test_center_seed_and_normalized_perturbation_map_exactly(self) -> None:
        for center in (self.centers[0], self.centers[-1]):
            denominator = next(
                value
                for value in range(3163, 3200)
                if gcd(value, center.denominator_root) == 1
            )
            seed_numerator = center.numerator * denominator**2
            seed_root = center.denominator_root * denominator
            value = homogeneous_square_value(seed_numerator, seed_root)
            expected_root = (2 * center.point[1] + center.point[0]) * seed_root**3
            self.assertEqual(expected_root.denominator, 1)
            self.assertEqual(value, expected_root.numerator**2)
            mapped = map_square_abscissa(
                seed_numerator, seed_root, abs(expected_root.numerator)
            )
            self.assertEqual(set(mapped), {center.point, point_negate(center.point)})

            offset = next(
                value for value in range(1, 50) if gcd(value, denominator) == 1
            )
            numerator, root = normalized_abscissa(center, denominator, offset)
            self.assertEqual(gcd(numerator, root), 1)
            self.assertEqual(
                Q(numerator, root**2) - center.x,
                Q(offset, denominator**2),
            )

    def test_generic_modular_mask_matches_exact_homogeneous_values(self) -> None:
        center = self.centers[7]
        radius = 9
        prime = 17
        denominator = next(
            value
            for value in range(3163, 3200)
            if gcd(value, center.denominator_root) == 1
        )
        residue_masks = build_offset_residue_masks(radius, (prime,))
        mask = allowed_offset_mask(center, denominator, prime, residue_masks)
        for offset in range(-radius, radius + 1):
            numerator = (
                center.numerator * denominator**2
                + offset * center.denominator_root**2
            )
            root = center.denominator_root * denominator
            value = homogeneous_square_value(numerator, root) % prime
            residue = any(square * square % prime == value for square in range(prime))
            self.assertEqual(bool(mask & (1 << (offset + radius))), residue)

    def test_generated_artifact_is_complete_disjoint_and_self_pinned(self) -> None:
        if not ARTIFACT.exists():
            self.skipTest("the companion-center sieve artifact has not been generated")
        data = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        self.assertEqual(
            data["reproduction"]["script_sha256"],
            sha256(SCRIPT.read_bytes()).hexdigest(),
        )
        self.assertEqual(data["center_manifest"]["count"], 32)
        self.assertTrue(data["center_manifest"]["public_x_excluded"])
        self.assertTrue(
            data["exact_nonoverlap_proof"]["all_exact_nonoverlap_checks_passed"]
        )
        result = data["search_result"]
        self.assertTrue(result["one_pass_no_retry"])
        self.assertTrue(result["search_complete"])
        self.assertFalse(result["wall_cap_reached"])
        self.assertEqual(
            result["processed_primitive_candidate_count"],
            result["declared_primitive_candidate_count"],
        )
        self.assertEqual(result["exact_square_abscissa_count"], 0)
        self.assertFalse(result["rank30_target_hit"])
        self.assertEqual(
            data["status"], "bounded_search_no_certified_30th_point"
        )


if __name__ == "__main__":
    unittest.main()
