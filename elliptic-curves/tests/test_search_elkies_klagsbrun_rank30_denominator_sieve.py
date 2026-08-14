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
    point_add,
    point_negate,
)
from search_elkies_klagsbrun_rank30_denominator_sieve import (  # noqa: E402
    PREVIOUS_DEEP_X_OFFSET_HEIGHT,
    allowed_offset_mask,
    anchor_data,
    build_offset_residue_masks,
    homogeneous_square_value,
    map_square_abscissa,
    normalized_abscissa,
    positive_coprime_offset_count,
    small_companion_lookup,
)


Q = Fraction
SCRIPT = CAS / "search_elkies_klagsbrun_rank30_denominator_sieve.py"
ARTIFACT = (
    GENERATED / "elliptic_elkies_klagsbrun_rank30_denominator_sieve.json"
)


class ElkiesKlagsbrunDenominatorSieveTests(unittest.TestCase):
    def test_all_public_x_denominators_are_squares(self) -> None:
        data = [anchor_data(index) for index in range(len(PUBLISHED_POINTS))]
        self.assertEqual(len(data), 29)
        self.assertEqual(sum(root == 1 for _, root in data), 28)
        self.assertEqual(sum(root == 9 for _, root in data), 1)
        for index, (numerator, root) in enumerate(data):
            self.assertEqual(PUBLISHED_POINTS[index][0], Q(numerator, root * root))

    def test_normalized_box_has_exact_new_square_denominator(self) -> None:
        denominator = 3163
        offset = 17
        numerator, root = normalized_abscissa(0, denominator, offset)
        self.assertEqual(gcd(numerator, root), 1)
        x_value = Q(numerator, root * root)
        self.assertEqual(x_value - PUBLISHED_POINTS[0][0], Q(offset, denominator**2))
        self.assertEqual(
            (x_value - PUBLISHED_POINTS[0][0]).denominator,
            denominator**2,
        )
        self.assertGreater(denominator**2, PREVIOUS_DEEP_X_OFFSET_HEIGHT)

    def test_homogeneous_identity_maps_the_public_seed_exactly(self) -> None:
        # k=0 is deliberately excluded from the search, but it calibrates the
        # homogeneous equation and inverse map without PARI or floating point.
        for index, denominator in ((0, 3163), (22, 3164)):
            u_value, anchor_root = anchor_data(index)
            self.assertEqual(gcd(denominator, anchor_root), 1)
            numerator = u_value * denominator**2
            root_denominator = anchor_root * denominator
            value = homogeneous_square_value(numerator, root_denominator)
            point = PUBLISHED_POINTS[index]
            expected_root = (2 * point[1] + point[0]) * root_denominator**3
            self.assertEqual(expected_root.denominator, 1)
            self.assertEqual(value, expected_root.numerator**2)
            mapped = map_square_abscissa(
                numerator, root_denominator, abs(expected_root.numerator)
            )
            self.assertEqual(set(mapped), {point, point_negate(point)})

    def test_bitset_residue_sieve_matches_direct_modular_values(self) -> None:
        radius = 11
        prime = 13
        anchor_index = 22
        denominator = 3164
        residue_masks = build_offset_residue_masks(radius, (prime,))
        mask = allowed_offset_mask(
            anchor_index, denominator, prime, residue_masks
        )
        u_value, anchor_root = anchor_data(anchor_index)
        for offset in range(-radius, radius + 1):
            numerator = u_value * denominator**2 + offset * anchor_root**2
            root_denominator = anchor_root * denominator
            value = homogeneous_square_value(numerator, root_denominator) % prime
            is_residue = any(square * square % prime == value for square in range(prime))
            self.assertEqual(
                bool(mask & (1 << (offset + radius))),
                is_residue,
                (offset, value),
            )

    def test_coprime_count_and_small_companion_replay(self) -> None:
        for denominator in (1, 2, 6, 35, 3163):
            expected = sum(gcd(offset, denominator) == 1 for offset in range(1, 41))
            self.assertEqual(
                positive_coprime_offset_count(40, denominator), expected
            )
        companions = small_companion_lookup()
        self.assertIn(PUBLISHED_POINTS[0], companions)
        doubled = point_add(PUBLISHED_POINTS[0], PUBLISHED_POINTS[0])
        difference = point_add(PUBLISHED_POINTS[0], point_negate(PUBLISHED_POINTS[1]))
        self.assertEqual(companions[doubled][0], 2)
        self.assertEqual(companions[difference][0], 1)
        self.assertEqual(companions[difference][1], -1)

    def test_generated_artifact_is_complete_and_self_pinned(self) -> None:
        if not ARTIFACT.exists():
            self.skipTest("the denominator-sieve artifact has not been generated")
        data = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        self.assertEqual(
            data["reproduction"]["script_sha256"],
            sha256(SCRIPT.read_bytes()).hexdigest(),
        )
        self.assertEqual(
            data["status"], "bounded_search_no_certified_30th_point"
        )
        result = data["search_result"]
        self.assertTrue(result["search_complete"])
        self.assertFalse(result["wall_cap_reached"])
        self.assertFalse(result["stopped_at_overlap_gate"])
        self.assertEqual(
            result["processed_primitive_candidate_count"],
            result["declared_primitive_candidate_count"],
        )
        self.assertFalse(result["rank30_target_hit"])
        self.assertEqual(result["certified_independent_30th_point_count"], 0)
        self.assertLess(
            data["overlap_calibration"]["overlap_fraction"],
            data["parameters"]["overlap_stop_fraction"],
        )


if __name__ == "__main__":
    unittest.main()
