#!/usr/bin/env python3

import sys
from pathlib import Path
import unittest


CAS = Path(__file__).resolve().parents[1] / "cas"
sys.path.insert(0, str(CAS))

import icarm_curve90 as curve90
from mod2_reduction_independence import (
    combined_mod2_rank,
    find_mod2_reduction_certificate,
)
import search_icarm_curve273_rank31 as engine


class Curve90SaturationTest(unittest.TestCase):
    def test_all_recorded_points_lie_on_curve(self):
        pools = (
            curve90.POINTS,
            curve90.SATURATED_POINTS,
            curve90.REDUCED_SATURATED_POINTS,
        )
        self.assertTrue(all(curve90.on_curve(point) for pool in pools for point in pool))

    def test_half_point_relation_replays_exactly(self):
        engine.load_curve_data(90)
        doubled = engine.point_multiply(curve90.HALF_POINT, 2)
        relation = None
        for coefficient, point in zip(
            curve90.HALF_POINT_RELATION, curve90.POINTS, strict=True
        ):
            relation = engine.point_add(
                relation, engine.point_multiply(point, coefficient)
            )
        self.assertEqual(doubled, relation)
        self.assertEqual(curve90.HALF_POINT_RELATION[1], 1)

    def test_reduced_basis_has_exact_rank_19_certificate(self):
        short_points = tuple(
            curve90.to_short(point) for point in curve90.REDUCED_SATURATED_POINTS
        )
        signatures = find_mod2_reduction_certificate(
            curve90.short_coefficients(), short_points, prime_bound=3000
        )
        self.assertEqual(combined_mod2_rank(signatures, len(short_points)), 19)
        self.assertEqual(
            tuple(signature.prime for signature in signatures),
            (11, 13, 17, 23, 29, 31, 43, 47, 59, 67, 71, 73, 103, 127, 139),
        )


if __name__ == "__main__":
    unittest.main()
