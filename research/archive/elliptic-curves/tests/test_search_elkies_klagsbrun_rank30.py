#!/usr/bin/env python3

from fractions import Fraction
from pathlib import Path
import sys
import unittest


CAS = Path(__file__).resolve().parents[1] / "cas"
sys.path.insert(0, str(CAS))

from elkies_klagsbrun_rank29 import (  # noqa: E402
    PUBLISHED_POINTS,
    point_on_general_curve,
)
from search_elkies_klagsbrun_rank30 import (  # noqa: E402
    affine_substitute,
    build_charts,
    exact_linear_combination,
    point_add,
    point_multiply,
    point_negate,
    poly_evaluate,
    secant_slope,
    slope_discriminant,
    slope_to_points,
    x_polynomial,
)


Q = Fraction


class ElkiesKlagsbrunRank30SearchTests(unittest.TestCase):
    def test_slope_quartic_maps_a_known_secant_exactly(self) -> None:
        base = PUBLISHED_POINTS[0]
        target = PUBLISHED_POINTS[15]
        slope = secant_slope(0, 15)
        discriminant = slope_discriminant(base)
        value = poly_evaluate(discriminant, slope)
        # The quotient roots are the target and -(base+target).
        expected_other = point_negate(point_add(base, target))
        square_root = target[0] - expected_other[0]
        self.assertEqual(value, square_root * square_root)
        self.assertEqual(
            set(slope_to_points(base, slope, square_root)),
            {target, expected_other},
        )

    def test_exact_group_law(self) -> None:
        point = PUBLISHED_POINTS[0]
        self.assertIsNone(point_add(point, point_negate(point)))
        doubled = point_add(point, point)
        self.assertEqual(doubled, point_multiply(point, 2))
        self.assertTrue(point_on_general_curve(doubled))
        relation = [0] * 29
        relation[0] = 2
        self.assertEqual(exact_linear_combination(relation), doubled)

    def test_x_affine_chart_seeds(self) -> None:
        left = PUBLISHED_POINTS[0]
        right = PUBLISHED_POINTS[1]
        polynomial = affine_substitute(
            x_polynomial(), left[0], right[0] - left[0]
        )
        for parameter, point in ((Q(0), left), (Q(1), right)):
            square_root = 2 * point[1] + point[0]
            self.assertEqual(
                poly_evaluate(polynomial, parameter), square_root * square_root
            )

    def test_default_chart_counts(self) -> None:
        charts = build_charts(
            x_pair_height=10,
            x_offset_height=10,
            slope_offset_height=10,
            slope_pair_height=10,
            slope_pair_count=40,
        )
        counts = {
            kind: sum(chart.kind == kind for chart in charts)
            for kind in {chart.kind for chart in charts}
        }
        self.assertEqual(counts["x_pair_affine"], 406)
        self.assertEqual(counts["x_integer_offset"], 29)
        self.assertEqual(counts["slope_integer_offset"], 812)
        self.assertEqual(counts["slope_pair_affine"], 40)


if __name__ == "__main__":
    unittest.main()
