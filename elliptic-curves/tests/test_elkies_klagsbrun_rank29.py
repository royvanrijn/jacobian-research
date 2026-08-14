#!/usr/bin/env python3

from fractions import Fraction
from pathlib import Path
import sys
import unittest


CAS = Path(__file__).resolve().parents[1] / "cas"
sys.path.insert(0, str(CAS))

from elkies_klagsbrun_rank29 import (  # noqa: E402
    PUBLISHED_POINTS,
    curve_discriminant,
    from_short_point,
    point_on_general_curve,
    point_on_short_curve,
    published_short_points,
    short_weierstrass_coefficients,
    to_short_point,
)
from verify_elkies_klagsbrun_rank29 import build_artifact  # noqa: E402


Q = Fraction


class ElkiesKlagsbrunRank29Tests(unittest.TestCase):
    def test_all_published_points_are_exact(self) -> None:
        self.assertEqual(len(PUBLISHED_POINTS), 29)
        self.assertTrue(all(point_on_general_curve(point) for point in PUBLISHED_POINTS))

    def test_short_transport_is_exact(self) -> None:
        points = published_short_points()
        self.assertEqual(len(points), 29)
        self.assertEqual(points[0], to_short_point(PUBLISHED_POINTS[0]))
        self.assertEqual(from_short_point(points[0]), PUBLISHED_POINTS[0])
        self.assertTrue(all(point_on_short_curve(point) for point in points))

    def test_short_model_is_integral_and_nonsingular(self) -> None:
        coefficients = short_weierstrass_coefficients()
        self.assertEqual(coefficients[:3], (Q(0), Q(0), Q(0)))
        self.assertTrue(all(value.denominator == 1 for value in coefficients))
        self.assertNotEqual(curve_discriminant(), 0)

    def test_exact_rank29_certificate(self) -> None:
        artifact = build_artifact(certificate_prime_bound=500)
        certificate = artifact["finite_reduction_certificate"]
        self.assertEqual(certificate["combined_exact_rank_over_F2"], 29)
        self.assertEqual(certificate["two_torsion_certificate_prime"], 67)
        self.assertFalse(artifact["claim"]["target_rank_30_achieved"])


if __name__ == "__main__":
    unittest.main()
