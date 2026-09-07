#!/usr/bin/env python3
"""Exact completed-square and ratpoints parser tests for q12o5867 search."""

from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import os
import subprocess
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC_ROOT = ROOT / "elliptic-curves"
CAS = ELLIPTIC_ROOT / "cas"
sys.path.insert(0, str(ELLIPTIC_ROOT))
sys.path.insert(0, str(CAS))

from ecsearch.q12o5867_point_search import (  # noqa: E402
    affine_substitute_polynomial,
    completed_square_coefficients,
    integral_square_scaled_coefficients,
    parse_ratpoints_abscissae,
    points_from_completed_square_abscissa,
)
from elliptic_candidate_record import is_on_weierstrass_curve  # noqa: E402


RATPOINTS = ROOT / "tmp/ratpoints/root/usr/bin/ratpoints"
RATPOINTS_LIBRARY = ROOT / "tmp/ratpoints/root/usr/lib/x86_64-linux-gnu"


class Q12O5867PointSearchTests(unittest.TestCase):
    MODEL = (0, 0, 1, -7, 6)

    def test_completed_square_and_exact_inverse_map(self) -> None:
        self.assertEqual(completed_square_coefficients(self.MODEL), (25, -28, 0, 4))
        points = points_from_completed_square_abscissa(self.MODEL, Fraction(-3))
        self.assertEqual(points, ((Fraction(-3), Fraction(0)), (Fraction(-3), Fraction(-1))))
        self.assertTrue(all(is_on_weierstrass_curve(self.MODEL, point) for point in points))

    def test_quiet_abscissa_parser(self) -> None:
        output = "(1 : 0)\n(-3 : 1)\n(-2 : 1)\n(-3 : 1)\n"
        self.assertEqual(
            parse_ratpoints_abscissae(output),
            (Fraction(-3), Fraction(-2)),
        )
        with self.assertRaises(ValueError):
            parse_ratpoints_abscissae("diagnostic noise\n")

    def test_exact_section_normalized_affine_chart(self) -> None:
        original = completed_square_coefficients(self.MODEL)
        transformed = affine_substitute_polynomial(
            original, Fraction(-3), Fraction(1)
        )
        for value in (Fraction(-5, 7), Fraction(0), Fraction(1), Fraction(9, 4)):
            left = sum(
                coefficient * value**index
                for index, coefficient in enumerate(transformed)
            )
            x_coordinate = Fraction(-3) + value
            right = sum(
                coefficient * x_coordinate**index
                for index, coefficient in enumerate(original)
            )
            self.assertEqual(left, right)
        integral, scale = integral_square_scaled_coefficients(transformed)
        self.assertEqual(scale, 1)
        self.assertEqual(tuple(Fraction(value) for value in integral), transformed)
        rational_integral, rational_scale = integral_square_scaled_coefficients(
            (Fraction(1, 2), Fraction(1, 3))
        )
        self.assertEqual(rational_scale, 6)
        self.assertEqual(rational_integral, (18, 12))

    @unittest.skipUnless(
        RATPOINTS.exists() and RATPOINTS_LIBRARY.exists(),
        "unpacked ratpoints fixture is unavailable",
    )
    def test_unpacked_ratpoints_known_point_fixture(self) -> None:
        environment = os.environ.copy()
        environment["LD_LIBRARY_PATH"] = str(RATPOINTS_LIBRARY)
        completed = subprocess.run(
            [str(RATPOINTS), "25 -28 0 4", "20", "-q", "-y"],
            text=True,
            capture_output=True,
            check=True,
            timeout=10,
            env=environment,
        )
        abscissae = parse_ratpoints_abscissae(completed.stdout)
        self.assertIn(Fraction(-3), abscissae)
        points = points_from_completed_square_abscissa(self.MODEL, Fraction(-3))
        self.assertIn((Fraction(-3), Fraction(0)), points)


if __name__ == "__main__":
    unittest.main()
