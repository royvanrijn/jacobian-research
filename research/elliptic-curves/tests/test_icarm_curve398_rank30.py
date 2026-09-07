from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(CAS))

from icarm_curve398 import (  # noqa: E402
    POINTS,
    SHORT_POINTS,
    on_curve,
    short_coefficients,
)
from mod2_reduction_independence import (  # noqa: E402
    combined_mod2_rank,
    find_mod2_reduction_certificate,
    find_two_torsion_certificate_prime,
)


class IcarmCurve398Rank30Tests(unittest.TestCase):
    def test_all_public_points_are_exact(self) -> None:
        self.assertEqual(len(POINTS), 30)
        self.assertTrue(all(on_curve(point) for point in POINTS))
        self.assertEqual(len(set(POINTS)), 30)

    def test_short_transport_is_exact(self) -> None:
        coefficients = short_coefficients()
        coefficient_a = Fraction(coefficients[3])
        coefficient_b = Fraction(coefficients[4])
        for x_value, y_value in SHORT_POINTS:
            self.assertEqual(
                y_value * y_value,
                x_value**3 + coefficient_a * x_value + coefficient_b,
            )

    def test_finite_reduction_certificate_has_rank_30(self) -> None:
        coefficients = short_coefficients()
        torsion_prime = find_two_torsion_certificate_prime(
            coefficients, prime_bound=500
        )
        self.assertIsNotNone(torsion_prime)
        signatures = find_mod2_reduction_certificate(
            coefficients,
            SHORT_POINTS,
            prime_bound=2000,
        )
        self.assertEqual(combined_mod2_rank(signatures, 30), 30)


if __name__ == "__main__":
    unittest.main()
