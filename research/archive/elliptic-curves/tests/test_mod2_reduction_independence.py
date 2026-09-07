from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(CAS))

from mod2_reduction_independence import (  # noqa: E402
    combined_mod2_rank,
    finite_add,
    finite_curve_points,
    finite_multiply,
    find_mod2_reduction_certificate,
    find_two_torsion_certificate_prime,
    gf2_rank,
    mod2_reduction_signature,
    short_curve_has_no_rational_2_torsion_modular_certificate,
)
from nagao_1994 import rank13_base_changed_short_jacobian_coefficients  # noqa: E402
from verify_nagao_u42_rank17 import (  # noqa: E402
    CERTIFICATE_PRIMES,
    TWO_TORSION_CERTIFICATE_PRIME,
    load_basis,
)


Q = Fraction


class Mod2ReductionIndependenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.coefficients = rank13_base_changed_short_jacobian_coefficients(Q(42))
        cls.points = load_basis(
            ROOT / "archive/elliptic-curves/artifacts/generated-results/elliptic_nagao_u42_height_10000000.json"
        )

    def test_finite_group_enumeration_and_doubling(self) -> None:
        # y^2=x^3-x over F_5 has eight points and full rational 2-torsion.
        points = finite_curve_points(-1, 0, 5)
        self.assertEqual(len(points), 8)
        for point in points:
            doubled = finite_multiply(point, 2, -1, 5)
            self.assertEqual(doubled, finite_add(point, point, -1, 5))

    def test_binary_rank(self) -> None:
        self.assertEqual(gf2_rank(((1, 0, 1), (0, 1, 1), (1, 1, 0)), 3), 2)
        self.assertEqual(gf2_rank(((1, 0, 0), (0, 1, 0), (0, 0, 1)), 3), 3)

    def test_composite_modulus_and_off_curve_input_are_rejected(self) -> None:
        with self.assertRaises(ValueError):
            mod2_reduction_signature(self.coefficients, self.points, 15)
        tampered = self.points[:-1] + ((self.points[-1][0], self.points[-1][1] + 1),)
        with self.assertRaises(ValueError):
            mod2_reduction_signature(self.coefficients, tampered, 11)

    def test_u42_certificate_has_exact_full_column_rank(self) -> None:
        signatures = tuple(
            mod2_reduction_signature(self.coefficients, self.points, prime)
            for prime in CERTIFICATE_PRIMES
        )
        self.assertEqual(combined_mod2_rank(signatures, 17), 17)
        self.assertEqual(sum(signature.quotient_dimension for signature in signatures), 18)

    def test_automatic_certificate_search_reaches_full_rank(self) -> None:
        signatures = find_mod2_reduction_certificate(
            self.coefficients, self.points, prime_bound=137
        )
        self.assertEqual(tuple(item.prime for item in signatures), CERTIFICATE_PRIMES)
        self.assertEqual(combined_mod2_rank(signatures, 17), 17)

    def test_duplicate_point_destroys_full_rank(self) -> None:
        tampered = self.points[:-1] + (self.points[0],)
        signatures = tuple(
            mod2_reduction_signature(self.coefficients, tampered, prime)
            for prime in CERTIFICATE_PRIMES
        )
        self.assertLess(combined_mod2_rank(signatures, 17), 17)

    def test_mod31_proves_no_rational_two_torsion(self) -> None:
        self.assertTrue(
            short_curve_has_no_rational_2_torsion_modular_certificate(
                self.coefficients, TWO_TORSION_CERTIFICATE_PRIME
            )
        )
        self.assertEqual(find_two_torsion_certificate_prime(self.coefficients), 31)


if __name__ == "__main__":
    unittest.main()
