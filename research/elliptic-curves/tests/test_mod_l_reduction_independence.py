from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(CAS))

from mod_l_reduction_independence import (  # noqa: E402
    combined_mod_l_rank,
    find_mod_l_reduction_certificate,
    find_no_rational_l_torsion_prime,
    gf_l_rank,
    mod_l_reduction_signature,
    no_rational_l_torsion_reduction_certificate,
)
from nagao_1994 import (  # noqa: E402
    primitive_visible_points,
    quartic_point_to_short_jacobian,
    short_jacobian_coefficients,
)
from nagao_1994_section7 import (  # noqa: E402
    SECTION7_CONSTRUCTION,
    SECTION7_LINEAR_COMPANION_SECTIONS,
)
from verify_nagao_section7_linear_sections import (  # noqa: E402
    EXPECTED_CERTIFICATE_PRIMES,
    EXPECTED_TORSION_CERTIFICATE_PRIME,
)


Q = Fraction


class ModLReductionIndependenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        parameter = Q(1)
        visible = tuple(
            quartic_point_to_short_jacobian(
                SECTION7_CONSTRUCTION, parameter, point
            )
            for point in primitive_visible_points(
                SECTION7_CONSTRUCTION, parameter
            )
        )
        cls.points = visible[:11] + (
            SECTION7_LINEAR_COMPANION_SECTIONS[0].jacobian_point(parameter),
        )
        cls.coefficients = short_jacobian_coefficients(
            SECTION7_CONSTRUCTION, parameter
        )

    def test_rank_over_prime_fields(self) -> None:
        rows = ((1, 2, 0), (2, 1, 0), (0, 0, 2))
        self.assertEqual(gf_l_rank(rows, 3, 3), 2)
        self.assertEqual(gf_l_rank(((1, 0), (0, 1)), 2, 5), 2)
        with self.assertRaises(ValueError):
            gf_l_rank(rows, 3, 4)
        with self.assertRaises(ValueError):
            gf_l_rank(((1, 2),), 3, 3)

    def test_trivial_mod3_quotient(self) -> None:
        # y^2=x^3-x over F_5 has order eight, so multiplication by three is
        # bijective and E(F_5)/3E(F_5) is trivial.
        signature = mod_l_reduction_signature(
            (Q(0), Q(0), Q(0), Q(-1), Q(0)), (), 5, 3
        )
        self.assertEqual(signature.group_order, 8)
        self.assertEqual(signature.multiple_subgroup_order, 8)
        self.assertEqual(signature.quotient_dimension, 0)
        self.assertEqual(signature.rows, ())

        # y^2=x^3+2 over F_7 has group (Z/3Z)^2.  This exercises the
        # multi-dimensional quotient branch rather than only cyclic quotients.
        full = mod_l_reduction_signature(
            (Q(0), Q(0), Q(0), Q(0), Q(2)), (), 7, 3
        )
        self.assertEqual(full.group_order, 9)
        self.assertEqual(full.multiple_subgroup_order, 1)
        self.assertEqual(full.quotient_dimension, 2)

    def test_section7_mod3_certificate_has_full_rank(self) -> None:
        signatures = find_mod_l_reduction_certificate(
            self.coefficients, self.points, modulus=3, prime_bound=200
        )
        self.assertEqual(
            tuple(signature.prime for signature in signatures),
            EXPECTED_CERTIFICATE_PRIMES,
        )
        self.assertEqual(combined_mod_l_rank(signatures, 12, 3), 12)
        self.assertEqual(
            find_no_rational_l_torsion_prime(
                self.coefficients, modulus=3, prime_bound=200
            ),
            EXPECTED_TORSION_CERTIFICATE_PRIME,
        )
        self.assertTrue(
            no_rational_l_torsion_reduction_certificate(
                self.coefficients, EXPECTED_TORSION_CERTIFICATE_PRIME, 3
            )
        )

    def test_invalid_moduli_and_tampered_points_are_rejected(self) -> None:
        with self.assertRaises(ValueError):
            mod_l_reduction_signature(self.coefficients, self.points, 11, 4)
        tampered = self.points[:-1] + (
            (self.points[-1][0], self.points[-1][1] + 1),
        )
        with self.assertRaises(ValueError):
            mod_l_reduction_signature(self.coefficients, tampered, 11, 3)


if __name__ == "__main__":
    unittest.main()
