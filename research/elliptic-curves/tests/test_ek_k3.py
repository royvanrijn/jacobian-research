from __future__ import annotations

from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

CAS_DIRECTORY = Path(__file__).resolve().parents[1] / "cas"
if str(CAS_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(CAS_DIRECTORY))

from ek_k3 import (
    EKK3Family,
    fraction_mod,
    hensel_lift_simple_root,
    legendre_symbol,
    polynomial_eval,
    rational_square_root,
    valuation,
)


class ElementaryArithmeticTests(unittest.TestCase):
    def test_fraction_mod_and_valuation(self) -> None:
        self.assertEqual(fraction_mod(Q(3, 2), 5), 4)
        self.assertEqual(valuation(Q(3 * 11**4, 5 * 11), 11), 3)
        with self.assertRaises(ValueError):
            fraction_mod(Q(1, 5), 5)
        with self.assertRaises(ValueError):
            valuation(Q(0), 5)

    def test_symbols_and_square_roots(self) -> None:
        self.assertEqual(legendre_symbol(0, 7), 0)
        self.assertEqual(legendre_symbol(2, 7), 1)
        self.assertEqual(legendre_symbol(3, 7), -1)
        self.assertEqual(rational_square_root(Q(49, 121)), Q(7, 11))
        self.assertIsNone(rational_square_root(Q(2)))
        self.assertIsNone(rational_square_root(Q(-1)))

    def test_generic_hensel_lift(self) -> None:
        # x^2 - 2 has the root 3 modulo 7.
        root = hensel_lift_simple_root((-2, 0, 1), 3, 7, 4)
        self.assertEqual(polynomial_eval((-2, 0, 1), root, 7**4), 0)
        self.assertEqual(root % 7, 3)
        with self.assertRaises(ValueError):
            hensel_lift_simple_root((-2, 0, 1), 1, 7, 2)
        with self.assertRaises(ValueError):
            hensel_lift_simple_root((0, 0, 1), 0, 7, 2)


class EKK3FamilyTests(unittest.TestCase):
    FAMILY = EKK3Family(Q(2, 5), Q(2))
    PARAMETER = Q(-1468, 21)

    def test_shimura_parameter_guard(self) -> None:
        with self.assertRaises(ValueError):
            EKK3Family(Q(2, 5), Q(3))

    def test_factor_formula_and_invariants(self) -> None:
        expected = (
            (Q(-186, 25), Q(8, 5)),
            (Q(-206, 25), Q(-12, 5)),
            (Q(-206, 25), Q(-236, 25)),
            (Q(-186, 25), Q(256, 25)),
            (Q(134, 25), Q(56, 25)),
            (Q(274, 25), Q(-36, 25)),
            (Q(254, 25), Q(84, 25)),
            (Q(154, 25), Q(-24, 25)),
        )
        self.assertEqual(self.FAMILY.factor_coefficients(), expected)
        a = self.FAMILY.a(self.PARAMETER)
        b = self.FAMILY.b(self.PARAMETER)
        invariants = self.FAMILY.invariants(self.PARAMETER)
        self.assertEqual(invariants["c4"], 16 * (4 * a**2 - 3 * b))
        self.assertEqual(invariants["c6"], 64 * a * (9 * b - 8 * a**2))
        self.assertEqual(invariants["discriminant"], 64 * b**2 * (a**2 - b))

    def test_power_roots_and_exact_valuations(self) -> None:
        expected = {
            11: (5, 5),
            17: (109, 8),
            19: (102, 7),
        }
        for prime, (residue, factor_index) in expected.items():
            roots = self.FAMILY.power_roots(prime, 2)
            match = [
                root for root in roots
                if root.residue == residue and root.factor_index == factor_index
            ]
            self.assertEqual(len(match), 1)
            self.assertEqual(
                (self.PARAMETER.numerator
                 - residue * self.PARAMETER.denominator) % prime**2,
                0,
            )
            factor_valuations = [
                valuation(value, prime)
                for value in self.FAMILY.b_factors(self.PARAMETER)
            ]
            self.assertEqual(factor_valuations[factor_index - 1], 2)
            self.assertEqual(sum(value != 0 for value in factor_valuations), 1)
            self.assertEqual(valuation(self.FAMILY.a(self.PARAMETER), prime), 0)
            self.assertEqual(
                valuation(self.FAMILY.invariants(self.PARAMETER)["discriminant"], prime),
                4,
            )

    def test_bad_and_good_local_data_are_distinct(self) -> None:
        bad = self.FAMILY.local_data(5, 11)
        self.assertFalse(bad.good_reduction)
        self.assertEqual(bad.vanishing_factors, (5,))
        self.assertFalse(bad.split_multiplicative)
        self.assertIsNone(bad.trace)

        good = self.FAMILY.local_data(0, 11)
        self.assertTrue(good.good_reduction)
        self.assertEqual((good.point_count, good.trace), (18, -6))
        self.assertEqual(good.vanishing_factors, ())
        self.assertIsNone(good.split_multiplicative)

    def test_published_points_satisfy_curve(self) -> None:
        points = self.FAMILY.known_points(self.PARAMETER)
        self.assertEqual(len(points), 9)
        self.assertEqual(len({x for x, _ in points}), 9)
        a = self.FAMILY.a(self.PARAMETER)
        b = self.FAMILY.b(self.PARAMETER)
        for x, y in points:
            self.assertEqual(y**2, x**3 + 2 * a * x**2 + b * x)

    def test_singular_exact_linear_root_and_valid_alternative(self) -> None:
        singular = Q(-28, 67)
        alternative = Q(155, 152)
        self.assertEqual(self.FAMILY.b_factors(singular)[4], 0)
        self.assertFalse(self.FAMILY.is_nonsingular(singular))
        self.assertTrue(self.FAMILY.is_nonsingular(alternative))
        self.assertEqual(valuation(self.FAMILY.b_factors(alternative)[4], 11), 4)


if __name__ == "__main__":
    unittest.main()
