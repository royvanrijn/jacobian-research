#!/usr/bin/env python3
"""Regression tests for the degree-five quantum-residue audit."""

import unittest

from sympy.polys.domains import GF, QQ

from explore_degree_five_a2_subprincipal import degree_five_sample
from explore_degree_five_quantum_residue import (
    EXCEPTIONAL_PROFILE,
    degree_five_family,
    degree_five_exceptional_family,
    first_order_audit,
    fifth_order_audit,
    laurent_monomials,
    third_order_family,
)
from analyze_degree_five_cubic_fifth_order import cubic_number_field


def _qq(poly):
    return {
        monomial: QQ(value.numerator, value.denominator)
        for monomial, value in poly.items()
    }


class QuantumResidueRegression(unittest.TestCase):
    def test_reconstructed_cubic_fifth_order_branch_coordinates(self):
        field, a, tau = cubic_number_field()
        self.assertEqual(
            94 * a**3 + 335 * a**2 + 400 * a + 160,
            field.zero,
        )
        self.assertEqual(
            8 * tau + 658 * a**2 + 1593 * a + 976,
            field.zero,
        )
        u34 = -(
            field(101668771215) * a**2 / field(2097152)
            + field(487549466415) * a / field(4194304)
            + field(18474132105) / field(262144)
        )
        u19 = (
            field(243)
            * (
                field(122116896574) * a**2
                + field(292049112895) * a
                + field(176583624080)
            )
            / field(10485760)
        )
        self.assertEqual(
            field(1099511627776) * u34**2
            + (
                field(106607433445539840) * a**2
                + field(255616334647787520) * a
                + field(154972252369059840)
            )
            * u34
            + field(30500510715227703150) * a**2
            + field(72079857867156799275) * a
            + field(43297155944526166650),
            field.zero,
        )
        self.assertEqual(
            field(94371840) * u19
            + (
                field(8766095360) * a
                + field(19587399680)
            )
            * u34
            + field(186954838846362) * a**2
            + field(447513347427885) * a
            + field(270836273654640),
            field.zero,
        )
        hbar7_period = (
            field(2189187)
            * (
                field(587583566) * a**2
                + field(1388701707) * a
                + field(831388850)
            )
            / field(83886080)
        )
        self.assertNotEqual(hbar7_period, field.zero)

    def test_sparse_family_matches_existing_rational_sample(self):
        expected_s, expected_t = degree_five_sample()
        got_s, got_t = degree_five_family(QQ, QQ(-1, 2), QQ.one)
        self.assertEqual(got_s, _qq(expected_s))
        self.assertEqual(got_t, _qq(expected_t))

    def test_laurent_band_contains_polynomials_and_only_declared_poles(self):
        polynomial = set(laurent_monomials(17, 0, 0))
        one_pole = set(laurent_monomials(17, 0, 1))
        self.assertLess(polynomial, one_pole)
        self.assertEqual(
            min(x_degree for x_degree, _, _ in one_pole),
            -1,
        )
        self.assertTrue(
            all(
                x_degree + q_degree + 3 * z_degree <= 17
                for x_degree, q_degree, z_degree in one_pole
            )
        )

    def test_known_sample_obstruction_survives_first_laurent_band(self):
        field = GF(32003)
        a = -field.one / field(2)
        S, T = degree_five_family(field, a, field.one)
        first = first_order_audit(S, T, field)
        self.assertEqual(first.kernel_dimension, 57)
        self.assertEqual(first.gauge_rank, 19)
        self.assertEqual(first.quotient_dimension, 38)

        polynomial_third = third_order_family(S, T, field, 0)
        self.assertEqual(polynomial_third.column_count, 1569)
        self.assertEqual(polynomial_third.operator_rank, 1527)
        self.assertEqual(len(polynomial_third.kernel), 42)
        polynomial_fifth = fifth_order_audit(
            S, T, polynomial_third, field, 0
        )
        self.assertEqual(polynomial_fifth.correction_rank, 594)
        self.assertTrue(polynomial_fifth.constant_outside_span)

        laurent_third = third_order_family(S, T, field, 1)
        self.assertEqual(laurent_third.column_count, 1719)
        self.assertEqual(laurent_third.operator_rank, 1675)
        self.assertEqual(len(laurent_third.kernel), 44)
        laurent_fifth = fifth_order_audit(S, T, laurent_third, field, 2)
        self.assertEqual(laurent_fifth.correction_rank, 721)
        self.assertTrue(laurent_fifth.constant_outside_span)

    def test_exceptional_chart_has_its_own_exact_obstruction_stratum(self):
        field = GF(32003)
        S, T = degree_five_exceptional_family(field, field.one)
        first = first_order_audit(S, T, field, EXCEPTIONAL_PROFILE)
        self.assertEqual(first.kernel_dimension, 50)
        self.assertEqual(first.gauge_rank, 16)
        self.assertEqual(first.quotient_dimension, 34)
        third = third_order_family(
            S,
            T,
            field,
            0,
            EXCEPTIONAL_PROFILE,
        )
        self.assertEqual(third.column_count, 1853)
        self.assertEqual(third.operator_rank, 1817)
        self.assertEqual(len(third.kernel), 36)
        fifth = fifth_order_audit(
            S,
            T,
            third,
            field,
            0,
            EXCEPTIONAL_PROFILE,
        )
        self.assertEqual(fifth.correction_rank, 769)
        self.assertTrue(fifth.constant_outside_span)


if __name__ == "__main__":
    unittest.main()
