#!/usr/bin/env python3
"""Regression checks for the compact Fermigier local continuation."""

import unittest
from fractions import Fraction

from probe_mestre_fermigier_two_section_local_continuation import (
    BASE_U,
    BASE_V,
    SECOND_SECTION,
    base_coordinates,
    formal_branch,
)


Q = Fraction


class MestreFermigierTwoSectionLocalContinuationTest(unittest.TestCase):
    def test_base_coordinates_match_the_screened_affine_lines(self):
        coordinates = base_coordinates()
        self.assertEqual(coordinates[4:6], (Q(-58, 11), Q(-17, 11)))
        self.assertEqual(coordinates[6:], SECOND_SECTION)

    def test_first_order_branch_has_the_exact_tangent(self):
        result = formal_branch(order=1)
        series = result["coordinate_series"]
        self.assertEqual(series["u"][:2], [str(BASE_U), "1"])
        self.assertEqual(series["v"][1], str(Q(7, 9)))
        self.assertEqual(series["second_intercept"][1], str(Q(4011, 7744)))
        self.assertEqual(series["second_slope"][1], str(Q(-7, 121)))
        self.assertEqual(
            result["reconstructed_rational_model"]["matches_recursive_formal_lift_through_order"],
            1,
        )


if __name__ == "__main__":
    unittest.main()
