#!/usr/bin/env python3
"""Fast checks for the exact Fermigier two-section component certificate."""

import unittest
from fractions import Fraction

from probe_mestre_fermigier_two_section_local_continuation import (
    BASE_U,
    SECOND_SECTION,
)
from verify_mestre_fermigier_two_section_component import (
    component_coordinates,
    leading_invariant,
    leading_square,
)


Q = Fraction


class MestreFermigierTwoSectionComponentTest(unittest.TestCase):
    def test_local_seed_is_recovered(self):
        coordinates = component_coordinates(BASE_U)
        self.assertEqual(coordinates[4:6], (Q(-58, 11), Q(-17, 11)))
        self.assertEqual(coordinates[6:], SECOND_SECTION)
        self.assertEqual(leading_invariant(coordinates), leading_square(BASE_U) ** 2)


if __name__ == "__main__":
    unittest.main()
