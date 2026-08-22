#!/usr/bin/env python3
"""Regression test for the finite signed intersection audit."""

import unittest

from audit_mestre_fermigier_two_section_intersection import replay


class MestreFermigierTwoSectionIntersectionTest(unittest.TestCase):
    def test_seed_intersection_and_quotient_non_escape(self):
        result = replay()
        self.assertEqual(result["finite_intersection"], {"T": "-479/56", "x": "445/56", "y": "-1141635/28"})
        self.assertEqual(result["finite_reduction"]["visible_mod3_rank"], 9)
        self.assertEqual(result["finite_reduction"]["augmented_mod3_rank"], 9)


if __name__ == "__main__":
    unittest.main()
