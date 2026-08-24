#!/usr/bin/env python3
"""Regression test for the bounded component quotient screen."""

import unittest

from screen_mestre_fermigier_two_section_component_escape import replay


class MestreFermigierTwoSectionComponentEscapeTest(unittest.TestCase):
    def test_small_grid_has_no_second_section_escape(self):
        result = replay(height=2, t_bound=2, prime_bound=101)
        self.assertEqual(result["second_section_escape_count"], 0)
        self.assertEqual(result["second_section_mod2_escape_count"], 0)
        self.assertGreater(result["record_count"], 0)


if __name__ == "__main__":
    unittest.main()
