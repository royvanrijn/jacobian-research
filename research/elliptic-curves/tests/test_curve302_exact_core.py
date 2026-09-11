"""Exact regressions for the distinction between parity, Q and Z containment."""
import sys
import json
from pathlib import Path
import unittest

from sage.all import ZZ, matrix

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"cas"))
import run_curve302_exact_core as core


class ExactCoreTests(unittest.TestCase):
    def test_parity_core_need_not_lift_to_rational_core(self):
        a = matrix(ZZ, [[1, 2]])
        b = matrix(ZZ, [[1, 0]])
        self.assertTrue(core.lab.in_span(1, [1], 2))
        self.assertIsNone(core.word_witness(a, [1, 0]))
        result = core.common_core([("a", a), ("b", b)], 2)
        self.assertEqual(result["common_mod2_dimension"], 1)
        self.assertEqual(result["common_rational_dimension"], 0)
        self.assertEqual(result["common_integral_rank"], 0)

    def test_rational_containment_has_minimal_odd_multiple(self):
        g = matrix(ZZ, [[3, 0], [0, 5]])
        w = core.word_witness(g, [1, 1])
        self.assertEqual(w["minimum_positive_multiple"], 15)
        self.assertEqual(w["integer_coefficients_for_multiple"], [5, 3])
        self.assertFalse(w["integral"])

    def test_integral_axis_reconstructed_from_mixed_generators(self):
        g = matrix(ZZ, [[1, 2], [1, 3]])
        w = core.word_witness(g, [1, 0])
        self.assertTrue(w["integral"])
        self.assertEqual(w["integer_coefficients_for_multiple"], [3, -2])

    def test_full_integral_intersection_preserves_indices(self):
        a = matrix(ZZ, [[3, 0], [0, 1]])
        b = matrix(ZZ, [[1, 0], [0, 5]])
        r = core.common_core([("a", a), ("b", b)], 2)
        self.assertEqual(r["common_rational_dimension"], 2)
        self.assertEqual(r["common_integral_basis_hnf"], [[3, 0], [0, 5]])
        self.assertEqual(r["witnesses_in_each_run"]["a"]["rational_basis"][0]["minimum_positive_multiple"], 3)
        self.assertEqual(json.loads(json.dumps(r))["common_integral_rank"], 2)

    def test_exact_equality_accepts_large_integer_words(self):
        k = 2**60+1
        g = matrix(ZZ, [[1, k], [0, 1]])
        w = core.word_witness(g, [1, 0])
        self.assertEqual(w["integer_coefficients_for_multiple"], [1, -k])


if __name__ == "__main__":
    unittest.main()
