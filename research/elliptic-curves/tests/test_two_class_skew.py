import unittest
from math import prod
import two_class_relation_core as core
import two_class_skew_core as skew


class SkewTests(unittest.TestCase):
    def test_exact_triangle_minimum(self):
        c = [13, -7, 3, 1]
        box = skew.rectangle(c, 8)
        self.assertEqual(box['A']*box['B'], 256)
        expected = min(sum(abs(x)*2**(i*k+(3-i)*(8-k)) for i, x in enumerate(c)) for k in range(9))
        self.assertEqual(int(box['norm_triangle_upper_bound']), expected)

    def test_sieve_matches_gcd_both_axes_and_degenerate_lines(self):
        for c in [[1, -1, 0, 1], [6, 5, -3, 2]]:
            support = [2, 3, 5, 7, 11, 13]
            roots = [(p, [x for x in range(p) if core.form_value(c, x, 1) % p == 0],
                         [x for x in range(p) if core.form_value(c, 1, x) % p == 0]) for p in support]
            for A, B in [(16, 2), (2, 16)]:
                points = [v for row in skew.sieve_lines(c, {'A': A, 'B': B}, roots) for v in row]
                self.assertEqual(len({(a, b) for a, b, _, _ in points}), len(points))
                for a, b, n, r in points:
                    self.assertEqual(n, core.form_value(c, a, b))
                    self.assertEqual(r, core.strip_support(n, prod(support)))


if __name__ == '__main__':
    unittest.main()
