"""Exact finite-allocation checks; no new bounds or point searches."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'cas'))
from height_portfolio_cost import ceil_fourth_root, solve
from verify_height_portfolio_cost import assignment_minima


class CostAllocation(unittest.TestCase):
    def test_rational_thresholds_and_boundary(self):
        for value in (F(1,1000), F(16), F(16)-F(1,100), F(16)+F(1,100),
                      F(125000**4), F(10**180,7)):
            H = ceil_fourth_root(value)
            self.assertLess(F((H-1)**4), value)
            self.assertGreaterEqual(F(H**4), value)

    def test_unequal_boxes_and_model_prices(self):
        table = [[1,256], [81,16]]
        result = solve(table, 1, [0,0], [1,1])
        self.assertEqual(result['optimum']['heights'], [1,2])
        self.assertEqual(F(result['optimum']['objective']), 5)
        self.assertEqual(F(result['best_uniform_by_exact_model_count']['2']['objective']), 8)
        self.assertEqual(solve(table,1,[0,0],[1,4])['optimum']['heights'], [3,0])
        self.assertEqual(solve(table,1,[10,10],[1,1])['optimum']['heights'], [3,0])
        paid = solve(table,1,[0,0],[1,1],100)
        self.assertEqual(paid['optimum']['heights'], [1,2])
        self.assertEqual(F(paid['optimum']['objective']), 105)

    def test_independent_exhaustion_of_small_boxes(self):
        # All assignments and all integer boxes are different completeness checks.
        for a,b,c,d in product((1,2), repeat=4):
            table = [[a**4,b**4],[c**4,d**4]]
            result = solve(table,1,[F(1,3),F(2,3)],[2,1],7)
            costs = []
            for h,k in product(range(3), repeat=2):
                if (h >= a or k >= b) and (h >= c or k >= d):
                    costs.append(F(7)+(F(1,3)+2*h*h if h else 0)+(F(2,3)+k*k if k else 0))
            expected = min(costs)
            self.assertEqual(F(result['optimum']['objective']), expected)
            self.assertEqual(assignment_minima(result['thresholds'],[F(1,3),F(2,3)],[2,1],F(7))[2], expected)

    def test_incomplete_costs_and_cap_fail_closed(self):
        for table,fixed,rates in (([],[],[]), ([[1],[1,2]],[0],[1]),
                                  ([[1]],[0],[-1]), ([[1]],[-1],[1]),
                                  ([[1.0]],[0],[1]), ([[1]],[],[1])):
            with self.assertRaises(ValueError):
                solve(table,1,fixed,rates)
        with self.assertRaises(ValueError):
            solve([[1,2],[3,4]],1,[0,0],[1,1],node_cap=1)


if __name__ == '__main__':
    unittest.main()
