import itertools
import random
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'cas'))
from visibility_lattice_fast import IntegerExactParity
from visibility_lattice_v2 import ExactParity


class IntegerCVPTests(unittest.TestCase):
    def test_all_small_parities_and_node_accounting(self):
        rng = random.Random(302)
        for n in range(1, 5):
            for _ in range(4):
                a = [[rng.randrange(-3, 4) for j in range(n)] for i in range(n)]
                g = [[sum(a[k][i]*a[k][j] for k in range(n)) + (i == j)
                      for j in range(n)] for i in range(n)]
                old, new = ExactParity(g), IntegerExactParity(g)
                for p in itertools.product((0, 1), repeat=n):
                    expected = old.solve(p, p)
                    self.assertEqual(expected, new.solve(p, p))
                    if expected['nodes'] > 1:
                        with self.assertRaises(RuntimeError):
                            new.solve(p, p, expected['nodes'] - 1)

    def test_ties_and_zero(self):
        solver = IntegerExactParity([[2, 0], [0, 2]])
        self.assertEqual(len(solver.solve([1, 1], [3, -3])['minima']), 4)
        self.assertEqual(solver.solve([0, 0], [2, 2])['minima'], [(0, 0)])

    def test_independent_small_ellipsoid(self):
        # q(x,y) >= 2(x*x+y*y); seeds have q <= 9, so [-5,5]^2
        # contains every candidate that could match or improve the seed.
        solver = IntegerExactParity([[4, 1], [1, 3]])
        for p in itertools.product((0, 1), repeat=2):
            rows = [(4*x*x + 2*x*y + 3*y*y, (x, y))
                    for x in range(-5, 6) for y in range(-5, 6)
                    if (x % 2, y % 2) == p]
            best = min(q for q, _ in rows)
            result = solver.solve(p, p)
            self.assertEqual(result['norm'], best)
            self.assertEqual(result['minima'], sorted(w for q, w in rows if q == best))

    def test_invalid_input(self):
        with self.assertRaises(ArithmeticError):
            IntegerExactParity([[1, 2], [2, 1]])
        solver = IntegerExactParity([[2]])
        with self.assertRaises(ArithmeticError):
            solver.solve([1], [0])
        with self.assertRaises(ArithmeticError):
            solver.solve([1, 0], [1, 0])


if __name__ == '__main__':
    unittest.main()
