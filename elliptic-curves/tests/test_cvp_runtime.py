"""Lazy CVP compared with exhaustive small controls, plus wide lazy requests."""
from fractions import Fraction
from itertools import product
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'cas'))
from research_runtime.cvp import VoronoiIterator


class CvpTests(unittest.TestCase):
    def test_exact_correlated_order_and_all_coset_minima(self):
        for gram, target in [([[5, 2], [2, 1]], [0, 0]),
                             ([[4, 1, 0], [1, 3, 1], [0, 1, 2]], ['1/3', '-1/5', '1/7'])]:
            n = len(gram); expected = {}
            target = list(map(Fraction, target))
            for z in product(range(-5, 6), repeat=n):
                mask = sum((c % 2) << i for i, c in enumerate(z))
                d = [Fraction(c, 2)-t for c, t in zip(z, target)]
                norm = sum(d[i]*gram[i][j]*d[j] for i in range(n) for j in range(n))
                expected[mask] = min(norm, expected.get(mask, norm))
            actual = list(VoronoiIterator(gram, target=target, include_zero=True))
            self.assertEqual([h.squared_distance for h in actual], sorted(expected.values()))
            self.assertEqual({h.mask: h.squared_distance for h in actual}, expected)

    def test_checkpoint_diversity_and_budget_preserve_unseen(self):
        gram = [[1 if i == j else 0 for j in range(4)] for i in range(4)]
        iterator = VoronoiIterator(gram, binding='state-1/policy-16')
        first = iterator.next_holes(3, diversity_window=2)
        restored = VoronoiIterator.resume(iterator.checkpoint(), binding='state-1/policy-16')
        self.assertEqual([h.record() for h in iterator], [h.record() for h in restored])
        self.assertEqual(len(first), 3)
        with self.assertRaises(ValueError):
            VoronoiIterator.resume(iterator.checkpoint(), binding='changed-basis')
        iterator = VoronoiIterator(gram)
        with self.assertRaises(TimeoutError):
            iterator.next_holes(10, node_budget=2)
        self.assertEqual(len(list(iterator)), 15)
        with self.assertRaises(ValueError):
            VoronoiIterator([[1, 2], [2, 1]])

    def test_high_dimension_is_lazy_and_pruning_is_queried(self):
        rank = 60
        gram = [[i+1 if i == j else 0 for j in range(rank)] for i in range(rank)]
        iterator = VoronoiIterator(gram)
        holes = iterator.next_holes(8, node_budget=5000)
        self.assertEqual(len({h.mask for h in holes}), 8)
        self.assertLess(len(iterator.heap), 5000)
        # This did not construct a parity list of size 2^60.
        self.assertLess(len(iterator.seen), 20)
        iterator = VoronoiIterator([[1, 0], [0, 1]])
        queued = iterator.next_holes(1, diversity_window=3)
        blocked = {h.mask for h in iterator.pending}
        self.assertEqual(iterator.next_holes(2, allowed=lambda h: h.mask not in blocked), [])


if __name__ == '__main__':
    unittest.main()
