from fractions import Fraction as F
from pathlib import Path
import sys
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'cas'))
from research_runtime.pointed_orbit_compression import compress
from research_runtime.search_state import raw_state
from research_runtime.finite_reduction import ReductionCache
from research_runtime.memory_store import MemoryFactStore


class PointedOrbitCompressionTests(unittest.TestCase):
    def test_pair_witness_and_same_discovered_basis(self):
        model = (0, 0, 0, -7, 10)
        basis = [(1, 2)]
        points = [(2, 2), (13, 46)]
        compressed = compress(model, basis, [1], points)
        self.assertEqual(compressed['kept_indices'], [0])
        self.assertEqual(compressed['skipped'][0]['partner_index'], 0)
        caches = [ReductionCache(MemoryFactStore()) for _ in range(2)]
        states = [raw_state(model, basis, cache=c, prime_bound=31) for c in caches]
        for i, roster in enumerate((points, [points[j] for j in compressed['kept_indices']])):
            for p in roster:
                states[i] = states[i].adjoin(p, cache=caches[i], extra_primes=(5, 7, 11, 13, 17, 19, 23, 29, 31))
        self.assertEqual(states[0].basis, states[1].basis)
        self.assertEqual(states[0].reductions, states[1].reductions)
        self.assertEqual(states[1].rank, 2)

    def test_fixed_point_and_invalid_subgroup_words(self):
        model = (0, 0, 0, -7, 10)
        # C=2Q fixes Q under the involution; retain Q.
        self.assertEqual(compress(model, [(2, 2)], [2], [(2, 2)])['kept_indices'], [0])
        for word in ([0], [F(1, 2)], []):
            with self.assertRaises(ValueError):
                compress(model, [(1, 2)], word, [(2, 2)])
        with self.assertRaises(ValueError):
            compress(model, [(1, 2)], [1], [(2, 3)])


if __name__ == '__main__':
    unittest.main()
