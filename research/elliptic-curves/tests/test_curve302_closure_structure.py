from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import unittest

CAS = Path(__file__).resolve().parents[1] / "cas"
sys.path.insert(0, str(CAS))
import curve302_closure_structure_core as c


def synthetic_states():
    dimension = 4
    rows = []
    for mask in range(1 << dimension):
        costs = []
        for target in range(dimension):
            if mask >> target & 1:
                costs.append(None)
            elif target == 0:
                costs.append(5)
            elif target == 1:
                costs.append(4 if mask & 1 else 12)
            elif target == 2:
                costs.append(3 if (mask & 3) == 3 else 20)
            else:
                costs.append(2 if (mask & 7) == 7 else 25)
        rows.append({"state_mask":mask,"retained_numerators":costs})
    return rows


class ClosureStructureTests(unittest.TestCase):
    def test_cascade_and_waves(self):
        final, waves = c.closure_waves(synthetic_states(), 1, 5, 4)
        self.assertEqual(final, 15)
        self.assertEqual(waves, [[1],[2],[3]])

    def test_minimal_trigger_hyperedges(self):
        rows = synthetic_states()
        self.assertEqual(c.minimal_triggers(rows, 2, 5, 4), [3])
        self.assertEqual(c.minimal_triggers(rows, 3, 5, 4), [7])

    def test_arrival_minimax(self):
        arrivals, full = c.arrival_thresholds(synthetic_states(), 1, 4)
        self.assertEqual(arrivals, [0,4,4,4])
        self.assertEqual(full, 4)

    def test_closure_axioms(self):
        rows = synthetic_states()
        closures = [c.closure_mask(rows, mask, 5, 4) for mask in range(16)]
        self.assertTrue(all(c.closure_law_audit(closures, 4)["closure_axioms"].values()))

    def test_binary_subspace_signature_is_order_independent(self):
        self.assertEqual(c.rref_binary([0b0011,0b0101,0b0110],4), (5,6))
        self.assertEqual(c.binary_rank([3,5,6],4), 2)
        self.assertEqual(c.rref_binary([5,3],4), c.rref_binary([3,5],4))

    def test_primitive_vector(self):
        self.assertEqual(c.primitive_vector([-4,0,2]), (2,0,-1))
        self.assertEqual(c.primitive_vector([0,-3,3]), (0,1,-1))

    def test_strict_membership(self):
        strict = c.rref_binary([1,2],4)
        self.assertTrue(c.is_strict_mod2(3,strict))
        self.assertFalse(c.is_strict_mod2(4,strict))

    def test_spearman(self):
        self.assertAlmostEqual(c.spearman([1,2,3],[10,20,30]), 1.0)
        self.assertAlmostEqual(c.spearman([1,2,3],[30,20,10]), -1.0)


if __name__ == "__main__":
    unittest.main()
