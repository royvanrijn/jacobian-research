import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'cas'))
from visibility_generic_bank import enumerate_bank
from v3_transfer_orbits import shell_bank


class GenericBankTests(unittest.TestCase):
    def test_matches_reference_exact_ellipsoid(self):
        for g, u in [([[4, 1], [1, 3]], [[1, 0], [0, 1]]),
                     ([[2, 0], [0, 2]], [[1, 1], [0, 1]])]:
            old = shell_bank(g, u, bound=10, node_limit=10000)
            new = enumerate_bank(g, u, bound=10, shells=(8, 10), node_limit=10000)
            self.assertEqual(new['norm_counts'], old['norm_counts'])
            self.assertEqual(new['represented_parity_minima'], old['all_represented_parity_minima'])
            self.assertEqual(new['rows'], old['rows'])

    def test_budget_stops_without_bank(self):
        with self.assertRaises(RuntimeError):
            enumerate_bank([[2]], [[1]], bound=10, shells=(2,), node_limit=1)


if __name__ == '__main__':
    unittest.main()
