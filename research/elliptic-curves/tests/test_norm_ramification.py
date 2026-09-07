import sys
from pathlib import Path
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'cas'))
from research_runtime.norm_ramification import preflight, cubic_data, isolated_remainder, Q


class NormRamificationTests(unittest.TestCase):
    f = [-1, -1, 0, 1]

    def test_shared_support_cancellation_is_preserved(self):
        alpha = [2, 1, 0]  # norm seven; pi(alpha)^2 is a global square.
        self.assertEqual(preflight(self.f, [alpha])['forced_zero_count'], 1)
        result = preflight(self.f, [alpha, alpha])
        self.assertEqual(result['forced_zero_count'], 0)
        self.assertEqual(result['additional_strict_classes'], 'UNKNOWN')
        self.assertEqual(result['whole_curve_rank_decision'], 'UNKNOWN')

    def test_peeling_rounds_have_valid_dependency_order(self):
        # First norm 7*61 forces coordinate zero. Only then is 7 isolated.
        result = preflight(self.f, [[8, 6, 1], [2, 1, 0]])
        self.assertEqual(result['norms'], ['427', '7'])
        self.assertEqual([[r['index'] for r in wave] for wave in result['rounds']], [[0], [1]])
        self.assertEqual(result['unresolved_indices'], [])

    def test_full_prime_support_is_removed(self):
        self.assertEqual(isolated_remainder(Q(7**9 * 61), 7, []), 61)
        self.assertEqual(isolated_remainder(Q(7**9 * 61), 1, [Q(7)]), 61)

    def test_unknown_is_not_an_acceptance_certificate(self):
        result = preflight(self.f, [[3, 1, 0]])  # norm 25, no isolated odd valuation.
        self.assertEqual(result['norms'], ['25'])
        self.assertEqual(result['unresolved_indices'], [0])
        self.assertEqual(result['additional_strict_classes'], 'UNKNOWN')

    def test_rational_scaling_and_norm(self):
        norms, _ = cubic_data(self.f, [['2/3', '1/3', 0]])
        self.assertEqual(norms, [Q(7, 27)])
        self.assertEqual(preflight(self.f, [['2/3', '1/3', 0]])['forced_zero_count'], 1)
        self.assertEqual(preflight(self.f, [[2, 0, 0]])['forced_zero_count'], 0)

    def test_reject_invalid_arithmetic(self):
        for f, elements in [([0, 0, 0, 1], [[1, 0, 0]]),
                            (self.f, [[0, 0, 0]]), (self.f, [[2.0, 1, 0]]),
                            (self.f, [[True, 1, 0]]), ([0, -1, 0, 1], [[0, 1, 0]])]:
            with self.assertRaises(ValueError):
                preflight(f, elements)


if __name__ == '__main__':
    unittest.main()
