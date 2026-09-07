import unittest
from itertools import combinations, product

import block_rank_theory as b


def alternating_from_bits(n, bits):
    matrix = [[0]*n for _ in range(n)]
    for bit, (i, j) in enumerate(combinations(range(n), 2)):
        matrix[i][j] = matrix[j][i] = bits >> bit & 1
    return matrix


class BlockRankTheoryTests(unittest.TestCase):
    def test_saturation_is_not_rank(self):
        # A=Z, M=2Z, P=1. The verified rational relation P in M tensor Q.
        result = b.signature_sandwich(1, 0, [[0]], [[1]], [[1]])
        self.assertEqual(result['observed_signature_increment'], 1)
        self.assertEqual(result['exact_quotient_rank'], 0)

    def test_torsion_is_not_rank(self):
        result = b.signature_sandwich(1, 1, [[1, 0]], [[0, 1]], [[2]])
        self.assertEqual(result['exact_quotient_rank'], 0)

    def test_even_new_point_can_be_invisible(self):
        result = b.signature_sandwich(1, 0, [[1, 0]], [[0, 0]], [])
        self.assertEqual((result['quotient_lower_bound'], result['quotient_upper_bound']), (0, 1))
        self.assertIsNone(result['exact_quotient_rank'])

    def test_signatures_and_relations_close_a_block(self):
        result = b.signature_sandwich(1, 0, [[1, 0, 0]],
                                     [[0, 1, 0], [0, 0, 1], [0, 1, 1]], [[1, 1, -1]])
        self.assertEqual(result['exact_quotient_rank'], 2)
        with self.assertRaises(ValueError):
            b.signature_sandwich(0, 0, [], [[1]], [[1]])
        with self.assertRaises(ValueError):
            b.mod_rank([[1]], 4)

    def test_modular_rank_cannot_pool_different_primes(self):
        points = [[2, 0], [0, 3]]
        self.assertEqual(b.mod_rank(points, 2), 1)
        self.assertEqual(b.mod_rank(points, 3), 1)
        self.assertEqual(b.qrank(points), 2)
        # Conversely, repeated information on one vector cannot be added.
        self.assertEqual(b.mod_rank([[1]], 2) + b.mod_rank([[1]], 3), 2)
        self.assertEqual(b.qrank([[1]]), 1)

    def test_overlap_and_diminishing_increment(self):
        generic = [[1, 0, 0, 0]]
        first = [[0, 1, 0, 0], [0, 0, 1, 0]]
        second = [[0, 0, 1, 0], [0, 0, 0, 1]]
        result = b.block_overlap(generic, first, second)
        self.assertEqual((result['union_rank'], result['intersection_rank']), (3, 1))
        self.assertEqual(result['second_increment_after_first'], 1)
        smaller = b.block_overlap(generic, first[:1], second)
        self.assertEqual(smaller['second_increment_after_first'], 2)

    def test_all_alternating_matrices_through_five(self):
        count = 0
        for n in range(6):
            for bits in range(1 << (n*(n-1)//2)):
                full = alternating_from_bits(n, bits)
                full_rank = b.mod_rank(full, 2, n)
                for d in range(n+1):
                    old = [row[:d] for row in full[:d]]
                    cross = [row[d:] for row in full[:d]]
                    result = b.radical_partner_bound(old, cross, n-d)
                    self.assertLessEqual(result['certified_pairing_rank_lower_bound'], full_rank)
                    if n-d == 1:
                        self.assertEqual(result['certified_pairing_rank_lower_bound'], full_rank)
                count += 1
        self.assertEqual(count, 1100)

    def test_unknown_partner_pairings_cannot_undo_the_bound(self):
        # Old H plus a three-dimensional radical, with two independent hits.
        old = alternating_from_bits(5, 1)
        cross = [[1, 0, 1], [0, 1, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]]
        result = b.radical_partner_bound(old, cross)
        self.assertEqual(result['certified_pairing_rank_lower_bound'], 6)
        for bits in range(8):
            new = alternating_from_bits(3, bits)
            full = [old[i]+cross[i] for i in range(5)]
            full += [[cross[i][j] for i in range(5)]+new[j] for j in range(3)]
            self.assertGreaterEqual(b.mod_rank(full, 2), 6)

    def test_zero_radical_column_does_not_mean_full_rank_known(self):
        result = b.radical_partner_bound([[0]], [[0, 0]])
        full = [[0, 0, 0], [0, 0, 1], [0, 1, 0]]
        self.assertEqual(result['certified_pairing_rank_lower_bound'], 0)
        self.assertEqual(b.mod_rank(full, 2), 2)
        with self.assertRaises(ValueError):
            b.radical_partner_bound([[1]], [[0]])

    def test_exclusion_requires_upper_and_strict_threshold(self):
        self.assertEqual(b.rank_exclusion(None, 0, 16, 20)['rank_upper_bound'], None)
        self.assertEqual(b.rank_exclusion(22, 0, 2, 20)['status'], 'NOT_EXCLUDED')
        result = b.rank_exclusion(22, 0, 4, 20)
        self.assertEqual(result['status'], 'TARGET_EXCLUDED')
        self.assertEqual(result['even_pairing_rank_sufficient_for_exclusion'], 4)
        with self.assertRaises(ValueError):
            b.rank_exclusion(18, 0, 16, 20, known_rank_lower=3)
        with self.assertRaises(ValueError):
            b.rank_exclusion(22, 0, 3, 20)


if __name__ == '__main__':
    unittest.main()
