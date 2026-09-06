import unittest
from fractions import Fraction as F
from itertools import combinations, product
from math import gcd, isqrt, prod

import mechanism_theory as t


class MechanismTheoryTests(unittest.TestCase):
    def test_quotient_accounting_without_assuming_preserved_basis(self):
        self.assertEqual(t.quotient_dimension([[2, 0, 0], [4, 0, 0]],
                                             [[1, 1, 0], [3, -1, 0], [0, 0, 2]]), 2)
        self.assertEqual(t.quotient_dimension([], [[0, 2], [0, 4]]), 1)
        self.assertEqual(t.qrank([[1, F(1, 3)], [3, 1]]), 1)
        with self.assertRaises(ValueError):
            t.qrank([[1], [1, 2]])

    def test_signed_cycle_requires_characteristic_zero(self):
        balanced = [(0, 1, 1, -1), (1, 2, 1, -1), (0, 2, 1, -1)]
        self.assertEqual(t.signed_relation_bound(3, balanced)['quotient_upper_bound'], 1)
        balanced[-1] = (0, 2, 1, 1)
        self.assertEqual(t.signed_relation_bound(3, balanced)['quotient_upper_bound'], 0)
        self.assertEqual(t.signed_relation_bound(4, balanced)['quotient_upper_bound'], 1)
        with self.assertRaises(ValueError):
            t.signed_relation_bound(2, [(0, 1, 2, 1)])

    def test_all_small_signed_graphs(self):
        # Every absent/positive/negative edge pattern on four vertices.
        pairs = list(combinations(range(4), 2))
        for choices in product((0, 1, -1), repeat=len(pairs)):
            edges = [(i, j, 1, -s) for (i, j), s in zip(pairs, choices) if s]
            bound = t.signed_relation_bound(4, edges)
            self.assertEqual(bound['relation_rank'] + bound['quotient_upper_bound'], 4)

    def test_all_small_collision_graphs_against_mask_hitting(self):
        pairs = [sum(1 << i for i in pair) for pair in combinations(range(4), 2)]
        for included in range(1 << len(pairs)):
            clusters = [edge for i, edge in enumerate(pairs) if included >> i & 1]
            result = t.minimum_square_tests(4, clusters)
            feasible = [i for i in range(16) if all(i & edge for edge in clusters)]
            smallest = min(i.bit_count() for i in feasible)
            self.assertEqual(result['minimum_tests'], smallest)
            self.assertEqual(sorted(result['index_masks']), sorted(i for i in feasible if i.bit_count() == smallest))
        self.assertEqual(t.minimum_square_tests(4, [15])['minimum_tests'], 3)
        self.assertEqual(t.minimum_square_tests(4, [3, 5, 9])['minimum_tests'], 1)
        with self.assertRaises(ValueError):
            t.minimum_square_tests(17, [])

    def test_resultant_exact_formula_and_nonsquare_scaling(self):
        for a, b in product(range(1, 5), repeat=2):
            self.assertEqual(t.resultant([a, 0, 1], [b, 0, 1]), (a-b)**2)
        forms = [[1, 0, 1], [1, 1, 2]]
        self.assertEqual(t.lift_at_parameter(forms, 0, 1)['status'], 'NATIVE_LIFT')
        # Scaling both by 2 preserves a square product but destroys the lift.
        self.assertEqual(t.lift_at_parameter([[2*x for x in q] for q in forms], 0, 1)['status'], 'NONTRIVIAL_NATIVE_CLASS')

    def test_signs_infinity_branch_and_nonproduct(self):
        forms = [[1, 0, 1], [1, 1, 2]]
        negatives = [[-x for x in q] for q in forms]
        self.assertEqual(t.lift_at_parameter(negatives, 0, 1)['status'], 'NONTRIVIAL_NATIVE_CLASS')
        self.assertEqual(t.lift_at_parameter(forms, 1, 1)['status'], 'NOT_PRODUCT_POINT')
        self.assertEqual(t.lift_at_parameter([[1, 0, 1], [2, 1, 1]], 1, 0)['status'], 'NATIVE_LIFT')
        self.assertEqual(t.lift_at_parameter([[1, 0, 2], [2, 1, 2]], 1, 0)['status'], 'NONTRIVIAL_NATIVE_CLASS')
        self.assertEqual(t.lift_at_parameter([[0, 1, 1], [1, 0, 1]], 0, 1)['status'], 'BRANCH_REQUIRES_SEPARATE_CHART')
        with self.assertRaises(ValueError):
            t.lift_at_parameter(forms, 2, 2)
        with self.assertRaises(ValueError):
            t.quadratic_support([[1, 2, 1], [1, 0, 1]])
        with self.assertRaises(ValueError):
            t.quadratic_support([forms[0], forms[0]])

    def test_support_without_factoring_and_small_product_points(self):
        self.assertEqual(t.strip_support(2**7*3**4*5**2, 6), (2**7*3**4, 25))
        self.assertEqual(t.strip_support(49, 1), (1, 49))
        cases = 0
        for a, b in product(range(1, 5), repeat=2):
            if a == b:
                continue
            forms = [[a, 0, 1], [b, 0, 1]]
            for T, Z in product(range(-6, 7), range(0, 7)):
                if gcd(T, Z) != 1:
                    continue
                values = [a*Z*Z+T*T, b*Z*Z+T*T]
                if isqrt(prod(values))**2 != prod(values):
                    continue
                result = t.lift_at_parameter(forms, T, Z)
                self.assertEqual(result['status'] == 'NATIVE_LIFT', all(isqrt(v)**2 == v for v in values))
                cases += 1
        self.assertEqual(cases, 26)

    def test_simple_collision_character_count(self):
        for p in (5, 7, 11, 13):
            for a, b in combinations(range(p), 2):
                for leading in range(1, p):
                    chi = 1 if pow(leading, (p-1)//2, p) == 1 else -1
                    count = sum(pow(leading*(s-a)*(s-b) % p, (p-1)//2, p) == 1 for s in range(p))
                    self.assertEqual(count, (p-2-chi)//2)

    def test_trace_coset_strict_boundary(self):
        vectors = list(product(range(-3, 4), repeat=2))
        norm = lambda v: sum(x*x for x in v)
        for v, w in combinations(vectors, 2):
            if v == tuple(-x for x in w) or any((x-y) % 2 for x, y in zip(v, w)):
                continue
            self.assertGreaterEqual(norm(v)+norm(w), 4)
        self.assertEqual(norm((1, 1)) + norm((1, -1)), 4)


if __name__ == '__main__':
    unittest.main()
