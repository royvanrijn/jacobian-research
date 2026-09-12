"""Ordinary fail-closed tests; no expensive arithmetic or class groups."""
from fractions import Fraction
from itertools import islice
from math import gcd
import unittest

import audit_wide_arithmetic_censoring as audit
import two_class_relation_core as core


class RelationTests(unittest.TestCase):
    def test_duplicate_and_zero_relations_do_not_increase_rank(self):
        stats = core.matrix_stats([3, 3, 0, 5], 3, 1)
        self.assertEqual(stats['rank_mod2'], 2)
        self.assertEqual(stats['deficiency'], 1)
        self.assertEqual(stats['rank_gain_beyond_canonical'], 1)

    def test_missing_columns_are_explicit(self):
        stats = core.matrix_stats([1], 5, 1)
        self.assertEqual(stats['deficiency'], 4)
        self.assertEqual(stats['zero_columns'], 4)
        self.assertEqual(stats['singleton_columns'], 1)

    def test_row_outside_factor_base_fails(self):
        for row in [-1, 8]:
            with self.assertRaises(ValueError):
                core.matrix_stats([row], 3)

    def test_deficiency_monotone_for_fixed_factor_base(self):
        values = [core.matrix_stats([3, 5, 6, 7][:i], 3)['deficiency'] for i in range(5)]
        self.assertEqual(values, sorted(values, reverse=True))

    def test_signed_exponents_have_correct_parity(self):
        self.assertEqual(sum((e%2)*(1<<i) for i, e in [[0,-3],[1,-2],[2,1]]), 5)

    def test_pair_budget_prefix_and_primitivity(self):
        a, b = list(core.primitive_pairs(100)), list(core.primitive_pairs(1000))
        self.assertEqual(a, b[:100])
        self.assertEqual(len(set(b)), 1000)
        self.assertTrue(all(y>0 and gcd(x,y)==1 for x,y in b))
        heights = [max(abs(x),y) for x,y in b]
        self.assertEqual(heights, sorted(heights))

    def test_support_strips_powers_and_keeps_unknown_cofactor(self):
        self.assertEqual(core.strip_support(-2**40*3**11*17, 2*3*5), 17)
        self.assertEqual(core.strip_support(2**40*3**11, 2*3*5), 1)
        with self.assertRaises(ValueError):
            core.strip_support(0, 30)

    def test_binary_form_transformation_exact(self):
        c = [7, -3, 2, 5]
        for matrix in ([1,7,0,1], [0,-1,1,0], [1,0,-2,1], [1,0,0,-1]):
            transformed = core.transform(c, matrix)
            p,q,r,s = matrix
            for x,y in [(2,3), (-7,2), (0,1)]:
                lhs = sum(transformed[i]*x**i*y**(3-i) for i in range(4))
                rhs = sum(c[i]*(p*x+q*y)**i*(r*x+s*y)**(3-i) for i in range(4))
                self.assertEqual(lhs,rhs)

    def test_nonunimodular_rebase_rejected(self):
        with self.assertRaises(ValueError):
            core.transform([1,2,3,4], [2,0,0,1])

    def test_reduction_transport_up_to_global_sign(self):
        c = [-2, -2, 0, -1]
        reduced, matrix, steps, capped = core.reduce_form(c)
        transported = core.transform(c, matrix)
        self.assertIn(reduced, [transported, [-v for v in transported]])
        self.assertFalse(capped)
        self.assertLessEqual(max(map(abs,reduced)), max(map(abs,c)))


def row(i, known=True, ram=10, size=400, rank=17):
    r = {'curve_key':str(i),'cohort':'prospective_broad_2080','selection_mode':'frozen_prospective',
         'family':'11952','local_status':'PASS' if known else 'UNKNOWN_TIMEOUT',
         'log2_abs_minimal_discriminant':size,'final_rank_lower_bound':rank,
         'rank_bucket':'17-18' if rank<=18 else '23'}
    if known:
        r['field_ramified_prime_count']=ram
    return r


class CensoringTests(unittest.TestCase):
    def test_unknown_bounds_not_zero_imputation(self):
        stats = audit.stats([row(1), row(2,False),row(3,ram=3),row(4,False)])
        self.assertEqual(stats['ramified_ge9_population_fraction_bounds'],[0.25,0.75])
        self.assertEqual(stats['ramification_known'],2)

    def test_all_unknown_retained(self):
        stats = audit.stats([row(1,False)])
        self.assertIsNone(stats['median_ramified_primes'])
        self.assertEqual(stats['ramified_ge9_population_fraction_bounds'],[0,1])

    def test_historical_excluded(self):
        r=row(1)
        r['cohort']='historical_external'
        with self.assertRaises(ValueError):
            audit.analyze([r])

    def test_size_cuts_blind_to_rank_and_completion(self):
        rows=[row(i, known=i%2==0, size=300+i, rank=17+i%2) for i in range(50)]
        first=audit.analyze(rows)
        for r in rows:
            r['final_rank_lower_bound']=23
            r['local_status']='UNKNOWN_TIMEOUT'
            r.pop('field_ramified_prime_count',None)
        self.assertEqual(first['base_size_quintile_cuts'],audit.analyze(rows)['base_size_quintile_cuts'])

    def test_markdown_headers_escape_pipes(self):
        self.assertIn(r'\|Delta\|',audit.text_table(['|Delta|'],[[3]]))


if __name__=='__main__':
    unittest.main()
