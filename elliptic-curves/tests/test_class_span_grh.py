"""Small arithmetic controls and deliberate invalid-certificate mutations."""
from copy import deepcopy
from pathlib import Path
import random
import sys
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'cas'))
from sage.all import GF, matrix
from class_span_grh import RelationSpan, FieldAudit, quadratic_margin, sparse_row, verify_document
from class_span_fixtures import CASES, fixture


class ClassSpanTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fixtures = {name:fixture(poly) for name,poly in CASES}

    def test_seven_fields_against_certified_class_group_oracles(self):
        for name, (document, oracle) in self.fixtures.items():
            with self.subTest(field=name):
                result = verify_document(document)
                self.assertEqual(result['status'], 'CERTIFIED_UNDER_GRH')
                self.assertEqual(result['class_two_rank_upper_bound_under_grh'], oracle['unconditional_class_two_rank'])
                self.assertIsNone(result['exact_class_two_rank'])
                self.assertIsNone(result['elliptic_curve_rank'])

    def test_proper_span_is_unknown(self):
        for name in ['imaginary_c2','imaginary_c2_squared']:
            doc = deepcopy(self.fixtures[name][0])
            doc['anchors'] = []
            result = verify_document(doc)
            self.assertEqual(result['status'], 'UNKNOWN')
            self.assertIsNone(result['class_two_rank_upper_bound_under_grh'])

    def test_analytic_completion_with_unresolved_formal_coordinate(self):
        doc = deepcopy(self.fixtures['imaginary_trivial'][0])
        doc['relations'] = [r for r in doc['relations']
                            if all(doc['columns'][c]['p'] != 97 for c,e in r['factorization'])]
        result = verify_document(doc)
        self.assertGreater(result['formal_quotient_dimension'], 0)
        self.assertEqual(result['status'], 'CERTIFIED_UNDER_GRH')
        self.assertEqual(result['class_two_rank_upper_bound_under_grh'], 0)

    def test_missing_prime_coverage_rejected(self):
        doc = deepcopy(self.fixtures['imaginary_trivial'][0])
        doc['columns'] = []
        doc['relations'] = []
        doc['anchors'] = []
        with self.assertRaisesRegex(ValueError, 'missing prime'):
            verify_document(doc)

    def test_duplicate_and_forged_prime_columns_rejected(self):
        original = self.fixtures['imaginary_trivial'][0]
        doc = deepcopy(original)
        doc['columns'].append(doc['columns'][0])
        with self.assertRaisesRegex(ValueError, 'duplicate'):
            verify_document(doc)
        doc = deepcopy(original)
        doc['columns'][0]['f'] += 1
        with self.assertRaisesRegex(ValueError, 'claimed prime'):
            verify_document(doc)

    def test_wrong_field_or_integral_basis_rejected(self):
        for key in ['discriminant','integral_basis','signature']:
            doc = deepcopy(self.fixtures['imaginary_trivial'][0])
            if key == 'discriminant':doc['field'][key] = '-8'
            if key == 'signature':doc['field'][key] = [2,0]
            if key == 'integral_basis':doc['field'][key][0][0] = '2'
            with self.assertRaisesRegex(ValueError, 'field discriminant'):
                verify_document(doc)

    def test_forged_and_incomplete_relations_rejected(self):
        doc = deepcopy(self.fixtures['imaginary_c2'][0])
        relation = next(r for r in doc['relations'] if r['factorization'])
        relation['element'] = ['99991','0']
        with self.assertRaisesRegex(ArithmeticError, 'principal ideal'):
            verify_document(doc)
        doc = deepcopy(self.fixtures['imaginary_c2'][0])
        relation = next(r for r in doc['relations'] if r['factorization'])
        relation['factorization'] = []
        with self.assertRaisesRegex(ArithmeticError, 'principal ideal'):
            verify_document(doc)

    def test_relation_with_fractional_generator(self):
        doc = deepcopy(self.fixtures['imaginary_trivial'][0])
        factors = [[i,-c['e']] for i,c in enumerate(doc['columns']) if c['p'] == 2]
        doc['relations'].append({'element':['1/2','0'], 'factorization':factors})
        self.assertEqual(verify_document(doc)['status'], 'CERTIFIED_UNDER_GRH')

    def test_outside_columns_require_actual_cancellation(self):
        span = RelationSpan(3)
        span.add_relation(0b101)
        self.assertNotIn(0, span.analyze([])['known_coordinates'])
        span.add_relation(0b100)
        self.assertIn(0, span.analyze([])['known_coordinates'])

    def test_dependent_anchors_only_give_formal_upper_dimension(self):
        span = RelationSpan(3)
        span.add_relation(0b011)
        result = span.analyze([1,2,1])
        self.assertEqual(result['anchor_count'], 3)
        self.assertEqual(result['anchor_image_dimension_upper_bound'], 1)

    def test_random_row_spaces_against_sage_matrices(self):
        rng = random.Random(61723)
        for width in range(1,13):
            rows = [rng.randrange(1 << width) for _ in range(width)]
            anchors = [rng.randrange(1 << width) for _ in range(3)]
            def rank(vectors):
                return matrix(GF(2), [[(v >> c)&1 for c in range(width)] for v in vectors]).rank()
            span = RelationSpan(width)
            for row in rows:span.add_relation(row)
            result = span.analyze(anchors)
            self.assertEqual(result['relation_rank'], rank(rows))
            self.assertEqual(result['anchor_image_dimension_upper_bound'], rank(rows+anchors)-rank(rows))
            for c in range(width):
                self.assertEqual(c in result['known_coordinates'], rank(rows+anchors+[1 << c]) == rank(rows+anchors))

    def test_reducible_and_inexact_inputs_rejected(self):
        with self.assertRaises(ValueError):FieldAudit([-1,0,1])
        with self.assertRaises(ValueError):FieldAudit([1.0,0,1])
        with self.assertRaises(ValueError):sparse_row([[0,1],[0,2]], 1)
        with self.assertRaises(ValueError):quadratic_margin(-20,[2,0],[2],[],100)

    def test_unknown_even_powers_and_strict_cutoff(self):
        data = quadratic_margin(-20,[0,1],[2,9,100],[],100)
        self.assertEqual(data['prime_ideals_contributing'], 2)
        self.assertEqual(data['prime_power_terms'], 8)  # 2^1..2^6; 9^1..9^2
        self.assertEqual(data['unknown_odd_prime_power_terms'], 4)


if __name__ == '__main__':
    unittest.main()
