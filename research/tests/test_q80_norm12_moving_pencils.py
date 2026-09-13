import copy
import importlib.util
from pathlib import Path
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('q80_norm12_test', ROOT / 'elkies-k3/scripts/verify_q80_norm12_moving_pencils.py')
V = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(V)


class MovingPencilReplay(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.packet = V.read(V.DEFAULT / 'input.json.gz')
        cls.ctx = V.trace_context(cls.packet, 521)
        cls.records = V.read(V.DEFAULT / 'frames-521.json.gz')['records']

    def test_valid_other_trace_cannot_be_relabelled(self):
        record = copy.deepcopy(self.records[1])
        record['index'] = 0
        V.validate_norm12([record], self.ctx, identify=False)
        with self.assertRaisesRegex(ValueError, 'height-bounded norm12 trace identification'):
            V.validate_norm12([record], self.ctx)

    def test_interpolation_degree_bound_is_mandatory(self):
        record = copy.deepcopy(self.records[0])
        record['Ny'] += [0] * (20-len(record['Ny']))
        with self.assertRaisesRegex(ValueError, 'norm12 degree bounds'):
            V.validate_norm12([record], self.ctx)

    def test_monic_normalization_is_required_for_leading_cancellation(self):
        record = copy.deepcopy(self.records[0])
        record['h'][-1] = 2
        with self.assertRaisesRegex(ValueError, 'norm12 degree bounds'):
            V.validate_norm12([record], self.ctx)

    def test_fixed_member_cannot_replace_the_moving_pencil(self):
        record = copy.deepcopy(self.records[0])
        record['branch_matrix'] = [[row[4], 0, 0, 0, row[4]] for row in record['branch_matrix']]
        with self.assertRaisesRegex(ValueError, 'degree28 moving branch identity'):
            V.validate_norm12([record], self.ctx)

    def test_23_trace_comparisons_are_required(self):
        ctx = self.ctx | {'identification_sites': self.ctx['identification_sites'][:22],
                          'identified_word_values': self.ctx['identified_word_values'][:22]}
        with self.assertRaisesRegex(ValueError, '23 affine trace comparisons required'):
            V.validate_norm12([self.records[0]], ctx)

    def test_every_member_of_a_shared_projective_bucket_is_retained(self):
        identity = np.eye(5, dtype=int).tolist()
        old = [{'index': 10, 'branch_matrix': identity}]
        deep = [{'index': i, 'branch_matrix': identity} for i in (0, 1)]
        self.assertEqual(V.cross_pairs(old, deep, 5), {(10, 0), (10, 1)})

    def test_first_prime_alone_cannot_certify_completion(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)
            for prefix in ('prime', 'frames'):
                (path / (prefix+'-521.json.gz')).touch()
            with self.assertRaisesRegex(ValueError, 'prime-523'):
                V.require_completed_stages(path)


if __name__ == '__main__':
    unittest.main()
