"""Certificate corruption controls and geometric boundary regressions."""
import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    'common_quartic_replay', ROOT / 'elkies-k3/scripts/verify_common_quartic_singularity_locus.py')
C = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(C)


class CommonQuarticTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.packet = json.loads((C.DEFAULT / 'input.json').read_text())
        cls.result = json.loads((C.DEFAULT / 'result.json').read_text())
        cls.parents = json.loads((ROOT / C.SOURCE).read_text())['parents']

    def test_complete_replay(self):
        answer = C.verify(C.DEFAULT)
        self.assertEqual([c['genus'] for c in answer['controls']], [1, 0])
        self.assertEqual(len(answer['parents']), 4)
        self.assertFalse(answer['mw17_endpoint_complete'])

    def test_critical_norm_tamper(self):
        row = copy.deepcopy(self.result['parents'][0])
        row['critical_value_norm_monic_mod_p'][0] ^= 1
        with self.assertRaisesRegex(ValueError, 'critical value norm identity'):
            C.verify_parent(self.parents[0], row)

    def test_changed_parent_coefficients(self):
        parent = copy.deepcopy(self.parents[1])
        parent['A'][0] = str(C.Q(parent['A'][0])+1)
        with self.assertRaisesRegex(ValueError, 'critical polynomial identity'):
            C.verify_parent(parent, self.result['parents'][1])

    def test_missing_parent_is_not_an_exclusion(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory)
            (p / 'input.json').write_bytes((C.DEFAULT / 'input.json').read_bytes())
            result = copy.deepcopy(self.result)
            result['parents'].pop()
            (p / 'result.json').write_text(json.dumps(result))
            with self.assertRaisesRegex(ValueError, 'complete parent coverage'):
                C.verify(p)

    def test_constant_twist_and_duplicate_point_rejected(self):
        row = copy.deepcopy(self.result['controls'][0])
        row['d'] = [str(2*C.Q(v)) for v in row['d']]
        with self.assertRaisesRegex(ValueError, 'control input projection'):
            C.verify_control(self.packet['controls'][0], row)
        row = copy.deepcopy(self.result['controls'][0])
        row['points_and_sum'][1] = copy.deepcopy(row['points_and_sum'][0])
        with self.assertRaisesRegex(ValueError, 'point or sum attachment'):
            C.verify_control(self.packet['controls'][0], row)

    def test_infinity_and_higher_singularities(self):
        self.assertEqual(C.genus_from_multiplicities([1]*10),
                         {'infinity_multiplicity': 2, 'delta': 1, 'branch_degree': 10, 'genus': 4})
        self.assertEqual(C.genus_from_multiplicities([2]*4+[1]),
                         {'infinity_multiplicity': 3, 'delta': 5, 'branch_degree': 2, 'genus': 0})
        self.assertEqual(C.genus_from_multiplicities([3, 3, 2, 2, 1, 1])['genus'], 1)
        with self.assertRaisesRegex(ValueError, 'split or constant'):
            C.genus_from_multiplicities([2]*6)

    def test_false_base_rank_and_parent_rank_rejected(self):
        row = copy.deepcopy(self.result['controls'][0])
        row['rational_base_certificate']['three_times_point'] = ['0', '1']
        with self.assertRaisesRegex(ValueError, 'base point multiples'):
            C.verify_base_certificate(row)
        row = copy.deepcopy(self.result['controls'][1])
        row['inherited_parent_rank'] = 17
        with self.assertRaisesRegex(ValueError, 'parent rank scope'):
            C.verify_control(self.packet['controls'][1], row)

    def test_resultant_sign_leading_coefficient_and_common_root(self):
        p = 1009
        # f=2(t-1)(t-2), g=3(t-4); Res(f,g)=2*9*6=108.
        self.assertEqual(C.resultant([4, -6 % p, 2], [-12 % p, 3], p), 108)
        self.assertEqual(C.resultant([-12 % p, 3], [4, -6 % p, 2], p), 108)
        self.assertEqual(C.resultant([-1 % p, 1], [-2 % p, 1], p), -1 % p)
        self.assertEqual(C.resultant([4, -6 % p, 2], [-1 % p, 1], p), 0)


if __name__ == '__main__':
    unittest.main()
