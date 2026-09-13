import copy
import importlib.util
from pathlib import Path
import tempfile
import unittest

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
SPEC=importlib.util.spec_from_file_location('q80_pairs_test',ROOT/'elkies-k3/scripts/verify_q80_complete_genus_one_pairs.py')
V=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(V)


class CompletePencilReplay(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.packet=V.read(V.DEFAULT/'input.json.gz');cls.ctx=V.context(cls.packet,521)
        cls.records=V.read(V.DEFAULT/'prime-521/000000.json.gz')['records']

    def test_another_valid_trace_cannot_be_relabelled(self):
        record=copy.deepcopy(self.records[1]);record['index']=self.records[0]['index']
        V.validate_records([record],self.ctx,identify=False)
        with self.assertRaisesRegex(ValueError,'height-bounded trace identification'):
            V.validate_records([record],self.ctx)

    def test_interpolation_degree_bound_is_required(self):
        record=copy.deepcopy(self.records[0]);record['Nx'] += [0]*(10-len(record['Nx']))
        with self.assertRaisesRegex(ValueError,'frame degree bounds'):V.validate_records([record],self.ctx)

    def test_rational_functions_are_scaled_jointly(self):
        row={'numerator_coefficients_low_to_high':['1/5','1/5'],'denominator_coefficients_low_to_high':['1/5']}
        self.assertEqual(V.primitive_rf(row,5),([1,1],[1]))

    def test_projective_scaling_and_infinity_are_retained(self):
        identity=np.eye(5,dtype=int).tolist();scaled=(2*np.eye(5,dtype=int)).tolist()
        points=V.images([{'index':0,'branch_matrix':identity},{'index':1,'branch_matrix':scaled}],5)
        self.assertEqual(V.Counter(V.buckets(points)),V.Counter({(0,1):6}))

    def test_zero_projective_infinity_is_rejected(self):
        bad=np.eye(5,dtype=int);bad[:,4]=0
        with self.assertRaisesRegex(ValueError,'zero projective branch image'):
            V.images([{'index':0,'branch_matrix':bad.tolist()}],5)

    def test_adaptive_coverage_retains_untested_vertices(self):
        stages=[(15,{0:3,1:3,2:12,3:12})]
        self.assertEqual(V.remaining_pairs(4,stages),[[0,1],[2,3]])
        stages.append((7,{}));self.assertEqual(V.remaining_pairs(4,stages),[[2,3]])
        stages.append((12,{}));self.assertEqual(V.remaining_pairs(4,stages),[])

    def test_last_prime_is_mandatory(self):
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)
            for name in ['prime-521.json.gz','prime-523.json.gz','embedding-input.json','embedding-result.json.gz']:(path/name).touch()
            with self.assertRaisesRegex(ValueError,'prime-541'):V.require_completed_stages(path)

    def test_bad_surface_prime_is_rejected(self):
        with self.assertRaisesRegex(ValueError,'bad model reduction'):V.context(self.packet,509)


if __name__=='__main__':unittest.main()
