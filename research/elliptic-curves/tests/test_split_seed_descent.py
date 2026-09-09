"""Run with sage -python -m unittest discover ...; no research point search."""
import copy
from pathlib import Path
import sys
import unittest
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cas'))
from sage.all import QQ,EllipticCurve
from split_seed_descent import build_frame,Classifier,point_record

class DescentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.E=EllipticCurve(QQ,[-4,1]);cls.P=cls.E([0,1]);cls.Q=cls.E([2,1])
        cls.frame=build_frame(cls.E,[cls.P],101);cls.engine=Classifier(cls.frame)

    def test_integral_dependence(self):
        r=self.engine.classify(-3*self.P,8)
        self.assertEqual(r['status'],'INHERITED_RATIONAL_SPAN')
        self.assertEqual(int(r['relation_multiplier'])*(-3*self.P),int(r['relation_word'][0])*self.P)

    def test_odd_saturation_cycle(self):
        e=Classifier(build_frame(self.E,[3*self.P],101));r=e.classify(self.P,8)
        self.assertEqual(r['reason'],'EXACT_HALVING_CYCLE')
        self.assertEqual(r['relation_multiplier'],'3')
        self.assertEqual(r['relation_word'],['1'])

    def test_even_hidden_independent_direction(self):
        r=self.engine.classify(4*self.Q,8)
        self.assertEqual(r['status'],'NEW_INDEPENDENT_DIRECTION')
        self.assertEqual(r['reason'],'GLOBAL_NONHALVING_ESCAPE')
        self.assertEqual(r['steps'],2)

    def test_no_candidate_reads_or_point_search(self):
        # The core API has only explicit mathematical inputs; neither it nor
        # its already imported finite kernel needs to read any artifact.
        with patch.object(Path,'open',side_effect=AssertionError('unexpected artifact read')):
            r=self.engine.classify(self.Q,8)
        self.assertEqual(r['reason'],'GLOBAL_NONHALVING_ESCAPE')

    def test_caps_remain_unknown(self):
        self.assertEqual(self.engine.classify(self.P,0)['status'],'UNKNOWN_STEP_CAP')

    def test_frame_tampering_fails(self):
        f=copy.deepcopy(self.frame);f['finite_rows'][0][0]^=1
        with self.assertRaisesRegex(ArithmeticError,'frame rows'):Classifier(f)

    def test_basis_sign_preserves_classification(self):
        other=Classifier(build_frame(self.E,[-self.P],101))
        for P in [self.Q,4*self.Q,-3*self.P]:
            self.assertEqual(self.engine.classify(P,8)['status'],other.classify(P,8)['status'])

    def test_exact_halves(self):
        halves,cert=self.engine.halves(2*self.Q)
        self.assertEqual(halves,[self.Q])
        self.assertEqual(cert['rational_halves'],[point_record(self.Q)])

if __name__=='__main__':unittest.main()
