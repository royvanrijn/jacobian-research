"""Adversarial checks of the portable height-bound certificate."""
from copy import deepcopy
from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'elliptic-curves/cas'))
from verify_pointed_height_bounds import verify
from verify_pointed_minimal_neighbours import verify as verify_neighbours
from verify_pointed_height_portfolio import run as verify_portfolio

FOLDER=ROOT/'artifacts/generated-results/elliptic-curves/pointed_height_bounds_first_chart_v1'


class HeightBoundCertificates(unittest.TestCase):
    def packet(self,kind='factor_free'):
        return json.loads((FOLDER/f'{kind}-bounds.json').read_text())

    def test_both_models_replay(self):
        for kind in ('factor_free','preconditioned_full'):
            self.assertEqual(verify(self.packet(kind),ROOT)['status'],'PASS_FRACTION_AND_STURM_REPLAY')

    def test_changed_covering_map_rejected(self):
        p=self.packet();p['numerator'][0]=str(int(p['numerator'][0])+1)
        with self.assertRaises(ArithmeticError):verify(p,ROOT)

    def test_missing_projective_chart_rejected(self):
        p=self.packet();p['bezout'].pop()
        with self.assertRaises(ArithmeticError):verify(p,ROOT)

    def test_omitted_local_residue_rejected(self):
        p=self.packet()
        def first_split(n):
            if n.get('roots'):return n
            for c in n.get('children',[]):
                found=first_split(c['node'])
                if found:return found
        n=first_split(p['local_trees'][0]['affine'])
        self.assertIsNotNone(n)
        n['roots'].pop();n['children'].pop()
        with self.assertRaises(ArithmeticError):verify(p,ROOT)

    def test_understated_cancellation_rejected(self):
        p=self.packet();p['local_trees'][0]['upper']-=1
        with self.assertRaises(ArithmeticError):verify(p,ROOT)

    def test_false_real_enclosure_rejected(self):
        p=self.packet();p['real']['lower']=str(Fraction(p['real']['lower'])*Fraction(1000001,1000000))
        with self.assertRaises(ArithmeticError):verify(p,ROOT)

    def test_minimal_models_and_joint_bound_replay(self):
        self.assertEqual(verify_neighbours(FOLDER,ROOT)['status'],'PASS_FRACTION_MINIMALITY_AND_TRANSPORT_REPLAY')
        result=verify_portfolio(FOLDER,ROOT)
        self.assertEqual(result['selected_by_address_work_proxy']['models'],[0,2])
        self.assertLess(result['display']['address_work_ratio_to_best_single'],1)


if __name__=='__main__':unittest.main()
