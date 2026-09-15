"""Failure controls for the exact k1 norm gate; no repeated cubic census."""
import copy
import importlib.util
import os
from pathlib import Path
import resource
import unittest
ROOT=Path(__file__).resolve().parents[1]
PATH=ROOT/'elkies-k3/scripts/verify_q80_genus_one_k1_reciprocity.py'
spec=importlib.util.spec_from_file_location('q80_k1',PATH);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
OUT=Path(os.environ.get('Q80_K1_TEST_DIR',str(m.DEFAULT)))

class CubicReciprocityControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        resource.setrlimit(resource.RLIMIT_CPU,(20,25));resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
        cls.rat=m.read(ROOT/m.DEGREE/'rational-fibres.json')['rows']
        cls.quad=m.read(ROOT/m.DEGREE/'quadratic-fibres.json')['rows']
        cls.agreement,_=m.rational_data(cls.rat)
        cls.cubic=m.read(OUT/'cubic-census.json')

    def test_complete_composite_and_special_cubic_rows(self):
        m.check_composite(m.read(OUT/'composite-norms.json'),self.rat,self.quad)
        m.check_cubic_result(self.cubic,self.agreement)

    def test_partial_cubic_census_is_not_complete(self):
        row=copy.deepcopy(self.cubic);row['orbits']-=1
        with self.assertRaisesRegex(ValueError,'complete cubic orbit'):
            m.check_cubic_result(row,self.agreement)

    def test_nodal_row_cannot_be_omitted(self):
        row=copy.deepcopy(self.cubic);row['zero_norm_rows']=[r for r in row['zero_norm_rows'] if r['smooth']]
        with self.assertRaisesRegex(ValueError,'zero-norm rows'):
            m.check_cubic_result(row,self.agreement)

    def test_unused_zero_alone_is_not_an_agreement(self):
        self.assertIsNone(m.permitted_agreement([5,5,0],{6:(2,11)},[]))
        self.assertEqual(m.permitted_agreement([5,5,0],{5:(2,11)},[]),(2,11))

    def test_agreement_and_disagreement_residue_collision_is_rejected(self):
        self.assertIsNone(m.permitted_agreement([5,5,0],{5:(2,11)},[2]))

    def test_nonzero_unused_norm_is_rejected(self):
        self.assertIsNone(m.permitted_agreement([5,6,3],{5:(2,11)},[]))

    def test_character_product_must_hold(self):
        with self.assertRaisesRegex(ValueError,'character product'):
            m.permitted_agreement([5,6,0],{5:(2,11)},[])

    def test_rational_character_roster_must_be_injective(self):
        row=copy.deepcopy(self.rat);indices=[i for i,r in enumerate(row) if r['rational_codes']]
        row[indices[1]]['rational_codes'][0]=row[indices[0]]['rational_codes'][0]
        with self.assertRaisesRegex(ValueError,'distinct nonzero'):
            m.rational_data(row)

    def test_matching_agreement_code_cannot_be_hidden(self):
        row=copy.deepcopy(self.cubic);r=row['zero_norm_rows'][0];code=next(iter(self.agreement))
        r['norm_codes']=[0,code,code]
        with self.assertRaisesRegex(ValueError,'does not supply agreement'):
            m.check_cubic_result(row,self.agreement)

if __name__=='__main__':unittest.main()
