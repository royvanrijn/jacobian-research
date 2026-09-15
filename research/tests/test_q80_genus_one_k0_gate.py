"""Positive root controls and failure cases for the exact Q80 k0 gates."""
import copy
import importlib.util
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1]
SPEC=importlib.util.spec_from_file_location('q80_k0_gate',ROOT/'elkies-k3/scripts/verify_q80_genus_one_k0_gate.py')
M=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(M)
class GateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp=tempfile.TemporaryDirectory(prefix='q80-k0-controls-');cls.executables=[]
        for kind in ['certify','replay']:
            exe=Path(cls.tmp.name)/kind
            subprocess.run(['g++','-O3','-std=c++17',str(ROOT/f'elkies-k3/scripts/{kind}_q80_quartic_norm_census.cpp'),'-o',str(exe)],check=True,capture_output=True,timeout=30)
            cls.executables.append(exe)
        cls.points=M.read(ROOT/M.NORM4/'norm4-sections.json')['records']
        cls.contacts=M.read(M.DEFAULT/'preflight/square-branch-pairs-preview.json')
    @classmethod
    def tearDownClass(cls):cls.tmp.cleanup()
    def constant_input(self,a,b):
        # Synthetic constant-family controls; no actual MW17 point assertion.
        return f'1 {a}\n1 {b}\n17\n'+'1 0\n1 1\n'*17
    def census(self,exe,a,b):
        p=subprocess.run([str(exe),'0','1'],input=self.constant_input(a,b),text=True,capture_output=True,check=True,timeout=10)
        return json.loads(p.stdout)
    def test_tower_finds_split_positive_and_regularized_zero(self):
        r=self.census(self.executables[0],130,0)
        self.assertEqual((r['orbits'],r['smooth_split'],len(r['candidate_rows'])),(17161,17161,17161))
        self.assertTrue(all(sorted(row[1:])==[0,1,130] for row in r['candidate_rows']))
    def test_independent_finds_split_positive(self):
        r=self.census(self.executables[1],130,0)
        self.assertEqual(r['first_failed_character'],[0]*17+[17161])
    def test_tower_rejects_irreducible_cubic(self):
        # X^3+X+3 has no F_131 root, and degree three stays irreducible in degree four.
        self.assertTrue(all((x**3+x+3)%131 for x in range(131)))
        r=self.census(self.executables[0],1,3);self.assertEqual(r['smooth_split'],0)
    def test_independent_rejects_irreducible_cubic(self):
        r=self.census(self.executables[1],1,3);self.assertEqual(r['smooth_split'],0)
    def test_repeated_quadratic_and_infinity_divisors(self):
        q=[88,62,1];f=M.mul(M.mul(q,q),[130,1])
        self.assertEqual(sorted(M.low_factors(f)),sorted([(q,2),([130,1],1)]))
        divisors=M.quadratic_divisors(f)
        self.assertIn(tuple(q),divisors);self.assertIn((130,1,0),divisors)
        self.assertNotIn((1,0,0),divisors)
    def test_omitted_contact_rejected(self):
        r=copy.deepcopy(self.contacts);r['distinct_root_contact_pairs'].pop()
        with self.assertRaisesRegex(ValueError,'pair roster'):M.check_contacts(r,self.points)
    def test_higher_coefficient_lift_changes_error(self):
        row=M.read(M.DEFAULT/'preflight/square-lifts-preview.json')['rows'][0]
        model=M.read(ROOT/M.SOURCE)['weierstrass_model'];before=M.first_lift(row['pair'],row['H'],self.points,model)
        model=copy.deepcopy(model);from fractions import Fraction
        key='B_coefficients_low_to_high';model[key][0]=str(Fraction(model[key][0])+131)
        after=M.first_lift(row['pair'],row['H'],self.points,model)
        self.assertNotEqual(before[2],after[2]);self.assertEqual(before[:2],(24,25))
    def test_nodal_quadratic_retains_multiplicity(self):
        rows=M.read(ROOT/M.DEGREE/'quadratic-fibres.json')['rows'];r=next(x for x in rows if x['t']==4161)
        self.assertEqual(sorted(r['norm_codes']),[0,13412,13412])
        self.assertEqual(len(r['roots']),3);self.assertEqual(len(set(r['roots'])),2)
if __name__=='__main__':unittest.main()
