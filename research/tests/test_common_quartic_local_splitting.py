"""Rejection controls for three distinct local roots, including bad primes."""
from fractions import Fraction as F
import importlib.util
from pathlib import Path
import unittest

PATH=Path(__file__).resolve().parents[1]/'elkies-k3/scripts/verify_common_quartic_local_splitting.py'
SPEC=importlib.util.spec_from_file_location('local_splitting',PATH)
V=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(V)

def ball(a,b,p,z):
 z=F(z)
 return {'centre':str(z),'f_valuation':V.valuation(z**3+a*z+b,p),
         'derivative_valuation':V.valuation(3*z*z+a,p)}

class LocalSplittingTests(unittest.TestCase):
 def test_three_exact_roots(self):
  V.check_balls(F(-1),F(0),5,[ball(-1,0,5,z) for z in [-1,0,1]])

 def test_three_copies_of_one_root_are_rejected(self):
  with self.assertRaisesRegex(ValueError,'duplicate root'):
   V.check_balls(F(-1),F(0),5,[ball(-1,0,5,1)]*3)

 def test_distinct_centres_approximating_one_root_are_rejected(self):
  with self.assertRaisesRegex(ValueError,'overlapping root balls'):
   V.check_balls(F(-1),F(0),5,[ball(-1,0,5,1+5**i) for i in [1,2,3]])

 def test_hensel_inequality_must_be_strict(self):
  with self.assertRaisesRegex(ValueError,'Hensel inequality'):
   V.check_balls(F(-3),F(9),3,[ball(-3,9,3,0)]*3)

 def test_nonintegral_cubic_is_rejected(self):
  with self.assertRaisesRegex(ValueError,'integral cubic'):
   V.check_balls(F(1,5),F(1),5,[])

 def test_singular_fibre_is_rejected(self):
  with self.assertRaisesRegex(ValueError,'singular fibre'):
   V.check_balls(F(0),F(0),5,[ball(0,0,5,0)]*3)


if __name__=='__main__':unittest.main()
