"""The finite character witness must be a smooth, nonvanishing nonsquare."""
import importlib.util
from pathlib import Path
import unittest

PATH=Path(__file__).resolve().parents[1]/'elkies-k3/scripts/verify_q80_single_branch_reciprocity.py'
SPEC=importlib.util.spec_from_file_location('single_branch',PATH)
V=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(V)

class SingleBranchCharacterTests(unittest.TestCase):
 def test_valid_nonsquare_character(self):
  self.assertEqual(V.check_character(-1,0,0,2,1,5),2)

 def test_square_character_does_not_obstruct(self):
  with self.assertRaisesRegex(ValueError,'not a nonzero nonsquare'):
   V.check_character(-1,0,1,2,1,5)

 def test_wrong_root_is_rejected(self):
  with self.assertRaisesRegex(ValueError,'not a cubic root'):
   V.check_character(-1,0,2,2,1,5)

 def test_wrong_section_is_rejected(self):
  with self.assertRaisesRegex(ValueError,'section point identity'):
   V.check_character(-1,0,0,2,2,5)

 def test_zero_ordinate_requires_a_different_witness(self):
  with self.assertRaisesRegex(ValueError,'meets two-torsion'):
   V.check_character(-1,0,4,0,0,5)

 def test_singular_fibre_is_rejected(self):
  with self.assertRaisesRegex(ValueError,'singular parent fibre'):
   V.check_character(0,0,0,1,1,5)


if __name__=='__main__':unittest.main()
