"""The valuation gate must not promote a nonminimal integral equation."""
import sys,unittest
from fractions import Fraction as Q
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cas'))
from certify_next24_high_rank_minimal_v4 import minimality
class LocalMinimality(unittest.TestCase):
    def test_scaled_model_is_rejected(self):
        with self.assertRaisesRegex(ArithmeticError,'smaller integral model exists'):
            minimality(tuple(map(Q,[0,0,0,-16,64])))
    def test_normalized_obstruction_at_two(self):
        model=tuple(map(Q,[0,1,0,-24299287133098241035429438341240578512260214633,1445678882372225264547035837781374479248971363000092369735413098335463]))
        proof=minimality(model);p2=next(r for r in proof['local_exclusion_rows'] if r['prime']==2)
        self.assertEqual(p2['valuations_c4_c6_discriminant'],[6,9,13])
        self.assertEqual(len(p2['p_denominator_obstruction_indices']),12)
        self.assertTrue(all(p2['p_denominator_obstruction_indices']))
if __name__=='__main__':unittest.main()
