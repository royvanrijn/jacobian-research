import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cas'))
from pointed_quartic_search import PointedQuarticSearch
from pari_pointed_backend import witnesses

class PariWitnessTests(unittest.TestCase):
    def setUp(self):
        self.search=PointedQuarticSearch(curve=[0,-2],subgroup=[(3,5)],centre={'coefficients':[1]},coordinate_policy={'kind':'raw','matrix':['-27','10','10','0']})
    def test_infinity_is_a_separate_exact_point(self):
        hits,points,ms=witnesses(self.search,'MS|0\nDONE|0\n','bounded_search_complete',5)
        self.assertTrue(any(d==0 for n,d,r in hits))
        self.assertIn((3,-5),points)
        self.assertEqual(ms,0)
        # Timeout output has no asserted affine coverage, but infinity is exact.
        self.assertEqual(witnesses(self.search,'X|garbled','bounded_search_timeout',5)[:2],(hits,points))
    def test_rejects_partial_framing_and_out_of_box(self):
        with self.assertRaises(ArithmeticError):witnesses(self.search,'MS|0\nX|0\n','bounded_search_complete',5)
        with self.assertRaises(ArithmeticError):witnesses(self.search,'MS|0\nX|6\nDONE|1\n','bounded_search_complete',5)
if __name__=='__main__':unittest.main()
