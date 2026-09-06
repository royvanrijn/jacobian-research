import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cas'))
import select_retained_mw16_score_strata as selector


class MatchingTests(unittest.TestCase):
    def test_pairwise_caliper_rejects_opposite_extremes(self):
        p={'j_height_ratio_maximum':4,'parameter_height_ratio_maximum':2}
        rows=[{'j_height':str(h),'parameter_height':100} for h in [100,25,400]]
        self.assertTrue(selector.compatible(rows[:2],p))
        self.assertTrue(selector.compatible([rows[0],rows[2]],p))
        self.assertFalse(selector.compatible(rows,p))

    def test_caliper_boundary_and_parameter_height(self):
        p={'j_height_ratio_maximum':4,'parameter_height_ratio_maximum':2}
        rows=[{'j_height':'100','parameter_height':100},{'j_height':'400','parameter_height':200}]
        self.assertTrue(selector.compatible(rows,p))
        rows[1]['parameter_height']=201
        self.assertFalse(selector.compatible(rows,p))

    def test_arithmetic_height_ignores_homogeneous_coordinate_scale(self):
        f={'A_coefficients_low_to_high':['-1','2'], 'B_coefficients_low_to_high':['3','4']}
        a,ja,ha=selector.model(f,{'numerator':7,'denominator':5})
        b,jb,hb=selector.model(f,{'numerator':14,'denominator':10})
        self.assertEqual((ja,ha),(jb,hb))
        self.assertEqual(b[3],a[3]*2**8)
        self.assertEqual(b[4],a[4]*2**12)

if __name__=='__main__':unittest.main()
