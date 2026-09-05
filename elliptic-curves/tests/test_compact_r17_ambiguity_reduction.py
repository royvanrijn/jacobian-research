from fractions import Fraction as F
from pathlib import Path
import sys
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cas'))
from audit_compact_r17_ambiguous import signature,pivots
from mod_l_reduction_independence import mod_l_reduction_signature


class ProjectiveAmbiguityReductionTests(unittest.TestCase):
    def test_nonintegral_point_reduces_to_identity_without_losing_other_columns(self):
        model=(0,0,0,-1,1);q=(F(1),F(1))
        # 4*(0,1) on y^2=x^3-x+1, with good reduction at 7.
        r=(F(-223,784),F(24655,21952))
        self.assertEqual(r[1]**2,r[0]**3-r[0]+1)
        old=mod_l_reduction_signature(model,[q],7,3)
        new=signature(model,[q,r],7,3)
        self.assertEqual(new.group_order,old.group_order)
        self.assertEqual(new.rows,tuple((row[0],0) for row in old.rows))
        self.assertEqual(len(pivots(new.rows,3)),1)

    def test_entire_nonintegral_sublist_has_zero_quotient_image(self):
        model=(0,0,0,-1,1);r=(F(-223,784),F(24655,21952))
        row=signature(model,[r,(r[0],-r[1])],7,3)
        self.assertTrue(row.quotient_dimension)
        self.assertEqual(pivots(row.rows,3),[])


if __name__=='__main__':unittest.main()
