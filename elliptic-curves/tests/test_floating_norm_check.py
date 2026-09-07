import sys
import unittest
from math import inf, nan, ulp
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cas'))
from research_runtime.floating_norm_check import checked_distance_error


class FloatingNormCheck(unittest.TestCase):
    def test_retained_large_norm_roundoff(self):
        norm=2183933217
        error=checked_distance_error(norm,(norm+3*ulp(float(norm)))/4,2)
        self.assertEqual(error,1.430511474609375e-6)

    def test_scale_and_material_disagreement(self):
        for norm in (0,1,10**6,10**9,10**14):
            self.assertEqual(checked_distance_error(norm,norm/4,2),0)
            with self.assertRaises(ArithmeticError):
                checked_distance_error(norm,(norm+max(1,1024*ulp(float(norm))))/4,2)

    def test_nonfinite_and_negative_reports_fail_closed(self):
        for value in (nan,inf,-inf,-1):
            with self.assertRaises(ArithmeticError):checked_distance_error(1,value,2)

    def test_exact_inputs_required(self):
        for norm,degree in ((1.0,2),(True,2),(-1,2),(1,0),(1,True)):
            with self.assertRaises(ValueError):checked_distance_error(norm,1,degree)


if __name__=='__main__':unittest.main()
