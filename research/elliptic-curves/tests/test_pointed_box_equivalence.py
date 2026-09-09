"""Height-preserving projective transformations, including zero and infinity."""
from fractions import Fraction as F
from math import gcd
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'cas'))
from pointed_box_equivalence import box_key


class BoxEquivalenceTest(unittest.TestCase):
    def test_signed_permutations_and_scalar(self):
        m = (F(2,3),F(-7,5),F(11,13),F(17,19))
        a,b,c,d = m
        for n in (m,(-a,b,-c,d),(b,a,d,c),(-b,a,-d,c)):
            self.assertEqual(box_key(m), box_key(n))
            self.assertEqual(box_key(m), box_key([F(-23,29)*x for x in n]))

    def test_full_projective_height_boxes_preserved(self):
        def key(a,b):
            if b == 0:
                return None
            return F(a,b)
        for height in (1,2,3,7):
            pairs = [(a,b) for a in range(-height,height+1)
                     for b in range(-height,height+1) if gcd(a,b)==1]
            original = {key(a,b) for a,b in pairs}
            self.assertIn(None, original)
            self.assertIn(F(0), original)
            for operation in (lambda a,b:(-a,b),lambda a,b:(b,a),lambda a,b:(-b,a)):
                self.assertEqual(original,{key(*operation(a,b)) for a,b in pairs})

    def test_translation_not_merged_and_singular_rejected(self):
        self.assertNotEqual(box_key((1,0,0,1)),box_key((1,1,0,1)))
        with self.assertRaisesRegex(ArithmeticError,'singular'):
            box_key((1,2,2,4))


if __name__ == '__main__':
    unittest.main()
