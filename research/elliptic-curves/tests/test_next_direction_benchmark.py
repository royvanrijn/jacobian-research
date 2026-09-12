"""Exact transport and certificate-serialization regressions; no point search."""
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cas'))
from next_direction_benchmark import homogeneous, change_coordinates, SHEARS
from search_observability import multiply
from pointed_box_equivalence import box_key
from research_runtime.store import digest
from memory_rank_certificate import checked_rank


def evaluate(coefficients, n, d):
    degree = len(coefficients)-1
    return sum(F(v)*n**i*d**(degree-i) for i,v in enumerate(coefficients))


class NextDirectionTests(unittest.TestCase):
    def test_projective_transport_including_poles_and_infinity(self):
        for coefficients in ([3,-2,5], [7,3,-2,0,5]):
            for W in SHEARS:
                changed = homogeneous(coefficients,W)
                a,b,c,d=W
                for n,q in [(0,1),(1,0),(1,1),(-1,1),(2,3),(3,2)]:
                    self.assertEqual(evaluate(changed,n,q),evaluate(coefficients,a*n+b*q,c*n+d*q))

    def test_cover_discriminant_and_composed_map(self):
        P,Q = [1,2,3,4,5],[3,-2,1]
        disc = [4*P[i]+sum(Q[j]*Q[i-j] for j in range(3) if 0<=i-j<3) for i in range(5)]
        mapping = {'first_matrix':['2','1','1','1'],'second_matrix':['1','3','0','1'],
            'matrix':['2','7','1','4'],'reduced_P':P,'reduced_Q':Q,'discriminant_quartic':disc}
        for W in SHEARS:
            changed = change_coordinates(mapping,W)
            for n,d in [(1,0),(0,1),(2,3),(-1,1)]:
                self.assertEqual(evaluate(changed['discriminant_quartic'],n,d),
                    4*evaluate(changed['reduced_P'],n,d)+evaluate(changed['reduced_Q'],n,d)**2)
            self.assertEqual(tuple(map(F,changed['matrix'])),multiply(tuple(map(F,changed['first_matrix'])),tuple(map(F,changed['second_matrix']))))

    def test_height_symmetries_are_deduplicated_but_shears_are_not(self):
        M=(F(2),F(3),F(1),F(2))
        for W in [(1,0,0,1),(-1,0,0,1),(0,1,1,0),(0,-1,1,0)]:
            self.assertEqual(box_key(M),box_key(multiply(M,W)))
        self.assertNotEqual(box_key(M),box_key(multiply(M,SHEARS[1])))

    def test_real_finite_certificate_survives_json_roundtrip(self):
        proof=checked_rank((F(0),F(0),F(0),F(-1),F(1)),[(F(0),F(1))],[3,5,7,11,13,17,19],3)
        self.assertEqual(proof['rank_lower_bound'],1)
        self.assertEqual(digest(proof),digest(json.loads(json.dumps(proof))))


if __name__ == '__main__':
    unittest.main()
