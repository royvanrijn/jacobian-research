"""Tests of new subspace, local-character and reciprocity arguments."""
import itertools
from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parent))
import retrospective as r
import local_collision as lc
import reciprocity as rec


class LocalCollisionTests(unittest.TestCase):
    def test_orthogonal_and_intersection_against_exhaustive_vectors(self):
        for n in range(1,5):
            for a in range(1<<n):
                for b in range(1<<n):
                    rows=[a,b];ker=lc.orthogonal(rows,n)
                    brute={v for v in range(1<<n) if all((v&w).bit_count()%2==0 for w in rows)}
                    self.assertEqual({lc.lift(v,ker) for v in range(1<<len(ker))},brute)
                    self.assertEqual(lc.orthogonal(ker,n),lc.canonical(rows))
                    span_a=lc.canonical([a]);span_b=lc.canonical([b]);common=lc.intersection(span_a,span_b,n=n)
                    self.assertEqual({lc.lift(v,common) for v in range(1<<len(common))},{0,a}&{0,b})

    def test_coordinate_lifts_and_external_pairing_partner(self):
        base=[0b1101,0b1010,0b0100]
        for w in range(8):self.assertEqual(lc.coordinates(lc.lift(w,base),base),w)
        with self.assertRaises(ValueError):lc.coordinates(1,[2,4])
        # An isotropic line has zero self-pairing but is obstructed by an external partner.
        B=[[0,1],[1,0]]
        self.assertEqual(lc.pairing(1,1,B),0)
        self.assertEqual(lc.pairing(1,2,B),1)

    def test_nodal_finite_points_realize_split_unramified_kernel(self):
        # The local proof needs an s with s^2-d nonsquare for every d != 0.
        for p in (3,5,7,11,13,17,19,23,29,31):
            for d in range(1,p):
                chars=[]
                for s in range(p):
                    z=(s*s-d)%p
                    chars.append(0 if z==0 else 1 if pow(z,(p-1)//2,p)==1 else -1)
                self.assertEqual(sum(chars),-1)
                self.assertIn(-1,chars)

    def test_hilbert_bilinearity_square_invariance_and_known_values(self):
        for p in (2,3,5,7,'infinity'):
            for a,b,c in itertools.product((-3,-1,2,5),repeat=3):
                self.assertEqual(rec.hilbert(a*b,c,p),rec.hilbert(a,c,p)*rec.hilbert(b,c,p))
                self.assertEqual(rec.hilbert(a,b,p),rec.hilbert(b,a,p))
                self.assertEqual(rec.hilbert(a*c*c,b,p),rec.hilbert(a,b,p))
        self.assertEqual(rec.hilbert(-1,-1,2),-1)
        self.assertEqual(rec.hilbert(-1,-1,'infinity'),-1)
        for p in (3,5,7):self.assertEqual(rec.hilbert(p,p,p),-1 if p%4==3 else 1)

    def test_reciprocity_for_small_rationals(self):
        values=sorted({r.F(a,b) for a in range(-7,8) if a for b in range(1,8)})
        for a in values:
            for b in values:
                product=1
                for p in (2,3,5,7,'infinity'):product*=rec.hilbert(a,b,p)
                self.assertEqual(product,1)

    def test_against_pari_when_available(self):
        try:
            from sage.all import pari
        except ImportError:
            self.skipTest('optional independent PARI comparison requires sage -python')
        values=['-7','-5','-3','-2','-1','1','2','3','5','7','2/3','-5/7','7/4','-3/8']
        for a,b in itertools.product(values,repeat=2):
            for p in (2,3,5,7,11,'infinity'):
                q=0 if p=='infinity' else p
                self.assertEqual(rec.hilbert(a,b,p),int(pari.hilbert(pari(a),pari(b),q)))


if __name__=='__main__':unittest.main()
