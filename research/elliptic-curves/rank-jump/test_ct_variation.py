"""Independent checks for the retained simultaneous CT pencil."""
import unittest
import retrospective as r
import local_collision as lc
import ct_variation as cv


class CTVariationTests(unittest.TestCase):
    def test_full_two_torsion_action_has_only_scalar_centralizer(self):
        # S swaps the two vectors; R is the order-three action. Test all 16
        # linear endomorphisms, including singular ones.
        S=[2,1];R=[2,3];commuting=[]
        for packed in range(16):
            M=[packed&3,packed>>2]
            if all(cv.apply(M,cv.apply(G,v))==cv.apply(G,cv.apply(M,v)) for G in (S,R) for v in range(4)):
                commuting.append(M)
        self.assertEqual(commuting,[[0,0],[1,2]])

    def test_inconsistent_system_and_affine_solution(self):
        with self.assertRaises(ValueError):cv.solve([3,3],1,2)
        for b in range(8):
            A=[3,6,4];x=cv.solve(A,b,3)
            self.assertEqual(cv.apply(A,x),b)

    def test_chain_formula_against_original_matrices(self):
        inp=r.read(lc.INPUT);out=r.read(cv.OUTPUT)
        p=out['nondegenerate_difference_pencils'][0]
        E=p['nilpotent_chain_e_in_anchor_coordinates'];F=p['nilpotent_chain_f_in_anchor_coordinates']
        bases={int(x['parameter_u']):x['W_u_basis'] for x in inp['rows']}
        forms={x['u']:x['matrix'] for x in inp['ct']}
        for u in (-1,1):
            for i,e in enumerate(E):
                for j,f in enumerate(F):
                    actual=lc.pairing(lc.coordinates(e,bases[u]),lc.coordinates(f,bases[u]),forms[u])
                    expected=int(i+j==3)^(int(i+j==4) if u==1 else 0)
                    self.assertEqual(actual,expected)

    def test_independent_sage_linear_algebra_when_available(self):
        try:from sage.all import GF,matrix
        except ImportError:self.skipTest('optional Sage independent matrix check')
        out=r.read(cv.OUTPUT)
        for row in out['comparisons']:
            n=row['common_dimension']
            A,B=[matrix(GF(2),[[(a>>j)&1 for j in range(n)] for a in row[k]]) for k in ('first_form_packed_rows','second_form_packed_rows')]
            self.assertEqual((A+B).rank(),row['difference_rank'])
            if (row['u'],row['v'])==(-1,1):
                T=(A+B).inverse()*A;t=T.charpoly().parent().gen()
                self.assertEqual(T.charpoly(),t**10*(t+1)**2)
                self.assertEqual(T.minimal_polynomial(),t**5*(t+1))


if __name__=='__main__':unittest.main()
