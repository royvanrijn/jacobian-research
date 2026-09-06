from fractions import Fraction as F
import unittest
import retrospective as r
import oblique_split_cubic as block


class ObliqueBlockTests(unittest.TestCase):
    def test_independent_fingerprint_replay(self):
        anchor=r.read(block.OUTPUT)['anchor'];a1,a2,a3,a4,a6=map(F,anchor['integral_ainvariants'])
        joint=list(anchor['real_signatures']);offset=1
        for row in anchor['finite_fingerprints']:
            p=row['prime'];roots=[x for x in range(p) if (x**3+int(a2)*x*x+int(a4)*x+int(a6))%p==0]
            self.assertEqual(roots,row['roots']);observed=[]
            for i,P in enumerate(anchor['independent_points']):
                x=F(P[0]);bits=[]
                for root in roots:
                    value=(r.mod(x,p)-root)%p
                    if not value:value=(3*root*root+2*int(a2)*root+int(a4))%p
                    bits.append(int(pow(value,(p-1)//2,p)==p-1))
                sig=r.pack(bits);observed.append(sig);joint[i]|=sig<<offset
            self.assertEqual(observed,row['signatures']);offset+=3
        self.assertEqual(joint,anchor['joint_fingerprints']);self.assertEqual(r.rank(joint),3)

    def test_oblique_line_and_completed_ordinates(self):
        # The rational shear makes raw y zero; completed ordinates remain different.
        roots=[F(-1),F(2),F(1,2)];n=5
        for x in roots:
            self.assertEqual(x**3-F(3,2)*x*x-F(3,2)*x+1,0)
            self.assertEqual((n*x)**2,x**3+(n*n-F(3,2))*x*x-F(3,2)*x+1)
        self.assertEqual(len({abs(2*n*x) for x in roots}),3)

    def test_horizontal_norm_obstruction(self):
        self.assertEqual(F(2227,12),F(17*131,12))
        zeros=[(x,y) for x in range(17) for y in range(17) if (x*x+3*y*y)%17==0]
        self.assertEqual(zeros,[(0,0)])

    def test_universal_shared_cover_duplication_is_dependent(self):
        from sage.all import QQ,PolynomialRing
        R=PolynomialRing(QQ,names=('A','B','x'));A,B,x=R.gens();f=x**3+A*x+B
        N=x**4-2*A*x*x-8*B*x+A*A
        J=x**6+5*A*x**4+20*B*x**3-5*A*A*x*x-4*A*B*x-8*B*B-A**3
        self.assertEqual(N**3+16*A*N*f*f+64*B*f**3,J*J)
        # R2=N/(4f), H2=J/(8f^2), hence F(R2)=F(x)*H2^2.
        # The identity follows directly from doubling, so it cannot add a direction.
        K=R.fraction_field();R2=K(N)/(4*f)
        H2=-1+(3*x*x+A)*(x-R2)/(2*f)
        self.assertEqual(H2,K(J)/(8*f*f))

    def test_branch_specialization_has_only_torsion_in_the_new_block(self):
        from sage.all import QQ,EllipticCurve
        E=EllipticCurve(QQ,[0,-QQ(3)/2,0,-QQ(3)/2,1])
        for x in [-1,2,QQ(1)/2]:self.assertTrue((2*E(x,0)).is_zero())


if __name__=='__main__':unittest.main()
