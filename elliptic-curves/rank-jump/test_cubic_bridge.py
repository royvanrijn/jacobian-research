"""Checks of the cubic arithmetic and the distinction from the norm conic."""
import unittest
import retrospective as r
import local_collision as lc
import cubic_bridge as cb


class CubicBridgeTests(unittest.TestCase):
    def test_cubic_field_arithmetic(self):
        K=cb.Cubic(-1,1)
        for x in ((2,1,0),(0,1,1),(3,-2,1)):
            x=tuple(map(r.F,x))
            self.assertEqual(K.mul(x,K.inverse(x)),K.one)
            for y in ((1,2,3),(0,1,0)):
                y=tuple(map(r.F,y))
                self.assertEqual(K.norm(K.mul(x,y)),K.norm(x)*K.norm(y))
                self.assertEqual(K.mul(x,y),K.mul(y,x))
        self.assertEqual(K.mul(K.square(K.theta),K.theta),K.sub(K.theta,K.one))

    def test_obstructed_chain_classes_against_original_ct(self):
        inp=r.read(lc.INPUT)
        for u,partners in ((-1,[1,1,1]),(1,[163,260,260])):
            base=next(x['W_u_basis'] for x in inp['rows'] if int(x['parameter_u'])==u)
            B=next(x['matrix'] for x in inp['ct'] if x['u']==u)
            for m,p in zip((317529,491700,631775),partners):
                self.assertEqual(lc.pairing(lc.coordinates(m,base),lc.coordinates(p,base),B),1)

    def test_independent_sage_symbolic_isogeny_and_norm_witnesses(self):
        try:from sage.all import QQ,PolynomialRing,NumberField
        except ImportError:self.skipTest('optional Sage exact arithmetic check')
        R=PolynomialRing(QQ,'a,b,x');a,b,x=R.gens();F=R.fraction_field()
        xp=x+a+b/x
        self.assertEqual(F(x*(x*x+a*x+b)*(1-b/x**2)**2),
                         F(xp*(xp*xp-2*a*xp+a*a-4*b)))
        inp=r.read(lc.INPUT);out=r.read(cb.OUTPUT)
        A,B=map(QQ,inp['anchor']['short_model_ainvariants'][3:])
        S=PolynomialRing(QQ,'t');t=S.gen();K=NumberField(t**3+A*t+B,'th')
        def element(cs):return sum(QQ(c)*K.gen()**i for i,c in enumerate(cs))
        delta=element(out['constant_quadratic_discriminant'])
        for w in out['norm_witnesses']:
            beta,s,t=[element(w[k]) for k in ('beta','s','t')]
            self.assertEqual(s*s-delta*t*t,beta)
            self.assertTrue(beta.norm().is_square())
        for model in out['models']:
            self.assertEqual(element(model['quadratic_discriminant']),delta*element(model['gamma'])**2)


if __name__=='__main__':unittest.main()
