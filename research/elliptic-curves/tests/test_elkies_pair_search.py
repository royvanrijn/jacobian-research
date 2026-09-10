"""Pure-Python policy tests plus independent symbolic/finite reference checks."""
from fractions import Fraction as Q
import importlib.util
from math import gcd
from pathlib import Path
import sys
import unittest

CAS=Path(__file__).resolve().parents[1]/'cas'
sys.path.insert(0,str(CAS))
import elkies_pair_policy as p


class PolicyTests(unittest.TestCase):
    def test_rational_squares(self):
        self.assertEqual(p.square_root('81/100'),Q(9,10))
        self.assertIsNone(p.square_root('-1'))
        self.assertIsNone(p.square_root('2'))
        self.assertEqual(p.square_root(0),0)

    def test_squareclass_retains_scalar(self):
        self.assertEqual(p.same_quadratic_class([1,2,3],[4,8,12]),2)
        self.assertIsNone(p.same_quadratic_class([1,2,3],[2,4,6]))
        self.assertIsNone(p.same_quadratic_class([1,2,3],[-1,-2,-3]))
        self.assertIsNone(p.same_quadratic_class([1,2,3],[1,2,4]))
        self.assertEqual(p.same_quadratic_class([1,2,3,0],[1,2,3]),1)

    def test_same_cover_does_not_starve_pair_fallback(self):
        g=[[10,1,3],[1,10,1],[3,1,10]]
        rows=[{'mask':i,'word':[int(i==j) for j in range(3)],'q':q,'coefficient_bits':10}
              for i,q in enumerate(([1,0,1],[4,0,4],[2,1,1]))]
        order=p.pair_order(rows,g,1)
        self.assertEqual([r['kind'] for r in order],['same-cover','odd-intersection'])
        self.assertEqual(order[1]['intersection']%2,1)
        self.assertNotIn('rank_lower_bound',order[0])

    def test_even_intersection_is_not_odd_divisor_proof(self):
        rows=[{'mask':i,'word':[1-i,i],'q':q,'coefficient_bits':1}
              for i,q in enumerate(([1,0,1],[2,1,1]))]
        self.assertEqual(p.pair_order(rows,[[10,2],[2,10]],4),[])

    def test_invalid_lattice_rejected(self):
        with self.assertRaises(ValueError):
            p.pair_order([{'mask':1,'word':[1],'q':[1,0,1],'coefficient_bits':1}],[[4]],2)

    def test_controls_do_not_depend_on_scores(self):
        rows=[{'id':str(i),'shallow':{'score_units':i},'j_bits':10} for i in range(20)]
        a,_=p.split_controls(rows,4)
        for r in rows:r['shallow']['score_units']=-1000*int(r['id'])
        b,_=p.split_controls(list(reversed(rows)),4)
        self.assertEqual([r['id'] for r in a],[r['id'] for r in b])

    def test_controls_removed_before_shortlist(self):
        rows=[{'id':str(i),'shallow':{'score_units':i},'j_bits':10} for i in range(20)]
        controls,_=p.split_controls(rows,4)
        choices=p.shortlist(rows,controls,5)
        self.assertFalse({r['id'] for r in controls}&{r['id'] for r in choices})

    def test_deep_score_reorders_ranked_not_controls(self):
        controls=[{'id':'control','deep':{'score_units':-999},'j_bits':1}]
        rows=[{'id':str(i),'deep':{'score_units':-i},'j_bits':10} for i in range(4)]
        selected=p.final_selection(rows,controls,2)
        self.assertEqual([r['id'] for r in selected],['0','control','1'])
        self.assertEqual([r['control'] for r in selected],[False,True,False])

    def test_missing_traces_are_not_zeroes(self):
        a=p.mestre_score([(5,None),(7,0)])
        b=p.mestre_score([(5,0),(7,0)])
        self.assertEqual(a['missing_primes'],[5]);self.assertIsNone(a['rank_bound'])
        self.assertNotEqual(a['score_units'],b['score_units'])

    def test_bad_trace_rejected(self):
        with self.assertRaises(ValueError):p.mestre_score([(5,6)])
        with self.assertRaises(ValueError):p.mestre_score([(5,0),(5,1)])

    def test_twists_not_deduplicated_by_j_alone(self):
        self.assertTrue(p.short_isomorphic(-1,1,-16,64))
        self.assertFalse(p.short_isomorphic(-1,1,-4,8))
        self.assertFalse(p.short_isomorphic(-1,1,-1,-1))

    def test_special_j_isomorphisms(self):
        self.assertTrue(p.short_isomorphic(0,1,0,64))
        self.assertFalse(p.short_isomorphic(0,1,0,8))
        self.assertTrue(p.short_isomorphic(1,0,16,0))
        self.assertFalse(p.short_isomorphic(1,0,4,0))
        with self.assertRaises(ValueError):p.short_isomorphic(0,0,0,1)

    def test_height_and_roots(self):
        self.assertEqual(p.rational_height('3/1024'),11)
        self.assertEqual(p.nth_root(Q(81,16),4),Q(3,2))
        self.assertIsNone(p.nth_root(Q(4),4))
        self.assertEqual(p.nth_root(0,6),0)


@unittest.skipUnless(importlib.util.find_spec('sympy'),'independent symbolic checks require SymPy')
class FormulaTests(unittest.TestCase):
    def test_euclidean_trace_and_lift(self):
        import sympy as s
        h,m,b,k,q,z=s.symbols('h m b k q z')
        A=m*k-s.Rational(3,4)*b*b-h*h*q/4
        B=(m*m*q+b**3-2*m*b*k+h*h*(k*k-b*q))/4
        N=m*m-h*h*b;V=-m**3+s.Rational(3,2)*h*h*m*b-h**4*k/2
        self.assertEqual(s.expand(V*V-N**3-A*N*h**4-B*h**6),0)
        X=(b+h*z)/2;Y=-(h*k+m*z)/2
        self.assertEqual(s.rem(s.expand(Y*Y-X**3-A*X-B),z*z-q,z),0)

    def test_pointed_quartic_maps_are_birational(self):
        import sympy as s
        u,z,z0,c1,c2,c3,c4=s.symbols('u z z0 c1 c2 c3 c4')
        f=z0*z0+c1*u+c2*u*u+c3*u**3+c4*u**4
        x=(2*z0*(z+z0)+c1*u)/(u*u)
        y=((x*x-4*z0*z0*c4)*u-c1*x-2*z0*z0*c3)/(2*z0)
        a4=c1*c3-4*z0*z0*c4;a6=z0*z0*c3*c3+c1*c1*c4-4*z0*z0*c2*c4
        numerator=s.fraction(s.cancel(y*y-x**3-c2*x*x-a4*x-a6))[0]
        self.assertEqual(s.rem(numerator,z*z-f,z),0)
        ui=(2*z0*y+c1*x+2*z0*z0*c3)/(x*x-4*z0*z0*c4)
        zi=x*ui*ui/(2*z0)-z0-c1*ui/(2*z0)
        self.assertEqual(s.cancel(ui-u),0);self.assertEqual(s.cancel(zi-z),0)

    def test_branch_point_maps(self):
        import sympy as s
        u,z,c1,c2,c3,c4=s.symbols('u z c1 c2 c3 c4')
        f=c1*u+c2*u*u+c3*u**3+c4*u**4;x=c1/u;y=c1*z/u**2
        num=s.fraction(s.cancel(y*y-x**3-c2*x*x-c1*c3*x-c1*c1*c4))[0]
        self.assertEqual(s.rem(num,z*z-f,z),0)

    def test_odd_quartic_divisor_reduces_to_a_rational_point(self):
        import sympy as s
        u=s.symbols('u');g=u**3+u+1;v=2*u*u+u+1
        f=s.expand(v*v-g*(3*u+2))
        self.assertEqual(s.degree(f,u),4)
        # Leading coefficient is 1, a rational infinity case; change coefficients
        # to 3*u+2 -> 2*u+2, giving nonsquare leading coefficient 2 instead.
        f=s.expand(v*v-g*(2*u+2));q,r=s.div(v*v-f,g,u)
        self.assertEqual(r,0);root=-1
        self.assertEqual(f.subs(u,root),v.subs(u,root)**2)

    def test_elkies_reference_pair_has_a_nontorsion_deck_point(self):
        import sympy as s
        u=s.symbols('u');n=289444-u*u;d=130*u-38636
        f=s.Poly(s.expand(54756*n*n-3269604*n*d+22473889*d*d),u)
        coefficients=[Q(int(f.nth(4-i))) for i in range(5)]
        z0=Q(234);self.assertEqual(coefficients[0],z0*z0)
        _,c1,a2,c3,c4=coefficients
        a4=c1*c3-4*z0*z0*c4;a6=z0*z0*c3*c3+c1*c1*c4-4*z0*z0*a2*c4
        x=c1*c1/(4*z0*z0)-a2;y=-(c1*x+2*z0*z0*c3)/(2*z0)
        P=x,y;self.assertEqual(y*y,x**3+a2*x*x+a4*x+a6)
        def add(P,R):
            if P is None:return R
            if R is None:return P
            x,y=P;v,w=R
            if x==v and y==-w:return None
            slope=(3*x*x+2*a2*x+a4)/(2*y) if P==R else (w-y)/(v-x)
            xx=slope*slope-a2-x-v
            return xx,-y+slope*(x-xx)
        def multiply(n,P):
            out=None
            while n:
                if n%2:out=add(out,P)
                P=add(P,P);n//=2
            return out
        # Cubic discriminant for y^2=x^3+a2*x^2+a4*x+a6.
        disc=a2*a2*a4*a4-4*a4**3-4*a2**3*a6-27*a6*a6+18*a2*a4*a6
        orders=[];bound=0
        for prime in [5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71]:
            if any(a.denominator%prime==0 for a in (a2,a4,a6,disc)):continue
            mod=lambda a:a.numerator*pow(a.denominator,-1,prime)%prime
            if mod(disc)==0:continue
            a,b,c=map(mod,(a2,a4,a6))
            count=1
            for xx in range(prime):
                yy=(xx**3+a*xx*xx+b*xx+c)%prime
                count+=1 if yy==0 else (2 if pow(yy,(prime-1)//2,prime)==1 else 0)
            orders.append([prime,count]);bound=gcd(bound,count)
            if len(orders)==3:break
        self.assertGreaterEqual(len(orders),2)
        self.assertIsNotNone(multiply(bound,P))
        # This is X948's published reference pair, NOT an X1092 search result.


if __name__=='__main__':unittest.main()
