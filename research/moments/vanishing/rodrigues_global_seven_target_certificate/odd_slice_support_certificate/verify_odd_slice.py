"""Exact first-six-moment support closure on the even-zero, s1-zero slice.

Uses rational sphere moments and a small mod-101 coprimality certificate.
No msolve or global multivariate Groebner computation is used.
"""
from fractions import Fraction
from functools import lru_cache
from math import factorial
from pathlib import Path
import json
import time
import sympy as s

started=time.monotonic()
z,t,w,R,A,B,C,x,T=s.symbols('z t w R A B C x T')
F0=w+z*t-z**3/2
F=F0+R*w*w+A*t*w*w+B*(t**3-3*z*t*w)+C*z*z*t

# Full cyclic-frame adapter, including the oriented determinant and all
# symmetric tensor slots. Here Bmat is the quadratic operator; B is a scalar.
metric=s.Matrix([[0,0,1],[0,1,0],[1,0,0]])
v=s.Matrix([1,0,0])
Bmat=s.Matrix([[0,0,R],[s.Rational(1,2),0,0],[0,s.Rational(1,2),0]])
frame=s.Matrix.hstack(v,Bmat*v,Bmat*Bmat*v)
assert frame.det()==s.Rational(1,8)
assert Bmat.T*metric==metric*Bmat
assert s.trace(Bmat)==0 and s.trace(Bmat*Bmat)==0
assert Bmat.det()==R/4
G=frame.T*metric*frame
assert G==s.Matrix([[0,0,s.Rational(1,4)],[0,s.Rational(1,4),0],[s.Rational(1,4),0,0]])
coords=s.Matrix([z,t,w])
L=list(frame.T*metric*coords)
U={(0,0,1):2*A/3,(0,1,2):-4*B,(1,1,1):8*B,(1,2,2):32*C/3,(2,2,2):-32}
tensor=lambda i,j,k:U.get(tuple(sorted((i,j,k))),s.S.Zero)
for k in range(3):
    assert s.expand(sum(G[i,j]*tensor(i,j,k) for i in range(3) for j in range(3)))==0
H3=sum(tensor(i,j,k)*L[i]*L[j]*L[k]
       for i in range(3) for j in range(3) for k in range(3))
adapted=(v.T*metric*coords)[0]+(coords.T*metric*Bmat*coords)[0]+H3
assert s.expand(adapted-F)==0
assert s.diag(4,s.Rational(3,2),-s.Rational(1,4),s.Rational(3,32)).det()==-s.Rational(9,64)
print('Oriented cyclic adapter, all ten cubic slots and trace equations: PASS',flush=True)

@lru_cache(None)
def df(n):
    return 1 if n<2 else n*df(n-2)

def weak_tuples(total,count):
    if count==0:
        yield ()
    else:
        for k in range(total+1):
            for rest in weak_tuples(total-k,count-1):
                yield (k,)+rest

def exact_moment(n):
    """Multinomial enumeration; torus balance solves the first two counts."""
    out={}
    for k2,k3,k4,k5,k6,k7 in weak_tuples(n,6):
        left=n-k2-k3-k4-k5-k6-k7
        twice=left-3*k2+2*k3+2*k4-2*k7
        if twice<0 or twice%2:continue
        k1=twice//2;k0=left-k1
        if k0<0:continue
        degree_t=k1+k4+3*k5+k6+k7
        if degree_t%2:continue
        degree_z=k1+3*k2+k6+2*k7
        numerator=factorial(n)*(-1)**k2*(-3)**k6*factorial(degree_z)*df(degree_t-1)
        denominator=2**k2*df(2*degree_z+degree_t+1)
        for k in (k0,k1,k2,k3,k4,k5,k6,k7):denominator*=factorial(k)
        exponent=(k3,k4,k5+k6,k7)
        out[exponent]=out.get(exponent,Fraction(0))+Fraction(numerator,denominator)
    return s.Poly.from_dict({e:s.Rational(c.numerator,c.denominator)
                            for e,c in out.items() if c},(R,A,B,C))

def direct_moment(n):
    """Independent direct spatial expansion, followed by the sphere rule."""
    out=s.S.Zero
    for (iz,it,iw),coefficient in s.Poly(s.expand(F**n),z,t,w).terms():
        if iz==iw and it%2==0:
            out+=coefficient*s.Rational(factorial(iz)*df(it-1),df(2*iz+it+1))
    return s.Poly(s.expand(out),R,A,B,C)

raw={n:exact_moment(n) for n in range(2,7)}
for n,p in raw.items():
    assert p==direct_moment(n)
moments={n:p.clear_denoms()[1].primitive()[1] for n,p in raw.items()}
assert moments[2].as_expr()==4*A*C+15*B*B
assert moments[3].as_expr()==-4*B*C*R-9*B-6*R
assert moments[5].as_expr().subs({B:0,C:0,R:0})==-78*A
assert moments[5].as_expr().subs({B:0,A:0,R:0})==286*C
print('Moments 2 through 6 agree with independent direct expansion: PASS',flush=True)

# First-nine scheme closure: first-six support is origin-only, and these four
# odd moments have independent linear parts there. Independently check those
# derivatives by expanding only F0^(n-1) times a fixed spatial marker.
odd_rows=[]
markers=[w*w,t*w*w,t**3-3*z*t*w,z*z*t]
for n in (3,5,7,9):
    p=exact_moment(n)
    for variable,marker in zip((R,A,B,C),markers):
        direct=s.S.Zero
        for (iz,it,iw),coefficient in s.Poly(s.expand(n*F0**(n-1)*marker),z,t,w).terms():
            if iz==iw and it%2==0:
                direct+=coefficient*s.Rational(factorial(iz)*df(it-1),df(2*iz+it+1))
        assert s.factor(direct-p.coeff_monomial(variable))==0
    primitive=p.clear_denoms()[1].primitive()[1]
    odd_rows.append([primitive.coeff_monomial(variable) for variable in (R,A,B,C)])
odd_jacobian=s.Matrix(odd_rows)
assert odd_jacobian==s.Matrix([[-6,0,-9,0],[0,-78,0,286],[646,0,6460,0],[0,5244,0,0]])
assert odd_jacobian.det()==49411883664
print('Odd moments 3,5,7,9 have Jacobian determinant 49411883664: PASS',flush=True)

# All possible nonzero points have B*C*(2*B*C+3) != 0.
sub={A:-15*B**3/(4*x),R:-9*B/(2*(2*x+3)),C:x/B}
reduced={}
for n in (4,5,6):
    numerator,denominator=s.fraction(s.factor(moments[n].as_expr().subs(sub,simultaneous=True)))
    p=s.Poly(numerator,B)
    assert all(exponents[0]%4==0 for exponents in p.monoms())
    reduced[n]=s.expand(sum(c*T**(e[0]//4) for e,c in p.terms()))

h=128*x**3+1684*x**2+6528*x+6435
g=176*x**2+843*x+990
assert s.expand(reduced[4]+2*(3*h*T+26*x*x*g))==0
T_value=-26*x*x*g/(3*h)
univariate=[]
for n in (5,6):
    numerator,denominator=s.fraction(s.factor(reduced[n].subs(T,T_value)))
    p=s.Poly(numerator,x).primitive()[1]
    while p.eval(0)==0:p=p.exquo(s.Poly(x,x))
    if p.LC()<0:p=-p
    univariate.append(p)
assert [p.degree() for p in univariate]==[8,10]

# Full-degree mod-101 coprimality implies coprimality over Q, hence C.
# Export concrete Bezout polynomials, verified by plain polynomial identity.
bezout={}
for name,first,second in [('h_g',s.Poly(h,x),s.Poly(g,x)),
                          ('P8_Q10',univariate[0],univariate[1])]:
    first_mod=s.Poly(first,x,modulus=101)
    second_mod=s.Poly(second,x,modulus=101)
    assert first_mod.degree()==first.degree()
    assert second_mod.degree()==second.degree()
    left,right,gcd=s.gcdex(first_mod,second_mod)
    assert gcd==s.Poly(1,x,modulus=101)
    assert left*first_mod+right*second_mod==gcd
    assert s.gcd(first,second).degree()==0
    left_Q,right_Q,gcd_Q=s.gcdex(first,second)
    assert gcd_Q==s.Poly(1,x,domain=s.QQ)
    assert left_Q*first+right_Q*second==gcd_Q
    bezout[name]={'modulus':101,'first':str(first.as_expr()),'second':str(second.as_expr()),
                 'left':str(left.as_expr()),'right':str(right.as_expr()),
                 'left_over_Q':str(left_Q.as_expr()),'right_over_Q':str(right_Q.as_expr())}
    print(name,'Bezout identity mod101:',left.as_expr(),';',right.as_expr(),flush=True)

data={'phase':str(F),'primitive_moments':{str(n):str(p.as_expr()) for n,p in moments.items()},
      'substitution':{str(k):str(v) for k,v in sub.items()},
      'reduced_equations':{str(n):str(p) for n,p in reduced.items()},
      'T_value':str(T_value),'bezout':bezout,
      'odd_jacobian':[[int(v) for v in row] for row in odd_jacobian.tolist()],
      'odd_jacobian_determinant':int(odd_jacobian.det()),
      'conclusion':'The first six moments have only R=A=B=C=0 as a common complex zero in this family.'}
target=Path(__file__).with_name('odd_slice_support_data.json')
target.write_text(json.dumps(data,indent=2)+'\n')
print('ALL ASSERTIONS PASSED: first-six support is the Rodrigues origin.',flush=True)
print('Elapsed seconds:',time.monotonic()-started,flush=True)
