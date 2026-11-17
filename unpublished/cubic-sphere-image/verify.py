#!/usr/bin/env python3
"""Exact identities for the interval-cubic and constant-profile propositions.

Requires SymPy 1.14.0. No certificate file, search, floating point arithmetic,
Singular, or generated Groebner basis is used. The four-variable moment
polynomials are derived from the displayed multinomial sphere formula.
"""
from functools import lru_cache
from math import factorial,gcd,lcm
import sympy as S
from sympy.polys.domains import QQ
from sympy.polys.rings import ring

if not __debug__:
    raise RuntimeError('Run without -O: this verifier requires assertions')

assert S.__version__=='1.14.0'


def average(poly,t):
    return S.expand(sum(c/S.Integer(k+1) for (k,),c in S.Poly(S.expand(poly),t).terms() if k%2==0))


def interval_cubic():
    t,h,v,c,z=S.symbols('t h v c z')
    u=h*t**3+v*(t*t-S.Rational(1,3))+c*t
    M2,M3,M4=[average(u**m,t) for m in (2,3,4)]
    assert S.expand(M2-(105*c*c+126*c*h+45*h*h+28*v*v)/315)==0
    assert S.expand(S.rem(M3,M2,v)-8*v*(56*c*c+105*c*h+45*h*h)/2205)==0
    R=28028*c**4+84084*c**3*h+93275*c*c*h*h+46410*c*h**3+8955*h**4
    assert S.expand(S.rem(M4,M2,v)+16*R/735735)==0
    assert S.expand(S.rem(M3.subs(h,0),M2.subs(h,0),c)+S.Rational(256,4725)*v**3)==0
    R=R.subs({h:1,c:z});P=56*z*z+105*z+45;B=105*z*z+126*z+45
    assert S.rem(R,P,z)==-S.Rational(315,128)*(299*z+159)
    assert S.rem(R,B,z)==-S.Rational(84,125)*(117*z-110)
    assert P.subs(z,-S.Rational(159,299))==S.Rational(446976,89401)
    assert B.subs(z,S.Rational(110,117))==S.Rational(1169375,4563)


@lru_cache(None)
def df(n):
    return 1 if n<=0 else n*df(n-2)


def det_mod(matrix,prime):
    matrix=[[int(c)%prime for c in row] for row in matrix]
    determinant=1
    for column in range(len(matrix)):
        pivot=next((i for i in range(column,len(matrix)) if matrix[i][column]),None)
        if pivot is None:return 0
        if pivot!=column:
            matrix[column],matrix[pivot]=matrix[pivot],matrix[column]
            determinant=-determinant
        entry=matrix[column][column];determinant=determinant*entry%prime
        for row in range(column+1,len(matrix)):
            factor=matrix[row][column]*pow(entry,-1,prime)%prime
            matrix[row]=[(a-factor*b)%prime for a,b in zip(matrix[row],matrix[column])]
    return determinant%prime


def constant_profile():
    x,y,t,a,c,k=S.symbols('x y t a c k')
    v0,v1,d,e=S.symbols('v0 v1 d e');rr=S.symbols('r0:4')
    original=a*y+(v0+v1*x)*t*t+(c+d*x+e*x*x)*t+sum(rr[j]*x**j for j in range(4))
    rotation={x:x,t:t+k*x,y:y-k*t-k*k*x/2}
    assert S.expand((t*t+2*x*y).subs(rotation,simultaneous=True)-(t*t+2*x*y))==0
    transformed=S.Poly(S.expand(original.subs(rotation,simultaneous=True)),x,y,t)
    assert transformed.coeff_monomial(y)==a
    assert transformed.coeff_monomial(t)==c-a*k
    assert transformed.coeff_monomial(t*t)==v0 and transformed.coeff_monomial(x*t*t)==v1
    assert all(sum(powers)<=3 and powers[1]<=1 for powers,coefficient in transformed.terms())
    R,v0,v1,d,e,r0,r1,r2,r3=ring('v0,v1,d,e,r0,r1,r2,r3',QQ)
    variables=(v0,v1,d,e,r0,r1,r2,r3)
    weights=(1,2,2,3,1,2,3,4)
    phase=(0,1,1,2,0,1,2,3)
    tpower=(2,2,1,1,0,0,0,0)
    def vectors(n,index=0,prefix=()):
        if index==8:
            if n==0:yield prefix
            return
        for count in range(n//weights[index]+1):
            yield from vectors(n-count*weights[index],index+1,prefix+(count,))
    def moment(n):
        out={}
        for powers in vectors(n):
            k=sum(a*b for a,b in zip(powers,phase));j=sum(a*b for a,b in zip(powers,tpower))
            if j%2:continue
            denominator=df(2*k+j+1)
            for count in powers:denominator*=factorial(count)
            out[powers]=QQ(factorial(n)*df(j-1),denominator)
        return R.from_dict(out)
    assert moment(0)==1
    assert moment(1)==r0+v0/3
    assert moment(2)==r0*r0+2*r0*v0/3+v0*v0/5+2*r1/3+2*v1/15
    rr=(r0,r1,r2,r3)
    solutions=(-v0/3,-(2*v0*v0+3*v1)/15,-4*v0*(31*v0*v0+54*v1)/945,-(14175*d*d+4624*v0**4+10152*v0*v0*v1+3024*v1*v1)/28350)
    coefficients=(QQ(1),QQ(2,3),QQ(2,5),QQ(8,35))
    for n in range(1,5):
        value=moment(n)
        for variable,solution in zip(rr[:n-1],solutions[:n-1]):value=value.compose(variable,solution)
        assert value==coefficients[n-1]*(rr[n-1]-solutions[n-1])
    reduced=[];N=[]
    for n in range(5,9):
        value=moment(n)
        for variable,solution in zip(rr,solutions):value=value.compose(variable,solution)
        reduced.append(value)
        denominator=lcm(*(int(c.denominator) for c in value.values()))
        integers={powers:int(c*denominator) for powers,c in value.items()}
        content=gcd(*(abs(c) for c in integers.values()))
        if integers[max(integers)]<0:content=-content
        N.append(R.from_dict({powers:QQ(c//content) for powers,c in integers.items()}))
    assert reduced[3].compose(v0,R.zero)==-QQ(65536,41023125)*v1**4
    assert reduced[1].compose(v0,R.zero).compose(v1,R.zero)==QQ(8,231)*e*e
    BR,b,D,E,F=ring('b,D,E,F',QQ)
    rows=[]
    for n,polynomial in enumerate(N,5):
        row=BR.zero
        for powers,c in polynomial.items():
            i,j,h,k=powers[:4]
            assert not any(powers[4:]) and i+2*j+2*h+3*k==n
            row+=c*b**j*{(0,0):BR.one,(2,0):D,(1,1):E,(0,2):F}[(h,k)]
        rows.append(S.expand(row.as_expr()))
    bs,Ds,Es,Fs=S.symbols('b D E F')
    matrix=S.Matrix([[row.coeff(x) for x in (Ds,Es,Fs)]+[row.subs({Ds:0,Es:0,Fs:0})] for row in rows])
    expected=S.Matrix([
      [155925,155925,0,35504+94992*bs+52704*bs**2],
      [77395500-116093250*bs,77395500,-58046625,30907712+63343584*bs+3500496*bs**2-12573792*bs**3],
      [226800-297675*bs,226800+127575*bs,-212625,53504+17328*bs-211968*bs**2-108864*bs**3],
      [4819500-216877500*bs,4819500-65063250*bs,-75907125,7986496-27544704*bs-90351504*bs**2-16632864*bs**3+9144576*bs**4]])
    assert matrix==expected
    assert matrix[:,0]-matrix[:,1]-2*bs*matrix[:,2]==S.zeros(4,1)
    f=2485161*bs**3+6647562*bs**2+4908168*bs+1097984
    g=10216206*bs**4+24356619*bs**3-25372278*bs**2-41158800*bs-11711872
    assert S.Poly(matrix[[0,1,2],[0,2,3]].det(),bs).primitive()[1].as_expr()==f
    assert S.Poly(matrix[[0,1,3],[0,2,3]].det(),bs).primitive()[1].as_expr()==-g
    fc=list(map(int,S.Poly(f,bs).all_coeffs()));gc=list(map(int,S.Poly(g,bs).all_coeffs()))
    sylvester=[]
    for shift in range(4):sylvester.append([0]*shift+fc+[0]*(3-shift))
    for shift in range(3):sylvester.append([0]*shift+gc+[0]*(2-shift))
    assert det_mod(sylvester,19)==2
    corrupted=[row[:] for row in sylvester];corrupted[0]=[0]*7
    assert det_mod(corrupted,19)==0


if __name__=='__main__':
    interval_cubic();constant_profile()
    print('PASS interval cubics: four moments force zero')
    print('PASS constant-profile family with no t^3 term: all eight-moment elimination identities and the mod-19 Sylvester determinant')
    print('PASS corrupted determinant rejected; no external polynomial certificate is needed')
