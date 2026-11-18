"""Uniform odd-shear polynomial Keller maps over Q(i).

For odd m>=3, p=(m-1)//2, coefficients c=m(m-i)/(m^2-1),
 d=(1-c)^(-m). Compact A,B in Lambda,w have a polynomial lift
 F=(A(1+xy,xz-1)/x, B(1+xy,xz-1)/x, x(1+xy)).
Polynomiality, determinant and monodromy are proved in the companion note.
"""
from __future__ import annotations
from math import comb
import sympy as s


def binomial_terms(m):
    if m < 3 or m % 2 != 1:
        raise ValueError('m must be odd and at least three')
    p=(m-1)//2
    # (rational coefficient, c exponent, d exponent, Lambda exponent,
    #  C1 exponent, C2 exponent). Denominators are always nonzero.
    terms=[(s.Rational(1),0,0,0,2,0)]
    terms += [(s.Rational(-2*comb(m,j)*(-1)**j,(p+1)*j-p),1,j,
               p*(j+1),1+m*j,m-j)for j in range(m+1)]
    terms += [(s.Rational(comb(2*m,j)*(-1)**j,(p+1)*j-m),2,j,
               p*(j+2),m*j,2*m-j)for j in range(2*m+1)]
    return terms


def coefficient_pair(m):
    if m < 3 or m % 2 != 1:
        raise ValueError('m must be odd and at least three')
    c=s.Rational(m*m,m*m-1)-s.I*s.Rational(m,m*m-1)
    field=s.QQ_I
    d=field.to_sympy((field.one-field.convert(c))**(-m))
    return c,d


def expanded_pair(m=3):
    L,w=s.symbols('Lambda w');p=(m-1)//2;c,d=coefficient_pair(m)
    C1=s.Poly(1+c*L**p*w**m,L,w,domain=s.QQ_I)
    C2=s.Poly(w,L,w,domain=s.QQ_I)+s.Poly(d*L**p,L,w,domain=s.QQ_I)*C1**m
    terms=binomial_terms(m)
    C1powers=[s.Poly(1,L,w,domain=s.QQ_I)]
    C2powers=[C1powers[0]]
    for _ in range(max(row[4]for row in terms)):C1powers.append(C1powers[-1]*C1)
    for _ in range(max(row[5]for row in terms)):C2powers.append(C2powers[-1]*C2)
    B=s.Poly(0,L,w,domain=s.QQ_I)
    for a,b,e,f,g,h in terms:
        B += s.Poly(a*c**b*d**e*L**f,L,w,domain=s.QQ_I)*C1powers[g]*C2powers[h]
    return {'m':m,'p':p,'c':c,'d':d,'Lambda':L,'w':w,'A':C1*C2,'B':B}


def expanded_map(m=3):
    pair=expanded_pair(m);A,B=pair['A'],pair['B']
    x,y,z=s.symbols('x y z');field=s.QQ_I
    Lsub=s.Poly(1+x*y,x,y,z,domain=field)
    wsub=s.Poly(x*z-1,x,y,z,domain=field)
    Lpowers=[s.Poly(1,x,y,z,domain=field)];wpowers=[Lpowers[0]]
    for _ in range(max(i for(i,j),c in B.terms())):Lpowers.append(Lpowers[-1]*Lsub)
    for _ in range(max(j for(i,j),c in B.terms())):wpowers.append(wpowers[-1]*wsub)
    xpoly=s.Poly(x,x,y,z,domain=field)
    def lift(poly):
        out=s.Poly(0,x,y,z,domain=field)
        for(i,j),co in poly.terms():out+=Lpowers[i]*wpowers[j]*co
        return out.exquo(xpoly)
    return (x,y,z),[lift(A),lift(B),xpoly*Lsub]


def laurent_primitive(m):
    if m < 3 or m % 2 != 1:
        raise ValueError('m must be odd and at least three')
    p=(m-1)//2;X,Y,T,c,d=s.symbols('X Y T c d')
    P=X
    for j in range(m+1):
        e=(p+1)*j-p
        P +=s.Rational(-2*comb(m,j)*(-1)**j,e)*c*d**j*Y**(m-j)*X**e
    for j in range(2*m+1):
        e=(p+1)*j-m
        P +=s.Rational(comb(2*m,j)*(-1)**j,e)*c*c*d**j*Y**(2*m-j)*X**e
    Phi=X**(p+1)-c*(Y-d*X**(p+1))**m
    return {'m':m,'p':p,'X':X,'Y':Y,'T':T,'c':c,'d':d,'P':P,'Phi':Phi}