"""Rational recurrences for the gap-three Wronskian and its cyclic sections."""
from functools import lru_cache
import sympy as s

@lru_cache(maxsize=None)
def centered_coefficients(n,beta):
    beta=s.sympify(beta)
    if n==0:return (s.Integer(1),)
    hi=centered_coefficients(n-1,beta+1);lo=centered_coefficients(n-1,beta)
    a=4*n+beta-2
    at=lambda row,k:row[k] if 0<=k<len(row) else s.Integer(0)
    return tuple(s.expand(a*(a-1-k)*at(hi,k)+a*at(hi,k-1)+at(lo,k-2)) for k in range(2*n+1))

def centered(n,beta,z):
    return s.Poly(sum(c*z**k for k,c in enumerate(centered_coefficients(n,s.sympify(beta)))),z)

def wronskian(n,gamma,z):
    return centered(n,s.sympify(gamma)-n+1,z)

def next_fixed_shape(U,V,n,gamma,z):
    a=3*n+gamma+2
    return s.Poly(s.expand(a*((a-1+z)*U-z*s.diff(U,z))+z*z*V),z)

def cyclic_sections(W,z,u):
    return [s.Poly(sum(W.nth(k)*u**((k-r)//3) for k in range(r,W.degree()+1,3)),u) for r in range(3)]

def fifth_section_coefficients(g):
    A=5*(3*g*g+33*g+98)
    B=(g+10)*(51*g**4+1050*g**3+7945*g*g+26210*g+33264)
    C=5*(3*g**8+168*g**7+3962*g**6+51156*g**5+393407*g**4+1833972*g**3+5029028*g*g+7384944*g+4484480)
    return A,B,C

FIFTH_DISCRIMINANT_DESCENDING=(1701,98820,2542590,38257920,374339853,2503498860,11681645760,37903208880,82533243696,109328526720,66701465600)
