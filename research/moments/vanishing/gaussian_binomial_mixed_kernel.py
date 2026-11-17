"""Exact spatial normal forms for Phi_m[A]=E[A exp(tW(1+Z)^mY)].

The companion proof establishes all-order kernel completeness and
independence; this module implements its finite polynomial descent.
Operator-image sums are vector-space sums, not polynomial ideals.
"""
from __future__ import annotations
import sympy as sp

Z,W,Y,Q=sp.symbols("Z W Y q")


def validate_m(m):
    value=sp.sympify(m)
    if not value.is_Integer or value<1:
        raise ValueError("m must be a positive integer")
    return int(value)


def delta_one(amplitude):
    value=sp.sympify(amplitude)
    return sp.expand((Z*W-Y**2)*value-W*sp.diff(value,W)+Y*sp.diff(value,Y))


def delta_two(amplitude,m):
    m=validate_m(m)
    value=sp.sympify(amplitude)
    return sp.expand((W-(m-1)*Z*W+(m-1))*value
                     -(1+Z)*sp.diff(value,Z)+m*W*sp.diff(value,W))


def normal_form(amplitude,m):
    """Return N,B1,B2 with A=N+delta_one(B1)+delta_two(B2,m).

    N is the unique canonical representative: no ZW or Y^2 when m>=2;
    no W or Y^2 when m=1. Total degree never increases.
    """
    m=validate_m(m)
    value=sp.Poly(sp.expand(amplitude),Z,W,Y).as_expr()
    first=sp.Integer(0)
    second=sp.Integer(0)
    while True:
        bad=[(powers,coefficient) for powers,coefficient in sp.Poly(value,Z,W,Y).terms()
             if powers[2]>=2 or (powers[1]>=1 if m==1 else powers[0]>=1 and powers[1]>=1)]
        if not bad:
            return sp.expand(value),sp.expand(first),sp.expand(second)
        (a,b,j),coefficient=max(bad,key=lambda item:(sum(item[0]),item[0]))
        if (m==1 and b>=1) or (m>=2 and a>=1 and b>=1):
            if m==1:
                add_first=sp.Integer(0)
                add_second=coefficient*Z**a*W**(b-1)*Y**j
            else:
                add_first=sp.Integer(0)
                add_second=-coefficient*Z**(a-1)*W**(b-1)*Y**j/sp.Integer(m-1)
        else:
            quotient=coefficient*Z**a*W**b*Y**(j-2)
            add_first=-quotient
            add_second=Z*quotient if m==1 else -quotient/sp.Integer(m-1)
        value=sp.expand(value-delta_one(add_first)-delta_two(add_second,m))
        first+=add_first
        second+=add_second


def coefficient_classes(amplitude,m):
    """Return exact rational R(q),S(q) for the literal even/odd coefficients."""
    m=validate_m(m)
    normal,_,_=normal_form(amplitude,m)
    answer=[sp.Integer(0),sp.Integer(0)]
    for (a,b,j),coefficient in sp.Poly(normal,Z,W,Y).terms():
        n=2*Q+j
        if b:
            ratio=sp.ff((m-1)*n,b)
        else:
            ratio=sp.ff(n,a)/sp.rf((m-1)*n+1,a)
        answer[j]+=coefficient*ratio
    return tuple(sp.cancel(value) for value in answer)


def is_in_kernel(amplitude,m):
    """Exact original spatial-kernel decision, with no period truncation."""
    return normal_form(amplitude,m)[0]==0


def canonical_monomials(degree,m):
    m=validate_m(m)
    degree=int(degree)
    if degree<0:
        return []
    answer=[Z**a for a in range(degree+1)]
    if m>=2:
        answer.extend(W**b for b in range(1,degree+1))
    answer.extend(Y*Z**a for a in range(degree))
    if m>=2:
        answer.extend(Y*W**b for b in range(1,degree))
    return answer
