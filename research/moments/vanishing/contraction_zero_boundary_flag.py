"""Exact inputs on the complete top-two-jet flag under contraction zero.

The simple-root chart is not classified by this generator. The other
chart concerns a multiple quadratic root with a nonzero central part.
Moment formulas use the uniform t integral and a circular constant term.
"""
from functools import lru_cache
import math

import sympy as S

a, b, c, k, ell, d, g = S.symbols('a b c k ell d g')
f = S.Symbol('f')
T = S.Symbol('T')


def average_t(polynomial):
    return S.expand(sum(coefficient/S.Integer(power[0]+1)
                        for power, coefficient in S.Poly(polynomial, T).terms()
                        if power[0] % 2 == 0))


@lru_cache(maxsize=None)
def circular_moment(n, neg, central, pos, pos2, pos3):
    """h((w*neg(t)+central(t)+z*pos(t)+z^2*pos2(t)+z^3*pos3)^n).

    q=2zw+t^2=1. The uniform t expectation is (1/2) integral[-1,1].
    No Gaussian radial transfer is used.
    """
    result = S.Integer(0)
    product = (1-T*T)/2
    for p3 in range(n//4+1):
        for p2 in range((n-4*p3)//3+1):
            for p1 in range((n-4*p3-3*p2)//2+1):
                nzero = n-4*p3-3*p2-2*p1
                nneg = p1+2*p2+3*p3
                scalar = math.factorial(n)//math.prod(math.factorial(j)
                                  for j in (nneg,nzero,p1,p2,p3))
                polynomial = (neg*product)**nneg*central**nzero*pos**p1*pos2**p2
                result += scalar*pos3**p3*average_t(polynomial)
    return S.expand(result)


@lru_cache(maxsize=None)
def simple_moment(n):
    # H2=w*t+k*z*t+ell*z^2;
    # H3=d*(hb3[-1]-k*hb3[1]+ell*hb3[2]/2)+g*z^3.
    neg = c-d/4+T+5*d*T*T/4
    pos = a+k*d/4+k*T-5*k*d*T*T/4
    return circular_moment(n,neg,b*T,pos,ell*(1+d*T/2),g)


@lru_cache(maxsize=None)
def double_raw_moment(n):
    # H2=t^2-z*w+k*z^2; H3=f*z^2*t+g*z^3.
    return circular_moment(n,c,b*T+(3*T*T-1)/2,a,k+f*T,g)


@lru_cache(maxsize=None)
def double_reduced_moment(n):
    # c1 is legitimate only after the first3 equations prove c!=0.
    # Those equations then solve a,k without any extra localization.
    return circular_moment(n,S.Integer(1),b*T+(3*T*T-1)/2,
                           -b*b/2-S.Rational(3,10),
                           -3*b*b/2-S.Rational(31,70)+f*T,g)


def primitive(expression, variables):
    return S.Poly(expression,*variables,domain=S.QQ).clear_denoms()[1].primitive()[1].as_expr()
