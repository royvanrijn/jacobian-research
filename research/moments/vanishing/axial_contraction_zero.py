"""Complete semisimple axial quadratic family under Delta^2(H2*H3)=0.

All three linear coefficients are arbitrary. For H2=t^2-z*w the
contraction condition kills exactly cubic charges -1,0,1, leaving the
four displayed cubic coefficients. No exclusion is asserted by this
input generator alone.
"""
from fractions import Fraction
from functools import lru_cache
import math

import sympy as S

from boundary_rank1_quadratic import compositions
from boundary_axial_quadratic import spatial_average

A, B, C, D, E, F, G = S.symbols('a b c d e f g')
VARIABLES = (A, B, C, D, E, F, G)
CHARGES = (1, 0, -1, -3, -2, 2, 3)
U = S.Symbol('u')


@lru_cache(maxsize=None)
def exact_moment(order):
    result = {}
    for total in range(order+1):
        quadratic_count = order-total
        for powers in compositions(total, 7):
            if sum(k*j for k, j in zip(powers, CHARGES)):
                continue
            a, b, c, d, e, f, g = powers
            if (b+e+f) % 2:
                continue
            zbase, wbase = a+2*f+3*g, c+3*d+2*e
            assert zbase == wbase
            average = spatial_average(zbase, b+e+f, quadratic_count, 0, 0)
            if average:
                result[powers] = math.factorial(order)//(
                    math.factorial(quadratic_count)*math.prod(math.factorial(k) for k in powers))*average
    return result


def raw_expression(order):
    return S.expand(sum(S.Rational(c.numerator, c.denominator)
                        *math.prod(v**e for v, e in zip(VARIABLES, pp))
                        for pp, c in exact_moment(order).items()))


def primitive(expression, variables):
    return S.Poly(expression, *variables, domain=S.QQ).clear_denoms()[1].primitive()[1].as_expr()


def amplitude_moment(order, central_linear):
    """mu_n/lambda^(n mod2), with c1 and b0 or1, u=lambda^2.

    Dividing odd moments uses the declared nonzero axial quadratic
    amplitude. No conclusion is made at amplitude zero.
    """
    assert central_linear in (0, 1)
    result = S.Integer(0)
    for pp, coefficient in exact_moment(order).items():
        if central_linear == 0 and pp[1]:
            continue
        quadratic_count = order-sum(pp)
        assert quadratic_count >= order % 2
        assert (quadratic_count-order % 2) % 2 == 0
        result += (S.Rational(coefficient.numerator, coefficient.denominator)
                   *A**pp[0]*D**pp[3]*E**pp[4]*F**pp[5]*G**pp[6]
                   *U**((quadratic_count-order % 2)//2))
    return S.expand(result)


def linear_normalized_expression(name, order):
    if name == 'mixed_linear':
        value = amplitude_moment(order, 1).subs(
            {A: (112-54*D*G-4*E*F)/217,
             U: (8*D*G-4*E*F-105)/31}, simultaneous=True)
        return primitive(S.expand(value), (D, E, F, G))
    if name == 'transverse_nonisotropic':
        value = amplitude_moment(order, 0).subs(
            {A: 1, U: 7+2*D*G}, simultaneous=True)
        return primitive(S.expand(value), (D, E, F, G))
    if name == 'transverse_isotropic':
        value = raw_expression(order).subs(
            {A: 0, B: 0, C: 1, F: -S.Rational(27, 4)/E,
             G: S.Rational(1, 2)/D}, simultaneous=True)
        numerator = S.fraction(S.cancel(value))[0]
        poly = S.Poly(numerator, D, E)
        # Both d,e are units by the already imposed product equations.
        min_d = min(pp[0] for pp, c in poly.terms())
        min_e = min(pp[1] for pp, c in poly.terms())
        return primitive(S.cancel(numerator/D**min_d/E**min_e), (D, E))
    raise ValueError(name)


def central_linear_expression(order):
    selected = S.Poly(raw_expression(order).subs({A: 0, C: 0}), B, D, E, F, G)
    value = S.Integer(0)
    for (b, d, e, f, g), coefficient in selected.terms():
        assert d == g and e == f  # valid for the two used moments4,5
        value += (coefficient*B**b*((1+7*B*B)/2)**d
                  *(-(27+77*B*B)/4)**e)
    return primitive(S.expand(value), (B,))
