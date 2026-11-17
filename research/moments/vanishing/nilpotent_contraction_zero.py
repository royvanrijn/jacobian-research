"""Complete regular-nilpotent family after Delta^2(H2*H3)=0.

H2=z*t. The covariant condition kills the cubic charges -2,-1,0,
without restricting any linear coefficient. The recoverable linear-w
power is restored from the existing unrestricted moment generator.
"""
from collections import defaultdict
from functools import lru_cache
from fractions import Fraction
import math

import sympy as S

from nilpotent_cyclic_slice import CHARGES, exact_moment as cyclic_moment

A, B, C, D, E, F, G = S.symbols('a b c d e f g')
VARIABLES = (A, B, C, D, E, F, G)


@lru_cache(maxsize=None)
def exact_moment(order):
    result = {}
    for powers, coefficient in cyclic_moment(order).items():
        if powers[3] or powers[4] or powers[5]:
            continue
        cubic_count = sum(powers[2:])
        charge = sum(j*k for j, k in zip(CHARGES, powers[2:]))
        numerator = order-cubic_count-2*powers[0]-powers[1]-charge
        assert numerator % 2 == 0
        quadratic_count = numerator//2
        linear_w = charge+quadratic_count+powers[0]
        assert quadratic_count >= 0 and linear_w >= 0
        target = (powers[0], powers[1], linear_w, powers[2],
                  powers[6], powers[7], powers[8])
        assert target not in result
        result[target] = coefficient
    return result


def raw_expression(order):
    return S.expand(sum(S.Rational(c.numerator, c.denominator)
                        *math.prod(v**e for v, e in zip(VARIABLES, pp))
                        for pp, c in exact_moment(order).items()))


def primitive(expression, variables):
    return S.Poly(expression, *variables, domain=S.QQ).clear_denoms()[1].primitive()[1].as_expr()


def d_nonzero(order):
    value = raw_expression(order).subs({D: 1, F: -S.Rational(21, 2)*B*C,
                                      G: -S.Rational(35, 12)*(2*A*C+B*B)}, simultaneous=True)
    return primitive(S.expand(value), (A, B, C, E))


def d_zero_c_nonzero(order):
    value = raw_expression(order).subs({A: 0, B: 0, C: 1, D: 0}, simultaneous=True)
    return primitive(S.expand(value), (E, F, G))
