"""Exact full regular-nilpotent cyclic slice, normalized to ell_w=1.

F=w+ell_z*z+ell_t*t+z*t+arbitrary harmonic cubic. Every coefficient
shown is free. Global pure exhaustion is OPEN until exact identities
or a verified characteristic-zero computation establish it.
"""
from fractions import Fraction
from functools import lru_cache
import math

import sympy as S

from boundary_axial_quadratic import spatial_average
from boundary_rank1_quadratic import compositions

CHARGES = tuple(range(-3, 4))
LINEAR = S.symbols('ell_z ell_t')
CUBIC = S.symbols('gminus3 gminus2 gminus1 g0 g1 g2 g3')
VARIABLES = (*LINEAR, *CUBIC)
RODRIGUES = {v: (S.Rational(-1, 2) if v == CUBIC[-1] else S.Integer(0)) for v in VARIABLES}


@lru_cache(maxsize=None)
def exact_moment(order):
    result = {}
    factorials = [math.factorial(i) for i in range(3*order+2)]
    for cubic_count in range(order+1):
        for counts in compositions(cubic_count, 7):
            charge = sum(j*k for j, k in zip(CHARGES, counts))
            low3, low2, low1, middle, high1, high2, high3 = counts
            low_k = max(0, -charge)
            high_k = (order-cubic_count-charge)//2
            cubic_factorial = math.prod(factorials[k] for k in counts)
            for k in range(low_k, high_k+1):
                linear_t = order-cubic_count-charge-2*k
                linear_w = charge+k
                zbase = k+high1+2*high2+3*high3
                assert zbase == linear_w+3*low3+2*low2+low1
                for quadratic_count in range(order % 2, k+1, 2):
                    linear_z = k-quadratic_count
                    tbase = quadratic_count+linear_t+low2+middle+high2
                    assert tbase % 2 == 0
                    average = spatial_average(zbase, tbase, 0, low1+high1, middle)
                    if not average:
                        continue
                    multinomial = factorials[order]//(
                        factorials[quadratic_count]*factorials[linear_z]*factorials[linear_t]
                        *factorials[linear_w]*cubic_factorial)
                    powers = (linear_z, linear_t, *counts)
                    assert powers not in result
                    result[powers] = multinomial*average
    return result


def raw_expression(order):
    return S.expand(sum(S.Rational(c.numerator, c.denominator)
                        *math.prod(v**k for v, k in zip(VARIABLES, powers))
                        for powers, c in exact_moment(order).items()))


def primitive_expression(order):
    return S.Poly(raw_expression(order), *VARIABLES).clear_denoms()[1].primitive()[1].as_expr()


def rodrigues_linearization(order):
    """Differentiate the sparse moment at Rodrigues without expanding it."""
    result = [Fraction(0)]*len(VARIABLES)
    for powers, value in exact_moment(order).items():
        for index, exponent in enumerate(powers):
            if not exponent:
                continue
            other = list(powers)
            other[index] -= 1
            if any(other[:-1]):
                continue
            result[index] += value*exponent*Fraction(-1, 2)**other[-1]
    return [S.Rational(x.numerator, x.denominator) for x in result]
