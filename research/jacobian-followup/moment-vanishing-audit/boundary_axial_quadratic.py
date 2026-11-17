"""Exact input for H1=0, H2=t^2-z*w, with unrestricted harmonic H3.

This is the entire semisimple double-eigenvalue quadratic stratum, up
to scalar and orthogonal normalization. It does not include a nontrivial
Jordan block at the repeated eigenvalue. No purity classification is
asserted by this input generator.
"""
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
import math

import sympy as S

from boundary_rank1_quadratic import compositions, odd_factorial

CHARGES = tuple(range(-3, 4))
COEFFICIENTS = S.symbols('gminus3 gminus2 gminus1 g0 g1 g2 g3')


@lru_cache(maxsize=None)
def spatial_average(zbase, tbase, quadratic_count, height_power, middle):
    # Each factor is t^2-c*z*w. Convolve once, then integrate monomials.
    polynomial = [Fraction(1)]
    for power, scalar in [(quadratic_count, 1), (height_power, Fraction(1, 2)), (middle, 3)]:
        factor = [Fraction(math.comb(power, k))*(-scalar)**k for k in range(power+1)]
        product = [Fraction(0)]*(len(polynomial)+power)
        for i, a in enumerate(polynomial):
            for j, b in enumerate(factor):
                product[i+j] += a*b
        polynomial = product
    total = quadratic_count+height_power+middle
    denominator = odd_factorial(2*zbase+tbase+2*total+1)
    return sum((coefficient*math.factorial(zbase+k)
                *odd_factorial(tbase+2*(total-k)-1)/denominator
                for k, coefficient in enumerate(polynomial)), Fraction(0))


@lru_cache(maxsize=None)
def exact_moment(order):
    result = defaultdict(Fraction)
    for cubic_count in range(0, order+1, 2):
        quadratic_count = order-cubic_count
        for counts in compositions(cubic_count, 7):
            if sum(j*k for j, k in zip(CHARGES, counts)):
                continue
            low3, low2, low1, middle, high1, high2, high3 = counts
            zbase = high1+2*high2+3*high3
            assert zbase == 3*low3+2*low2+low1
            tbase = low2+middle+high2
            if tbase % 2:
                continue
            average = spatial_average(zbase, tbase, quadratic_count, low1+high1, middle)
            multinomial = math.factorial(order)//(
                math.factorial(quadratic_count)*math.prod(math.factorial(k) for k in counts))
            if average:
                result[counts] += multinomial*average
    return dict(result)


def raw_expression(order):
    return S.expand(sum(S.Rational(c.numerator, c.denominator)
                        *math.prod(v**k for v, k in zip(COEFFICIENTS, powers))
                        for powers, c in exact_moment(order).items()))


def primitive_expression(order):
    return S.Poly(raw_expression(order), *COEFFICIENTS).clear_denoms()[1].primitive()[1].as_expr()
