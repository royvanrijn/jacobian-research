"""Exact input for H1=0, H2=z^2+w^2, with all seven harmonic cubics.

This represents the complete non-nilpotent rank-two quadratic stratum
after scalar and oriented orthogonal normalization. Exclusion and purity
are unproved until an exact certificate is retained. No support is dropped.
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
def exact_moment(order):
    result = defaultdict(Fraction)
    for cubic_count in range(0, order+1, 2):
        quadratic_count = order-cubic_count
        for counts in compositions(cubic_count, 7):
            charge = sum(j*k for j, k in zip(CHARGES, counts))
            if (2*quadratic_count-charge) % 4:
                continue
            zquadratic = (2*quadratic_count-charge)//4
            if zquadratic < 0 or zquadratic > quadratic_count:
                continue
            wquadratic = quadratic_count-zquadratic
            low3, low2, low1, middle, high1, high2, high3 = counts
            zbase = 2*zquadratic+high1+2*high2+3*high3
            wbase = 2*wquadratic+3*low3+2*low2+low1
            tbase = low2+middle+high2
            assert zbase == wbase
            if tbase % 2:
                continue
            height_power = low1+high1
            average = Fraction(0)
            for i in range(height_power+1):
                for j in range(middle+1):
                    power_z = zbase+i+j
                    power_t = tbase+2*(height_power+middle-i-j)
                    scalar = Fraction(math.comb(height_power, i)*math.comb(middle, j)
                                      *(-1)**i*(-3)**j, 2**i)
                    average += scalar*Fraction(math.factorial(power_z)*odd_factorial(power_t-1),
                                               odd_factorial(2*power_z+power_t+1))
            multinomial = math.factorial(order)//(
                math.factorial(zquadratic)*math.factorial(wquadratic)
                *math.prod(math.factorial(k) for k in counts))
            if average:
                result[counts] += multinomial*average
    return dict(result)


def raw_expression(order):
    return S.expand(sum(S.Rational(c.numerator, c.denominator)
                        *math.prod(v**k for v, k in zip(COEFFICIENTS, powers))
                        for powers, c in exact_moment(order).items()))


def primitive_expression(order):
    return S.Poly(raw_expression(order), *COEFFICIENTS).clear_denoms()[1].primitive()[1].as_expr()
