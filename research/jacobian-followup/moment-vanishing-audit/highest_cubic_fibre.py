"""Exact weighted moments on the five-parameter-base zero fibre.

The seven remaining coefficient variables have weights3,4,5. This is a
homogeneous fibre used to test finiteness of the COMPLETE cubic chart
over its five retained base parameters, not an arbitrary support slice.
"""
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
import math

import sympy as S

from boundary_axial_quadratic import spatial_average
from highest_cubic_chart import DEPENDENT, DEPENDENT_WEIGHTS


@lru_cache(maxsize=None)
def weighted_monomials(total, weights):
    if len(weights) == 1:
        return ((total//weights[0],),) if total >= 0 and total % weights[0] == 0 else ()
    result = []
    for first in range(total//weights[0]+1):
        result.extend((first, *rest) for rest in weighted_monomials(total-first*weights[0], weights[1:]))
    return tuple(result)


@lru_cache(maxsize=None)
def exact_moment(order):
    result = defaultdict(Fraction)
    weights = (*DEPENDENT_WEIGHTS, 6)
    for powers in weighted_monomials(3*order, weights):
        base = order-sum(powers)
        assert base >= 0
        lw, lt, bm2, bm1, b0, cm2, cm1, cm3 = powers
        zbase = 3*base
        wbase = lw+2*bm2+bm1+2*cm2+cm1+3*cm3
        assert zbase == wbase
        tbase = lt+bm1+cm2
        if tbase % 2:
            continue
        average = spatial_average(zbase, tbase, b0, cm1, 0)
        if not average:
            continue
        scalar = math.factorial(order)//(math.factorial(base)*math.prod(math.factorial(k) for k in powers))
        for k in range(cm3+1):
            target = list(powers[:-1])
            target[1] += 2*k
            target[4] += 2*(cm3-k)
            result[tuple(target)] += scalar*average*math.comb(cm3, k)*Fraction(-35, 12)**k*Fraction(-7, 4)**(cm3-k)
    return {k: v for k, v in result.items() if v}


def raw_expression(order):
    return S.expand(sum(S.Rational(c.numerator, c.denominator)
                        *math.prod(v**k for v, k in zip(DEPENDENT, powers))
                        for powers, c in exact_moment(order).items()))
