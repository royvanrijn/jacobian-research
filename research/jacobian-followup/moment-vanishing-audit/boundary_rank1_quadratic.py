"""Exact moment input on H1=0, H2=w^2, with arbitrary harmonic cubic.

This is an exhaustive rank-one quadratic boundary architecture, not a
coefficient-support census. No exclusion or safety conclusion is asserted
until a retained exact certificate supplies it.
"""
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
import math

import sympy as S

CHARGES = (-2, -1, 0, 1, 2, 3)
COEFFICIENTS = S.symbols('cminus2 cminus1 c0 c1 c2 c3')
LOWEST = S.Symbol('cminus3')


def compositions(total, count):
    if count == 1:
        yield (total,)
    else:
        for first in range(total+1):
            for rest in compositions(total-first, count-1):
                yield (first, *rest)


@lru_cache(maxsize=None)
def odd_factorial(value):
    return math.prod(range(value, 0, -2))


@lru_cache(maxsize=None)
def exact_moment(order, include_lowest=False):
    """Coefficient dictionary after the w^3 cubic coefficient is removed."""
    result = defaultdict(Fraction)
    for cubic_count in range(0, order+1, 2):
        quadratic_count = order-cubic_count
        charges = (-3, *CHARGES) if include_lowest else CHARGES
        for counts in compositions(cubic_count, len(charges)):
            if sum(charge*count for charge, count in zip(charges, counts)) != 2*quadratic_count:
                continue
            if include_lowest:
                low3, low2, low1, middle, high1, high2, high3 = counts
            else:
                low3 = 0
                low2, low1, middle, high1, high2, high3 = counts
            zbase = high1+2*high2+3*high3
            wbase = 3*low3+2*low2+low1+2*quadratic_count
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
                math.factorial(quadratic_count)*math.prod(math.factorial(count) for count in counts))
            if average:
                result[counts] += multinomial*average
    return dict(result)


def expression(order, top_charge=None):
    if top_charge is None:
        variables = COEFFICIENTS
        terms = exact_moment(order)
    else:
        top = CHARGES.index(top_charge)
        variables = COEFFICIENTS[:top]
        accumulated = defaultdict(Fraction)
        for exponents, coefficient in exact_moment(order).items():
            if any(exponents[index] for index in range(top+1, 6)):
                continue
            accumulated[exponents[:top]] += coefficient
        terms = dict(accumulated)
    return S.expand(sum(S.Rational(value.numerator, value.denominator)
                        *math.prod(variable**exponent for variable, exponent in zip(variables, exponents))
                        for exponents, value in terms.items()))


def primitive(order, top_charge):
    variables = COEFFICIENTS[:CHARGES.index(top_charge)]
    poly = S.Poly(expression(order, top_charge), *variables)
    return poly.clear_denoms()[1].primitive()[1]


def next_coefficient_chart(order, top_charge):
    """Top cubic coefficient=1, next lower coefficient=0 via its unipotent orbit."""
    all_variables = (LOWEST, *COEFFICIENTS)
    top = top_charge+3
    variables = all_variables[:top-1]
    terms = defaultdict(Fraction)
    for powers, value in exact_moment(order, True).items():
        if powers[top-1] or any(powers[i] for i in range(top+1, 7)):
            continue
        terms[powers[:top-1]] += value
    expression = S.expand(sum(S.Rational(value.numerator, value.denominator)
                              *math.prod(variable**exponent for variable, exponent in zip(variables, powers))
                              for powers, value in terms.items()))
    return S.Poly(expression, *variables).clear_denoms()[1].primitive()[1]


def raw_expression(order):
    variables = (LOWEST, *COEFFICIENTS)
    return S.expand(sum(S.Rational(value.numerator, value.denominator)
                        *math.prod(variable**power for variable, power in zip(variables, powers))
                        for powers, value in exact_moment(order, True).items()))
