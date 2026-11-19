#!/usr/bin/env python3
"""Exact main-term algebra, not a test of the sixth-moment remainder.

Checks finite rational swap cancellations and polynomial generating-function
identities used to specify the ternary-divisor main. No numerical zeta moment
or unbounded correlation estimate is certified here.
"""

from fractions import Fraction as Q
from itertools import combinations
from math import factorial, prod
import json


def add(*polynomials):
    result = {}
    for polynomial in polynomials:
        for power, value in polynomial.items():
            result[power] = result.get(power, Q(0)) + value
    return {power: value for power, value in result.items() if value}


def scale(polynomial, scalar):
    return {power: scalar * value for power, value in polynomial.items()
            if scalar * value}


def multiply(left, right):
    result = {}
    for (a, j), x in left.items():
        for (b, k), y in right.items():
            power = (a + b, j + k)
            result[power] = result.get(power, Q(0)) + x * y
    return {power: value for power, value in result.items() if value}


def falling(polynomial, count):
    result = {(0, 0): Q(1)}
    for offset in range(count):
        result = multiply(result, add(polynomial, {(0, 0): Q(-offset)}))
    return scale(result, Q(1, factorial(count)))


def check_generating_functions():
    one = {(0, 0): Q(1)}
    alpha = {(1, 0): Q(1)}
    j = {(0, 1): Q(1)}
    # Coefficient of q^j in g_3(q,p^alpha)/(1-q)^3.
    b0 = scale(multiply(add(alpha, one), add(alpha, scale(one, 2))), Q(1, 2))
    b1 = scale(multiply(alpha, add(alpha, scale(one, 2))), Q(-1))
    b2 = scale(multiply(alpha, add(alpha, one)), Q(1, 2))
    numerator_coeff = add(
        multiply(b0, falling(add(j, scale(one, 2)), 2)),
        multiply(b1, falling(add(j, one), 2)),
        multiply(b2, falling(j, 2)),
    )
    expected = scale(multiply(add(alpha, j, one),
                              add(alpha, j, scale(one, 2))), Q(1, 2))
    assert numerator_coeff == expected

    # Sum d_3(p^j)^2 q^j = (1+4q+q^2)/(1-q)^5.
    square_coeff = add(falling(add(j, scale(one, 4)), 4),
                       scale(falling(add(j, scale(one, 3)), 4), Q(4)),
                       falling(add(j, scale(one, 2)), 4))
    d3 = scale(multiply(add(j, one), add(j, scale(one, 2))), Q(1, 2))
    assert square_coeff == multiply(d3, d3)
    return 2


def check_swaps(alpha, beta):
    terms = []
    for count in range(4):
        for chosen_a in combinations(alpha, count):
            for chosen_b in combinations(beta, count):
                a = tuple(x for x in alpha if x not in chosen_a)
                a += tuple(-x for x in chosen_b)
                b = tuple(x for x in beta if x not in chosen_b)
                b += tuple(-x for x in chosen_a)
                denominator = prod(x + y for x in a for y in b)
                assert denominator != 0
                terms.append((sum(chosen_a) + sum(chosen_b), denominator))
    assert len(terms) == 20
    for exponent in range(9):
        coefficient = sum((Q((-c) ** exponent, d) for c, d in terms), Q(0))
        assert coefficient == 0, (alpha, beta, exponent, coefficient)
    coefficient = sum((Q((-c) ** 9, d) for c, d in terms), Q(0))
    assert coefficient == 42, (alpha, beta, coefficient)
    return 20, 9


def main():
    identities = check_generating_functions()
    rays = [((1, 2, 3), (4, 5, 6)),
            ((1, 3, 7), (2, 5, 11)),
            ((2, 7, 13), (3, 11, 19)),
            ((3, 8, 21), (5, 13, 34))]
    swaps = polar = 0
    for alpha, beta in rays:
        count, cancelled = check_swaps(alpha, beta)
        swaps += count
        polar += cancelled
    diagonal = 2 * Q(3, 2) ** 9
    offdiagonal = Q(42) - diagonal
    assert diagonal == Q(19683, 256)
    assert offdiagonal == Q(-8931, 256)
    print(json.dumps({
        "status": "PASS",
        "arithmetic": "exact rational polynomial identities and Fractions",
        "generic_rays_checked": len(rays),
        "swap_terms_checked": swaps,
        "polar_coefficients_cancelled": polar,
        "all_order_generating_function_identities": identities,
        "total_coefficient_scaled_by_9_factorial": 42,
        "diagonal_coefficient_scaled_by_9_factorial": str(diagonal),
        "offdiagonal_coefficient_scaled_by_9_factorial": str(offdiagonal),
        "sixth_moment_remainder_proved": False,
    }, indent=2))


if __name__ == "__main__":
    main()
