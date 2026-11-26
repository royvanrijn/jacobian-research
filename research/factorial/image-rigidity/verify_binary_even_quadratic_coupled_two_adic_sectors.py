#!/usr/bin/env python3
"""Exact checks for the coupled two-adic sectors of a binary quadratic."""

from collections import Counter
from fractions import Fraction
from math import comb, factorial
import json


def valuation(value):
    value = Fraction(value)
    if not value:
        return float("inf")
    numerator, denominator = value.numerator, value.denominator
    return (numerator & -numerator).bit_length() - (denominator & -denominator).bit_length()


def univariate(order, coefficient):
    return sum(comb(order, index)*factorial(2*index)*coefficient**index
               for index in range(order+1))


def binary(order, a, b, c):
    result = Fraction()
    for i in range(order+1):
        for ell in range(order-i+1):
            for s in range(order-i-ell+1):
                multinomial = comb(order, i)*comb(order-i, ell)*comb(order-i-ell, s)
                result += multinomial*factorial(2*i+ell)*factorial(2*s+ell)*a**i*b**ell*c**s
    return result


def lower_bound(order, degree, lam):
    return (2+lam)*degree - degree.bit_count() - order.bit_count() + (order-degree).bit_count()


def check_beta_weights(counts):
    for h in range(22):
        for s in range(22):
            direct = sum(
                comb(h, ell)*2**ell*factorial(2*h-ell)*factorial(2*s+ell)
                for ell in range(h+1)
            )
            assert direct == 4**h*factorial(h)*factorial(h+s)*factorial(2*s)//factorial(s)
            counts["plane_beta_weights"] += 1


def check_carries(counts):
    threshold = Fraction(-4, 3)
    for j in range(4, 165, 4):
        for d in range(1, j+1):
            value = lower_bound(j, d, threshold)
            if d == 3 and j % 8 == 4:
                assert value == 0
            else:
                assert value > 0
            counts["multiple_of_four_carries"] += 1
    for j in range(1, 162, 8):
        for d in range(2, j+1):
            assert lower_bound(j, d, threshold) > threshold
            counts["one_mod_eight_carries"] += 1
    for j in range(6, 163, 8):
        for d in range(3, j+1):
            assert lower_bound(j, d, threshold) > 2+2*threshold
            counts["six_mod_eight_carries"] += 1

    for tau in range(3, 8):
        L = 1 << tau
        for w in (1, 3, 5, 9):
            j = L*w - 1
            def grouped(e):
                return (tau-3)*e + e.bit_count()+j.bit_count()-(j-e).bit_count()
            for e in range(j+1):
                h, b = divmod(e, L)
                assert grouped(e)-grouped(L*h) == (tau-3)*b+2*b.bit_count()
                counts["deep_digit_blocks"] += 1


def check_perturbations(counts):
    for a, c, epsilon in (
        (Fraction(1, 4), Fraction(1, 2), Fraction(2)),
        (Fraction(3, 8), Fraction(3, 2), Fraction(8, 3)),
        (Fraction(-5, 8), Fraction(1, 3), Fraction(4)),
        (Fraction(-3, 4), Fraction(7, 5), Fraction(32)),
        (Fraction(1, 16), Fraction(7, 2), Fraction(8)),
    ):
        b = 2*a*(1+epsilon)
        mom = [binary(j, a, b, c) for j in range(16)]
        plane = [binary(j, a, 2*a, c) for j in range(16)]
        r = valuation(4*a)
        for j in range(1, len(mom)):
            assert valuation(mom[j]-plane[j]) >= j*r-j.bit_count()+valuation(epsilon)
            counts["actual_perturbation_bounds"] += 1

        delta = a*c-b*b/4
        D = a+c-b
        for j in range(1, len(mom)-1):
            alpha = c-b/2
            beta = a-b/2
            boundary = alpha*univariate(j+1, c)+beta*univariate(j+1, a)
            assert D*mom[j+1] == 4*delta*(j+1)**2*mom[j]-4*delta*j*(j+1)*mom[j-1]+boundary
            counts["coupled_boundary_identities"] += 1

    for t in (Fraction(1, 2), Fraction(3, 2), Fraction(1, 3), Fraction(8, 7)):
        for j in range(1, 42):
            value = univariate(j, t)
            assert valuation(value-1-j*2*t) >= 1
            if j % 2 == 0:
                assert not valuation(value)
                assert valuation(value - 1 - j*2*t - 6*comb(j, 2)*(2*t)**2) >= 2
            counts["integral_boundary_congruences"] += 1


def run():
    counts = Counter()
    check_beta_weights(counts)
    check_carries(counts)
    check_perturbations(counts)
    print(json.dumps(counts, indent=2, sort_keys=True))


if __name__ == "__main__":
    run()
