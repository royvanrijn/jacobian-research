#!/usr/bin/env python3
"""Exact transcriptions for the prime-terminal proof, with no degree sweep premise."""

import json
from fractions import Fraction
from math import comb, factorial, prod


def add(a, b):
    result = [0] * max(len(a), len(b))
    for terms in (a, b):
        for i, term in enumerate(terms):
            result[i] += term
    return result


def multiply(a, b):
    result = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i + j] += x * y
    return result


def trim(a):
    while a and a[-1] == 0:
        a.pop()
    return a


def evaluate(coefficients, x):
    result = 0
    for coefficient in reversed(coefficients):
        result = result * x + coefficient
    return result


def homogeneous(j, a, b, c):
    return sum(
        comb(j, i)
        * comb(j - i, ell)
        * factorial(2 * i + ell)
        * factorial(2 * (j - i - ell) + ell)
        * a**i
        * b**ell
        * c ** (j - i - ell)
        for i in range(j + 1)
        for ell in range(j - i + 1)
    )


def beta_moment(j, c, d):
    return factorial(2 * j + 1) * sum(
        Fraction(
            comb(j, i) * comb(j - i, ell) * factorial(i + ell) * factorial(i),
            factorial(2 * i + ell + 1),
        )
        * (-d) ** i
        * (c - 1) ** ell
        for i in range(j + 1)
        for ell in range(j - i + 1)
    )


def rising_sum(r, i, t):
    return sum(prod(range(ell + 1, ell + 2 * i + 1)) * t**ell for ell in range(r + 1))


def binomial_moment(j, a, b, c, k):
    return sum(comb(j, i) * homogeneous(i, a, b, c) * k ** (j - i) for i in range(j + 1))


def perturbation_moment(j, t, d, k):
    return sum(
        comb(j, count)
        * comb(count, i)
        * k ** (count - i)
        * d**i
        * factorial(2 * (j - count))
        * rising_sum(2 * (j - count), i, t)
        for count in range(j + 1)
        for i in range(count + 1)
    )


def quartic_certificate():
    f = [1, 3, 3, 3, 1]
    double = [1, 0, 3, 0, 3, 0, 3, 0, 1]
    a = [-1, -36, 3, -27, -5, -30, -3, -12]
    b = [23, 39, 39, 12]
    assert trim(add(multiply(a, f), multiply(b, double))) == [22]
    for degree in range(1, 5):
        assert pow(11, degree, 23) != 1


def main():
    primes = (5, 7, 11, 13, 17, 23, 31)
    perturbations = betas = rising = transitions = 0
    for t, d, k in ((2, 3, 1), (-1, 4, 2), (3, -2, -1)):
        c = t**2 + d
        total_gap = 1 + c - 2 * t
        for j in range(10):
            assert perturbation_moment(j, t, d, k) == binomial_moment(j, 1, 2 * t, c, k)
            assert beta_moment(j, c, total_gap) == homogeneous(j, 1, 2 * t, c)
            perturbations += 1
            betas += 1

    for p in primes:
        h = (p - 1) // 2
        for i in range(h + 1):
            for t in range(p):
                rhs = factorial(2 * i) * (t - 1) ** (p - 1 - 2 * i)
                assert (rising_sum(p - 1, i, t) - rhs) % p == 0
                if i:
                    assert (rising_sum(p - 3, i, t) - rhs) % p == 0
                rising += 1
        for s in (Fraction(1, p), Fraction(2, p), Fraction(1, h), Fraction(2, h)):
            low = [int(i > h) + (p - 1 - i) * s for i in range(p)]
            top = [p * s] + [1 + int(i > h) + (p - i) * s for i in range(1, p)] + [2]
            expected_low = [h] if s < Fraction(1, h) else [h, p - 1] if s == Fraction(1, h) else [p - 1]
            expected_top = [0] if s < Fraction(2, p) else [0, p] if s == Fraction(2, p) else [p]
            assert [i for i, value in enumerate(low) if value == min(low)] == expected_low
            assert [i for i, value in enumerate(top) if value == min(top)] == expected_top
            upper = [
                int(i not in (0, 1, p, p + 1)) + (2 * i) // p + (p + 1 - i) * s
                for i in range(p + 2)
            ]
            expected_upper = [1] if s < Fraction(2, p) else [1, p + 1] if s == Fraction(2, p) else [p + 1]
            assert [i for i, value in enumerate(upper) if value == min(upper)] == expected_upper
            transitions += 3
        residue_k = pow(h * (-1) ** h, -1, p)
        rank = factorial(2 * p - 3) // p
        low = comb(p - 2, h) * residue_k * factorial(p - 2)
        assert (rank + low + pow(2, -1, p)) % p == 0

    # The initial off-diagonal residue conditions alone can actually cancel.
    for m in (4, 5):
        g = [comb(m, i) * factorial(2 * i) for i in range(m + 1)]
        assert evaluate(g, 10) % 13 == 0
    quartic_certificate()
    print(json.dumps({
        "prime_transition_faces": transitions,
        "rising_sum_residues": rising,
        "exact_beta_expansions": betas,
        "exact_perturbation_expansions": perturbations,
        "quartic_certificates": 2,
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
