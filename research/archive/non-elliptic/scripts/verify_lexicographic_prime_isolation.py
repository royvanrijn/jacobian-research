#!/usr/bin/env python3
"""Exact arithmetic checks for multi-radial lexicographic prime isolation."""

from __future__ import annotations

import math


def valuation(integer: int, prime: int) -> int:
    answer = 0
    while integer and integer % prime == 0:
        integer //= prime
        answer += 1
    return answer


def factorial_valuation(index: int, prime: int) -> int:
    answer = 0
    while index:
        index //= prime
        answer += index
    return answer


def radial_factorial(exponents: tuple[int, ...]) -> int:
    return math.prod(math.factorial(exponent) for exponent in exponents)


def odd_double_factorial(index: int) -> int:
    return math.prod(range(1, index + 1, 2))


def main() -> None:
    # A mixed-radix order separates these vectors, but p-adic Gaussian
    # factorial order sees only their common total degree.
    for mixed_radix in (2, 3, 7, 20):
        assert 1 + mixed_radix < 2 * mixed_radix

    for prime in (5, 7, 11, 17):
        left = (prime, prime)
        right = (0, 2 * prime)
        assert sum(factorial_valuation(n, prime) for n in left) == 2
        assert sum(factorial_valuation(n, prime) for n in right) == 2

    # Legendre's formula gives total radial degree at a prime dilation.
    vectors = ((0, 0), (1, 1), (0, 2), (3, 1), (1, 4, 2))
    for prime in (7, 11, 17):
        for vector in vectors:
            if max(vector, default=0) < prime:
                dilated = tuple(prime * n for n in vector)
                order = sum(factorial_valuation(n, prime) for n in dilated)
                assert order == sum(vector)

    # The normalized unit is (-1)^|n| product(n_i!) modulo p.
    for prime in (7, 11, 17):
        for vector in vectors:
            if max(vector, default=0) < prime:
                total = sum(vector)
                dilated = tuple(prime * n for n in vector)
                unit = radial_factorial(dilated) // (prime**total)
                expected = ((-1) ** total) * radial_factorial(vector)
                assert (unit - expected) % prime == 0

    # Long's first zero-weight radial face is invisible at every dilation.
    for exponent in (1, 2, 5, 6, 35):
        u1 = radial_factorial((exponent, 0))
        u2 = radial_factorial((0, exponent))
        assert -u1 + u2 == 0

    # At the real d=3 boundary, U=ZW and V=T^2 have factorial and odd
    # double-factorial weights.  Prime dilation still sees n+k.
    for prime in (7, 11, 17):
        for n, k in ((1, 0), (0, 1), (1, 2), (3, 1)):
            if max(n, 2 * k) < prime:
                hybrid = math.factorial(prime * n) * odd_double_factorial(
                    2 * prime * k - 1
                )
                assert valuation(hybrid, prime) == n + k

    # Lowest total degree does not control the Frobenius cross terms.
    # In (U1^3+U2^3+U3^3)^7, multiplicities (2,2,3) give coefficient 210
    # and exponent vector (6,6,9).  Its combined score is below that of
    # a pure Frobenius monomial U1^21.
    cross_coefficient = math.factorial(7) // (
        math.factorial(2) * math.factorial(2) * math.factorial(3)
    )
    cross_score = valuation(cross_coefficient, 7) + sum(
        factorial_valuation(n, 7) for n in (6, 6, 9)
    )
    frobenius_score = factorial_valuation(21, 7)
    assert cross_score == 2 < 3 == frobenius_score

    # Distinct residue characteristics cannot be successive in Z.
    assert math.gcd(5, 7) == 1

    # On a stretched lattice the p-adic order becomes the weighted degree.
    stretch = (1, 5, 25)
    for prime in (101, 103):
        for vector in ((1, 2, 0), (0, 1, 1), (3, 0, 2)):
            dilated = tuple(prime * a * n for a, n in zip(stretch, vector))
            order = sum(factorial_valuation(n, prime) for n in dilated)
            weighted_degree = sum(a * n for a, n in zip(stretch, vector))
            assert order == weighted_degree

    print("PASS ordinary prime dilation sees total radial degree")
    print("PASS mixed-radix order is not the Gaussian factorial order")
    print("PASS the real d=3 hybrid radial functional also sees total degree")
    print("PASS coordinate permutations remain exactly indistinguishable")
    print("PASS total degree alone misses lower-score Frobenius cross terms")
    print("PASS stretched pure Frobenius monomials realize weighted degree")


if __name__ == "__main__":
    main()
