#!/usr/bin/env python3
"""Exact regressions for the factorially weighted multitorus theorem.

The all-order proof is written in
extended-geometry/FACTORIALLY_WEIGHTED_MULTITORUS_THEOREM.md.  This script
checks representative finite identities; it is not a proof of the theorem.
"""

from __future__ import annotations

import math
from collections.abc import Iterable


Weight = tuple[int, ...]
Monomial = tuple[int, Weight]
Polynomial = dict[Monomial, int]


def add_term(answer: Polynomial, monomial: Monomial, coefficient: int) -> None:
    answer[monomial] = answer.get(monomial, 0) + coefficient
    if answer[monomial] == 0:
        del answer[monomial]


def multiply(left: Polynomial, right: Polynomial) -> Polynomial:
    answer: Polynomial = {}
    for (left_radial, left_weight), left_coefficient in left.items():
        for (right_radial, right_weight), right_coefficient in right.items():
            assert len(left_weight) == len(right_weight)
            monomial = (
                left_radial + right_radial,
                tuple(a + b for a, b in zip(left_weight, right_weight)),
            )
            add_term(
                answer,
                monomial,
                left_coefficient * right_coefficient,
            )
    return answer


def power(polynomial: Polynomial, exponent: int) -> Polynomial:
    assert exponent >= 0
    rank = len(next(iter(polynomial))[1])
    answer: Polynomial = {(0, (0,) * rank): 1}
    base = polynomial
    current = exponent
    while current:
        if current & 1:
            answer = multiply(answer, base)
        base = multiply(base, base)
        current //= 2
    return answer


def gamma_functional(polynomial: Polynomial) -> int:
    return sum(
        coefficient * math.factorial(radial)
        for (radial, weight), coefficient in polynomial.items()
        if all(entry == 0 for entry in weight)
    )


def coefficient(
    polynomial: Polynomial,
    radial: int,
    weight: Iterable[int],
) -> int:
    return polynomial.get((radial, tuple(weight)), 0)


def check_balanced_rank_two_prime_isolation() -> None:
    # The three angular weights (1,0), (0,1), and (-1,-1) balance only in
    # equal proportions.  Their radial degrees 1,2,0 have average one.
    polynomial: Polynomial = {
        (1, (1, 0)): 1,
        (2, (0, 1)): 1,
        (0, (-1, -1)): 1,
    }
    face_power = power(polynomial, 3)
    lowest_radial = 3
    lowest_coefficient = coefficient(face_power, lowest_radial, (0, 0))
    assert lowest_coefficient == 6
    assert all(
        radial >= lowest_radial
        for radial, weight in face_power
        if weight == (0, 0)
    )

    for prime in (5, 7):
        dilated = power(face_power, prime)
        moment = gamma_functional(dilated)
        divisor = math.factorial(lowest_radial * prime)
        assert moment != 0
        assert moment % divisor == 0
        assert (
            moment // divisor - pow(lowest_coefficient, prime, prime)
        ) % prime == 0


def check_strict_separation_and_mixed_cutoff() -> None:
    polynomial: Polynomial = {
        (0, (1, 0)): 2,
        (1, (2, 1)): -3,
    }
    multiplier: Polynomial = {
        (0, (-3, -1)): 5,
        (2, (-1, 0)): 7,
    }
    for exponent in range(1, 9):
        assert gamma_functional(power(polynomial, exponent)) == 0
    for exponent in range(4, 9):
        assert gamma_functional(
            multiply(multiplier, power(polynomial, exponent))
        ) == 0


def check_gaussian_embedding() -> None:
    # Under Z -> z and W -> u z^(-1), P=Z+W becomes z+u/z.
    embedded: Polynomial = {
        (0, (1,)): 1,
        (1, (-1,)): 1,
    }
    for exponent in range(1, 11):
        moment = gamma_functional(power(embedded, exponent))
        if exponent % 2:
            expected = 0
        else:
            half = exponent // 2
            expected = math.factorial(exponent) // math.factorial(half)
        assert moment == expected


def main() -> None:
    check_balanced_rank_two_prime_isolation()
    check_strict_separation_and_mixed_cutoff()
    check_gaussian_embedding()
    print("PASS multitorus: rank-two exposed coefficient and prime isolation")
    print("PASS multitorus: strict separation gives the mixed cutoff")
    print("PASS multitorus: the circular Gaussian embedding intertwines moments")


if __name__ == "__main__":
    main()
