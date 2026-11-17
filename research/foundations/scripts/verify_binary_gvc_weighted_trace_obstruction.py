#!/usr/bin/env python3
"""Exact regressions for the weighted finite-trace frontier.

The repeated-digit argument has a weighted extension.  If

    sum_i w_i * CT(f_i**n) = 0

for every positive ``n``, then for every tuple of positive integers
``(m_0, ..., m_s)`` one also has

    sum_i w_i * product_j CT(f_i**m_j) = 0.

Consequently the total weight of every nonzero class of Laurent
polynomials having the same complete constant-term moment sequence is
zero.  Equal weights cannot cancel in characteristic zero, which
recovers finite-trace separation.  Character weights can cancel,
however.

This script checks three exact warnings needed at the binary GVC
frontier:

* mixed base-p digit factorization for Laurent constant terms;
* an affine C2 character projection with all pure rows zero but a
  persistent mixed row, using the isoperiodic dilation pair
  z+z^-1 and z^2+z^-2;
* failure of the analogous digit factorization after inserting binary
  radial factorials.

The last example uses

    C(x,y) = y^2 + 4*x*y + 2*x^2.

For p=11 and N=1+11, its factorial moment has 11-adic valuation three,
whereas the no-carry prediction has valuation two and nonzero residue.
Thus Laurent finite-trace separation does not by itself expose a
factorial Hall/carry shell.
"""

from __future__ import annotations

from math import factorial


Laurent = dict[int, int]


def multiply(
    left: Laurent,
    right: Laurent,
    modulus: int | None = None,
) -> Laurent:
    answer: Laurent = {}
    for left_exponent, left_value in left.items():
        for right_exponent, right_value in right.items():
            exponent = left_exponent + right_exponent
            value = answer.get(exponent, 0) + left_value * right_value
            if modulus is not None:
                value %= modulus
            if value:
                answer[exponent] = value
            elif exponent in answer:
                del answer[exponent]
    return answer


def power(
    polynomial: Laurent,
    exponent: int,
    modulus: int | None = None,
) -> Laurent:
    answer: Laurent = {0: 1}
    base = polynomial
    remaining = exponent
    while remaining:
        if remaining & 1:
            answer = multiply(answer, base, modulus)
        remaining //= 2
        if remaining:
            base = multiply(base, base, modulus)
    return answer


def constant_term(polynomial: Laurent) -> int:
    return polynomial.get(0, 0)


def coefficient(polynomial: Laurent, exponent: int) -> int:
    return polynomial.get(exponent, 0)


def digit_index(digits: tuple[int, ...], prime: int) -> int:
    return sum(
        digit * prime**position
        for position, digit in enumerate(digits)
    )


def verify_mixed_digit_factorization() -> None:
    prime = 17
    digits = (2, 3, 1)
    polynomials = (
        {-2: 1, -1: -2, 0: 3, 1: 1},
        {-1: 2, 0: -1, 2: 1},
    )
    moment = digit_index(digits, prime)

    for polynomial in polynomials:
        largest_exponent = max(abs(value) for value in polynomial)
        assert prime > largest_exponent * max(digits)
        direct = constant_term(
            power(polynomial, moment, prime)
        ) % prime
        predicted = 1
        for digit in digits:
            predicted *= constant_term(
                power(polynomial, digit, prime)
            )
            predicted %= prime
        assert direct == predicted


def verify_affine_character_obstruction(depth: int = 12) -> None:
    first = {-1: 1, 1: 1}
    second = {-2: 1, 2: 1}
    first_power: Laurent = {0: 1}
    second_power: Laurent = {0: 1}

    for moment in range(1, depth + 1):
        first_power = multiply(first_power, first)
        second_power = multiply(second_power, second)
        # Dilation z -> z^2 preserves every constant-term power.
        assert constant_term(first_power) == constant_term(second_power)

        # The affine C2 character weights (1,-1) therefore cancel.
        affine_row = (
            constant_term(first_power)
            - constant_term(second_power)
        )
        assert affine_row == 0

    # A fixed Laurent multiplier distinguishes the two components at
    # every odd moment.  This is a trace obstruction, not a GVC pair.
    for moment in (1, 3, 5, 7):
        first_power = power(first, moment)
        second_power = power(second, moment)
        mixed_row = (
            coefficient(first_power, 1)
            - coefficient(second_power, 1)
        )
        assert mixed_row != 0


def homogeneous_factorial_moment(
    coefficients: tuple[int, ...],
    degree: int,
    moment: int,
) -> int:
    """Return L(C**moment), where L(x**a*y**b)=a!*b!."""

    coefficient_row = [1]
    for _ in range(moment):
        next_row = [0] * (len(coefficient_row) + degree)
        for current_exponent, current_value in enumerate(coefficient_row):
            for exponent, value in enumerate(coefficients):
                next_row[current_exponent + exponent] += (
                    current_value * value
                )
        coefficient_row = next_row

    total_degree = degree * moment
    return sum(
        value
        * factorial(x_exponent)
        * factorial(total_degree - x_exponent)
        for x_exponent, value in enumerate(coefficient_row)
    )


def valuation(value: int, prime: int) -> int:
    answer = 0
    while value and value % prime == 0:
        value //= prime
        answer += 1
    return answer


def verify_factorial_digit_failure() -> None:
    # Coefficient index is the x-exponent, so this is
    # y^2 + 4*x*y + 2*x^2.
    coefficients = (1, 4, 2)
    degree = 2
    prime = 11
    base_moment = homogeneous_factorial_moment(
        coefficients,
        degree,
        1,
    )
    assert base_moment == 10

    repeated_digit_index = 1 + prime
    large_moment = homogeneous_factorial_moment(
        coefficients,
        degree,
        repeated_digit_index,
    )
    no_carry_valuation = degree
    assert valuation(large_moment, prime) == 3
    normalized_residue = (
        large_moment // prime**no_carry_valuation
    ) % prime
    no_carry_prediction = base_moment**2 % prime
    assert normalized_residue == 0
    assert no_carry_prediction == 1
    assert normalized_residue != no_carry_prediction


def verify() -> None:
    verify_mixed_digit_factorization()
    verify_affine_character_obstruction()
    verify_factorial_digit_failure()
    print("PASS weighted mixed-digit factorization")
    print(
        "PASS affine C2 isoperiodic cancellation with a persistent "
        "mixed row"
    )
    print(
        "PASS explicit failure of factorial repeated-digit "
        "factorization"
    )
    print(
        "STATUS: exact obstruction to promoting finite Laurent traces "
        "to Hall/carry shells"
    )


if __name__ == "__main__":
    verify()
