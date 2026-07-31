#!/usr/bin/env python3
"""Exact torsion--torus trace identities for the binary GVC frontier.

Let ``F=C_q`` and let ``u`` lie in the group algebra

    QQ[Z x F] = QQ[z,z^-1][s]/(s^q-1).

The coefficient of the identity in both factors is related to the
regular representation by

    CT_(Z x F)(u^N) = (1/q) CT_z Tr(Reg(u)^N).

If ``D(t,z)=det(I-t Reg(u))``, Newton's trace identity gives

    sum_(N>=1) CT_(Z x F)(u^N) t^N
      = -(1/q) CT_z t*d/dt log D(t,z).

These formulas turn a finite-character cancellation into one
log-determinant identity.  For primes ``p=1 mod q``, Frobenius also
fixes the torsion coordinate:

    u^p = Frobenius(u)  in F_p[Z x F].

The script verifies all three identities exactly for cyclic groups of
orders two and three.  It additionally checks that the ``C_2`` regular
representation diagonalizes the coordinate-reversed pair
``z^-r G(z), z^r G(z)``.

This is an exact reduction, not the torsion--torus trace lemma itself
and not a proof of GVC(2).
"""

from __future__ import annotations

from fractions import Fraction

import sympy as sp


GroupElement = dict[tuple[int, int], Fraction]


def add_term(
    element: GroupElement,
    free_exponent: int,
    torsion_exponent: int,
    coefficient: Fraction,
    torsion_order: int,
) -> None:
    key = (free_exponent, torsion_exponent % torsion_order)
    element[key] = element.get(key, Fraction(0)) + coefficient
    if element[key] == 0:
        del element[key]


def multiply(
    left: GroupElement,
    right: GroupElement,
    torsion_order: int,
) -> GroupElement:
    answer: GroupElement = {}
    for (left_free, left_torsion), left_value in left.items():
        for (right_free, right_torsion), right_value in right.items():
            add_term(
                answer,
                left_free + right_free,
                left_torsion + right_torsion,
                left_value * right_value,
                torsion_order,
            )
    return answer


def power(
    element: GroupElement,
    exponent: int,
    torsion_order: int,
) -> GroupElement:
    answer: GroupElement = {(0, 0): Fraction(1)}
    base = element
    remaining = exponent
    while remaining:
        if remaining & 1:
            answer = multiply(answer, base, torsion_order)
        base = multiply(base, base, torsion_order)
        remaining //= 2
    return answer


def identity_coefficient(element: GroupElement) -> Fraction:
    return element.get((0, 0), Fraction(0))


def laurent_component(
    element: GroupElement,
    torsion_exponent: int,
    z: sp.Symbol,
    torsion_order: int,
) -> sp.Expr:
    return sp.expand(
        sum(
            sp.Rational(value.numerator, value.denominator) * z**free
            for (free, torsion), value in element.items()
            if torsion == torsion_exponent % torsion_order
        )
    )


def regular_matrix(
    element: GroupElement,
    torsion_order: int,
    z: sp.Symbol,
) -> sp.Matrix:
    components = tuple(
        laurent_component(element, torsion, z, torsion_order)
        for torsion in range(torsion_order)
    )
    return sp.Matrix(
        torsion_order,
        torsion_order,
        lambda row, column: components[(row - column) % torsion_order],
    )


def constant_term(expression: sp.Expr, z: sp.Symbol) -> sp.Expr:
    expanded = sp.expand(expression)
    answer = sp.Integer(0)
    for term in sp.Add.make_args(expanded):
        powers = term.as_powers_dict()
        if powers.get(z, 0) == 0:
            answer += term
    return sp.simplify(answer)


def verify_regular_trace(
    element: GroupElement,
    torsion_order: int,
    depth: int,
) -> tuple[Fraction, ...]:
    z = sp.Symbol("z")
    matrix = regular_matrix(element, torsion_order, z)
    direct_rows = []
    matrix_power = sp.eye(torsion_order)

    for moment in range(1, depth + 1):
        direct = identity_coefficient(
            power(element, moment, torsion_order)
        )
        matrix_power = matrix_power * matrix
        traced = constant_term(sp.trace(matrix_power), z) / torsion_order
        assert traced == sp.Rational(direct.numerator, direct.denominator)
        direct_rows.append(direct)

    return tuple(direct_rows)


def verify_log_determinant(
    element: GroupElement,
    torsion_order: int,
    depth: int,
) -> sp.Expr:
    z, t = sp.symbols("z t")
    matrix = regular_matrix(element, torsion_order, z)
    determinant = sp.factor(
        (sp.eye(torsion_order) - t * matrix).det()
    )
    logarithmic_derivative = sp.series(
        -t * sp.diff(determinant, t) / determinant,
        t,
        0,
        depth + 1,
    ).removeO()

    for moment in range(1, depth + 1):
        coefficient = sp.expand(logarithmic_derivative).coeff(t, moment)
        traced = constant_term(coefficient, z) / torsion_order
        direct = identity_coefficient(
            power(element, moment, torsion_order)
        )
        assert traced == sp.Rational(direct.numerator, direct.denominator)

    return determinant


def reduce_mod_prime(
    element: GroupElement,
    prime: int,
) -> dict[tuple[int, int], int]:
    answer = {}
    for key, value in element.items():
        denominator = value.denominator % prime
        assert denominator
        answer[key] = (
            value.numerator
            * pow(denominator, -1, prime)
        ) % prime
    return {key: value for key, value in answer.items() if value}


def multiply_mod_prime(
    left: dict[tuple[int, int], int],
    right: dict[tuple[int, int], int],
    torsion_order: int,
    prime: int,
) -> dict[tuple[int, int], int]:
    answer: dict[tuple[int, int], int] = {}
    for (left_free, left_torsion), left_value in left.items():
        for (right_free, right_torsion), right_value in right.items():
            key = (
                left_free + right_free,
                (left_torsion + right_torsion) % torsion_order,
            )
            answer[key] = (
                answer.get(key, 0) + left_value * right_value
            ) % prime
    return {key: value for key, value in answer.items() if value}


def power_mod_prime(
    element: dict[tuple[int, int], int],
    exponent: int,
    torsion_order: int,
    prime: int,
) -> dict[tuple[int, int], int]:
    answer = {(0, 0): 1}
    for _ in range(exponent):
        answer = multiply_mod_prime(
            answer,
            element,
            torsion_order,
            prime,
        )
    return answer


def verify_frobenius(
    element: GroupElement,
    torsion_order: int,
    prime: int,
) -> None:
    assert prime % torsion_order == 1
    reduced = reduce_mod_prime(element, prime)
    actual = power_mod_prime(
        reduced,
        prime,
        torsion_order,
        prime,
    )
    expected = {
        (prime * free, torsion): pow(value, prime, prime)
        for (free, torsion), value in reduced.items()
    }
    assert actual == expected


def reversal_element(
    coefficients: dict[int, Fraction],
    target_slope: int,
) -> GroupElement:
    """Return the C2 element with character components z^-r G and z^r G."""

    answer: GroupElement = {}
    for exponent, coefficient in coefficients.items():
        # e_+=(1+s)/2 carries z^-r G; e_-=(1-s)/2 carries z^r G.
        for shifted_exponent, identity_sign, torsion_sign in (
            (exponent - target_slope, 1, 1),
            (exponent + target_slope, 1, -1),
        ):
            add_term(
                answer,
                shifted_exponent,
                0,
                coefficient * identity_sign / 2,
                2,
            )
            add_term(
                answer,
                shifted_exponent,
                1,
                coefficient * torsion_sign / 2,
                2,
            )
    return answer


def verify_reversal_diagonalization(depth: int) -> tuple[Fraction, ...]:
    coefficients = {
        -2: Fraction(1),
        -1: Fraction(-1),
        0: Fraction(2),
        1: Fraction(3),
        2: Fraction(1),
    }
    element = reversal_element(coefficients, 1)
    rows = verify_regular_trace(element, 2, depth)

    z = sp.Symbol("z")
    polynomial = sum(
        sp.Rational(value.numerator, value.denominator) * z**exponent
        for exponent, value in coefficients.items()
    )
    for moment, row in enumerate(rows, 1):
        power_expression = sp.expand(polynomial**moment)
        expected = (
            power_expression.coeff(z, moment)
            + power_expression.coeff(z, -moment)
        ) / 2
        assert expected == sp.Rational(row.numerator, row.denominator)
    return rows


def verify() -> None:
    examples = {
        2: {
            (-2, 0): Fraction(2),
            (-1, 1): Fraction(-1),
            (0, 0): Fraction(3),
            (1, 1): Fraction(2),
            (2, 0): Fraction(1),
        },
        3: {
            (-2, 0): Fraction(1),
            (-1, 1): Fraction(2),
            (0, 2): Fraction(-1),
            (1, 0): Fraction(3),
            (2, 1): Fraction(1),
        },
    }
    depths = {2: 6, 3: 5}
    primes = {2: 5, 3: 7}

    for torsion_order, element in examples.items():
        rows = verify_regular_trace(
            element,
            torsion_order,
            depths[torsion_order],
        )
        determinant = verify_log_determinant(
            element,
            torsion_order,
            depths[torsion_order],
        )
        verify_frobenius(
            element,
            torsion_order,
            primes[torsion_order],
        )
        print(
            f"C_{torsion_order}: rows={rows}, "
            f"det(I-t Reg(u))={determinant}"
        )

    reversal_rows = verify_reversal_diagonalization(6)
    print(f"C_2 reversal trace rows: {reversal_rows}")
    print(
        "PASS torsion--torus regular trace, log determinant, "
        "and p=1 mod exp(F) Frobenius identities"
    )
    print(
        "STATUS: exact reduction established; finite-character "
        "separation is verified by the repeated-digit checker"
    )


if __name__ == "__main__":
    verify()
