#!/usr/bin/env python3
"""Independent stdlib replay of the 12-variable degree-three Keller map.

This file deliberately imports neither SymPy nor repository algebra helpers.
Polynomials are sparse dictionaries over fractions.Fraction.
"""

from __future__ import annotations

from fractions import Fraction


N = 12
ZERO = (0,) * N


def const(value: int | Fraction) -> dict[tuple[int, ...], Fraction]:
    coefficient = Fraction(value)
    return {ZERO: coefficient} if coefficient else {}


def var(index: int) -> dict[tuple[int, ...], Fraction]:
    exponents = [0] * N
    exponents[index] = 1
    return {tuple(exponents): Fraction(1)}


def monomial(
    coefficient: int | Fraction,
    *powers: tuple[int, int],
) -> dict[tuple[int, ...], Fraction]:
    exponents = [0] * N
    for index, exponent in powers:
        exponents[index - 1] += exponent
    value = Fraction(coefficient)
    return {tuple(exponents): value} if value else {}


def add(*polynomials: dict[tuple[int, ...], Fraction]) -> dict[tuple[int, ...], Fraction]:
    result: dict[tuple[int, ...], Fraction] = {}
    for polynomial in polynomials:
        for exponents, coefficient in polynomial.items():
            updated = result.get(exponents, Fraction(0)) + coefficient
            if updated:
                result[exponents] = updated
            else:
                result.pop(exponents, None)
    return result


def scale(
    polynomial: dict[tuple[int, ...], Fraction],
    coefficient: int | Fraction,
) -> dict[tuple[int, ...], Fraction]:
    value = Fraction(coefficient)
    return {
        exponents: value * old_coefficient
        for exponents, old_coefficient in polynomial.items()
        if value * old_coefficient
    }


def multiply(
    left: dict[tuple[int, ...], Fraction],
    right: dict[tuple[int, ...], Fraction],
) -> dict[tuple[int, ...], Fraction]:
    result: dict[tuple[int, ...], Fraction] = {}
    for left_exponents, left_coefficient in left.items():
        for right_exponents, right_coefficient in right.items():
            exponents = tuple(
                left_value + right_value
                for left_value, right_value in zip(left_exponents, right_exponents)
            )
            result[exponents] = (
                result.get(exponents, Fraction(0))
                + left_coefficient * right_coefficient
            )
    return {
        exponents: coefficient
        for exponents, coefficient in result.items()
        if coefficient
    }


def derivative(
    polynomial: dict[tuple[int, ...], Fraction],
    index: int,
) -> dict[tuple[int, ...], Fraction]:
    result: dict[tuple[int, ...], Fraction] = {}
    for exponents, coefficient in polynomial.items():
        exponent = exponents[index]
        if exponent:
            reduced = list(exponents)
            reduced[index] -= 1
            result[tuple(reduced)] = coefficient * exponent
    return result


def evaluate(
    polynomial: dict[tuple[int, ...], Fraction],
    point: tuple[Fraction, ...],
) -> Fraction:
    result = Fraction(0)
    for exponents, coefficient in polynomial.items():
        term = coefficient
        for value, exponent in zip(point, exponents):
            term *= value**exponent
        result += term
    return result


def determinant(
    matrix: list[list[dict[tuple[int, ...], Fraction]]],
) -> tuple[dict[tuple[int, ...], Fraction], int]:
    row_order = sorted(
        range(N),
        key=lambda row: sum(bool(entry) for entry in matrix[row]),
    )
    ordered = [matrix[row] for row in row_order]
    inversions = sum(
        row_order[left] > row_order[right]
        for left in range(N)
        for right in range(left + 1, N)
    )
    permutation_sign = -1 if inversions % 2 else 1
    memo: dict[int, dict[tuple[int, ...], Fraction]] = {}

    def recurse(mask: int, row: int) -> dict[tuple[int, ...], Fraction]:
        if row == N:
            return const(1)
        if mask in memo:
            return memo[mask]
        result: dict[tuple[int, ...], Fraction] = {}
        position = 0
        for column in range(N):
            if mask & (1 << column):
                if ordered[row][column]:
                    cofactor = multiply(
                        ordered[row][column],
                        recurse(mask ^ (1 << column), row + 1),
                    )
                    result = add(
                        result,
                        scale(cofactor, -1 if position % 2 else 1),
                    )
                position += 1
        memo[mask] = result
        return result

    value = scale(recurse((1 << N) - 1, 0), permutation_sign)
    return value, len(memo)


def main() -> None:
    z = [var(index) for index in range(N)]
    half = Fraction(1, 2)
    k = [
        add(
            z[0],
            monomial(half, (1, 2), (11, 1)),
            monomial(-Fraction(3, 2), (1, 2), (2, 1)),
            monomial(-1, (1, 1), (12, 1), (3, 1)),
            monomial(-1, (11, 1), (12, 1)),
        ),
        add(
            z[1],
            monomial(12, (1, 1), (2, 2)),
            monomial(-1, (1, 1), (2, 1), (9, 1)),
            monomial(-6, (1, 1), (3, 1), (8, 1)),
            monomial(3, (1, 1), (3, 1)),
            monomial(3, (1, 1), (5, 1), (8, 1)),
            monomial(-9, (2, 2), (6, 1)),
            monomial(3, (3, 1), (6, 1), (8, 1)),
            monomial(-3, (5, 1), (6, 1)),
            monomial(-1, (8, 1), (9, 1)),
        ),
        add(
            z[2],
            monomial(-1, (1, 1), (10, 1), (2, 1)),
            monomial(3, (1, 1), (2, 1), (3, 1)),
            monomial(1, (1, 1), (7, 1), (8, 1)),
            monomial(-1, (10, 1), (8, 1)),
            monomial(-3, (2, 2), (4, 1)),
            monomial(-7, (2, 2), (8, 1)),
            monomial(4, (2, 2)),
            monomial(-3, (2, 1), (3, 1), (6, 1)),
            monomial(1, (2, 1), (5, 1), (6, 1)),
            monomial(1, (3, 1), (4, 1), (8, 1)),
            monomial(-1, (4, 1), (5, 1)),
            monomial(-1, (6, 1), (7, 1)),
        ),
        add(
            z[3],
            monomial(-2, (1, 1), (2, 1), (8, 1)),
            monomial(-1, (8, 2)),
        ),
        add(
            z[4],
            monomial(1, (1, 1), (2, 1), (3, 1)),
            monomial(3, (2, 2)),
        ),
        add(z[5], monomial(1, (1, 2), (2, 1))),
        add(
            z[6],
            monomial(3, (2, 1), (3, 1)),
            monomial(-1, (2, 1), (5, 1)),
        ),
        add(z[7], monomial(1, (1, 1), (2, 1))),
        add(
            z[8],
            monomial(6, (1, 1), (3, 1)),
            monomial(-3, (1, 1), (5, 1)),
            monomial(-3, (3, 1), (6, 1)),
        ),
        add(
            z[9],
            monomial(-1, (1, 1), (7, 1)),
            monomial(7, (2, 2)),
            monomial(-1, (3, 1), (4, 1)),
        ),
        add(z[10], monomial(1, (1, 1), (3, 1))),
        add(z[11], monomial(-half, (1, 2))),
    ]

    nonlinear_degrees = {
        sum(exponents)
        for index, component in enumerate(k)
        for exponents, coefficient in add(component, scale(z[index], -1)).items()
        if coefficient
    }
    assert nonlinear_degrees == {2, 3}
    assert max(max(map(sum, component), default=0) for component in k) == 3

    p = tuple(
        [Fraction(0), Fraction(0), Fraction(-1, 4)] + [Fraction(0)] * 9
    )
    q = (
        Fraction(1),
        Fraction(-3, 2),
        Fraction(13, 2),
        Fraction(-9, 4),
        Fraction(3),
        Fraction(3, 2),
        Fraction(99, 4),
        Fraction(3, 2),
        Fraction(-3, 4),
        Fraction(-45, 8),
        Fraction(-13, 2),
        Fraction(1, 2),
    )
    image_p = tuple(evaluate(component, p) for component in k)
    image_q = tuple(evaluate(component, q) for component in k)
    assert p != q and image_p == image_q == p

    jacobian = [
        [derivative(k[row], column) for column in range(N)]
        for row in range(N)
    ]
    determinant_value, memo_size = determinant(jacobian)
    assert determinant_value == {ZERO: Fraction(1)}

    print("PASS independent F12: explicit 12-component map has degree 3")
    print("PASS independent F12: two distinct rational points have image p")
    print(
        "PASS independent F12: exact sparse determinant is 1 "
        f"(memoized minors={memo_size})"
    )


if __name__ == "__main__":
    main()
