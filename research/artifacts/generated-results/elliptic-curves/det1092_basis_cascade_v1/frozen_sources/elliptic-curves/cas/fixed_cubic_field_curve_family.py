#!/usr/bin/env python3
"""Exact algebra for the fixed-cubic-field varying-curve experiment.

This module is deliberately Sage-free so that the defining identities and
the GF(2) kernel calculation remain part of the ordinary Python test suite.
The arithmetic runner is ``run_fixed_cubic_field_curve_family.sage``.
"""

from __future__ import annotations

from fractions import Fraction
from math import gcd
from typing import Iterable, Sequence


Q = Fraction


def _q(value: object) -> Fraction:
    """Coerce Python, string, or Sage-style rational values exactly."""

    if isinstance(value, Fraction):
        return value
    return Fraction(str(value))


def fixed_field_cubic_coefficients(
    A: int | Fraction, B: int | Fraction, u: int | Fraction
) -> tuple[Fraction, Fraction, Fraction, Fraction]:
    """Return ascending coefficients of ``Norm(x-(theta+u*theta^2))``.

    Here ``theta^3 + A*theta + B = 0``.
    """

    A, B, u = map(_q, (A, B, u))
    return (
        B + A * B * u**2 - B**2 * u**3,
        A + 3 * B * u + A**2 * u**2,
        2 * A * u,
        Q(1),
    )


def cubic_discriminant(
    coefficients: Sequence[int | Fraction],
) -> Fraction:
    """Discriminant of a monic cubic in ascending coefficient order."""

    if len(coefficients) != 4 or _q(coefficients[3]) != 1:
        raise ValueError("expected four ascending coefficients of a monic cubic")
    c, b, a, _ = map(_q, coefficients)
    return a * a * b * b - 4 * b**3 - 4 * a**3 * c - 27 * c * c + 18 * a * b * c


def discriminant_multiplier(
    A: int | Fraction, B: int | Fraction, u: int | Fraction
) -> Fraction:
    A, B, u = map(_q, (A, B, u))
    return 1 + A * u**2 + B * u**3


def field_multiply(
    left: Sequence[object],
    right: Sequence[object],
    A: object,
    B: object,
) -> tuple[object, object, object]:
    """Multiply power-basis rows modulo ``theta^3+A*theta+B``.

    The entries need only support ordinary ring arithmetic; this is also used
    by tests with exact rational values.
    """

    if len(left) != 3 or len(right) != 3:
        raise ValueError("cubic power-basis rows must have length three")
    raw = [0, 0, 0, 0, 0]
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            raw[i + j] += a * b
    # theta^4=-A*theta^2-B*theta and theta^3=-A*theta-B.
    raw[2] -= A * raw[4]
    raw[1] -= B * raw[4]
    raw[1] -= A * raw[3]
    raw[0] -= B * raw[3]
    return tuple(raw[:3])


def field_product(
    rows: Iterable[Sequence[object]], A: object, B: object
) -> tuple[object, object, object]:
    answer: tuple[object, object, object] = (1, 0, 0)
    for row in rows:
        answer = field_multiply(answer, row, A, B)
    return answer


def inverse_theta_coefficients(
    A: int | Fraction, B: int | Fraction, u: int | Fraction
) -> tuple[Fraction, Fraction, Fraction]:
    """Express theta in the basis ``1, alpha_u, alpha_u^2``."""

    A, B, u = map(_q, (A, B, u))
    D = discriminant_multiplier(A, B, u)
    if not D:
        raise ZeroDivisionError("alpha_u does not generate the cubic algebra")
    return -2 * B * u**2 / D, (1 - A * u**2) / D, -u / D


def covering_values(
    beta: Sequence[object],
    gamma: Sequence[object],
    homogenizing_coordinate: object,
    A: object,
    B: object,
    u: object,
) -> tuple[object, object, object]:
    """Evaluate the two covering quadrics and the recovered x numerator.

    For ``gamma=a+b*theta+c*theta^2`` and projective coordinate ``d``, the
    equations are the theta coefficient of ``beta*gamma^2`` plus ``d^2`` and
    the theta^2 coefficient plus ``u*d^2``.  If both vanish and ``d != 0``,
    the x-coordinate is the returned constant coefficient divided by ``d^2``.
    """

    square = field_multiply(gamma, gamma, A, B)
    product = field_multiply(beta, square, A, B)
    d2 = homogenizing_coordinate * homogenizing_coordinate
    return product[1] + d2, product[2] + u * d2, product[0]


def f2_rank(rows: Iterable[Sequence[int]]) -> int:
    pivots: dict[int, int] = {}
    for row in rows:
        value = sum((int(bit) & 1) << index for index, bit in enumerate(row))
        while value:
            pivot = value.bit_length() - 1
            if pivot in pivots:
                value ^= pivots[pivot]
            else:
                pivots[pivot] = value
                break
    return len(pivots)


def f2_kernel_masks(rows: Sequence[Sequence[int]]) -> list[int]:
    """Basis masks for dependencies among the supplied GF(2) row vectors."""

    pivots: dict[int, tuple[int, int]] = {}
    dependencies: list[int] = []
    for index, row in enumerate(rows):
        value = sum((int(bit) & 1) << column for column, bit in enumerate(row))
        provenance = 1 << index
        while value:
            pivot = value.bit_length() - 1
            if pivot in pivots:
                value ^= pivots[pivot][0]
                provenance ^= pivots[pivot][1]
            else:
                pivots[pivot] = (value, provenance)
                break
        if value == 0:
            dependencies.append(provenance)
    return dependencies


def mask_indices(mask: int, width: int) -> list[int]:
    return [index for index in range(width) if (mask >> index) & 1]


def bounded_integer_parameters(bound: int) -> tuple[Fraction, ...]:
    if bound < 0:
        raise ValueError("parameter bound must be nonnegative")
    return tuple(Q(value) for value in range(-bound, bound + 1))


def bounded_rational_parameters(height: int) -> tuple[Fraction, ...]:
    """All reduced rationals p/q with max(|p|,q) <= height."""

    if height < 1:
        raise ValueError("rational height must be positive")
    values = set()
    for denominator in range(1, height + 1):
        for numerator in range(-height, height + 1):
            if gcd(abs(numerator), denominator) == 1:
                values.add(Q(numerator, denominator))
    return tuple(sorted(values))


__all__ = [
    "bounded_integer_parameters",
    "bounded_rational_parameters",
    "covering_values",
    "cubic_discriminant",
    "discriminant_multiplier",
    "f2_kernel_masks",
    "f2_rank",
    "field_multiply",
    "field_product",
    "fixed_field_cubic_coefficients",
    "inverse_theta_coefficients",
    "mask_indices",
]
