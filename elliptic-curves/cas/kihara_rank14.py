#!/usr/bin/env python3
"""Exact replay of Kihara's explicit rank-at-least-14 family.

The primary source is Shoichi Kihara, *On an elliptic curve over
Q(t) of rank >= 14*, Proc. Japan Acad. Ser. A 77 (2001), 50--51,
doi:10.3792/pjaa.77.50.  The paper constructs a quartic ``y^2=r(x)``
with fifteen displayed points.  It takes ``P15`` as the origin and reports
that ``P1,...,P14`` are independent, using a numerical canonical-height
determinant after specialization at ``t=2``.

This module checks the printed polynomial construction and all printed
point abscissae over exact rationals.  It deliberately does not promote the
paper's numerical height computation to a new repository-local independence
certificate.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
import json
from math import isqrt
from typing import Iterable, Sequence


Q = Fraction

PRIMARY_SOURCE_DOI = "https://doi.org/10.3792/pjaa.77.50"
PRIMARY_SOURCE_PDF = (
    "https://projecteuclid.org/download/pdf_1/euclid.pja/1148393079"
)
PUBLISHED_RANK_LOWER_BOUND = 14
PUBLISHED_INDEPENDENCE_SPECIALIZATION = Q(2)
PUBLISHED_HEIGHT_DETERMINANT_APPROX = "221792776617402574.10"


def _multiply(
    left: Sequence[Fraction], right: Sequence[Fraction]
) -> tuple[Fraction, ...]:
    answer = [Q(0)] * (len(left) + len(right) - 1)
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            answer[left_index + right_index] += Q(left_value) * Q(right_value)
    return tuple(answer)


def _evaluate(coefficients: Sequence[Fraction], value: Fraction) -> Fraction:
    answer = Q(0)
    value = Q(value)
    for coefficient in reversed(coefficients):
        answer = answer * value + Q(coefficient)
    return answer


def _monic_polynomial_from_roots(
    roots: Iterable[Fraction],
) -> tuple[Fraction, ...]:
    answer = (Q(1),)
    for root in roots:
        answer = _multiply(answer, (-Q(root), Q(1)))
    return answer


def _square_approximant(
    monic_degree_twelve: Sequence[Fraction],
) -> tuple[Fraction, ...]:
    """Return the unique monic sextic matching degrees 12 through 6."""

    product = tuple(Q(value) for value in monic_degree_twelve)
    if len(product) != 13 or product[-1] != 1:
        raise ValueError("a monic degree-twelve polynomial is required")
    approximant = [Q(0)] * 7
    approximant[6] = Q(1)
    for index in range(5, -1, -1):
        square = _multiply(approximant, approximant)
        degree = 6 + index
        approximant[index] = (product[degree] - square[degree]) / 2
    return tuple(approximant)


def _rational_square_root(value: Fraction) -> Fraction:
    value = Q(value)
    if value < 0:
        raise ValueError("the rational value is negative")
    numerator = isqrt(value.numerator)
    denominator = isqrt(value.denominator)
    if numerator**2 != value.numerator or denominator**2 != value.denominator:
        raise ValueError("the rational value is not a square")
    return Q(numerator, denominator)


def specialized_parameters(
    parameter_t: Fraction,
) -> tuple[Fraction, Fraction, Fraction]:
    """Return Kihara's printed ``(p,q,u)`` base change."""

    parameter_t = Q(parameter_t)
    if parameter_t == 0:
        raise ValueError("Kihara's printed base change has a pole at t=0")
    t = parameter_t
    p = t**2 * (8 + 3 * t**2)
    q = -6 * (2 + t**2) * (4 + t**2)
    u = (
        4
        * (2 + t**2)
        * (2304 + 2400 * t**2 + 928 * t**4 + 150 * t**6 + 9 * t**8)
        * (1152 + 1632 * t**2 + 860 * t**4 + 201 * t**6 + 18 * t**8)
        / t
    )
    return p, q, u


def kihara_a_values(p: Fraction, q: Fraction) -> tuple[Fraction, ...]:
    """Return the six printed values ``a1,...,a6``."""

    p, q = Q(p), Q(q)
    return (
        Q(0),
        (2 * p**2 + p * q + 2 * q**2) ** 2,
        2 * (p + q) ** 2 * (2 * p**2 + p * q + q**2),
        q**2 * (4 * p**2 - p * q + 4 * q**2),
        p * (2 * p - q) * (2 * p**2 + 4 * p * q + 5 * q**2),
        4 * p**4 + 8 * p**3 * q + 9 * p**2 * q**2 - 2 * p * q**3 + 2 * q**4,
    )


@dataclass(frozen=True)
class KiharaSpecialization:
    parameter_t: Fraction
    p: Fraction
    q: Fraction
    u: Fraction
    a_values: tuple[Fraction, ...]
    b_values: tuple[Fraction, ...]
    product_coefficients: tuple[Fraction, ...]
    approximant_coefficients: tuple[Fraction, ...]
    quartic_coefficients: tuple[Fraction, ...]


def kihara_specialization(parameter_t: Fraction) -> KiharaSpecialization:
    """Construct the specialized quartic, checking ``F=G^2-r`` exactly."""

    parameter_t = Q(parameter_t)
    p, q, u = specialized_parameters(parameter_t)
    a_values = kihara_a_values(p, q)
    b_values = tuple(u + value for value in a_values) + tuple(
        -u + value for value in a_values
    )
    product = _monic_polynomial_from_roots(b_values)
    approximant = _square_approximant(product)
    square = _multiply(approximant, approximant)
    remainder = tuple(left - right for left, right in zip(square, product))
    if any(remainder[index] for index in range(5, len(remainder))):
        raise AssertionError("Kihara's remainder did not reduce to a quartic")
    quartic = remainder[:5]
    if not quartic[-1]:
        raise AssertionError("Kihara's remainder lost its quartic term")
    return KiharaSpecialization(
        parameter_t=parameter_t,
        p=p,
        q=q,
        u=u,
        a_values=a_values,
        b_values=b_values,
        product_coefficients=product,
        approximant_coefficients=approximant,
        quartic_coefficients=quartic,
    )


def published_extra_abscissae(
    specialization: KiharaSpecialization,
) -> tuple[Fraction, Fraction, Fraction]:
    """Return the printed abscissae of ``P13,P14,P15``."""

    p, q, u, t = (
        specialization.p,
        specialization.q,
        specialization.u,
        specialization.parameter_t,
    )
    denominator_13 = 2 * p**2 + 2 * p * q + 3 * q**2
    numerator_13 = 2 * p**2 + 4 * p * q + 5 * q**2
    degree_six = (
        8 * p**6
        + 28 * p**5 * q
        + 58 * p**4 * q**2
        + 69 * p**3 * q**3
        + 76 * p**2 * q**4
        + 40 * p * q**5
        + 22 * q**6
    )
    x_13 = numerator_13 * u / denominator_13 + degree_six / denominator_13

    common = 1152 + 1632 * t**2 + 860 * t**4 + 201 * t**6 + 18 * t**8
    polynomial_14 = (
        10616832
        - 18579456 * t
        + 33619968 * t**2
        - 51535872 * t**3
        + 45895680 * t**4
        - 61848576 * t**5
        + 35397888 * t**6
        - 41945856 * t**7
        + 16968640 * t**8
        - 17591104 * t**9
        + 5232272 * t**10
        - 4675248 * t**11
        + 1035180 * t**12
        - 769824 * t**13
        + 126252 * t**14
        - 71874 * t**15
        + 8559 * t**16
        - 2916 * t**17
        + 243 * t**18
    )
    denominator_14 = t * (
        2304 + 3168 * t**2 + 1580 * t**4 + 339 * t**6 + 27 * t**8
    )
    x_14 = -4 * common * polynomial_14 / denominator_14
    x_15 = (
        4
        * (-48 + 24 * t - 34 * t**2 + 16 * t**3 - 6 * t**4 + 3 * t**5)
        * (96 + 80 * t**2 + 4 * t**3 + 18 * t**4 + 3 * t**5)
        * common
        / t
    )
    return x_13, x_14, x_15


def known_quartic_points(
    parameter_t: Fraction,
) -> tuple[tuple[Fraction, Fraction], ...]:
    """Return exact specializations of the printed points ``P1,...,P15``."""

    specialization = kihara_specialization(parameter_t)
    visible = tuple(
        (x_value, _evaluate(specialization.approximant_coefficients, x_value))
        for x_value in specialization.b_values
    )
    extra = tuple(
        (
            x_value,
            _rational_square_root(
                _evaluate(specialization.quartic_coefficients, x_value)
            ),
        )
        for x_value in published_extra_abscissae(specialization)
    )
    points = visible + extra
    if any(
        y_value**2
        != _evaluate(specialization.quartic_coefficients, x_value)
        for x_value, y_value in points
    ):
        raise AssertionError("a printed Kihara point failed the quartic exactly")
    return points


def binary_invariants(
    coefficients: Sequence[Fraction],
) -> tuple[Fraction, Fraction]:
    """Return the classical ``I,J`` of an ascending binary quartic."""

    if len(coefficients) != 5:
        raise ValueError("a binary quartic has five coefficients")
    e, d, c, b, a = (Q(value) for value in coefficients)
    invariant_i = 12 * a * e - 3 * b * d + c**2
    invariant_j = (
        72 * a * c * e
        + 9 * b * c * d
        - 27 * a * d**2
        - 27 * b**2 * e
        - 2 * c**3
    )
    return invariant_i, invariant_j


def short_jacobian_coefficients(
    parameter_t: Fraction,
) -> tuple[Fraction, ...]:
    """Return ``[0,0,0,-27I,-27J]`` for the specialized quartic."""

    invariant_i, invariant_j = binary_invariants(
        kihara_specialization(parameter_t).quartic_coefficients
    )
    return Q(0), Q(0), Q(0), -27 * invariant_i, -27 * invariant_j


def verify_rational_specialization(parameter_t: Fraction) -> dict[str, object]:
    """Return a compact exact verification record for one rational ``t``."""

    specialization = kihara_specialization(parameter_t)
    points = known_quartic_points(parameter_t)
    invariant_i, invariant_j = binary_invariants(
        specialization.quartic_coefficients
    )
    discriminant = (4 * invariant_i**3 - invariant_j**2) / 27
    return {
        "parameter_t": str(Q(parameter_t)),
        "remainder_degree": 4,
        "displayed_point_count": len(points),
        "distinct_point_count": len(set(points)),
        "distinct_abscissa_count": len({point[0] for point in points}),
        "all_points_exact": all(
            point[1] ** 2
            == _evaluate(specialization.quartic_coefficients, point[0])
            for point in points
        ),
        "quartic_discriminant_nonzero": bool(discriminant),
        "quartic_coefficient_numerator_digits": [
            len(str(abs(value.numerator)))
            for value in specialization.quartic_coefficients
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--t", default="2", help="rational specialization")
    arguments = parser.parse_args()
    print(json.dumps(verify_rational_specialization(Q(arguments.t)), indent=2))


if __name__ == "__main__":
    main()
