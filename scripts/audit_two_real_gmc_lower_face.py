#!/usr/bin/env python3
"""Independent exact regression for the lower-face proof of GMC(2)."""

from __future__ import annotations

import math
from fractions import Fraction


LaurentRadial = dict[tuple[int, int], int]
Radial = dict[int, int]


def multiply(left: LaurentRadial, right: LaurentRadial) -> LaurentRadial:
    answer: LaurentRadial = {}
    for (tw, uw), a in left.items():
        for (tv, uv), b in right.items():
            key = (tw + tv, uw + uv)
            answer[key] = answer.get(key, 0) + a * b
    return {key: value for key, value in answer.items() if value}


def power(polynomial: LaurentRadial, exponent: int) -> LaurentRadial:
    answer: LaurentRadial = {(0, 0): 1}
    base = polynomial
    current = exponent
    while current:
        if current & 1:
            answer = multiply(answer, base)
        base = multiply(base, base)
        current //= 2
    return answer


def constant_term(polynomial: LaurentRadial) -> Radial:
    return {
        radial: coefficient
        for (weight, radial), coefficient in polynomial.items()
        if weight == 0
    }


def factorial_functional(polynomial: Radial) -> int:
    return sum(
        coefficient * math.factorial(degree)
        for degree, coefficient in polynomial.items()
    )


def lower_face(polynomial: LaurentRadial) -> tuple[Fraction, dict[int, int]]:
    leading: dict[int, tuple[int, int]] = {}
    for (weight, radial), coefficient in polynomial.items():
        if weight not in leading or radial < leading[weight][0]:
            leading[weight] = (radial, coefficient)

    candidates: list[tuple[Fraction, Fraction]] = []
    if 0 in leading:
        candidates.append((Fraction(leading[0][0]), Fraction(0)))
    for negative in (weight for weight in leading if weight < 0):
        for positive in (weight for weight in leading if weight > 0):
            left = leading[negative][0]
            right = leading[positive][0]
            theta = Fraction(right - left, positive - negative)
            rho = Fraction(left) - theta * negative
            candidates.append((rho, theta))
    rho, theta = min(candidates, key=lambda pair: pair[0])
    face = {
        weight: coefficient
        for weight, (radial, coefficient) in leading.items()
        if Fraction(radial) == rho + theta * weight
    }
    assert min(face) <= 0 <= max(face)
    return rho, face


def laurent_univariate_power_ct(
    polynomial: dict[int, int], exponent: int
) -> int:
    lifted = {(weight, 0): coefficient for weight, coefficient in polynomial.items()}
    return constant_term(power(lifted, exponent)).get(0, 0)


def first_nonzero_face_moment(face: dict[int, int]) -> tuple[int, int]:
    for exponent in range(1, 25):
        value = laurent_univariate_power_ct(face, exponent)
        if value:
            return exponent, value
    raise AssertionError("sample exceeded the constant-term search bound")


def check_sample(polynomial: LaurentRadial, prime: int) -> None:
    rho, face = lower_face(polynomial)
    exponent, face_ct = first_nonzero_face_moment(face)
    baseline = rho * exponent
    assert baseline.denominator == 1
    radial_degree = baseline.numerator

    first = constant_term(power(polynomial, exponent))
    assert min(first) == radial_degree
    assert first[radial_degree] == face_ct

    dilated = constant_term(power(polynomial, exponent * prime))
    assert min(dilated) >= radial_degree * prime

    divisor = math.factorial(radial_degree * prime)
    value = factorial_functional(dilated)
    assert value % divisor == 0
    assert value // divisor % prime == pow(face_ct, prime, prime)

    frobenius = {
        degree * prime: pow(coefficient, prime, prime)
        for degree, coefficient in first.items()
    }
    all_degrees = set(dilated) | set(frobenius)
    assert all(
        (dilated.get(degree, 0) - frobenius.get(degree, 0)) % prime == 0
        for degree in all_degrees
    )


def main() -> None:
    samples: tuple[LaurentRadial, ...] = (
        {
            (1, 0): 1,
            (1, 1): 2,
            (0, 2): 3,
            (-1, 1): 5,
            (-2, 2): 7,
        },
        {
            (5, 1): 1,
            (5, 2): 2,
            (-2, 1): 3,
            (-3, 1): 4,
        },
        {
            (2, 2): 1,
            (1, 1): 2,
            (-1, 1): 3,
            (-2, 2): 5,
        },
        {
            (3, 4): 2,
            (0, 1): 7,
            (-4, 3): 11,
        },
    )
    for sample in samples:
        check_sample(sample, 5)
        check_sample(sample, 7)

    print("PASS lower face: supporting-line minima and exposed coefficients")
    print("PASS lower face: Frobenius constant-term dilation")
    print("PASS lower face: normalized factorial isolation")


if __name__ == "__main__":
    main()
