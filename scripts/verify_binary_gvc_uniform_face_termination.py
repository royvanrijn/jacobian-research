#!/usr/bin/env python3
"""Exact regressions for uniform binary Hall and weighted-face termination.

The all-degree statements are proved in the accompanying note.  This script
checks the Hall formula on a broad finite window, exhausts the valuation
inequality on small unequal-weight lattice segments, and gives an exact
moment example in which a cancellation at the first moment is broken at a
prime dilation by the unique lower endpoint.
"""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
from math import comb, factorial


def hall_regression(limit: int = 12) -> None:
    for r in range(1, limit + 1):
        for d in range(r + 1, limit + 2):
            for mu in range(1, r + 1):
                for c in range(d + 1):
                    hall_fails = d - c < mu
                    formula = c >= d - mu + 1
                    assert hall_fails == formula


def line_points(u: int, v: int, weight: int) -> list[tuple[int, int]]:
    return [
        (a, b)
        for a in range(weight // u + 1)
        for b in range(weight // v + 1)
        if u * a + v * b == weight
    ]


def valuation_regression(limit_weight: int = 30) -> None:
    """Check the coefficient-independent inequalities in Theorem 3.1."""

    for u in range(1, 7):
        for v in range(1, 7):
            if u == v:
                continue
            for weight in range(1, limit_weight + 1):
                points = line_points(u, v, weight)
                if not points:
                    continue
                points.sort()
                for i in range(len(points)):
                    for j in range(i, len(points)):
                        for h in range(len(points)):
                            for k in range(h, len(points)):
                                lo = max(i, h)
                                hi = min(j, k)
                                if lo > hi:
                                    continue
                                intersection = points[lo : hi + 1]
                                alpha = min(intersection, key=lambda z: z[0] + z[1])
                                s = sum(alpha)
                                p = max(
                                    5,
                                    1
                                    + max(
                                        coordinate
                                        for point in points
                                        for coordinate in point
                                    ),
                                )
                                # It is enough to check all lattice points in
                                # the dilated intersection.  Non-p-multiples
                                # receive two Frobenius factors.
                                left_a, right_a = points[i], points[j]
                                left_b, right_b = points[h], points[k]
                                min_x = max(p * left_a[0], p * left_b[0])
                                max_x = min(p * right_a[0], p * right_b[0])
                                for x_exp in range(min_x, max_x + 1):
                                    remaining = p * weight - u * x_exp
                                    if remaining < 0 or remaining % v:
                                        continue
                                    y_exp = remaining // v
                                    if x_exp % p == 0 and y_exp % p == 0:
                                        beta = (x_exp // p, y_exp // p)
                                        lower_bound = sum(beta)
                                        if beta == alpha:
                                            assert lower_bound == s
                                        else:
                                            assert lower_bound >= s + 1
                                    else:
                                        factorial_value = x_exp // p + y_exp // p
                                        assert factorial_value + 2 >= s + 1


def multiply(
    left: dict[tuple[int, int], int],
    right: dict[tuple[int, int], int],
) -> dict[tuple[int, int], int]:
    answer: dict[tuple[int, int], int] = defaultdict(int)
    for (a, b), c in left.items():
        for (i, j), d in right.items():
            answer[a + i, b + j] += c * d
    return dict(answer)


def power(
    polynomial: dict[tuple[int, int], int], exponent: int
) -> dict[tuple[int, int], int]:
    answer = {(0, 0): 1}
    base = polynomial
    n = exponent
    while n:
        if n & 1:
            answer = multiply(answer, base)
        base = multiply(base, base)
        n //= 2
    return answer


def scalar_moment(
    operator: dict[tuple[int, int], int],
    polynomial: dict[tuple[int, int], int],
    exponent: int,
) -> int:
    op_power = power(operator, exponent)
    poly_power = power(polynomial, exponent)
    return sum(
        factorial(a) * factorial(b) * coefficient * poly_power.get((a, b), 0)
        for (a, b), coefficient in op_power.items()
    )


def translate(
    polynomial: dict[tuple[int, int], int], point: tuple[int, int]
) -> dict[tuple[int, int], int]:
    answer: dict[tuple[int, int], int] = defaultdict(int)
    zx, zy = point
    for (a, b), coefficient in polynomial.items():
        for i in range(a + 1):
            for j in range(b + 1):
                answer[i, j] += (
                    coefficient
                    * comb(a, i)
                    * comb(b, j)
                    * zx ** (a - i)
                    * zy ** (b - j)
                )
    return {exponent: coefficient for exponent, coefficient in answer.items() if coefficient}


def differential_value(
    operator: dict[tuple[int, int], int],
    polynomial: dict[tuple[int, int], int],
    exponent: int,
    point: tuple[int, int],
) -> int:
    op_power = power(operator, exponent)
    poly_power = power(polynomial, exponent)
    zx, zy = point
    answer = 0
    for (a, b), op_coefficient in op_power.items():
        for (i, j), poly_coefficient in poly_power.items():
            if i < a or j < b:
                continue
            falling_x = factorial(i) // factorial(i - a)
            falling_y = factorial(j) // factorial(j - b)
            answer += (
                op_coefficient
                * poly_coefficient
                * falling_x
                * falling_y
                * zx ** (i - a)
                * zy ** (j - b)
            )
    return answer


def ray_moment(
    operator: dict[tuple[int, int], int],
    shifted_polynomial: dict[tuple[int, int], int],
    delta: tuple[int, int],
    exponent: int,
) -> int:
    """Coefficient on the output ray exponent * delta.

    ``shifted_polynomial`` stores B after subtracting one copy of delta
    from every exponent.
    """

    op_power = power(operator, exponent)
    poly_power = power(shifted_polynomial, exponent)
    dx, dy = exponent * delta[0], exponent * delta[1]
    answer = 0
    for (a, b), coefficient in op_power.items():
        other = poly_power.get((a, b), 0)
        if not other:
            continue
        falling_x = factorial(a + dx) // factorial(dx)
        falling_y = factorial(b + dy) // factorial(dy)
        answer += falling_x * falling_y * coefficient * other
    return answer


def valuation(number: int, prime: int) -> int:
    assert number
    answer = 0
    while number % prime == 0:
        answer += 1
        number //= prime
    return answer


def prime_endpoint_example() -> None:
    # Weights (3,2), common weight 12.  The first moment cancels between
    # (4,0) and (2,3), but (4,0) is the unique least-ordinary-degree point
    # of the Newton-segment intersection.
    operator = {(0, 6): 1, (2, 3): 1, (4, 0): 1}
    polynomial = {(2, 3): -2, (4, 0): 1}
    assert scalar_moment(operator, polynomial, 1) == 0
    fifth = scalar_moment(operator, polynomial, 5)
    assert fifth != 0
    assert valuation(fifth, 5) == 4

    # The shifted-ray version has carrier delta=(1,2).  Its first ray
    # coefficient cancels, while the seventh again isolates (4,0).
    shifted_polynomial = {(2, 3): -1, (4, 0): 3}
    delta = (1, 2)
    assert ray_moment(operator, shifted_polynomial, delta, 1) == 0
    shifted_seventh = ray_moment(operator, shifted_polynomial, delta, 7)
    assert shifted_seventh != 0
    assert valuation(shifted_seventh, 7) == 4


def translated_multiradial_identity() -> None:
    operator = {(1, 0): 2, (0, 2): -1}
    polynomial = {(2, 1): 3, (0, 3): -2, (1, 0): 5}
    point = (2, -1)
    translated = translate(polynomial, point)
    for exponent in range(1, 5):
        assert scalar_moment(operator, translated, exponent) == differential_value(
            operator, polynomial, exponent, point
        )


def homogeneous_factorial_and_channel_regression() -> None:
    """Check the beta identity and the two compatibility warnings."""

    # C(U,V)=2 U^3-3 U^2 V+5 U V^2-7 V^3.
    polynomial = {(3, 0): 2, (2, 1): -3, (1, 2): 5, (0, 3): -7}
    for exponent in range(1, 6):
        expanded = power(polynomial, exponent)
        factorial_moment = sum(
            coefficient * factorial(a) * factorial(b)
            for (a, b), coefficient in expanded.items()
        )
        # Integral of U^a(1-U)^b is a!b!/(a+b+1)!.
        beta_moment = sum(
            Fraction(
                coefficient * factorial(a) * factorial(b),
                factorial(a + b + 1),
            )
            for (a, b), coefficient in expanded.items()
        )
        assert factorial_moment == factorial(3 * exponent + 1) * beta_moment

    # For G=S U+S^{-1}V, the S-constant term is zero in odd powers and
    # binom(2k,k) U^k V^k in power 2k.
    def channel_constant_term(exponent: int) -> dict[tuple[int, int], int]:
        if exponent % 2:
            return {}
        half = exponent // 2
        return {(half, half): comb(exponent, half)}

    first_channel = channel_constant_term(1)
    second_channel = channel_constant_term(2)
    assert not first_channel
    assert second_channel == {(1, 1): 2}
    assert second_channel != power(first_channel, 2)

    # The toric exponent map (a,b)->(a+b,b) changes factorial weights by
    # a vector-dependent factor even on one ordinary-degree line.
    def blowup_ratio(a: int, b: int) -> Fraction:
        old = factorial(a) * factorial(b)
        new = factorial(a + b) * factorial(b)
        return Fraction(new, old)

    assert blowup_ratio(2, 0) == 1
    assert blowup_ratio(1, 1) == 2


def main() -> None:
    hall_regression()
    valuation_regression()
    prime_endpoint_example()
    translated_multiradial_identity()
    homogeneous_factorial_and_channel_regression()
    print("PASS: uniform Hall, weighted-face, and beta-channel regressions")


if __name__ == "__main__":
    main()
