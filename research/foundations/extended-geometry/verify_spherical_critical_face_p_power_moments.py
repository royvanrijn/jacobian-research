#!/usr/bin/env python3
"""Replay a critical-face moment over the rationals."""

from collections import defaultdict
from fractions import Fraction


def multiply(left, right):
    product = defaultdict(int)
    for (left_angle, left_height), left_coefficient in left.items():
        for (right_angle, right_height), right_coefficient in right.items():
            exponent = (left_angle + right_angle, left_height + right_height)
            product[exponent] += left_coefficient * right_coefficient
    return {exponent: coefficient for exponent, coefficient in product.items()
            if coefficient}


def power(polynomial, exponent):
    result = {(0, 0): 1}
    for _ in range(exponent):
        result = multiply(result, polynomial)
    return result


def verify():
    d, B, n, s, prime, r = 2, 2, 37, 3, 5, 2
    b = d * B
    q = prime**r
    assert (s * q - 1) % (d * n) == 0 and s * b < n
    order = (s * q - 1) // n

    difference = {(d, 0): 1, (0, 1): -1}
    quadratic = power(difference, 2)
    quadratic[0, 0] = -1
    radial_and_extra_root = {(-b, B): -2, (d-b, B-1): 1}
    phase = multiply(multiply(power(difference, n - 2*b), power(quadratic, b)),
                     radial_and_extra_root)
    # These are x*t^36 and 2*y*t^32 after Laurent substitution.
    phase[1, 36] = 1
    phase[-1, 32] = 1
    phase[-1, 34] = -1
    for angle, height in phase:
        assert height >= 0
        assert n - 2*b <= height + Fraction(angle, d) <= n

    moment = sum(
        Fraction(coefficient, height + 1)
        for (angle, height), coefficient in power(phase, order).items()
        if angle == 0 and height % 2 == 0
    )
    normalized = s * q * moment
    assert normalized.denominator % prime != 0
    residue = normalized.numerator * pow(normalized.denominator, -1, prime) % prime
    assert order == 2 and residue == 1
    print(f"Order {order}: {s*q} times the moment is {residue} modulo {prime}")


if __name__ == "__main__":
    verify()
