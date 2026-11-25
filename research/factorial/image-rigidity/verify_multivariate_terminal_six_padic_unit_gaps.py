#!/usr/bin/env python3
"""Replay the finite residue tables for terminal six at primes five and eleven."""

from itertools import product
from math import comb, factorial


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def remainder(dividend, divisor, prime):
    coefficients = list(dividend)
    inverse = pow(divisor[-1], -1, prime)
    while len(coefficients) >= len(divisor):
        offset = len(coefficients) - len(divisor)
        quotient = coefficients[-1] * inverse
        for index, coefficient in enumerate(divisor):
            position = offset + index
            coefficients[position] = (coefficients[position] - quotient * coefficient) % prime
        while coefficients and coefficients[-1] == 0:
            coefficients.pop()
    return coefficients


def gcd_degree(first, second, prime):
    while second:
        first, second = second, remainder(first, second, prime)
    return len(first) - 1


def moment_coefficient(order, index, increasing, decreasing, modulus):
    return (
        comb(order, index)
        * pow(factorial(index), increasing, modulus)
        * pow(factorial(order - index), decreasing, modulus)
    ) % modulus


def divide_by_five(coefficient):
    require(coefficient % 5 == 0, ("five is not a divisor", coefficient))
    return coefficient // 5


def verify_five_faces():
    degree_by_class = {}
    for increasing, decreasing in product(range(4), repeat=2):
        lower = [
            divide_by_five(moment_coefficient(5, k, increasing, decreasing, 25))
            for k in range(1, 5)
        ]
        upper = [
            divide_by_five(moment_coefficient(6, k, increasing, decreasing, 25))
            for k in range(2, 5)
        ]
        degree_by_class[increasing, decreasing] = gcd_degree(lower, upper, 5)
        if (increasing, decreasing) == (2, 2):
            require(lower == [1, 3, 3, 1], ("quintic exceptional face", lower))
            require(upper == [2, 4, 2], ("sextic exceptional face", upper))
    require(
        {key: degree for key, degree in degree_by_class.items() if degree} == {(2, 2): 2},
        degree_by_class,
    )
    return len(degree_by_class)


def translated_jet(order, increasing, decreasing, top_index):
    coefficients = [
        divide_by_five(moment_coefficient(order, k, increasing, decreasing, 125))
        for k in range(order + 1)
    ]
    return [
        sum(
            (-1) ** (k - index) * comb(k, index) * coefficients[k]
            for k in range(index, order + 1)
        ) % 25
        for index in range(top_index + 1)
    ]


def verify_five_lifts():
    representatives = (2, 6, 10, 14, 18, 22)
    possible = []
    for increasing, decreasing in product(representatives, repeat=2):
        constant, linear, quadratic = translated_jet(6, increasing, decreasing, 2)
        require((constant, linear % 5, quadratic % 5) == (15, 0, 2),
                ("sextic lift", increasing, decreasing, constant, linear, quadratic))
        constant, linear, quadratic, cubic = translated_jet(5, increasing, decreasing, 3)
        require((constant % 5, linear % 5, quadratic % 5, cubic % 5) == (0, 0, 0, 4),
                ("quintic lift", increasing, decreasing, constant, linear, quadratic, cubic))
        if (constant, linear) == (0, 5):
            possible.append((increasing, decreasing))
    require(possible == [(18, 18)], possible)
    return len(representatives) ** 2


def verify_eleven_faces():
    degree_by_class = {}
    for increasing, decreasing in product(range(10), repeat=2):
        first, second = (
            [moment_coefficient(order, k, increasing, decreasing, 11)
             for k in range(order + 1)]
            for order in (5, 6)
        )
        degree_by_class[increasing, decreasing] = gcd_degree(first, second, 11)
    require(
        {key: degree for key, degree in degree_by_class.items() if degree}
        == {(0, 0): 5, (2, 5): 1, (5, 2): 1},
        degree_by_class,
    )
    return len(degree_by_class)


if __name__ == "__main__":
    print(f"Terminal six modular tables passed: {verify_five_faces()} initial 5-adic cases, "
          f"{verify_five_lifts()} lifted 5-adic cases, {verify_eleven_faces()} cases modulo 11.")
