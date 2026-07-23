#!/usr/bin/env python3
"""Exact checks for toric Cox--Keller permutation maps."""

from __future__ import annotations

from itertools import product
from math import gcd

import sympy as sp


# Symbolic separated determinant in rank three.
x1, x2, x3, z1, z2, z3 = sp.symbols(
    "x1 x2 x3 z1 z2 z3", nonzero=True
)
n1, n2, n3 = 3, 5, 7
separated = (
    x1**n1,
    x2**n2,
    x3**n3,
    z1 / (n1 * x1 ** (n1 - 1)),
    z2 / (n2 * x2 ** (n2 - 1)),
    z3 / (n3 * x3 ** (n3 - 1)),
)
variables = (x1, x2, x3, z1, z2, z3)
separated_jacobian = sp.factor(
    sp.Matrix(separated).jacobian(variables).det()
)
assert separated_jacobian == 1


# Compressed determinant for the same three independent unit directions.
z = sp.symbols("z")
horizontal_jacobian = (
    n1 * n2 * n3
    * x1 ** (n1 - 1)
    * x2 ** (n2 - 1)
    * x3 ** (n3 - 1)
)
compressed = (
    x1**n1,
    x2**n2,
    x3**n3,
    z / horizontal_jacobian,
)
compressed_jacobian = sp.factor(
    sp.Matrix(compressed).jacobian((x1, x2, x3, z)).det()
)
assert compressed_jacobian == 1
print("PASS: separated and compressed toric Cox maps have determinant one")


def toric_map(
    point: tuple[int, ...],
    exponents: tuple[int, ...],
    prime: int,
) -> tuple[int, ...]:
    rank = len(exponents)
    xs = point[:rank]
    zs = point[rank:]
    horizontal = tuple(
        pow(value, exponent, prime)
        for value, exponent in zip(xs, exponents)
    )
    vertical = tuple(
        (
            value
            * pow(
                (
                    exponent
                    * pow(x_value, exponent - 1, prime)
                )
                % prime,
                -1,
                prime,
            )
        )
        % prime
        for x_value, value, exponent in zip(xs, zs, exponents)
    )
    return (*horizontal, *vertical)


def enumerate_map(
    exponents: tuple[int, ...], prime: int
) -> tuple[int, int]:
    rank = len(exponents)
    source = product(
        *([range(1, prime)] * rank),
        *([range(prime)] * rank),
    )
    images = {
        toric_map(point, exponents, prime)
        for point in source
    }
    expected = (prime - 1) ** rank * prime**rank
    return len(images), expected


for exponents, prime in (
    ((3,), 5),
    ((3, 5), 17),
):
    assert all(gcd(exponent, prime - 1) == 1 for exponent in exponents)
    image_size, expected_size = enumerate_map(exponents, prime)
    assert image_size == expected_size

# Failure example: cubes are not a permutation of F_7^*.
assert gcd(3, 7 - 1) != 1
image_size, expected_size = enumerate_map((3,), 7)
assert image_size < expected_size
print("PASS: finite-field bijectivity is exactly the power-map gcd criterion")
print("PASS toric Cox permutation maps")
