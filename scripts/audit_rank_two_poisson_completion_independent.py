#!/usr/bin/env python3
"""Dependency-free exact audit of the rank-two Poisson completion.

This intentionally shares no SymPy code with
``verify_rank_two_poisson_completion.py``.  It rebuilds the compact formulas
in a small sparse polynomial ring over ``fractions.Fraction`` and checks all
six brackets, the Jacobian determinant, and the rational collision.
"""

from fractions import Fraction
from itertools import permutations


Exponent = tuple[int, int, int, int]
Polynomial = dict[Exponent, Fraction]
ZERO: Exponent = (0, 0, 0, 0)


def clean(poly: Polynomial) -> Polynomial:
    return {exponent: coefficient for exponent, coefficient in poly.items() if coefficient}


def constant(value: int | Fraction) -> Polynomial:
    coefficient = Fraction(value)
    return {} if coefficient == 0 else {ZERO: coefficient}


def variable(index: int) -> Polynomial:
    exponent = [0, 0, 0, 0]
    exponent[index] = 1
    return {tuple(exponent): Fraction(1)}


def add(*polynomials: Polynomial) -> Polynomial:
    result: Polynomial = {}
    for polynomial in polynomials:
        for exponent, coefficient in polynomial.items():
            result[exponent] = result.get(exponent, Fraction(0)) + coefficient
    return clean(result)


def scale(value: int | Fraction, polynomial: Polynomial) -> Polynomial:
    coefficient = Fraction(value)
    return clean({exponent: coefficient * entry for exponent, entry in polynomial.items()})


def subtract(left: Polynomial, right: Polynomial) -> Polynomial:
    return add(left, scale(-1, right))


def multiply(*polynomials: Polynomial) -> Polynomial:
    result = constant(1)
    for polynomial in polynomials:
        product: Polynomial = {}
        for left_exponent, left_coefficient in result.items():
            for right_exponent, right_coefficient in polynomial.items():
                exponent = tuple(
                    left + right
                    for left, right in zip(left_exponent, right_exponent, strict=True)
                )
                product[exponent] = (
                    product.get(exponent, Fraction(0))
                    + left_coefficient * right_coefficient
                )
        result = clean(product)
    return result


def power(polynomial: Polynomial, exponent: int) -> Polynomial:
    assert exponent >= 0
    result = constant(1)
    factor = polynomial
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = multiply(result, factor)
        factor = multiply(factor, factor)
        remaining //= 2
    return result


def derivative(polynomial: Polynomial, index: int) -> Polynomial:
    result: Polynomial = {}
    for exponent, coefficient in polynomial.items():
        degree = exponent[index]
        if degree:
            derived = list(exponent)
            derived[index] -= 1
            result[tuple(derived)] = coefficient * degree
    return result


def poisson(left: Polynomial, right: Polynomial) -> Polynomial:
    """{p,x}={z,q}=1 in variable order (x,q,p,z)."""

    return add(
        multiply(derivative(left, 2), derivative(right, 0)),
        scale(-1, multiply(derivative(left, 0), derivative(right, 2))),
        multiply(derivative(left, 3), derivative(right, 1)),
        scale(-1, multiply(derivative(left, 1), derivative(right, 3))),
    )


def determinant(matrix: list[list[Polynomial]]) -> Polynomial:
    result: Polynomial = {}
    for permutation in permutations(range(4)):
        inversions = sum(
            permutation[i] > permutation[j]
            for i in range(4)
            for j in range(i + 1, 4)
        )
        term = multiply(*(matrix[row][permutation[row]] for row in range(4)))
        result = add(result, scale(-1 if inversions % 2 else 1, term))
    return result


def evaluate(polynomial: Polynomial, point: tuple[Fraction, ...]) -> Fraction:
    return sum(
        coefficient
        * point[0] ** exponent[0]
        * point[1] ** exponent[1]
        * point[2] ** exponent[2]
        * point[3] ** exponent[3]
        for exponent, coefficient in polynomial.items()
    )


one = constant(1)
x, q, p, z = (variable(index) for index in range(4))
X = x
Q = q
Z = add(
    scale(3, multiply(power(x, 2), p)),
    multiply(add(constant(2), scale(-6, multiply(x, q))), z),
)
E = add(
    scale(Fraction(1, 2), multiply(add(one, scale(3, multiply(x, q))), p)),
    scale(-3, multiply(power(q, 2), z)),
)
W = add(Z, scale(-9, power(Q, 2)))
Y = add(Q, scale(Fraction(-1, 3), multiply(X, W)))
U = add(one, multiply(X, Y))

R = add(scale(2, X), scale(-3, multiply(power(X, 2), Y)), scale(-1, multiply(power(X, 3), W)))
S = scale(
    Fraction(1, 2),
    add(
        multiply(power(U, 3), W),
        multiply(power(Y, 2), U, add(constant(4), scale(3, multiply(X, Y)))),
    ),
)
T = add(
    Y,
    scale(3, multiply(X, power(U, 2), W)),
    scale(3, multiply(X, power(Y, 2), add(constant(4), scale(3, multiply(X, Y))))),
)

completion_inside = add(
    scale(10, multiply(power(W, 3), power(X, 2))),
    scale(90, multiply(power(W, 2), X, Y)),
    scale(20, power(W, 2)),
    scale(-18, multiply(W, power(X, 3), power(Y, 5))),
    scale(-90, multiply(W, power(X, 2), power(Y, 4))),
    scale(-180, multiply(W, X, power(Y, 3))),
    scale(90, multiply(W, power(Y, 2))),
    scale(-54, multiply(power(X, 2), power(Y, 6))),
    scale(-234, multiply(X, power(Y, 5))),
    scale(-375, power(Y, 4)),
)
D = add(E, scale(Fraction(-1, 60), completion_inside))

assert R == add(scale(2, x), scale(-3, multiply(power(x, 2), q)))
assert poisson(D, R) == one
assert poisson(S, T) == one
assert poisson(R, S) == {}
assert poisson(R, T) == {}
assert poisson(D, S) == {}
assert poisson(D, T) == {}

outputs = (R, T, D, S)
jacobian = [[derivative(output, index) for index in range(4)] for output in outputs]
assert determinant(jacobian) == one

collision_points = (
    (Fraction(0), Fraction(0), Fraction(1, 24), Fraction(-1, 8)),
    (Fraction(1), Fraction(2, 3), Fraction(247, 96), Fraction(-89, 64)),
    (Fraction(-1), Fraction(-2, 3), Fraction(247, 96), Fraction(-89, 64)),
)
target = (Fraction(0), Fraction(0), Fraction(0), Fraction(-1, 8))
for point in collision_points:
    assert tuple(evaluate(output, point) for output in outputs) == target
assert len(set(collision_points)) == 3

assert tuple(len(output) for output in outputs) == (2, 22, 139, 47)

print("PASS independent sparse audit: all six Poisson brackets and determinant one")
print("PASS independent sparse audit: three rational points map to (0,0,0,-1/8)")
print("PASS independent sparse audit: output term counts are (2,22,139,47)")
