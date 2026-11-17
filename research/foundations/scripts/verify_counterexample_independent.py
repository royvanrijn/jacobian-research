#!/usr/bin/env python3
"""Dependency-free exact verifier for the finite 3D certificate.

This deliberately does not import SymPy or any ``jcsearch`` module.  It uses a
small local sparse-polynomial implementation over the integers, plus
``fractions.Fraction`` for evaluating the collision points.  Its purpose is to
avoid sharing an algebra engine or helper code with ``verify_counterexample.py``.
"""

from fractions import Fraction


Exponent = tuple[int, int, int]
Polynomial = dict[Exponent, int]
ZERO = (0, 0, 0)


def clean(poly: Polynomial) -> Polynomial:
    return {exponent: coefficient for exponent, coefficient in poly.items()
            if coefficient}


def constant(value: int) -> Polynomial:
    return {} if value == 0 else {ZERO: value}


def monomial(exponent: Exponent, coefficient: int = 1) -> Polynomial:
    return {} if coefficient == 0 else {exponent: coefficient}


def add(*polynomials: Polynomial) -> Polynomial:
    result: Polynomial = {}
    for polynomial in polynomials:
        for exponent, coefficient in polynomial.items():
            result[exponent] = result.get(exponent, 0) + coefficient
    return clean(result)


def scale(value: int, polynomial: Polynomial) -> Polynomial:
    return clean({exponent: value * coefficient
                  for exponent, coefficient in polynomial.items()})


def subtract(left: Polynomial, right: Polynomial) -> Polynomial:
    return add(left, scale(-1, right))


def multiply(*polynomials: Polynomial) -> Polynomial:
    result = constant(1)
    for polynomial in polynomials:
        product: Polynomial = {}
        for left_exponent, left_coefficient in result.items():
            for right_exponent, right_coefficient in polynomial.items():
                exponent = tuple(a + b for a, b in
                                 zip(left_exponent, right_exponent))
                product[exponent] = (
                    product.get(exponent, 0)
                    + left_coefficient * right_coefficient
                )
        result = clean(product)
    return result


def power(polynomial: Polynomial, exponent: int) -> Polynomial:
    if exponent < 0:
        raise ValueError("polynomial exponent must be nonnegative")
    result = constant(1)
    factor = polynomial
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = multiply(result, factor)
        factor = multiply(factor, factor)
        remaining //= 2
    return result


def derivative(polynomial: Polynomial, variable: int) -> Polynomial:
    result: Polynomial = {}
    for exponent, coefficient in polynomial.items():
        degree = exponent[variable]
        if degree:
            derived_exponent = list(exponent)
            derived_exponent[variable] -= 1
            result[tuple(derived_exponent)] = coefficient * degree
    return result


def determinant_3x3(matrix: list[list[Polynomial]]) -> Polynomial:
    a, b, c = matrix[0]
    d, e, f = matrix[1]
    g, h, i = matrix[2]
    return add(
        multiply(a, subtract(multiply(e, i), multiply(f, h))),
        scale(-1, multiply(b, subtract(multiply(d, i), multiply(f, g)))),
        multiply(c, subtract(multiply(d, h), multiply(e, g))),
    )


def evaluate(polynomial: Polynomial,
             point: tuple[Fraction, Fraction, Fraction]) -> Fraction:
    return sum(
        Fraction(coefficient)
        * point[0] ** exponent[0]
        * point[1] ** exponent[1]
        * point[2] ** exponent[2]
        for exponent, coefficient in polynomial.items()
    )


def total_degree(polynomial: Polynomial) -> int:
    return max(sum(exponent) for exponent in polynomial)


one = constant(1)
x = monomial((1, 0, 0))
y = monomial((0, 1, 0))
z = monomial((0, 0, 1))
xy = multiply(x, y)
u = add(one, xy)
four_plus_3xy = add(constant(4), scale(3, xy))

coordinate_1 = add(
    multiply(power(u, 3), z),
    multiply(power(y, 2), u, four_plus_3xy),
)
coordinate_2 = add(
    y,
    scale(3, multiply(x, power(u, 2), z)),
    scale(3, multiply(x, power(y, 2), four_plus_3xy)),
)
coordinate_3 = add(
    scale(2, x),
    scale(-3, multiply(power(x, 2), y)),
    scale(-1, multiply(power(x, 3), z)),
)
mapping = (coordinate_1, coordinate_2, coordinate_3)

jacobian = [
    [derivative(coordinate, variable) for variable in range(3)]
    for coordinate in mapping
]
determinant = determinant_3x3(jacobian)
assert determinant == {ZERO: -2}, determinant

points = (
    (Fraction(0), Fraction(0), Fraction(-1, 4)),
    (Fraction(1), Fraction(-3, 2), Fraction(13, 2)),
    (Fraction(-1), Fraction(3, 2), Fraction(13, 2)),
)
assert len(set(points)) == 3
target = (Fraction(-1, 4), Fraction(0), Fraction(0))
for point in points:
    image = tuple(evaluate(coordinate, point) for coordinate in mapping)
    assert image == target, (point, image)

degrees = tuple(total_degree(coordinate) for coordinate in mapping)
assert degrees == (7, 6, 4), degrees

print("PASS: dependency-free sparse determinant is -2")
print("PASS: three distinct rational points have image", target)
print("PASS: dependency-free component total degrees =", degrees)
