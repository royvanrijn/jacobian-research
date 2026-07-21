#!/usr/bin/env python3
"""Clean-room weighted marked-root theorem audit using only Python's standard library."""
from __future__ import annotations

from fractions import Fraction
from itertools import permutations


Monomial = tuple[int, int, int]
Polynomial = dict[Monomial, Fraction]


def clean(poly: Polynomial) -> Polynomial:
    return {monomial: value for monomial, value in poly.items() if value}


def constant(value: int | Fraction) -> Polynomial:
    value = Fraction(value)
    return {} if not value else {(0, 0, 0): value}


def variable(index: int) -> Polynomial:
    exponent = [0, 0, 0]
    exponent[index] = 1
    return {tuple(exponent): Fraction(1)}


def add(*polynomials: Polynomial) -> Polynomial:
    result: Polynomial = {}
    for polynomial in polynomials:
        for monomial, coefficient in polynomial.items():
            result[monomial] = result.get(monomial, Fraction(0)) + coefficient
    return clean(result)


def scale(polynomial: Polynomial, scalar: int | Fraction) -> Polynomial:
    scalar = Fraction(scalar)
    return clean({monomial: scalar * value for monomial, value in polynomial.items()})


def multiply(left: Polynomial, right: Polynomial) -> Polynomial:
    result: Polynomial = {}
    for lm, lc in left.items():
        for rm, rc in right.items():
            monomial = tuple(a + b for a, b in zip(lm, rm))
            result[monomial] = result.get(monomial, Fraction(0)) + lc * rc
    return clean(result)


def power(polynomial: Polynomial, exponent: int) -> Polynomial:
    result = constant(1)
    base = polynomial
    while exponent:
        if exponent & 1:
            result = multiply(result, base)
        base = multiply(base, base)
        exponent //= 2
    return result


def derivative(polynomial: Polynomial, index: int) -> Polynomial:
    result: Polynomial = {}
    for monomial, coefficient in polynomial.items():
        if monomial[index]:
            exponent = list(monomial)
            factor = exponent[index]
            exponent[index] -= 1
            result[tuple(exponent)] = coefficient * factor
    return clean(result)


def divide_x_power(polynomial: Polynomial, exponent: int) -> Polynomial:
    assert all(monomial[0] >= exponent for monomial in polynomial)
    return {
        (monomial[0] - exponent, monomial[1], monomial[2]): coefficient
        for monomial, coefficient in polynomial.items()
    }


def compose_univariate(coefficients: list[Fraction], argument: Polynomial) -> Polynomial:
    result: Polynomial = {}
    current = constant(1)
    for coefficient in coefficients:
        result = add(result, scale(current, coefficient))
        current = multiply(current, argument)
    return result


def compose_after_gamma_division(
    coefficients: list[Fraction],
    u: Polynomial,
    gamma: Polynomial,
    divisor: int,
) -> Polynomial:
    result: Polynomial = {}
    for exponent, coefficient in enumerate(coefficients):
        if not coefficient:
            continue
        assert exponent >= divisor
        term = multiply(power(u, exponent), power(gamma, exponent - divisor))
        result = add(result, scale(term, coefficient))
    return result


def univariate_multiply(left: list[Fraction], right: list[Fraction]) -> list[Fraction]:
    result = [Fraction(0)] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] += a * b
    return result


def univariate_derivative(coefficients: list[Fraction]) -> list[Fraction]:
    return [Fraction(index) * value for index, value in enumerate(coefficients)][1:]


def univariate_value(coefficients: list[Fraction], value: int) -> Fraction:
    return sum(coefficient * Fraction(value) ** exponent for exponent, coefficient in enumerate(coefficients))


def determinant(matrix: list[list[Polynomial]]) -> Polynomial:
    result: Polynomial = {}
    for permutation in permutations(range(3)):
        inversions = sum(
            permutation[i] > permutation[j]
            for i in range(3)
            for j in range(i + 1, 3)
        )
        term = constant(-1 if inversions % 2 else 1)
        for row, column in enumerate(permutation):
            term = multiply(term, matrix[row][column])
        result = add(result, term)
    return result


def audit_primitive(label: str, primitive: list[Fraction]) -> tuple[str, int]:
    while primitive and primitive[-1] == 0:
        primitive.pop()
    degree = len(primitive) - 1
    seed = univariate_derivative(primitive)
    c = -univariate_value(seed, 1)
    assert degree >= 3 and c
    assert primitive[0] == 0 and univariate_value(primitive, 1) == 0
    assert seed[0] == 0
    second = univariate_derivative(seed)
    kappa = univariate_value(second, 1) / c
    assert kappa != -2
    a = -(1 + kappa) / (2 + kappa)

    x, y, z = variable(0), variable(1), variable(2)
    v = multiply(x, y)
    source_s = multiply(power(x, 2), z)
    u = add(constant(1), v)
    gamma = add(constant(1), scale(v, a), source_s)
    marked_root = multiply(u, gamma)

    # q=(W H'(W)-H(W))/c, formed independently coefficient by coefficient.
    shifted_seed = [Fraction(0)] + seed
    q_length = max(len(shifted_seed), len(primitive))
    q = [Fraction(0)] * q_length
    for index in range(q_length):
        q[index] = (
            (shifted_seed[index] if index < len(shifted_seed) else 0)
            - (primitive[index] if index < len(primitive) else 0)
        ) / c
    assert q[0] == q[1] == 0

    beta_numerator = add(
        constant(c), compose_after_gamma_division(seed, u, gamma, 1)
    )
    alpha_numerator = add(
        u, compose_after_gamma_division(q, u, gamma, 2)
    )
    target_b = divide_x_power(beta_numerator, 1)
    target_a = divide_x_power(alpha_numerator, 2)
    target_c = multiply(x, gamma)

    jacobian = [
        [derivative(component, index) for index in range(3)]
        for component in (target_a, target_b, target_c)
    ]
    assert determinant(jacobian) == constant(c)

    # Independently check the marked incidence and its derivative identity.
    incidence = add(
        compose_univariate(primitive, marked_root),
        scale(multiply(multiply(target_b, target_c), marked_root), -1),
        scale(multiply(target_a, power(target_c, 2)), c),
    )
    incidence_derivative = add(
        compose_univariate(seed, marked_root),
        scale(multiply(target_b, target_c), -1),
        scale(gamma, c),
    )
    assert not incidence
    assert not incidence_derivative

    # The displayed reconstruction formulas reduce to these denominator-free
    # identities, so they hold in the common fraction field.
    assert target_c == multiply(x, gamma)
    assert marked_root == multiply(u, gamma)
    assert add(multiply(target_b, multiply(x, gamma)), scale(gamma, -c), scale(compose_univariate(seed, marked_root), -1)) == {}
    assert add(multiply(target_a, multiply(power(x, 2), power(gamma, 2))), scale(multiply(u, power(gamma, 2)), -1), scale(compose_univariate(q, marked_root), -c)) == {}
    return label, degree


audited: list[tuple[str, int]] = []
for degree in range(3, 10):
    primitive = [Fraction(0)] * (degree + 2)
    primitive[degree] = 1
    primitive[degree + 1] = -1
    audited.append(audit_primitive(f"canonical degree {degree + 1}", primitive))

# Independent noncanonical admissible primitives H=W^d(1-W)(1+lambda*W(1-W)).
for degree, parameter in ((2, -3), (3, -3), (4, -3)):
    base = [Fraction(0)] * degree + [Fraction(1), Fraction(-1)]
    deformation = [Fraction(1), Fraction(parameter), Fraction(-parameter)]
    audited.append(
        audit_primitive(
            f"deformed base {degree}, lambda={parameter}",
            univariate_multiply(base, deformation),
        )
    )

assert [degree for _, degree in audited[:7]] == list(range(4, 11))
print("PASS: clean-room polynomiality and constant Jacobian in degrees 4--10")
print("PASS: clean-room marked-incidence and reconstruction identities")
print("PASS: canonical and noncanonical seeds audited without SymPy/project imports")
