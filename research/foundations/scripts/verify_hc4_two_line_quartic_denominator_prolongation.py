#!/usr/bin/env python3
"""Exact first-prolongation obstruction for the clean P=x^3*y packet."""

from __future__ import annotations

from itertools import permutations, product

import sympy as sp


x, y, z, t = sp.symbols("x y z t")
variables = (x, y, z, t)
A, B, Gamma, D, a, b = sp.symbols("A B Gamma D a b")


def ternary_form(degree: int, prefix: str) -> tuple[sp.Expr, tuple[sp.Symbol, ...]]:
    """Return a general homogeneous ternary form and its coefficients."""

    monomials = []
    for x_degree in range(degree, -1, -1):
        for y_degree in range(degree - x_degree, -1, -1):
            z_degree = degree - x_degree - y_degree
            monomials.append(x**x_degree * y**y_degree * z**z_degree)
    coefficients = sp.symbols(f"{prefix}0:{len(monomials)}")
    return sum(
        coefficient * monomial
        for coefficient, monomial in zip(coefficients, monomials)
    ), coefficients


def quaternary_quadratic() -> tuple[sp.Expr, tuple[sp.Symbol, ...]]:
    """Return a general homogeneous quadratic in x,y,z,t."""

    monomials = []
    for i in range(2, -1, -1):
        for j in range(2 - i, -1, -1):
            for k in range(2 - i - j, -1, -1):
                ell = 2 - i - j - k
                monomials.append(x**i * y**j * z**k * t**ell)
    coefficients = sp.symbols(f"h0:{len(monomials)}")
    return sum(
        coefficient * monomial
        for coefficient, monomial in zip(coefficients, monomials)
    ), coefficients


def permutation_sign(permutation: tuple[int, ...]) -> int:
    inversions = sum(
        permutation[i] > permutation[j]
        for i in range(4)
        for j in range(i + 1, 4)
    )
    return -1 if inversions % 2 else 1


def determinant_layer(
    matrices: dict[int, sp.Matrix], target_degree: int
) -> sp.Expr:
    """Extract one lambda coefficient by determinant multilinearity."""

    result = 0
    available_degrees = tuple(matrices)
    for permutation in permutations(range(4)):
        sign = permutation_sign(permutation)
        for degrees in product(available_degrees, repeat=4):
            if sum(degrees) != target_degree:
                continue
            term = sp.Integer(sign)
            for row, degree in enumerate(degrees):
                term *= matrices[degree][row, permutation[row]]
            result += term
    return sp.expand(result)


r4, r4_coefficients = ternary_form(4, "r")
Q2, q2_coefficients = ternary_form(2, "q")
r3, r3_coefficients = ternary_form(3, "p")
h2, h2_coefficients = quaternary_quadratic()

h5 = (
    A * x * y**4
    + x**4 * (B * y + Gamma * z) / 24
    + D * x**5 / 120
)
s3 = a * x * y**2 + b * x**3
h4 = t * s3 + r4
h3 = a**2 * x * t**2 / (6 * A) + t * Q2 + r3

matrices = {
    3: sp.hessian(h5, variables),
    2: sp.hessian(h4, variables),
    1: sp.hessian(h3, variables),
    0: sp.hessian(h2, variables),
}

# The forced xt^2 coefficient exactly cancels the Schur face.
degree_ten = determinant_layer(matrices, 10)
assert degree_ten == 0

# The next face has one repair-independent nonzero extreme coefficient.
degree_nine = determinant_layer(matrices, 9)
immutable = sp.factor(
    sp.Poly(degree_nine, x, y, z, t).coeff_monomial(x**8 * t)
)
assert immutable == -Gamma**2 * a**3 / (54 * A)

repair_coefficients = set(
    r4_coefficients + q2_coefficients + r3_coefficients + h2_coefficients
)
assert not (immutable.free_symbols & repair_coefficients)
assert not (immutable.free_symbols & {B, D, b})

# The only degree-nine contribution using h2 is det(Hess(h5))*h2_tt.
h2_only = determinant_layer({3: matrices[3], 0: matrices[0]}, 9)
assert sp.factor(h2_only) == sp.factor(matrices[3][:3, :3].det() * matrices[0][3, 3])
assert sp.Poly(h2_only, x, y, z, t).coeff_monomial(x**8 * t) == 0

print("PASS: complete degree-ten Schur face cancels")
print("PASS: retained every quartic, cubic, and quadratic repair coefficient")
print("PASS: immutable degree-nine coefficient is -Gamma^2*a^3/(54*A)")
print("THEOREM: the clean x^7*y^2 two-line packet has no HC4 prolongation")
print("SCOPE: collision equations and the other P=x^3*y incidences are not used")
