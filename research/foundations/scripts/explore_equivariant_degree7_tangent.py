#!/usr/bin/env python3
"""Exact tangent audit for the full degree-at-most-seven equivariant support.

This deliberately makes no exhaustion claim.  It enlarges the published
sixteen-monomial, z-linear foundational slice to every monomial of total
degree at most seven with source weights (1,-1,-2) and target weights
(-2,-1,1).  After fixing the linear part and the effective diagonal gauge,
it computes the Zariski tangent space of the Keller coefficient scheme at
the foundational point.
"""

from itertools import product

import sympy as sp


u, v, tau = sp.symbols("u v tau")


def invariant_support(output_weight: int) -> tuple[tuple[int, int], ...]:
    """Return invariant exponents after removing x**output_weight."""

    support: list[tuple[int, int]] = []
    for i, j, k in product(range(8), repeat=3):
        if i + j + k > 7 or i - j - 2 * k != output_weight:
            continue
        # x**(-output_weight) times the source monomial is u**j*v**k.
        support.append((j, k))
    return tuple(support)


A_support = invariant_support(-2)
B_support = invariant_support(-1)
C_support = invariant_support(1)

assert tuple(map(len, (A_support, B_support, C_support))) == (10, 9, 7)
assert sum(map(len, (A_support, B_support, C_support))) == 26


def coefficient_family(
    prefix: str, support: tuple[tuple[int, int], ...]
) -> tuple[sp.Symbol, ...]:
    return tuple(sp.Symbol(f"{prefix}{i}{j}") for i, j in support)


A_coefficients = coefficient_family("A", A_support)
B_coefficients = coefficient_family("B", B_support)
C_coefficients = coefficient_family("C", C_support)


def invariant_polynomial(
    coefficients: tuple[sp.Symbol, ...],
    support: tuple[tuple[int, int], ...],
) -> sp.Expr:
    return sp.Add(
        *(coefficient * u**i * v**j for coefficient, (i, j) in zip(coefficients, support))
    )


A = invariant_polynomial(A_coefficients, A_support)
B = invariant_polynomial(B_coefficients, B_support)
C = invariant_polynomial(C_coefficients, C_support)

weighted_matrix = sp.Matrix(
    [
        (-2 * A, sp.diff(A, u), sp.diff(A, v)),
        (-B, sp.diff(B, u), sp.diff(B, v)),
        (C, sp.diff(C, u), sp.diff(C, v)),
    ]
)


def coefficient(
    family: tuple[sp.Symbol, ...],
    support: tuple[tuple[int, int], ...],
    exponent: tuple[int, int],
) -> sp.Symbol:
    return family[support.index(exponent)]


# Fix the linear part (z,y,2x), then use the effective two-dimensional
# diagonal gauge to normalize the u- and v-coefficients of C.
fixed = {
    coefficient(A_coefficients, A_support, (0, 1)): 1,
    coefficient(B_coefficients, B_support, (1, 0)): 1,
    coefficient(C_coefficients, C_support, (0, 0)): 2,
    coefficient(C_coefficients, C_support, (1, 0)): -3,
    coefficient(C_coefficients, C_support, (0, 1)): -1,
}

variables = tuple(
    variable
    for variable in (*A_coefficients, *B_coefficients, *C_coefficients)
    if variable not in fixed
)
assert len(variables) == 21

keller_polynomial = sp.Poly(sp.expand(weighted_matrix.det() + 2), u, v)
keller_equations = tuple(
    sp.expand(value.subs(fixed)) for _monomial, value in keller_polynomial.terms()
)
assert len(keller_equations) == 37

# The foundational point in invariant coordinates.
foundational_point = {variable: sp.Integer(0) for variable in variables}
for exponent, value in {
    (2, 0): 4,
    (1, 1): 3,
    (3, 0): 7,
    (2, 1): 3,
    (4, 0): 3,
    (3, 1): 1,
}.items():
    foundational_point[coefficient(A_coefficients, A_support, exponent)] = value
for exponent, value in {
    (0, 1): 3,
    (2, 0): 12,
    (1, 1): 6,
    (3, 0): 9,
    (2, 1): 3,
}.items():
    foundational_point[coefficient(B_coefficients, B_support, exponent)] = value

assert all(equation.subs(foundational_point) == 0 for equation in keller_equations)

linearization = sp.Matrix(keller_equations).jacobian(variables).subs(foundational_point)
assert linearization.rank() == 20
tangent_basis = linearization.nullspace()
assert len(tangent_basis) == 1
tangent = tangent_basis[0]

# Every coefficient added beyond the published z-linear (7,6,4) slice
# vanishes in the unique tangent vector.
legacy_A = {(2, 0), (1, 1), (3, 0), (2, 1), (4, 0), (3, 1)}
legacy_B = {(0, 1), (2, 0), (1, 1), (3, 0), (2, 1)}
legacy_C: set[tuple[int, int]] = set()
legacy_variables = {
    coefficient(A_coefficients, A_support, exponent) for exponent in legacy_A
} | {
    coefficient(B_coefficients, B_support, exponent) for exponent in legacy_B
} | {
    coefficient(C_coefficients, C_support, exponent) for exponent in legacy_C
}
assert all(
    tangent[index] == 0
    for index, variable in enumerate(variables)
    if variable not in legacy_variables
)

# The unique tangent is obstructed on its literal affine line at order two.
line_substitution = {
    variable: foundational_point[variable] + tau * tangent[index]
    for index, variable in enumerate(variables)
}
line_determinant = sp.Poly(
    sp.expand(weighted_matrix.det().subs(fixed).subs(line_substitution) + 2),
    tau,
)
assert line_determinant.coeff_monomial(tau) == 0
assert line_determinant.coeff_monomial(tau**2) != 0
assert line_determinant.degree() == 2

nonzero_tangent_variables = tuple(
    str(variable)
    for index, variable in enumerate(variables)
    if tangent[index] != 0
)

print("PASS: full degree<=7 equivariant support has 10+9+7=26 monomials")
print("PASS: linear and diagonal normalization leaves 21 coefficients and 37 equations")
print("PASS: foundational tangent rank is 20, hence tangent dimension one")
print("PASS: all ten added support coefficients vanish in the unique tangent")
print("PASS: the unique tangent has a nonzero quadratic obstruction on its affine line")
print("tangent support:", ", ".join(nonzero_tangent_variables))
