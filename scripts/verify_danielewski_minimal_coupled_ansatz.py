#!/usr/bin/env python3
"""Exact audit of the minimal coupled Danielewski target ansatz."""

from __future__ import annotations

from itertools import combinations
from random import Random

import sympy as sp


a, b, c = sp.symbols("a b c")
x = b * c
w = b * (x**2 + 1)

basis = (a, c, x, w, a * c, a * x, a * w)
names = ("a", "c", "x", "w", "ac", "ax", "aw")
triples = tuple(combinations(range(len(basis)), 3))
basis_jacobian = sp.Matrix(basis).jacobian((a, b, c))


def coefficient(polynomial: sp.Expr, monomial: sp.Expr) -> sp.Expr:
    return sp.Poly(sp.expand(polynomial), a, b, c).coeff_monomial(monomial)


minors: dict[tuple[int, int, int], sp.Expr] = {
    triple: sp.expand(basis_jacobian[list(triple), :].det())
    for triple in triples
}

carriers = []
for triple, minor in minors.items():
    constant = coefficient(minor, 1)
    cubic_trace = coefficient(minor, b**2 * c**2)
    assert sp.expand(cubic_trace - 3 * constant) == 0
    if constant or cubic_trace:
        carriers.append((triple, minor))

assert carriers == [((0, 1, 3), -3 * b**2 * c**2 - 1)]
assert minors[(0, 1, 3)] == sp.Matrix((a, c, w)).jacobian(
    (a, b, c)
).det()
print("PASS: all 35 Pluecker-ledger minors obey the coefficient lock")
print("PASS: det(a,c,w)=-(1+3*b^2*c^2) is the unique carrier")


def target_jacobian(coefficient_matrix: sp.Matrix) -> sp.Expr:
    targets = coefficient_matrix * sp.Matrix(basis)
    return sp.expand(targets.jacobian((a, b, c)).det())


def cauchy_binet(coefficient_matrix: sp.Matrix) -> sp.Expr:
    answer = 0
    for triple in triples:
        pluecker = coefficient_matrix[:, list(triple)].det()
        answer += pluecker * minors[triple]
    return sp.expand(answer)


random = Random(20260723)
for _ in range(40):
    matrix = sp.Matrix(
        3,
        len(basis),
        [random.randint(-3, 3) for _ in range(3 * len(basis))],
    )
    direct = target_jacobian(matrix)
    ledger = cauchy_binet(matrix)
    assert sp.expand(direct - ledger) == 0
    assert coefficient(direct, b**2 * c**2) == 3 * coefficient(direct, 1)
print("PASS: deterministic target matrices satisfy Cauchy--Binet and the lock")


# The ansatz genuinely permits several a-active target coordinates.
fully_active = sp.Matrix(
    [
        [1, 0, 0, 0, 1, 0, 0],
        [0, 1, 0, 0, 0, 1, 0],
        [0, 0, 0, 1, 0, 0, 1],
    ]
)
active_targets = tuple(fully_active * sp.Matrix(basis))
assert all(sp.diff(target, a) != 0 for target in active_targets)
active_jacobian = target_jacobian(fully_active)
assert coefficient(active_jacobian, b**2 * c**2) == 3 * coefficient(
    active_jacobian,
    1,
)
print("PASS: the excluded space includes triples with three active a-directions")


# Adjoining a^2 cannot alter either coefficient: any new minor has a factor
# a or vanishes because the a and a^2 differential rows are proportional.
extended_basis = (*basis, a**2)
extended_jacobian = sp.Matrix(extended_basis).jacobian((a, b, c))
for triple in combinations(range(len(extended_basis)), 3):
    minor = sp.expand(extended_jacobian[list(triple), :].det())
    assert coefficient(minor, b**2 * c**2) == 3 * coefficient(minor, 1)
print("PASS: adjoining a^2 does not escape the coefficient lock")


# All ambient-quadratic monomials except cx and xw obey a stronger lock.
diagonal_basis = (
    a,
    c,
    x,
    w,
    a * c,
    a * x,
    a * w,
    a**2,
    c**2,
    x**2,
    w**2,
    c * w,
)
diagonal_jacobian = sp.Matrix(diagonal_basis).jacobian((a, b, c))
diagonal_triples = tuple(combinations(range(len(diagonal_basis)), 3))
assert len(diagonal_triples) == 220
for triple in diagonal_triples:
    minor = sp.expand(diagonal_jacobian[list(triple), :].det())
    locked = (
        coefficient(minor, x**2)
        - coefficient(minor, x**4)
        + coefficient(minor, x**6)
        - 3 * coefficient(minor, 1)
    )
    assert sp.expand(locked) == 0
print("PASS: all 220 diagonal-quadratic minors obey the stronger lock")


# Normalize the minimal escape chart.  z_(5*i+j) is the coefficient of
# N_j in target row i, for N=(ac,ax,aw,cx,xw) and leading rows (a,c,w).
escape_basis = (a * c, a * x, a * w, c * x, x * w)
z = sp.symbols("z0:15")
escape_targets = tuple(
    leading
    + sum(z[5 * row + column] * term for column, term in enumerate(escape_basis))
    for row, leading in enumerate((a, c, w))
)
escape_error = sp.Poly(
    sp.expand(
        sp.Matrix(escape_targets).jacobian((a, b, c)).det() + 1
    ),
    a,
    b,
    c,
)
escape_equations = tuple(escape_error.coeffs())
linear_equations = tuple(
    equation
    for equation in escape_equations
    if sp.Poly(equation, *z).total_degree() <= 1
)
linear_solution = next(iter(sp.linsolve(linear_equations, z)))
linear_substitution = {
    variable: answer
    for variable, answer in zip(z, linear_solution)
    if variable != answer
}
assert linear_substitution == {
    z[0]: 0,
    z[1]: -2 * z[14] - 2 * z[8],
    z[2]: 0,
    z[5]: -z[12],
    z[9]: 0,
    z[13]: 0,
}

r, s_value = z[14], z[8]
p, q = z[3] * z[6], z[11] * z[4]
equation_by_monomial = {
    monomial: sp.factor(equation.subs(linear_substitution))
    for monomial, equation in escape_error.terms()
}
E2 = 2 * q + 4 * r**2 + 5 * r * s_value + 2 * p + 4 * s_value**2 - 3
E4 = 4 * q + 8 * r**2 + 11 * r * s_value + 4 * p + 8 * s_value**2
K = q * s_value + 2 * r**2 * s_value + r * p + 2 * r * s_value**2
E3 = 3 * K + 2 * (r + s_value)
assert sp.expand(equation_by_monomial[(0, 2, 2)] - E2) == 0
assert sp.expand(equation_by_monomial[(0, 4, 4)] - E4) == 0
assert sp.expand(equation_by_monomial[(0, 5, 5)] - 5 * K) == 0
assert sp.expand(equation_by_monomial[(0, 3, 3)] - E3) == 0
assert equation_by_monomial[(0, 1, 3)] == z[11] * z[3]
assert sp.expand(E4 - 2 * E2 - (r * s_value + 6)) == 0
assert sp.expand(E3 - 3 * K - 2 * (r + s_value)) == 0

core = (E2, E4, K, E3, z[11] * z[3])
free_variables = (z[3], z[4], z[6], z[8], z[11], z[14])
core_basis = sp.groebner(core, *free_variables, order="lex")
assert core_basis.polys == [sp.Poly(1, *free_variables)]
print("PASS: the normalized cx/xw escape chart has a five-equation unit ideal")
print("PASS minimal coupled Danielewski ansatz")
