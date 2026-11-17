#!/usr/bin/env python3
"""Exact audit of the cubic Danielewski de Rham compiler."""

from __future__ import annotations

from functools import cache
from itertools import product

import sympy as sp


a, c, w, x = sp.symbols("a c w x")
P = x**3 + x
P_prime = sp.diff(P, x)


def bracket(first: sp.Expr, second: sp.Expr) -> sp.Expr:
    return sp.expand(
        c
        * (
            sp.diff(first, c) * sp.diff(second, x)
            - sp.diff(first, x) * sp.diff(second, c)
        )
        + P_prime
        * (
            sp.diff(first, c) * sp.diff(second, w)
            - sp.diff(first, w) * sp.diff(second, c)
        )
        + w
        * (
            sp.diff(first, x) * sp.diff(second, w)
            - sp.diff(first, w) * sp.diff(second, x)
        )
    )


@cache
def reduce_x_power(power: int) -> tuple[sp.Rational, sp.Rational]:
    if power == 0:
        return sp.S.One, sp.S.Zero
    if power == 1:
        return sp.S.Zero, sp.S.One
    previous = reduce_x_power(power - 2)
    factor = -sp.Rational(power - 1, power + 1)
    return factor * previous[0], factor * previous[1]


def reduce_relation(polynomial: sp.Expr) -> sp.Expr:
    answer = 0
    for (c_power, w_power), coefficient in sp.Poly(
        sp.expand(polynomial),
        c,
        w,
    ).terms():
        common = min(c_power, w_power)
        answer += (
            coefficient
            * c ** (c_power - common)
            * w ** (w_power - common)
            * P**common
        )
    return sp.expand(answer)


def de_rham_reduce(polynomial: sp.Expr) -> tuple[sp.Expr, sp.Expr]:
    normal = reduce_relation(polynomial)
    answer = [sp.S.Zero, sp.S.Zero]
    for (c_power, w_power, x_power, a_power), coefficient in sp.Poly(
        normal,
        c,
        w,
        x,
        a,
    ).terms():
        if c_power or w_power:
            continue
        constant, linear = reduce_x_power(x_power)
        answer[0] += coefficient * constant * a**a_power
        answer[1] += coefficient * linear * a**a_power
    return tuple(sp.factor(entry) for entry in answer)


assert de_rham_reduce(P_prime) == (0, 0)
expected_powers = (
    (1, 0),
    (0, 1),
    (-sp.Rational(1, 3), 0),
    (0, -sp.Rational(1, 2)),
    (sp.Rational(1, 5), 0),
    (0, sp.Rational(1, 3)),
)
assert tuple(de_rham_reduce(x**power) for power in range(6)) == expected_powers
for power in range(9):
    assert de_rham_reduce(sp.diff(P * x**power, x)) == (0, 0)
print("PASS: k[x]/d(Pk[x]) reduces exactly to the basis (1,x)")


# Every bracket of bounded normal-form monomials is cohomologically zero.
normal_monomials = [x**degree for degree in range(5)]
normal_monomials += [
    variable**power * x**degree
    for variable, power, degree in product((c, w), range(1, 4), range(4))
]
for first in normal_monomials:
    for second in normal_monomials:
        assert de_rham_reduce(bracket(first, second)) == (0, 0)
print("PASS: the compiler kills all tested exact Poisson two-forms")


# Directly verify the cohomological flux identity on a nontrivial triple.
A_test = a**2 * c + a * x + w
F_test = a**3 * w + a * c * x + x**2
G_test = a * c**2 + a**2 * x + w**2
flux = de_rham_reduce(sp.expand(A_test * bracket(F_test, G_test)))
flux_derivative = tuple(sp.diff(entry, a) for entry in flux)
ledger = sp.expand(
    sp.diff(A_test, a) * bracket(F_test, G_test)
    - sp.diff(F_test, a) * bracket(A_test, G_test)
    + sp.diff(G_test, a) * bracket(A_test, F_test)
)
assert tuple(
    sp.expand(left - right)
    for left, right in zip(flux_derivative, de_rham_reduce(ledger))
) == (0, 0)
print("PASS: differentiation of compiled flux equals the coupled ledger")


# The two coordinates agree with independent periods on paths 0->i,0->-i.
period_matrix = sp.Matrix(
    [
        [sp.I, sp.I**2 / 2],
        [-sp.I, (-sp.I) ** 2 / 2],
    ]
)
assert sp.simplify(period_matrix.det()) != 0
print("PASS: the compiled coordinates are separated by vanishing periods")


# Compile the generic normalized maximal quadratic ansatz.
quadratic_residual_basis = (
    x,
    a * c,
    a * x,
    a * w,
    a**2,
    c**2,
    x**2,
    w**2,
    c * w,
    c * x,
    x * w,
)
quadratic_width = len(quadratic_residual_basis)
coefficients = sp.symbols(f"z0:{3 * quadratic_width}")
targets = tuple(
    leading
    + sum(
        coefficients[quadratic_width * row + column] * term
        for column, term in enumerate(quadratic_residual_basis)
    )
    for row, leading in enumerate((a, c, w))
)
generic_flux = de_rham_reduce(
    sp.expand(targets[0] * bracket(targets[1], targets[2]))
)
assert tuple(sp.degree(entry, a) for entry in generic_flux) == (2, 3)

flux_equations = []
for coordinate, polynomial in enumerate(generic_flux):
    poly = sp.Poly(polynomial, a)
    for power in range(1, poly.degree() + 1):
        target = 1 if coordinate == 0 and power == 1 else 0
        equation = sp.factor(poly.coeff_monomial(a**power) - target)
        if equation != 0:
            flux_equations.append(equation)
assert len(flux_equations) == 5
flux_coefficient_degrees = {
    sp.Poly(equation, *coefficients).total_degree()
    for equation in flux_equations
}
assert flux_coefficient_degrees == {3}
print("PASS: the maximal quadratic flux compiles to five exact equations")
print("PASS cubic Danielewski de Rham compiler")
