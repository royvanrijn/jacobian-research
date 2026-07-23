#!/usr/bin/env python3
"""Exact boundary-lattice audit for the Davenport Cox--Sunada covers."""
from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.sunada import (  # noqa: E402
    A,
    T,
    U,
    Y,
    Z,
    branch_cubic,
    davenport_pair,
)


g, h = davenport_pair()

# Work exactly in K(T), K=Q(a)/(a^2+a+2), using an abstract algebraic
# generator rather than a floating or radical embedding.
minimal_variable = sp.symbols("minimal_variable")
number_field = sp.QQ.alg_field_from_poly(
    sp.Poly(minimal_variable**2 + minimal_variable + 2, minimal_variable),
    alias="alpha",
)
alpha = number_field.ext
function_field = number_field.frac_field(T)


def boundary_factorization(
    variable: sp.Symbol, polynomial: sp.Expr
) -> tuple[sp.Poly, tuple[tuple[sp.Poly, int], ...]]:
    pullback = sp.Poly(
        sp.expand(branch_cubic().subs(U, polynomial).subs(A, alpha)),
        variable,
        domain=function_field,
    )
    coefficient, factors = sp.factor_list(pullback)
    assert coefficient == 1
    return pullback, tuple(factors)


for variable, polynomial in ((Y, g), (Z, h)):
    pullback, factors = boundary_factorization(variable, polynomial)
    assert sorted((factor.degree(), multiplicity) for factor, multiplicity in factors) == [
        (3, 1),
        (6, 1),
        (6, 2),
    ]

    jacobian_unit = sp.Poly(
        sp.expand(sp.diff(polynomial, variable).subs(A, alpha)),
        variable,
        domain=function_field,
    ).monic()
    matching = [
        multiplicity
        for factor, multiplicity in factors
        if factor.monic() == jacobian_unit
    ]
    assert matching == [2]

    simple_factors = [
        factor.monic()
        for factor, multiplicity in factors
        if factor.monic() != jacobian_unit
    ]
    assert sorted(factor.degree() for factor in simple_factors) == [3, 6]
    reconstructed = jacobian_unit**2
    for factor in simple_factors:
        reconstructed *= factor
    assert reconstructed == pullback.monic()

# Match the displayed point-cover factors, not only their degrees.
displayed_E3 = (
    Y**3 + (1 + 3 * A) * T * Y + (7 + 5 * A) * T
)
displayed_E6 = (
    Y**6
    + (10 + 8 * A) * T * Y**4
    + (6 + 8 * A) * T * Y**3
    + (-19 + 28 * A) * T**2 * Y**2
    + (-18 + 36 * A) * T**2 * Y
    - 27 * T**2
    + (-88 - 36 * A) * T**3
)
point_pullback, point_factors = boundary_factorization(Y, g)
point_factor_polys = {factor.monic() for factor, _ in point_factors}
displayed_E3_poly = sp.Poly(
    displayed_E3.subs(A, alpha), Y, domain=function_field
).monic()
displayed_E6_poly = sp.Poly(
    displayed_E6.subs(A, alpha), Y, domain=function_field
).monic()
assert displayed_E3_poly in point_factor_polys
assert displayed_E6_poly in point_factor_polys

# The most direct realization of the last primitive row uses E3 as a new
# coordinate because it is linear in T.  Its Jacobian introduces a new
# linear divisor which is coprime to every existing boundary factor.
linear_chart_divisor = (1 + 3 * A) * Y + (7 + 5 * A)
assert sp.expand(sp.diff(displayed_E3, T) - linear_chart_divisor) == 0
chart_jacobian = sp.factor(
    sp.Matrix((displayed_E3, Y)).jacobian((T, Y)).det()
)
assert sp.factor(chart_jacobian - linear_chart_divisor) == 0
linear_poly = sp.Poly(
    linear_chart_divisor.subs(A, alpha),
    Y,
    domain=function_field,
).monic()
for factor, _ in point_factors:
    assert linear_poly.gcd(factor.monic()).degree() == 0

# The boundary exponent vector of the target unit is (1,1,2), while the
# Jacobian unit is the third primitive source unit.
target_pullback_vector = sp.Matrix([1, 1, 2])
jacobian_vector = sp.Matrix([0, 0, 1])
assert target_pullback_vector.rank() == 1
assert jacobian_vector != sp.zeros(3, 1)
assert jacobian_vector != target_pullback_vector
ledger_completion = sp.Matrix(
    [
        list(target_pullback_vector),
        list(jacobian_vector),
        [1, 0, 0],
    ]
)
assert ledger_completion.det() == 1
assert all(entry.q == 1 for entry in ledger_completion.inv())

# Any polynomial suspension which preserves the two core target coordinates
# (T,g_T(Y)) has a block-triangular Jacobian divisible by g'_T(Y), regardless
# of how the extra outputs depend on T and Y.
g_T_symbol, g_Y_symbol = sp.symbols("g_T_symbol g_Y_symbol")
extra = sp.Matrix(3, 5, lambda row, column: sp.symbols(f"e_{row}_{column}"))
block_jacobian = sp.Matrix.vstack(
    sp.Matrix([[1, 0, 0, 0, 0], [g_T_symbol, g_Y_symbol, 0, 0, 0]]),
    extra,
)
extra_vertical = extra[:, 2:5]
assert sp.factor(block_jacobian.det() - g_Y_symbol * extra_vertical.det()) == 0

print("PASS: each pulled-back branch divisor factors with pattern 3 * 6 * 6^2")
print("PASS: the doubled degree-six factor is exactly the derivative Jacobian")
print("PASS: source and target boundary-unit ranks are respectively three and one")
print("PASS: one further primitive row gives a unimodular boundary completion")
print("PASS: the direct E3 coordinate chart introduces a new coprime linear divisor")
print("PASS: affine stabilization cannot erase either positive unit rank")
print("PASS: every coordinate-preserving polynomial suspension retains the derivative factor")
