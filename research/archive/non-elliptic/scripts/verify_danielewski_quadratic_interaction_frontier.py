#!/usr/bin/env python3
"""Exact audit of the quadratic Danielewski interaction frontier."""

from __future__ import annotations

from itertools import combinations

import sympy as sp


a, b, c = sp.symbols("a b c")
x = b * c
w = b * (x**2 + 1)

names = ("ac", "ax", "aw", "a2", "c2", "x2", "w2", "cw", "cx", "xw")
nonlinear_basis = (
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
common_indices = {0, 1, 2}
diagonal_indices = {3, 4, 5, 6, 7}
escape_indices = {8, 9}
extra_indices = tuple(sorted(diagonal_indices | escape_indices))

width = len(nonlinear_basis)
coefficients = sp.symbols(f"z0:{3 * width}")
targets = tuple(
    leading
    + sum(
        coefficients[width * row + column] * term
        for column, term in enumerate(nonlinear_basis)
    )
    for row, leading in enumerate((a, c, w))
)

# The normalized determinant has constant term -1 at the source origin.
error = sp.Poly(
    sp.expand(sp.Matrix(targets).jacobian((a, b, c)).det() + 1),
    a,
    b,
    c,
)
generic_terms = error.terms()
assert len(generic_terms) == 102


def unit_ideal_for_chart(chart: tuple[int, ...]) -> bool:
    active = common_indices | set(chart)
    zero_substitution = {
        coefficients[width * row + column]: sp.S.Zero
        for row in range(3)
        for column in range(width)
        if column not in active
    }
    equations = tuple(
        specialized
        for _, equation in generic_terms
        if (specialized := equation.xreplace(zero_substitution)) != 0
    )
    variables = tuple(
        sorted(
            set().union(*(equation.free_symbols for equation in equations))
            & set(coefficients),
            key=lambda variable: int(str(variable)[1:]),
        )
    )

    linear_equations = tuple(
        equation
        for equation in equations
        if sp.Poly(equation, *variables).total_degree() <= 1
    )
    linear_solution = next(iter(sp.linsolve(linear_equations, variables)))
    linear_substitution = {
        variable: answer
        for variable, answer in zip(variables, linear_solution)
        if variable != answer
    }
    free_variables = tuple(
        sorted(
            set().union(*(answer.free_symbols for answer in linear_solution))
            & set(variables),
            key=lambda variable: int(str(variable)[1:]),
        )
    )

    reduced_equations: list[sp.Expr] = []
    for equation in equations:
        reduced = sp.factor(equation.xreplace(linear_substitution))
        if reduced != 0 and reduced not in reduced_equations:
            reduced_equations.append(reduced)

    basis = sp.groebner(
        reduced_equations,
        *free_variables,
        order="grevlex",
    )
    return basis.polys == [sp.Poly(1, *free_variables)]


audited_charts: list[tuple[int, ...]] = []
for size in range(1, 5):
    for chart in combinations(extra_indices, size):
        if not (set(chart) & escape_indices):
            continue
        assert unit_ideal_for_chart(chart)
        audited_charts.append(chart)

assert len(audited_charts) == 68
assert sum(len(chart) == 1 for chart in audited_charts) == 2
assert sum(len(chart) == 2 for chart in audited_charts) == 11
assert sum(len(chart) == 3 for chart in audited_charts) == 25
assert sum(len(chart) == 4 for chart in audited_charts) == 30

labels = tuple(tuple(names[index] for index in chart) for chart in audited_charts)
assert ("cx",) in labels
assert ("xw",) in labels
assert ("a2", "cx", "xw") in labels
assert ("cw", "cx", "xw") in labels
print("PASS: all 2 one-direction interaction charts have unit ideal")
print("PASS: all 11 two-direction interaction charts have unit ideal")
print("PASS: all 25 three-direction interaction charts have unit ideal")
print("PASS: all 30 four-direction interaction charts have unit ideal")
print("PASS quadratic Danielewski interaction frontier (68 exact charts)")
