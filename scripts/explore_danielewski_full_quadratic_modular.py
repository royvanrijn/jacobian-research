#!/usr/bin/env python3
"""Modular regressions for the maximal ambient-quadratic coefficient ideal.

The characteristic-zero unit ideal is certified by the separate full
quadratic obstruction checker. This script reproduces three reductions.
"""

from __future__ import annotations

import shutil
import subprocess

import sympy as sp


a, b, c = sp.symbols("a b c")
x = b * c
w = b * (x**2 + 1)
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
error = sp.Poly(
    sp.expand(sp.Matrix(targets).jacobian((a, b, c)).det() + 1),
    a,
    b,
    c,
)
assert len(error.terms()) == 102

raw_equations = error.coeffs()
linear_equations = tuple(
    equation
    for equation in raw_equations
    if sp.Poly(equation, *coefficients).total_degree() <= 1
)
linear_matrix, _ = sp.linear_eq_to_matrix(linear_equations, coefficients)
assert linear_matrix.rank() == 3
linear_solution = next(iter(sp.linsolve(linear_equations, coefficients)))
linear_substitution = {
    variable: answer
    for variable, answer in zip(coefficients, linear_solution)
    if variable != answer
}
free_variables = tuple(
    sorted(
        set().union(*(answer.free_symbols for answer in linear_solution))
        & set(coefficients),
        key=lambda variable: int(str(variable)[1:]),
    )
)
assert len(free_variables) == 27

reduced_equations: list[sp.Expr] = []
for _, equation in error.terms():
    reduced = sp.factor(equation.xreplace(linear_substitution))
    if reduced == 0:
        continue
    if all(
        sp.expand(reduced - previous) != 0
        for previous in reduced_equations
    ):
        reduced_equations.append(sp.expand(reduced))
assert len(reduced_equations) == 89

singular = shutil.which("Singular")
if singular is None:
    raise SystemExit("Singular is required for the modular exploration")

polynomials = tuple(
    sp.sstr(equation).replace("**", "^")
    for equation in reduced_equations
)
for prime in (32003, 32009, 32027):
    program = (
        f"ring r={prime},({','.join(map(str, free_variables))}),dp;\n"
        "option(redSB);\n"
        f"ideal I={','.join(polynomials)};\n"
        "ideal G=slimgb(I);\n"
        "print(size(G));\n"
        "print(G[1]);\n"
    )
    result = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
        timeout=300,
    )
    assert result.stdout.strip() == "1\n1"
    print(f"PASS: maximal quadratic ideal is unit modulo {prime}")

print("PASS modular regressions for the exact rational quadratic obstruction")
