#!/usr/bin/env python3
"""Exact characteristic-zero audit of the full quadratic obstruction."""

from __future__ import annotations

import shutil
import subprocess

import sympy as sp


a, b, c = sp.symbols("a b c")
x = b * c
w = b * (x**2 + 1)
residual_basis = (
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
width = len(residual_basis)
coefficients = sp.symbols(f"z0:{3 * width}")
targets = tuple(
    leading
    + sum(
        coefficients[width * row + column] * term
        for column, term in enumerate(residual_basis)
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

linear_equations = tuple(
    equation
    for equation in error.coeffs()
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
assert len(free_variables) == 30

reduced_equations: list[sp.Expr] = []
seen_equations: set[str] = set()
for _, equation in error.terms():
    reduced = sp.factor(equation.xreplace(linear_substitution))
    if reduced == 0:
        continue
    reduced = sp.expand(reduced)
    key = sp.srepr(reduced)
    if key not in seen_equations:
        seen_equations.add(key)
        reduced_equations.append(reduced)
assert len(reduced_equations) == 89
print("PASS: the corrected normalized ledger has 30 variables and 89 equations")


# Column-major order: quadratic monomial first, then target row.
column_major = tuple(
    sorted(
        free_variables,
        key=lambda variable: (
            int(str(variable)[1:]) % width,
            int(str(variable)[1:]) // width,
        ),
    )
)
singular = shutil.which("Singular")
if singular is None:
    raise SystemExit("Singular is required for the rational certificate")

polynomials = tuple(
    sp.sstr(equation).replace("**", "^")
    for equation in reduced_equations
)
program = (
    f"ring r=0,({','.join(map(str, column_major))}),dp;\n"
    "option(redSB);\n"
    f"ideal I={','.join(polynomials)};\n"
    'LIB "modstd.lib";\n'
    'ideal G=modGB("slimgb",I,1);\n'
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
print("PASS: exact modular reconstruction over Q gives Groebner basis {1}")
print("PASS full ambient-quadratic Danielewski obstruction")
