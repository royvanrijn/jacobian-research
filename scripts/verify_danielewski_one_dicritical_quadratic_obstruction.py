#!/usr/bin/env python3
"""Verify the coupled quadratic obstruction for cw=x(x+1).

The target surface is the first smooth Danielewski completion with one
missing boundary plane.  We normalize the source-linear jet to (a,c,w) and
allow every ambient polynomial of degree at most two modulo cw=x(x+1).
The resulting constant-Jacobian coefficient ideal is tested over several
large finite fields and reconstructed exactly in characteristic zero.
"""

from __future__ import annotations

import shutil
import subprocess

import sympy as sp


a, b, c = sp.symbols("a b c")
x = b * c
w = b * (x + 1)

# The relation cw=x^2+x removes cw from the ambient-quadratic basis.
residual_basis = (
    x,
    a * c,
    a * x,
    a * w,
    a**2,
    c**2,
    x**2,
    w**2,
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

linear_equations = tuple(
    equation
    for equation in error.coeffs()
    if sp.Poly(equation, *coefficients).total_degree() <= 1
)
linear_matrix, _ = sp.linear_eq_to_matrix(linear_equations, coefficients)
linear_rank = linear_matrix.rank()
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

reduced_equations: list[sp.Expr] = []
seen_equations: set[str] = set()
for _, equation in error.terms():
    reduced = sp.expand(equation.xreplace(linear_substitution))
    if reduced == 0:
        continue
    key = sp.srepr(reduced)
    if key not in seen_equations:
        seen_equations.add(key)
        reduced_equations.append(reduced)

print(
    "LEDGER",
    f"source_terms={len(error.terms())}",
    f"coefficients={len(coefficients)}",
    f"linear_rank={linear_rank}",
    f"free_variables={len(free_variables)}",
    f"reduced_equations={len(reduced_equations)}",
)
assert len(error.terms()) == 73
assert len(coefficients) == 30
assert linear_rank == 3
assert len(free_variables) == 27
assert len(reduced_equations) == 64

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
    raise SystemExit("Singular is required")

polynomials = tuple(
    sp.sstr(equation).replace("**", "^") for equation in reduced_equations
)
for prime in (32003, 32009, 32027):
    program = (
        f"ring r={prime},({','.join(map(str, column_major))}),dp;\n"
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
    lines = result.stdout.strip().splitlines()
    assert lines == ["1", "1"]
    print("MODULAR", prime, "size=" + lines[0], "first=" + lines[1])

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
    timeout=600,
)
assert result.stdout.strip() == "1\n1"
print("RATIONAL size=1 first=1")
print("PASS one-dicritical ambient-quadratic Danielewski obstruction")
