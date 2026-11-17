#!/usr/bin/env python3
"""Exact obstruction for the degenerate cubic boundary branch in chart 0010.

The boundary Schur chain has common pivot

    L = -det Hess(f-5*b^3/6).

On L=0, the two-variable zero-Hessian normal form writes

    f = 5*b^3/6 + P(alpha*b+beta*c)

up to affine terms.  For a potential of degree at most three, nonzero
boundary determinant forces P'' to be a nonzero constant.  The remaining
boundary equation contains the term 4*beta*b^3, while derivatives of
g=V_a have degree at most one.  Hence beta=0, and after rescaling the linear
form the complete relevant cubic Cauchy data are

    f = 5*b^3/6 + lambda*b^2/2,
    g = (mu+2)*c + g1*b + g2*b^2,
    h = b^2/4 + h1*b + h2*b^2,

with lambda*mu != 0.  Add the general degree-three terms containing at least
two normal variables a,d.  The boundary determinant is then
lambda*mu^2/2.

This script restricts the full graph determinant to Y=D=0, sets it equal to
that boundary constant, saturates by lambda*mu, and verifies over Q that the
49 slice coefficients generate the unit ideal.
"""

from __future__ import annotations

import runpy
import subprocess

import sympy as sp


graph = runpy.run_path("scripts/search_hc4_graph_polarizations.py")
h_variables = graph["h_variables"]
Q = list(graph["position_coordinates_h"])
M = list(graph["momentum_coordinates_h"])
X, Y, W, D = h_variables

mask = (0, 0, 1, 0)
q = sp.Matrix([M[index] if mask[index] else Q[index] for index in range(4)])
m = sp.Matrix(
    [-Q[index] if mask[index] else M[index] for index in range(4)]
)

u = sp.symbols("u0:4")
a, b, c = u[:3]
d = u[3] - u[1]

(
    lam,
    mu,
    g1,
    g2,
    h1,
    h2,
    A0,
    Ab,
    Ac,
    B0,
    Bb,
    Bc,
    C0,
    Cb,
    Cc,
    n0,
    n1,
    n2,
    n3,
    z_saturation,
) = parameters = sp.symbols(
    "lam mu g1 g2 h1 h2 A0 Ab Ac B0 Bb Bc "
    "C0 Cb Cc n0 n1 n2 n3 z_saturation"
)

potential = (
    sp.Rational(5, 6) * b**3
    + lam * b**2 / 2
    + a * ((mu + 2) * c + g1 * b + g2 * b**2)
    + d * (h1 * b + (sp.Rational(1, 4) + h2) * b**2)
    + a**2 * (A0 + Ab * b + Ac * c) / 2
    + a * d * (B0 + Bb * b + Bc * c)
    + d**2 * (C0 + Cb * b + Cc * c) / 2
    + n0 * a**3 / 6
    + n1 * a**2 * d / 2
    + n2 * a * d**2 / 2
    + n3 * d**3 / 6
)

hessian = sp.hessian(sp.expand(potential), u).subs(
    dict(zip(u, q, strict=True)), simultaneous=True
)
candidate_jacobian = (
    m.jacobian(h_variables) + hessian * q.jacobian(h_variables)
).subs({Y: 0, D: 0})

target_determinant = lam * mu**2 / 2
slice_difference = sp.Poly(
    sp.expand(
        candidate_jacobian.det(method="berkowitz") - target_determinant
    ),
    X,
    W,
)
coefficient_equations = [
    coefficient for _, coefficient in slice_difference.terms()
]
assert len(coefficient_equations) == 49
coefficient_equations.append(z_saturation * lam * mu - 1)


def singular_expression(expression: sp.Expr) -> str:
    cleared = sp.Poly(
        expression, *parameters
    ).clear_denoms()[1].as_expr()
    return str(cleared).replace("**", "^")


singular_script = (
    f"ring r=0,({','.join(map(str, parameters))}),dp;\n"
    f"ideal I={','.join(singular_expression(eq) for eq in coefficient_equations)};\n"
    "option(redSB);\n"
    "ideal G=slimgb(I);\n"
    "reduce(1,G);\n"
)
result = subprocess.run(
    ["Singular", "-q"],
    input=singular_script,
    text=True,
    capture_output=True,
    timeout=300,
    check=True,
)
remainder = result.stdout.strip()
assert remainder == "0"

print("PASS: the normalized L=0 cubic family has 49 exact slice equations")
print("PASS: saturation by lambda*mu is the unit ideal over Q")
print("PASS: chart 0010 has no constant determinant in this cubic branch")
print("SCOPE: the generic boundary branch L!=0 remains open")
