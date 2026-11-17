#!/usr/bin/env python3
"""Exact full-graph obstruction for the cubic boundary family in chart 1000.

The boundary classification reduces every degree-at-most-three potential
with nonzero constant determinant to

    V = 5*b^3/6 + lambda*b^2/2
      + c*((mu+1)*a/2 + g1*b + g2*b^2)
      + d*(h1*b + h2*b^2)
      + c^2*(C0+Ca*a+Cb*b)/2
      + c*d*(N0+Na*a+Nb*b)
      + d^2*(R0+Ra*a+Rb*b)/2
      + n0*c^3/6 + n1*c^2*d/2 + n2*c*d^2/2 + n3*d^3/6,

where u=(a,b,c,b+d), lambda*mu != 0, and the target determinant is
lambda*mu^2/2.

The full determinant identity is first restricted to Y=0, with X,W,D
left free.  Its 185 coefficients plus saturation have a 22-element exact
Groebner basis over Q.  Quotienting by that basis and adding determinant
evaluations at (X,Y,W,D)=(1,1,0,0) and (1,1,1,0) gives the unit ideal.
Therefore this complete cubic boundary family cannot satisfy the global
constant-determinant identity.
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

mask = (1, 0, 0, 0)
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
    C0,
    Ca,
    Cb,
    N0,
    Na,
    Nb,
    R0,
    Ra,
    Rb,
    n0,
    n1,
    n2,
    n3,
    z_saturation,
) = parameters = sp.symbols(
    "lam mu g1 g2 h1 h2 C0 Ca Cb N0 Na Nb "
    "R0 Ra Rb n0 n1 n2 n3 z_saturation"
)

potential = (
    sp.Rational(5, 6) * b**3
    + lam * b**2 / 2
    + c * ((mu + 1) * a / 2 + g1 * b + g2 * b**2)
    + d * (h1 * b + h2 * b**2)
    + c**2 * (C0 + Ca * a + Cb * b) / 2
    + c * d * (N0 + Na * a + Nb * b)
    + d**2 * (R0 + Ra * a + Rb * b) / 2
    + n0 * c**3 / 6
    + n1 * c**2 * d / 2
    + n2 * c * d**2 / 2
    + n3 * d**3 / 6
)

hessian = sp.hessian(sp.expand(potential), u)
jq = q.jacobian(h_variables)
jm = m.jacobian(h_variables)
target_determinant = lam * mu**2 / 2


def candidate_jacobian(substitution: dict[sp.Symbol, int]) -> sp.Matrix:
    q_value = [coordinate.subs(substitution) for coordinate in q]
    evaluated_hessian = hessian.subs(
        dict(zip(u, q_value, strict=True)), simultaneous=True
    )
    return jm.subs(substitution) + evaluated_hessian * jq.subs(substitution)


slice_jacobian = (
    jm
    + hessian.subs(dict(zip(u, q, strict=True)), simultaneous=True) * jq
).subs(Y, 0)
slice_difference = sp.Poly(
    sp.expand(
        slice_jacobian.det(method="berkowitz") - target_determinant
    ),
    X,
    W,
    D,
)
slice_equations = [
    coefficient for _, coefficient in slice_difference.terms()
]
assert len(slice_equations) == 185
slice_equations.append(z_saturation * lam * mu - 1)

evaluation_points = ((1, 1, 0, 0), (1, 1, 1, 0))
point_equations = []
for point in evaluation_points:
    substitution = dict(zip(h_variables, point, strict=True))
    point_equations.append(
        sp.expand(
            candidate_jacobian(substitution).det(method="berkowitz")
            - target_determinant
        )
    )


def singular_expression(expression: sp.Expr) -> str:
    cleared = sp.Poly(
        expression, *parameters
    ).clear_denoms()[1].as_expr()
    return str(cleared).replace("**", "^")


singular_script = (
    f"ring r=0,({','.join(map(str, parameters))}),dp;\n"
    f"ideal I={','.join(singular_expression(eq) for eq in slice_equations)};\n"
    "option(redSB);\n"
    "ideal G=slimgb(I);\n"
    f"ideal J=G,{','.join(singular_expression(eq) for eq in point_equations)};\n"
    "ideal H=slimgb(J);\n"
    "reduce(1,H);\n"
    "size(G);\n"
)
result = subprocess.run(
    ["Singular", "-q"],
    input=singular_script,
    text=True,
    capture_output=True,
    timeout=300,
    check=True,
)
output = [line.strip() for line in result.stdout.splitlines() if line.strip()]
assert output == ["0", "22"]

print("PASS: the normalized chart 1000 cubic family has 185 Y=0 equations")
print("PASS: the saturated Y=0 ideal has a 22-element Groebner basis")
print("PASS: two exact nonzero-Y graph points make the quotient unit over Q")
print("PASS: chart 1000 has no constant determinant in this cubic branch")
