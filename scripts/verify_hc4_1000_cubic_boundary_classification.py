#!/usr/bin/env python3
"""Classify all nonzero cubic boundary data in PC(2) chart 1000.

Parameterize the restrictions to the caustic image c=d=0 of an arbitrary
potential of degree at most three:

    f=V, g=V_c, h=V_d,
    V_cc, V_cd, V_dd.

Affine terms in f and constants in g,h are omitted because they do not
affect the Hessian.  The resulting 26 parameters capture all boundary
two-jet data; four pure-normal cubic coefficients occur only away from the
boundary.

Set the W^2 and W^1 boundary coefficients to zero and W^0 to a nonzero
constant kappa.  Exact Groebner reduction over Q proves that every field
point has

    f = 5*b^3/6 + f02*b^2/2,
    g = g10*a + g01*b + g02*b^2/2,
    h = h01*b + h02*b^2/2,
    f02*(2*g10-1)^2 = 2*kappa != 0.

Consequently every cubic boundary solution lies on L=0.  The companion
full-graph checker includes all pure-normal cubic terms and eliminates this
family using the Y=0 slice and two exact nonzero-Y graph points.
"""

from __future__ import annotations

import runpy
import subprocess

import sympy as sp


boundary = runpy.run_path(
    "scripts/verify_hc4_1000_boundary_schur_chain.py"
)
a = boundary["a"]
b = boundary["b"]

(
    f20,
    f11,
    f02,
    f30,
    f21,
    f12,
    f03,
    g10,
    g01,
    g20,
    g11,
    g02,
    h10,
    h01,
    h20,
    h11,
    h02,
    C0,
    Ca,
    Cb,
    N0,
    Na,
    Nb,
    R0,
    Ra,
    Rb,
    kappa,
    z_saturation,
) = parameters = sp.symbols(
    "f20 f11 f02 f30 f21 f12 f03 "
    "g10 g01 g20 g11 g02 h10 h01 h20 h11 h02 "
    "C0 Ca Cb N0 Na Nb R0 Ra Rb kappa z_saturation"
)

f = (
    f20 * a**2 / 2
    + f11 * a * b
    + f02 * b**2 / 2
    + f30 * a**3 / 6
    + f21 * a**2 * b / 2
    + f12 * a * b**2 / 2
    + f03 * b**3 / 6
)
g = (
    g10 * a
    + g01 * b
    + g20 * a**2 / 2
    + g11 * a * b
    + g02 * b**2 / 2
)
h = (
    h10 * a
    + h01 * b
    + h20 * a**2 / 2
    + h11 * a * b
    + h02 * b**2 / 2
)

boundary_substitution = {
    boundary["f_aa"]: sp.diff(f, a, 2),
    boundary["f_ab"]: sp.diff(f, a, b),
    boundary["f_bb"]: sp.diff(f, b, 2),
    boundary["g_a"]: sp.diff(g, a),
    boundary["g_b"]: sp.diff(g, b),
    boundary["h_a"]: sp.diff(h, a),
    boundary["h_b"]: sp.diff(h, b),
    boundary["V_cc"]: C0 + Ca * a + Cb * b,
    boundary["V_cd"]: N0 + Na * a + Nb * b,
    boundary["V_dd"]: R0 + Ra * a + Rb * b,
}

boundary_polynomials = [
    sp.Poly(
        sp.expand(
            boundary[name].subs(
                boundary_substitution, simultaneous=True
            )
            - (kappa if name == "coefficient_0" else 0)
        ),
        a,
        b,
    )
    for name in ("coefficient_2", "coefficient_1", "coefficient_0")
]
coefficient_equations = [
    coefficient
    for polynomial in boundary_polynomials
    for _, coefficient in polynomial.terms()
]
assert len(coefficient_equations) == 52
coefficient_equations.append(z_saturation * kappa - 1)

# These are selected elements of the exact Groebner basis.  The nilpotent
# equations suffice for statements about field-valued solutions.
forced_relations = (
    h11,
    h20,
    h10,
    g20,
    f12,
    f21,
    f30,
    f11,
    f20,
    g11**2,
    f03 * g11 - 5 * g11,
    2 * f03 * g10 + 4 * f02 * g11 - f03 - 10 * g10 + 5,
    (f03 - 5) ** 2,
    f02 * (2 * g10 - 1) ** 2 - 2 * kappa,
    4 * f02**2 * g10 * g11
    - 2 * f02**2 * g11
    + (f03 - 5) * kappa,
)


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
    f"ideal R={','.join(singular_expression(eq) for eq in forced_relations)};\n"
    "reduce(R,G);\n"
)
result = subprocess.run(
    ["Singular", "-q"],
    input=singular_script,
    text=True,
    capture_output=True,
    timeout=300,
    check=True,
)
remainders = [
    line.strip()
    for line in result.stdout.splitlines()
    if line.strip()
]
assert remainders == [f"_[{index}]=0" for index in range(1, 16)]

print("PASS: the chart 1000 cubic boundary ideal has 52 coefficients")
print("PASS: saturation by kappa forces the degenerate Schur normal form")
print("PASS: every field-valued cubic boundary solution has L=0")
print("SCOPE: combine with the full-graph cubic-family certificate")
