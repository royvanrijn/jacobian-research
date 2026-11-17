#!/usr/bin/env python3
"""Classify all nonzero cubic boundary data in PC(2) chart 0010.

Parameterize the restrictions to the caustic plane a=d=0 of an arbitrary
potential of degree at most three:

    f=V, g=V_a, h=V_d,
    A=V_aa, B=V_ad, C=V_dd.

Affine terms in f and constants in g,h do not affect the Hessian and are
omitted.  The resulting 26 parameters capture all boundary two-jet data;
the four pure-normal cubic coefficients do not occur until the normal
equations.

Set the W^2 and W^1 boundary coefficients to zero and W^0 to a nonzero
constant kappa.  Exact Groebner reduction over Q proves that every field
point lies on the degenerate Schur branch:

    f = 5*b^3/6 + f20*b^2/2,
    g = g10*b + g01*c + g20*b^2/2,
    h = h10*b + h20*b^2/2,
    f20*(g01-2)^2 = 2*kappa != 0.

The companion degenerate-branch checker includes the omitted pure-normal
cubic terms and proves that this family cannot satisfy the full determinant
identity even on the smaller graph slice Y=D=0.
"""

from __future__ import annotations

import runpy
import subprocess

import sympy as sp


boundary = runpy.run_path(
    "scripts/verify_hc4_0010_boundary_schur_chain.py"
)
b = boundary["Y"]
c = boundary["D"]

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
    A0,
    Ab,
    Ac,
    B0,
    Bb,
    Bc,
    C0,
    Cb,
    Cc,
    kappa,
    z_saturation,
) = parameters = sp.symbols(
    "f20 f11 f02 f30 f21 f12 f03 "
    "g10 g01 g20 g11 g02 h10 h01 h20 h11 h02 "
    "A0 Ab Ac B0 Bb Bc C0 Cb Cc kappa z_saturation"
)

f = (
    f20 * b**2 / 2
    + f11 * b * c
    + f02 * c**2 / 2
    + f30 * b**3 / 6
    + f21 * b**2 * c / 2
    + f12 * b * c**2 / 2
    + f03 * c**3 / 6
)
g = (
    g10 * b
    + g01 * c
    + g20 * b**2 / 2
    + g11 * b * c
    + g02 * c**2 / 2
)
h = (
    h10 * b
    + h01 * c
    + h20 * b**2 / 2
    + h11 * b * c
    + h02 * c**2 / 2
)

boundary_substitution = {
    boundary["f_bb"]: sp.diff(f, b, 2),
    boundary["f_bc"]: sp.diff(f, b, c),
    boundary["f_cc"]: sp.diff(f, c, 2),
    boundary["g_b"]: sp.diff(g, b),
    boundary["g_c"]: sp.diff(g, c),
    boundary["h_b"]: sp.diff(h, b),
    boundary["h_c"]: sp.diff(h, c),
    boundary["A"]: A0 + Ab * b + Ac * c,
    boundary["B"]: B0 + Bb * b + Bc * c,
    boundary["C"]: C0 + Cb * b + Cc * c,
}

boundary_polynomials = [
    sp.Poly(
        sp.expand(
            boundary[name].subs(
                boundary_substitution, simultaneous=True
            )
            - (kappa if name == "coefficient_0" else 0)
        ),
        b,
        c,
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

forced_relations = (
    h02,
    h11,
    h01,
    g02,
    f03,
    f12,
    f21,
    f02,
    f11,
    g11**2,
    (f30 - 5) ** 2,
    f20 * (g01 - 2) ** 2 - 2 * kappa,
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
assert remainders == [f"_[{index}]=0" for index in range(1, 13)]

print("PASS: the chart 0010 cubic boundary ideal has 52 exact coefficients")
print("PASS: saturation by kappa forces the degenerate Schur normal form")
print("PASS: every cubic boundary solution has L=0")
print("SCOPE: combine with the degenerate-branch slice certificate")
