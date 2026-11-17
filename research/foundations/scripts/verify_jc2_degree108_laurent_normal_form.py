#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "artifacts" / "generated-results" / "jc2_degree108_laurent_normal_form.json"
OUT.parent.mkdir(parents=True, exist_ok=True)

X, s = sp.symbols("X s", nonzero=True)
A = sp.Function("A")(X)
D = sp.Function("D")(X)
phi = sp.Function("phi")(X)
f = sp.Function("f")(X)
h = sp.Function("h")(X)

# Top equation is 2 A D' - 3 D A' = X^2.
D1 = (X**2 + 3 * D * sp.diff(A, X)) / (2 * A)
D2 = sp.simplify(sp.diff(D1, X).subs(sp.diff(D, X), D1))
D3 = sp.simplify(
    sp.diff(D2, X).subs({sp.diff(D, X, 2): D2, sp.diff(D, X): D1})
)
subs_top = {
    sp.diff(D, X, 3): D3,
    sp.diff(D, X, 2): D2,
    sp.diff(D, X): D1,
}

# First layer.
B = sp.expand(phi * (sp.diff(A, X) - A / X) - sp.diff(phi, X) * A / 2)
E = sp.expand(
    phi * (sp.diff(D, X) - sp.Rational(3, 2) * D / X)
    - sp.Rational(3, 4) * sp.diff(phi, X) * D
)
first_eq = sp.expand(
    -2 * A * sp.diff(E, X)
    - B * sp.diff(D, X)
    + 3 * D * sp.diff(B, X)
    + 2 * E * sp.diff(A, X)
)
assert sp.factor(sp.together(first_eq.subs(subs_top))) == 0
assert sp.expand(
    B + A**3 / (2 * X**2) * sp.diff(X**2 * phi / A**2, X)
) == 0

# Resonant first-layer mode is a target shear.
rho = A**2 / X**2
B_rho = sp.factor(B.subs(phi, rho).doit())
E_rho = sp.factor(sp.together(E.subs(phi, rho).doit().subs(subs_top)))
assert B_rho == 0
assert sp.expand(E_rho - A / 2) == 0

# Density-preserving first source field.
J = -X**2 * s**4
Vx = phi / s
Vs = -(phi / (2 * X) + sp.diff(phi, X) / 4)
assert sp.expand(sp.diff(J * Vx, X) + sp.diff(J * Vs, s)) == 0
assert sp.expand(Vx * sp.diff(A * s**2, X) + Vs * sp.diff(A * s**2, s) - B * s) == 0
assert sp.factor(
    sp.together(
        (Vx * sp.diff(D * s**3, X) + Vs * sp.diff(D * s**3, s) - E * s**2).subs(subs_top)
    )
) == 0

# Particular second layer supplied by 1/2 V^2.
def V(expr: sp.Expr) -> sp.Expr:
    return sp.expand(Vx * sp.diff(expr, X) + Vs * sp.diff(expr, s))

C0 = sp.expand(V(V(A * s**2)) / 2)
F0 = sp.expand((V(V(D * s**3)) / 2).coeff(s, 1))
second_particular = sp.expand(
    -2 * A * sp.diff(F0, X)
    - B * sp.diff(E, X)
    + 3 * D * sp.diff(C0, X)
    + 2 * E * sp.diff(B, X)
    + F0 * sp.diff(A, X)
)
assert sp.factor(sp.together(second_particular.subs(subs_top))) == 0

# Complete homogeneous second layer.
Chat = sp.expand(
    f * sp.diff(A, X)
    - (sp.Rational(2, 3) * sp.diff(f, X) + sp.Rational(4, 3) * f / X) * A
)
Fhat = sp.expand(
    f * sp.diff(D, X) - (sp.diff(f, X) + 2 * f / X) * D
)
second_homogeneous = sp.expand(
    -2 * A * sp.diff(Fhat, X) + 3 * D * sp.diff(Chat, X) + Fhat * sp.diff(A, X)
)
assert sp.factor(sp.together(second_homogeneous.subs(subs_top))) == 0
assert sp.expand(
    Fhat + D**2 / X**2 * sp.diff(X**2 * f / D, X)
) == 0

# Resonant second-layer mode is a target translation.
f_res = D / X**2
Chat_res = sp.factor(sp.together(Chat.subs(f, f_res).doit().subs(subs_top)))
Fhat_res = sp.factor(sp.together(Fhat.subs(f, f_res).doit().subs(subs_top)))
assert sp.expand(Chat_res + sp.Rational(1, 3)) == 0
assert Fhat_res == 0

# Density-preserving second and third source fields.
Wx = f / s**2
Ws = -(sp.diff(f, X) / 3 + 2 * f / (3 * X)) / s
assert sp.expand(sp.diff(J * Wx, X) + sp.diff(J * Ws, s)) == 0

Ux = h / s**3
Us = -(sp.diff(h, X) / 2 + h / X) / s**2
assert sp.expand(sp.diff(J * Ux, X) + sp.diff(J * Us, s)) == 0
P_minus_one = sp.factor(
    sp.expand(Ux * sp.diff(A * s**2, X) + Us * sp.diff(A * s**2, s)).coeff(s, -1)
)
assert sp.expand(P_minus_one - (h * sp.diff(A, X) - (sp.diff(h, X) + 2 * h / X) * A)) == 0

# Elimination of G from the last two equations.
B0, E0, C, F, G = sp.symbols("B0 E0 C F G")
Bp, Cp, Fp, Gp = sp.symbols("Bp Cp Fp Gp")
# s1 gives 2 A G' = -B F' + 2 E C' + F B'.
Gp_expr = (-B0 * Fp + 2 * E0 * Cp + F * Bp) / (2 * A)
residue = sp.factor(-B0 * Gp_expr + F * Cp)
expected_residue = (
    B0**2 * Fp - 2 * B0 * E0 * Cp - B0 * F * Bp + 2 * A * F * Cp
) / (2 * A)
assert sp.expand(residue - expected_residue) == 0

# Universal leading edge.
a8, d12, t = sp.symbols("a8 d12 t", nonzero=True)
b8 = sp.Rational(13, 2) * t * a8
e12 = sp.Rational(39, 4) * t * d12
c8, f12, g12 = sp.symbols("c8 f12 g12")
eq2_top = -2 * a8 * 12 * f12 - b8 * 12 * e12 + 3 * d12 * 8 * c8 + 2 * e12 * 8 * b8 + f12 * 8 * a8
eq1_top = -2 * a8 * 12 * g12 - b8 * 12 * f12 + 2 * e12 * 8 * c8 + f12 * 8 * b8
eq0_top = -b8 * 12 * g12 + f12 * 8 * c8
sol = sp.solve([eq2_top, eq1_top, eq0_top], [c8, f12, g12], dict=True)
assert sol == [{
    c8: sp.Rational(169, 16) * a8 * t**2,
    f12: sp.Rational(507, 16) * d12 * t**2,
    g12: sp.Rational(2197, 64) * d12 * t**3,
}]
shift = sp.Rational(13, 4) * t
assert sp.expand(
    a8 * (s + shift)**2
    - (a8 * s**2 + b8 * s + sol[0][c8])
) == 0
assert sp.expand(
    d12 * (s + shift)**3
    - (d12 * s**3 + e12 * s**2 + sol[0][f12] * s + sol[0][g12])
) == 0

result = {
    "scope": "degree-108 no-vertical-edge JC2 Laurent normal form",
    "status": "graded source-symmetry identities verified",
    "verified": [
        "complete first-layer formulas modulo top equation",
        "first-layer integrating factor and target-shear resonance",
        "density preservation of V_phi",
        "second-layer V_phi^2 particular solution",
        "homogeneous second-layer formulas and target-translation resonance",
        "density preservation of W_f and U_h",
        "absence equation for a regular level-3 action on P",
        "G-eliminated residue",
        "universal square/cube leading edge",
    ],
    "leading_shift": "13*t/4",
    "level3_residue": "B^2 F' - 2 B E C' - B F B' + 2 A F C'",
}
OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
print(json.dumps(result, indent=2))
