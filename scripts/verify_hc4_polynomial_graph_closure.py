#!/usr/bin/env python3
"""Verify the exact determinant identities in HC4RSD72.

This checker covers the moving ternary-scroll calculation.  The
three-variable quasi-translation classification and the closed/generative
polynomial theorem are external algebraic inputs documented in the note.
"""
from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "hc4_polynomial_graph_closure.json"
)
OUTPUT.parent.mkdir(parents=True, exist_ok=True)

t, r, s, z, H = sp.symbols("t r s z H")
a = sp.Function("a")(t)
b = sp.Function("b")(t)
c = sp.Function("c")(t)
q = sp.Function("q")(t)
P = sp.Function("P")
Q = sp.Function("Q")

h = r + q * s
P_base = P(t, H)
Q_base = Q(t, H)
P_h = sp.diff(P_base, H).subs(H, h)
P_t = sp.diff(P_base, t).subs(H, h)
Q_h = sp.diff(Q_base, H).subs(H, h)
Q_hh = sp.diff(Q_base, H, 2).subs(H, h)

L = s * P(t, h) + Q(t, h)
Hamiltonian = a * r + b * s + c
variables = (t, r, s)

LH = sp.hessian(L, variables)
HH = sp.hessian(Hamiltonian, variables)
poly = sp.Poly(sp.expand((LH - z * HH).det()), z)
coeff_z = poly.coeff_monomial(z)

# Moving-scroll relation b' = a' q and hence b'' = a''q + a'q'.
subs = {
    sp.diff(b, t): sp.diff(a, t) * q,
    sp.diff(b, t, 2): sp.diff(a, t, 2) * q + sp.diff(a, t) * sp.diff(q, t),
}
coeff_z = sp.factor(coeff_z.xreplace(subs))
expected = sp.factor(
    P_h
    * (
        (sp.diff(a, t, 2) * h + sp.diff(c, t, 2)) * P_h
        - 3 * s * sp.diff(a, t) * sp.diff(q, t) * P_h
        - 2 * sp.diff(a, t) * sp.diff(q, t) * Q_h
        - 2 * sp.diff(a, t) * P_t
    )
)
assert sp.simplify(coeff_z - expected) == 0

# After the moving coefficient forces P_h=0, write P=P0(t).
P0 = sp.Function("P0")(t)
L_reduced = s * P0 + Q(t, h)
det_reduced = sp.factor(sp.hessian(L_reduced, variables).det())
expected_det = sp.factor(
    -(sp.diff(P0, t) + sp.diff(q, t) * Q_h) ** 2 * Q_hh
)
assert sp.simplify(det_reduced - expected_det) == 0

result = {
    "scope": "HC4RSD72 polynomial-graph moving ternary scroll",
    "status": "exact determinant identities verified",
    "identities": {
        "linear_pencil_coefficient": "[z] det(L''-zH'') = P_h*((a''h+c'')P_h-3s a'q'P_h-2a'q'Q_h-2a'P_t)",
        "reduced_legendre_determinant": "det L''=-(P'+q'Q_h)^2 Q_hh",
    },
    "external_inputs": [
        "dimension-three quasi-translation classification",
        "existence of closed/generative polynomials",
    ],
}
OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
print(json.dumps(result, indent=2))
