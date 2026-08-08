#!/usr/bin/env python3
"""Verify the exact determinant identities in HC4RSD72 by generic 2-jets.

The determinant calculation only uses second derivatives, so arbitrary
quadratic jets at a generic point give an exact finite symbolic certificate.
The three-variable quasi-translation classification and the closed/generative
polynomial theorem are external inputs documented in the companion note.
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
a0, a1, a2, c0, c1, c2, q0, q1, q2 = sp.symbols(
    "a0 a1 a2 c0 c1 c2 q0 q1 q2"
)
b0 = sp.symbols("b0")
p0, pt, ph, ptt, pth, phh = sp.symbols("p0 pt ph ptt pth phh")
qv0, qt, qh, qtt, qth, qhh = sp.symbols("Q0 Qt Qh Qtt Qth Qhh")

# Generic second jets at t=0.
a = a0 + a1 * t + a2 * t**2 / 2
q = q0 + q1 * t + q2 * t**2 / 2
c = c0 + c1 * t + c2 * t**2 / 2
# b'=a'q and b''=a''q+a'q' at the chosen point.
b = b0 + a1 * q0 * t + (a2 * q0 + a1 * q1) * t**2 / 2
P = p0 + pt * t + ph * H + ptt * t**2 / 2 + pth * t * H + phh * H**2 / 2
Q = qv0 + qt * t + qh * H + qtt * t**2 / 2 + qth * t * H + qhh * H**2 / 2

h = r + q * s
L = sp.expand(s * P.subs(H, h) + Q.subs(H, h))
Hamiltonian = sp.expand(a * r + b * s + c)
variables = (t, r, s)

pencil = sp.expand((sp.hessian(L, variables) - z * sp.hessian(Hamiltonian, variables)).det())
coeff_z = sp.expand(sp.Poly(pencil, z).coeff_monomial(z).subs(t, 0))

H0 = r + q0 * s
P_h = ph + phh * H0
P_t = pt + pth * H0
Q_h = qh + qhh * H0
expected = sp.expand(
    P_h
    * (
        (a2 * H0 + c2) * P_h
        - 3 * s * a1 * q1 * P_h
        - 2 * a1 * q1 * Q_h
        - 2 * a1 * P_t
    )
)
assert sp.factor(coeff_z - expected) == 0

# After the moving coefficient forces P_h=0, write P=P0(t).
P_reduced = p0 + pt * t + ptt * t**2 / 2
L_reduced = sp.expand(s * P_reduced + Q.subs(H, h))
det_reduced = sp.expand(sp.hessian(L_reduced, variables).det().subs(t, 0))
expected_det = sp.expand(-(pt + q1 * Q_h) ** 2 * qhh)
assert sp.factor(det_reduced - expected_det) == 0

result = {
    "scope": "HC4RSD72 polynomial-graph moving ternary scroll",
    "status": "exact generic-2-jet determinant identities verified",
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
