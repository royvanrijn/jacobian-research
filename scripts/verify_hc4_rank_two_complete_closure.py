#!/usr/bin/env python3
"""Exact symbolic certificates for HC4RSD59--HC4RSD60."""
from __future__ import annotations

import json
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "artifacts" / "generated-results" / "hc4_rank_two_complete_closure.json"
OUT.parent.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# 1. Constant-coordinate cotangent determinant.
# ---------------------------------------------------------------------------
x, y, z, w = sp.symbols("x y z w")
P = sp.Function("P")(x, w)
Q = sp.Function("Q")(x, w)
R = sp.Function("R")(x, w)
psi = y * P + z * Q + R
H = sp.hessian(psi, (x, y, z, w))
J = sp.diff(P, x) * sp.diff(Q, w) - sp.diff(P, w) * sp.diff(Q, x)
assert sp.factor(H.det() - J**2) == 0

# ---------------------------------------------------------------------------
# 2. Exceptional moving-frame adjugate and passive-variable obstruction.
# ---------------------------------------------------------------------------
h, p, q = sp.symbols("h p q")
HY = sp.Matrix([[h, p, q], [p, 0, 0], [q, 0, 0]])
r = sp.Matrix([0, q, -p])
assert sp.simplify(HY.adjugate() + r * r.T) == sp.zeros(3)

# A generic connection Omega and generic gradient of a passive-linear
# potential.  Only r^T n is needed; the H*Omega*Y contribution vanishes after
# contraction with r because r^T H=0.
omega = sp.symbols("o11:14 o21:24 o31:34")
Omega = sp.Matrix(3, 3, omega)
v, t = sp.symbols("v t")
Pu, Qu, Ru, Pv, Qv = sp.symbols("Pu Qu Ru Pw Qw")
g = sp.Matrix([v * Pu + t * Qu + Ru, sp.symbols("P0"), sp.symbols("Q0")])
gw = sp.Matrix([sp.symbols("g1w"), Pv, Qv])
s = Omega * r
contracted = sp.expand((r.T * gw)[0] - (s.T * g)[0])
assert sp.expand(contracted).coeff(v) == -s[0] * Pu
assert sp.expand(contracted).coeff(t) == -s[0] * Qu

# If the kernel ratio q:p genuinely varies with the active coordinate, the
# identity s_1 = omega_12*q - omega_13*p = 0 forces both connection entries to
# vanish.  The checker records the exact linear expression used in the proof.
assert sp.expand(s[0] - (Omega[0, 1] * q - Omega[0, 2] * p)) == 0

result = {
    "scope": "HC4RSD59--60 moving rank-two [2,2] closure",
    "status": "exact symbolic identities verified",
    "identities": {
        "cotangent_determinant": "det Hess(yP+zQ+R)=J(P,Q)^2",
        "exceptional_adjugate": "adj([[*,p,q],[p,0,0],[q,0,0]])=-r r^T",
        "passive_coefficients": "coeff_v(r^T n)=-(Omega r)_1 P_u and coeff_t(r^T n)=-(Omega r)_1 Q_u",
    },
}
OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
print(json.dumps(result, indent=2))
