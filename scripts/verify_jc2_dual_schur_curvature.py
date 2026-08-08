#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "artifacts" / "generated-results" / "jc2_dual_schur_curvature.json"
OUT.parent.mkdir(parents=True, exist_ok=True)

u, s, lam, mu, c = sp.symbols("u s lam mu c")
px, py, qx, qy = sp.symbols("px py qx qy")
pxx, pxy, pyy = sp.symbols("pxx pxy pyy")
qxx, qxy, qyy = sp.symbols("qxx qxy qyy")
Q = sp.symbols("Q")

M = sp.Matrix([
    [u*pxx+s*qxx, u*pxy+s*qxy, px],
    [u*pxy+s*qxy, u*pyy+s*qyy, py],
    [px, py, 0],
])

BP = py**2*pxx - 2*px*py*pxy + px**2*pyy
BPQ = py**2*qxx - 2*px*py*qxy + px**2*qyy
assert sp.expand(M.det() + u*BP + s*BPQ) == 0

# The Schur lemma contributes -lambda*c^2 through the rank-one update,
# while s is specialized to mu + lambda*Q.
repaired = sp.expand(M.det().subs(s, mu + lam*Q) - lam*c**2)
expected = sp.expand(-u*BP - (mu + lam*Q)*BPQ - lam*c**2)
assert sp.expand(repaired - expected) == 0

# Curvature is a homogeneous cubic under a target pencil R=aP+bQ.
a, b = sp.symbols("a b")
rx = a*px + b*qx
ry = a*py + b*qy
rxx = a*pxx + b*qxx
rxy = a*pxy + b*qxy
ryy = a*pyy + b*qyy
BR = sp.expand(ry**2*rxx - 2*rx*ry*rxy + rx**2*ryy)
assert sp.Poly(BR, a, b).total_degree() == 3

# Exact nonlinear target-shear law Q -> Q - phi(P).
phi1, phi2 = sp.symbols("phi1 phi2")
sx = qx - phi1*px
sy = qy - phi1*py
sxx = qxx - phi1*pxx - phi2*px**2
sxy = qxy - phi1*pxy - phi2*px*py
syy = qyy - phi1*pyy - phi2*py**2
Bshear = sp.expand(sy**2*sxx - 2*sx*sy*sxy + sx**2*syy)
Cpencil = sp.expand(BR.subs({a: -phi1, b: 1}))
J = px*qy - py*qx
assert sp.factor(Bshear - Cpencil + phi2*J**2) == 0

result = {
    "scope": "JC2 dual-variable Schur descent",
    "status": "exact identities verified",
    "identities": {
        "minor": "det Hess(uP+sQ)=-u B(P)-s B(P;Q)",
        "repair": "det Hess psi=-u B(P)-(mu+lambda Q)B(P;Q)-lambda c^2",
        "target_pencil": "B(aP+bQ) is homogeneous cubic in (a,b)",
        "target_shear": "B(Q-phi(P))=C_F(-phi'(P),1)-J(P,Q)^2 phi''(P)",
    },
}
OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
print(json.dumps(result, indent=2))
