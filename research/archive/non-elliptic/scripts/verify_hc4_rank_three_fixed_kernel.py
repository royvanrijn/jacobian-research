#!/usr/bin/env python3
"""Verify exact matrix identities in HC4RSD65."""
from __future__ import annotations

import json
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "artifacts" / "generated-results" / "hc4_rank_three_fixed_kernel.json"
OUT.parent.mkdir(parents=True, exist_ok=True)

s, X = sp.symbols("s X")
t = sp.symbols("t")
b0, b1, b2 = sp.symbols("b0 b1 b2")
c00, c01, c02, c11, c12, c22 = sp.symbols("c00 c01 c02 c11 c12 c22")
h00, h01, h02, h11, h12, h22 = sp.symbols("h00 h01 h02 h11 h12 h22")

b = sp.Matrix([b0, b1, b2])
C = sp.Matrix([[c00,c01,c02],[c01,c11,c12],[c02,c12,c22]])
H = sp.Matrix([[h00,h01,h02],[h01,h11,h12],[h02,h12,h22]])
M = sp.Matrix.vstack(
    sp.Matrix.hstack(sp.Matrix([[t]]), b.T),
    sp.Matrix.hstack(b, C + s * H),
)
detM = sp.expand(M.det(method="domain-ge"))
assert sp.factor(sp.Poly(detM, s).coeff_monomial(s**3) - t * H.det()) == 0

# Once t=0, this is exactly the scalar bordered-Hessian block.
M0 = M.subs(t, 0)
assert sp.expand(M0.det(method="domain-ge") + (b.T * (C+s*H).adjugate() * b)[0]) == 0

# If the lower block is L + X*K, the X^2 coefficient is the UFE obstruction.
k00, k01, k02, k11, k12, k22 = sp.symbols("k00 k01 k02 k11 k12 k22")
K = sp.Matrix([[k00,k01,k02],[k01,k11,k12],[k02,k12,k22]])
L = C
MX = sp.Matrix.vstack(
    sp.Matrix.hstack(sp.zeros(1,1), b.T),
    sp.Matrix.hstack(b, L + X*K),
)
detX = sp.Poly(sp.expand(MX.det(method="domain-ge")), X)
assert sp.factor(detX.coeff_monomial(X**2) + (b.T * K.adjugate() * b)[0]) == 0

result = {
    "scope": "fixed top kernel in rank-three [4] HC4 stratum",
    "status": "exact block reduction verified",
    "identities": {
        "top_s_coefficient": "[s^3] det(S+sT)=psi_xx det(Hess A)",
        "scalar_block": "det[[0,p^T],[p,L]]=-p^T adj(L)p",
        "ufe_coefficient": "[x^2] det=-p^T adj(Hess P)p",
    },
}
OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
print(json.dumps(result, indent=2))
