#!/usr/bin/env python3
"""Exact symbolic certificates for HC4RSD63."""
from __future__ import annotations

import json
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "artifacts" / "generated-results" / "hc4_rank_two_length_three_closure.json"
OUT.parent.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# 1. Composition law for the binary bordered-Hessian expression.
# ---------------------------------------------------------------------------
hx, hw, hxx, hxw, hww, p1, p2 = sp.symbols(
    "hx hw hxx hxw hww p1 p2"
)
fx = p1 * hx
fw = p1 * hw
fxx = p2 * hx**2 + p1 * hxx
fxw = p2 * hx * hw + p1 * hxw
fww = p2 * hw**2 + p1 * hww
Bf = sp.expand(fw**2 * fxx - 2 * fx * fw * fxw + fx**2 * fww)
Bh = hw**2 * hxx - 2 * hx * hw * hxw + hx**2 * hww
assert sp.factor(Bf - p1**3 * Bh) == 0

# ---------------------------------------------------------------------------
# 2. Full determinant in passive-kernel coordinates.
# ---------------------------------------------------------------------------
# Coordinates are (x,w,y,z), while h=y-q(x,w)z.  Work with abstract jets of
# Phi(x,w,h), q(x,w), and a(x,w), using the Hessian transformation formula.
t = sp.symbols("z")
q, qx, qw, qxx, qxw, qww = sp.symbols("q qx qw qxx qxw qww")
a_x, a_w, a_xx, a_xw, a_ww = sp.symbols("ax aw axx axw aww")
Ph, Pxx, Pxw, Pww, Pxh, Pwh, Phh = sp.symbols(
    "Ph Pxx Pxw Pww Pxh Pwh Phh"
)

# Hessian in moving coordinates (x,w,h,t) for Phi(x,w,h)+a(x,w)t.
H_moving = sp.Matrix(
    [
        [Pxx + a_xx * t, Pxw + a_xw * t, Pxh, a_x],
        [Pxw + a_xw * t, Pww + a_ww * t, Pwh, a_w],
        [Pxh, Pwh, Phh, 0],
        [a_x, a_w, 0, 0],
    ]
)

# Jacobian d(x,w,h,t)/d(x,w,y,z), h=y-qz.
J = sp.Matrix(
    [
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [-qx * t, -qw * t, 1, -q],
        [0, 0, 0, 1],
    ]
)

# Extra Hessian term Ph * Hess(h) caused by the active-dependent shear.
H_h = sp.Matrix(
    [
        [-qxx * t, -qxw * t, 0, -qx],
        [-qxw * t, -qww * t, 0, -qw],
        [0, 0, 0, 0],
        [-qx, -qw, 0, 0],
    ]
)
H_original = sp.expand(J.T * H_moving * J + Ph * H_h)
determinant = sp.expand(H_original.det(method="domain-ge"))
poly_z = sp.Poly(determinant, t)
assert poly_z.degree() <= 1

# The z coefficient is -Phi_hh times the binary bordered-Hessian expression
# of b_lambda=a-lambda q with lambda=Phi_h treated as a scalar under active
# differentiation.
r = sp.Matrix([a_x - Ph * qx, a_w - Ph * qw])
H_b = sp.Matrix(
    [
        [a_xx - Ph * qxx, a_xw - Ph * qxw],
        [a_xw - Ph * qxw, a_ww - Ph * qww],
    ]
)
expected_z = -Phh * (r.T * H_b.adjugate() * r)[0]
assert sp.factor(poly_z.coeff_monomial(t) - expected_z) == 0

# ---------------------------------------------------------------------------
# 3. Once q and a depend on one active coordinate x, the complete determinant
#    factors into two units.
# ---------------------------------------------------------------------------
one_variable = {
    qw: 0,
    qxw: 0,
    qww: 0,
    a_w: 0,
    a_xw: 0,
    a_ww: 0,
}
det_one = sp.factor(determinant.subs(one_variable))
expected_one = -(
    a_x - Ph * qx
) ** 2 * (Pww * Phh - Pwh**2)
assert sp.factor(det_one - expected_one) == 0

result = {
    "scope": "HC4RSD63 complete moving [3,1] closure",
    "status": "exact symbolic identities verified",
    "identities": {
        "binary_composition": "B2(P(h))=P'(h)^3 B2(h)",
        "z_coefficient": "[z] det Hess(psi)=-Phi_hh B2(a-Phi_h q)",
        "one_characteristic_factorization": "det Hess(psi)=-(a'-q' Phi_h)^2 det Hess_{w,h}(Phi)",
    },
}
OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
print(json.dumps(result, indent=2))
