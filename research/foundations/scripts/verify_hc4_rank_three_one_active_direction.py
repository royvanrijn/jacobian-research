#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "artifacts" / "generated-results" / "hc4_rank_three_one_active_direction.json"
OUT.parent.mkdir(parents=True, exist_ok=True)

x,y,z,w,s,U = sp.symbols("x y z w s U")
F = sp.Function("F")
R = sp.Function("R")
P = sp.Function("P")
Q = sp.Function("Q")
H = sp.Function("H")

A = x*y + F(x)*z + R(x,w)
vars = (x,y,z,w)
T = sp.hessian(A, vars)
f1 = sp.diff(F(x), x)
u = y + f1*z

# Top equation D^2 psi=0 is integrated as psi=z*P(x,u,w)+Q(x,u,w).
psi = z*P(x,u,w) + Q(x,u,w)
S = sp.hessian(psi, vars)
poly = sp.Poly(sp.expand((S+s*T).det()), s)
c2 = sp.factor(sp.simplify(poly.coeff_monomial(s**2).subs(y, U-f1*z).doit()))
expected_c2 = (
    3*z*sp.diff(F(x),x,2)*sp.diff(R(x,w),w,2)*sp.diff(P(x,U,w),U)**2
    +2*sp.diff(F(x),x,2)*sp.diff(R(x,w),w,2)*sp.diff(P(x,U,w),U)*sp.diff(Q(x,U,w),U)
    -sp.diff(R(x,w),w,2)*sp.diff(R(x,w),x,2)*sp.diff(P(x,U,w),U)**2
    +sp.diff(R(x,w),x,w)**2*sp.diff(P(x,U,w),U)**2
    -2*sp.diff(R(x,w),x,w)*sp.diff(P(x,U,w),U)*sp.diff(P(x,U,w),w)
    +2*sp.diff(R(x,w),w,2)*sp.diff(P(x,U,w),U)*sp.diff(P(x,U,w),x)
    +sp.diff(P(x,U,w),w)**2
)
assert sp.expand(c2-expected_c2) == 0

# After c2 forces P=P(x), compute the final two identities.
psi2 = z*H(x) + Q(x,u,w)
S2 = sp.hessian(psi2, vars)
poly2 = sp.Poly(sp.expand((S2+s*T).det()), s)
c1 = sp.factor(sp.simplify(poly2.coeff_monomial(s).subs(y, U-f1*z).doit()))
c0 = sp.factor(sp.simplify(poly2.coeff_monomial(1).subs(y, U-f1*z).doit()))
expected_c1 = -(
    sp.diff(F(x),x,2)*sp.diff(Q(x,U,w),U) + sp.diff(H(x),x)
)**2 * sp.diff(Q(x,U,w),U,2) * sp.diff(R(x,w),w,2)
expected_c0 = -(
    sp.diff(F(x),x,2)*sp.diff(Q(x,U,w),U) + sp.diff(H(x),x)
)**2 * (
    sp.diff(Q(x,U,w),U,2)*sp.diff(Q(x,U,w),w,2)
    - sp.diff(Q(x,U,w),U,w)**2
)
assert sp.expand(c1-expected_c1) == 0
assert sp.expand(c0-expected_c0) == 0

result = {
    "scope": "rank-three [4] one-active-direction family A=xy+F(x)z+R(x,w)",
    "status": "moving kernel impossible",
    "conclusion": "F''=0",
    "identities": {
        "s2": "coefficient of z is 3 F'' R_ww P_u^2; then P_w^2",
        "s1": "-(F'' Q_u+H')^2 Q_uu R_ww",
        "detS": "-(F'' Q_u+H')^2 (Q_uu Q_ww-Q_uw^2)"
    }
}
OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
print(json.dumps(result, indent=2))
