#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "artifacts" / "generated-results" / "hc4_rank_three_wronskian_generic_jet.json"
OUT.parent.mkdir(parents=True, exist_ok=True)

x,w,y,z,s,U = sp.symbols("x w y z s U")
Axx, Axw, Aww = sp.symbols("Axx Axw Aww")
Rxx, Rxw, Rww = sp.symbols("Rxx Rxw Rww")
beta, gamma = sp.symbols("beta gamma")

# Universal normalized active two-jet: dh=(1,0) at the origin.
h = x + sp.Rational(1,2)*Axx*x**2 + Axw*x*w + sp.Rational(1,2)*Aww*w**2
# Normalize p'=1,q'=0.  beta is the first projective derivative/Wronskian.
F = sp.Rational(1,2)*beta*h**2 + sp.Rational(1,6)*gamma*h**3
R = sp.Rational(1,2)*Rxx*x**2 + Rxw*x*w + sp.Rational(1,2)*Rww*w**2
A = h*y + F*z + R
vars = (x,w,y,z)
T = sp.hessian(A, vars)

P = sp.Function("P")
Q = sp.Function("Q")
r = beta*h + sp.Rational(1,2)*gamma*h**2
u = y + r*z
psi = z*P(x,w,u) + Q(x,w,u)
S = sp.hessian(psi, vars)

# Freeze at a generic active point after differentiating, but retain passive U,z.
S0 = S.subs({x:0, w:0, y:U}).doit()
T0 = T.subs({x:0, w:0, y:U}).doit()
poly = sp.Poly(sp.expand((S0+s*T0).det(method="berkowitz")), s)
c2 = sp.factor(poly.coeff_monomial(s**2))
Pu = sp.diff(P(0,0,U), U)
Pw = sp.Subs(sp.diff(P(0,sp.Symbol('_xi_2'),U), sp.Symbol('_xi_2')), sp.Symbol('_xi_2'), 0)
# Avoid relying on SymPy's dummy-symbol spelling: identify the z coefficient directly.
assert sp.factor(sp.diff(c2, z) - 3*beta*(Aww*U+Rww)*Pu**2) == 0
# Once P_u=0, every mixed derivative containing u vanishes; the surviving
# transverse square is P_w^2.  Recompute with P=P(x) at the next stage below.

H = sp.Function("H")
psi2 = z*H(h) + Q(x,w,u)
S2 = sp.hessian(psi2, vars)
S20 = S2.subs({x:0,w:0,y:U}).doit()
poly2 = sp.Poly(sp.expand((S20+s*T0).det(method="berkowitz")), s)
c1 = sp.factor(poly2.coeff_monomial(s))
c0 = sp.factor(poly2.coeff_monomial(1))
Qu = sp.diff(Q(0,0,U), U)
Quu = sp.diff(Q(0,0,U), U, 2)
H1 = sp.diff(H(sp.Symbol('hh')), sp.Symbol('hh')).subs(sp.Symbol('hh'),0)
# SymPy may represent H'(0) with a Subs object; get it directly from c1 by
# constructing the same object through differentiation before substitution.
H1 = sp.diff(H(h), x).subs({x:0,w:0})
expected_c1 = -(Aww*U+Rww)*(beta*Qu+H1)**2*Quu
assert sp.simplify(c1-expected_c1) == 0

# Imposing Q_uu=0 leaves the universal constant-term square.
c0_reduced = sp.factor(c0.subs(Quu,0))
Quw = sp.diff(Q(x,w,U), U, w).subs({x:0,w:0})
expected_c0 = (beta*Qu+H1)**2*Quw**2
assert sp.simplify(c0_reduced-expected_c0) == 0

result = {
    "scope": "generic nonlinear-generator jet for moving rank-three [4]",
    "status": "Wronskian identities verified",
    "identities": {
        "s2_t_coefficient": "3 Theta Lambda C_u^2",
        "s1": "-Lambda E_uu (Theta E_u+B)^2",
        "det_after_Euu_zero": "(E_u_tau)^2 (Theta E_u+B)^2"
    },
    "parameters": {
        "Theta": "beta",
        "Lambda": "Aww*U+Rww"
    }
}
OUT.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
print(json.dumps(result, indent=2))
