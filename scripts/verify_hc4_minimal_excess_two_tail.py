#!/usr/bin/env python3
"""Verify the all-degree h=0 two-tail closure identities."""
from __future__ import annotations

import json
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT_DIR = ROOT / "artifacts" / "generated-results"
ARTIFACT_DIR.mkdir(parents=True, exist_ok=True)

x,y,z = sp.symbols("x y z")
r,k,alpha,beta,qc,t = sp.symbols("r k alpha beta qc t")
B = sp.Function("B")(x,y)
P = sp.Function("P")(x,y)


def bordered(c):
    grad = sp.Matrix([sp.diff(c,v) for v in (x,y,z)])
    hess = sp.hessian(c,(x,y,z))
    return sp.Poly(sp.expand((grad.T*hess.adjugate()*grad)[0]),z)

# ---------------------------------------------------------------------------
# 1. Strict strip 5o/3 < d < 2o.
# Use unnormalized z-tail coefficients; factorials only multiply equations by
# nonzero constants and do not change their zero sets.
# ---------------------------------------------------------------------------
R4 = x**r
R3 = x**(r+k-1)*(alpha*x+beta*y)
c = R4*z**4 + R3*z**3 + qc*B*z**2
J = bordered(c)

D = 16*k**2 - 32*k - 3*r**2 - 12*r + 16
c12 = sp.factor(J.coeff_monomial(z**12))
expected12 = -x**(2*r-2) * (
    4*qc*r*(r+4)*x**r*sp.diff(B,y,2)
    + D*beta**2*x**(2*k+2*r-2)
)
assert sp.factor(c12-expected12) == 0
assert sp.factor(D.subs(k,r+1) - r*(13*r-12)) == 0
assert sp.factor(D-r*(13*r-12)-16*(k-r-1)*(k+r-1)) == 0

# k>r>=1 => D>0.  Once beta=0 and qc=0, the highest coefficient linear in
# the next lower P*z is independent of alpha.
pp = sp.symbols("pp")
c0 = x**r*z**4 + alpha*x**(r+k)*z**3 + pp*P*z
expr = bordered(c0).as_expr()
linear = sp.Poly(sp.expand(sp.diff(expr,pp).subs(pp,0)),z)
assert sp.factor(
    linear.LC() + 4*r*(r+4)*x**(3*r-2)*sp.diff(P,y,2)
) == 0

# ---------------------------------------------------------------------------
# 2. Boundary d=5o/3, so o=3k and d=5k.
# ---------------------------------------------------------------------------
R4b = x**k
R3b = x**(2*k-1)*(alpha*x+beta*y)
cb = t*z**5 + R4b*z**4 + R3b*z**3 + qc*B*z**2
Jb = bordered(cb)

c14 = sp.factor(Jb.coeff_monomial(z**14))
expected14 = 25*t**2/x**2 * (
    x**k*k*(k-1)*qc*sp.diff(B,y,2)
    - (2*k-1)**2*beta**2*x**(4*k-2)
)
assert sp.factor(c14-expected14) == 0

# With qc=0, beta=0, the next lower P*z has a nonzero P_yy coefficient.
c0b = t*z**5 + x**k*z**4 + alpha*x**(2*k)*z**3 + pp*P*z
exprb = bordered(c0b).as_expr()
linearb = sp.Poly(sp.expand(sp.diff(exprb,pp).subs(pp,0)),z)
assert sp.factor(
    linearb.LC() - 25*k*(k-1)*t**2*x**(k-2)*sp.diff(P,y,2)
) == 0

result = {
    "scope": "h=0 two-tail scalar reverse-Schur closure",
    "closed_region": "d >= 5o/3",
    "strict_strip": "5o/3 < d < 2o",
    "boundary": "d = 5o/3",
    "recurrence_candidate": "deg_y R_(J-l) <= l",
}
(ARTIFACT_DIR / "hc4_minimal_excess_two_tail.json").write_text(
    json.dumps(result,indent=2)+"\n",encoding="utf-8"
)
print(json.dumps(result,indent=2))
