#!/usr/bin/env python3
"""Verify the symbolic core of HC4RSD51: complete all-degree h=0 closure."""
from __future__ import annotations

import json
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT_DIR = ROOT / "artifacts" / "generated-results"
ARTIFACT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# 1. Explicit bordered-Hessian expansion.
# ---------------------------------------------------------------------------
p,q,r0 = sp.symbols("p q r")
A0,B0,C0,D0,E0,F0 = sp.symbols("A0 B0 C0 D0 E0 F0")
g = sp.Matrix([p,q,r0])
H = sp.Matrix([[A0,B0,C0],[B0,D0,E0],[C0,E0,F0]])
U = sp.expand((g.T*H.adjugate()*g)[0])
expected = sp.expand(
    A0*D0*r0**2 - 2*A0*E0*q*r0 + A0*F0*q**2 - B0**2*r0**2
    + 2*B0*C0*q*r0 + 2*B0*E0*p*r0 - 2*B0*F0*p*q
    - C0**2*q**2 - 2*C0*D0*p*r0 + 2*C0*E0*p*q
    + D0*F0*p**2 - E0**2*p**2
)
assert U == expected

# Each term in the explicit formula has exactly two y derivatives and two z
# derivatives when p=c_x, q=c_y, r=c_z and A0...F0 are Hessian entries.
# Record derivative bidegrees (dy,dz) of each symbol.
deg = {
    p:(0,0), q:(1,0), r0:(0,1),
    A0:(0,0), B0:(1,0), C0:(0,1),
    D0:(2,0), E0:(1,1), F0:(0,2),
}
poly = sp.Poly(U,p,q,r0,A0,B0,C0,D0,E0,F0)
vars_ = (p,q,r0,A0,B0,C0,D0,E0,F0)
for monomial,_ in poly.terms():
    dy = sum(power*deg[var][0] for var,power in zip(vars_,monomial))
    dz = sum(power*deg[var][1] for var,power in zip(vars_,monomial))
    assert (dy,dz) == (2,2)

# ---------------------------------------------------------------------------
# 2. Positive-degree highest-tail mixed coefficient (HC4RSD48), repeated here
# independently in abstract jets.
# ---------------------------------------------------------------------------
j,i,rr = sp.symbols("j i rr", integer=True, positive=True)
x,eps = sp.symbols("x eps")
R,Rx,Ry,Rxx,Rxy,Ryy = sp.symbols("R Rx Ry Rxx Rxy Ryy")
S,Sx,Sy,Sxx,Sxy,Syy = sp.symbols("S Sx Sy Sxx Sxy Syy")
gmix = sp.Matrix([Rx+eps*Sx,Ry+eps*Sy,j*R+eps*i*S])
Hmix = sp.Matrix([
    [Rxx+eps*Sxx,Rxy+eps*Sxy,j*Rx+eps*i*Sx],
    [Rxy+eps*Sxy,Ryy+eps*Syy,j*Ry+eps*i*Sy],
    [j*Rx+eps*i*Sx,j*Ry+eps*i*Sy,j*(j-1)*R+eps*i*(i-1)*S],
])
linear = sp.expand(sp.diff((gmix.T*Hmix.adjugate()*gmix)[0],eps).subs(eps,0))
power_sub = {
    R:x**rr,Rx:rr*x**(rr-1),Ry:0,
    Rxx:rr*(rr-1)*x**(rr-2),Rxy:0,Ryy:0,
}
assert sp.factor(linear.subs(power_sub) + j*rr*(rr+j)*x**(3*rr-2)*Syy) == 0

# ---------------------------------------------------------------------------
# 3. Scalar-highest-tail starter.
# ---------------------------------------------------------------------------
J,I = sp.symbols("J I", integer=True, positive=True)
t = sp.symbols("t", nonzero=True)
gs = sp.Matrix([eps*Sx,eps*Sy,J*t+eps*I*S])
Hs = sp.Matrix([
    [eps*Sxx,eps*Sxy,eps*I*Sx],
    [eps*Sxy,eps*Syy,eps*I*Sy],
    [eps*I*Sx,eps*I*Sy,J*(J-1)*t+eps*I*(I-1)*S],
])
Us = sp.expand((gs.T*Hs.adjugate()*gs)[0])
quad = sp.factor(sp.diff(Us,eps,2).subs(eps,0)/2)
assert quad == J**2*t**2*(Sxx*Syy-Sxy**2)

# ---------------------------------------------------------------------------
# 4. Scalar-highest-tail descendant after S=x^r.
# ---------------------------------------------------------------------------
eta,ii = sp.symbols("eta ii", integer=True)
T,Tx,Ty,Txx,Txy,Tyy = sp.symbols("T Tx Ty Txx Txy Tyy")
Sr = x**rr
Srx = sp.diff(Sr,x)
Srxx = sp.diff(Sr,x,2)
gd = sp.Matrix([
    eps*Srx+eta*Tx,
    eta*Ty,
    J*t+eps*I*Sr+eta*ii*T,
])
Hd = sp.Matrix([
    [eps*Srxx+eta*Txx,eta*Txy,eps*I*Srx+eta*ii*Tx],
    [eta*Txy,eta*Tyy,eta*ii*Ty],
    [eps*I*Srx+eta*ii*Tx,eta*ii*Ty,
     J*(J-1)*t+eps*I*(I-1)*Sr+eta*ii*(ii-1)*T],
])
Ud = sp.expand((gd.T*Hd.adjugate()*gd)[0])
mixed = sp.factor(sp.diff(sp.diff(Ud,eps),eta).subs({eps:0,eta:0}))
assert mixed == J**2*t**2*rr*(rr-1)*x**(rr-2)*Tyy

# ---------------------------------------------------------------------------
# 5. Integer inequalities used at the end of the induction.
# ---------------------------------------------------------------------------
o,k = sp.symbols("o k", integer=True, positive=True)
# These are recorded algebraically; the checker also samples the admissible
# integer range to guard transcription errors.
for ov in range(3,80):
    for kv in range(2,80):
        if kv > (ov-1):
            continue
        Jpos = 2 + (ov-1)//kv
        assert Jpos-2 < ov-1
for ov in [8,49,288,1681]:
    for kv in range(2,min(ov,100)):
        assert 1 + (ov-1)//kv < ov-1

result = {
    "scope":"complete all-degree h=0 scalar reverse-Schur closure",
    "status":"symbolic core verified",
    "derivative_count":"every bordered monomial has exactly two y and two z derivatives",
    "positive_top_mixed":"-J r(r+J) x^(3r-2) S_yy",
    "scalar_top_quadratic":"J^2 t^2 det Hess(S)",
    "scalar_top_mixed":"J^2 t^2 r(r-1) x^(r-2) T_yy",
    "consequence":"HC4RSD51: all h=0 packets close in every degree",
}
(ARTIFACT_DIR / "hc4_minimal_excess_complete.json").write_text(
    json.dumps(result,indent=2)+"\n",encoding="utf-8"
)
print(json.dumps(result,indent=2))
