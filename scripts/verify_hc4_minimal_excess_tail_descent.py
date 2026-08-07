#!/usr/bin/env python3
"""Verify the all-degree minimal-excess tail-descent identities for HC4."""
from __future__ import annotations

import json
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT_DIR = ROOT / "artifacts" / "generated-results"
ARTIFACT_DIR.mkdir(parents=True, exist_ok=True)

# Abstract binary jets.
A,B,R = sp.symbols("A B R")
Ax,Ay,Axx,Axy,Ayy = sp.symbols("Ax Ay Axx Axy Ayy")
Bx,By,Bxx,Bxy,Byy = sp.symbols("Bx By Bxx Bxy Byy")
Rx,Ry,Rxx,Rxy,Ryy = sp.symbols("Rx Ry Rxx Rxy Ryy")
a,s,z = sp.symbols("a s z")

# ---------------------------------------------------------------------------
# 1. Universal highest-tail identity.
# ---------------------------------------------------------------------------
j = sp.symbols("j")
g = sp.Matrix([Rx, Ry, j*R])
H = sp.Matrix([
    [Rxx,Rxy,j*Rx],
    [Rxy,Ryy,j*Ry],
    [j*Rx,j*Ry,j*(j-1)*R],
])
highest = sp.factor((g.T*H.adjugate()*g)[0])
Q = Rx**2*Ryy - 2*Rx*Ry*Rxy + Ry**2*Rxx
expected = j*R*(j*R*(Rxx*Ryy-Rxy**2) - (j+1)*Q)
assert sp.expand(highest-expected) == 0

# ---------------------------------------------------------------------------
# 2. Universal highest mixed coefficient.  All terms linear in epsilon have
#    the same z-degree 3j+i-2, so z may be set to one after differentiation.
# ---------------------------------------------------------------------------
i,r,eps,X = sp.symbols("i r eps X")
S,Sx,Sy,Sxx,Sxy,Syy = sp.symbols("S Sx Sy Sxx Sxy Syy")
gmix = sp.Matrix([Rx+eps*Sx, Ry+eps*Sy, j*R+eps*i*S])
Hmix = sp.Matrix([
    [Rxx+eps*Sxx, Rxy+eps*Sxy, j*Rx+eps*i*Sx],
    [Rxy+eps*Sxy, Ryy+eps*Syy, j*Ry+eps*i*Sy],
    [j*Rx+eps*i*Sx, j*Ry+eps*i*Sy,
     j*(j-1)*R+eps*i*(i-1)*S],
])
linear = sp.expand(sp.diff((gmix.T*Hmix.adjugate()*gmix)[0],eps).subs(eps,0))
linear_power = sp.factor(linear.subs({
    R:X**r, Rx:r*X**(r-1), Ry:0,
    Rxx:r*(r-1)*X**(r-2), Rxy:0, Ryy:0,
}))
assert sp.factor(linear_power + j*r*(r+j)*X**(3*r-2)*Syy) == 0

# ---------------------------------------------------------------------------
# 3. Complete h=0 face with a z^3 tail R/6.
# ---------------------------------------------------------------------------
P = A**2 + a*z*A + s*a**2*z**2
cx = Bx*P + B*(2*A+a*z)*Ax + z**3*Rx/6
cy = By*P + B*(2*A+a*z)*Ay + z**3*Ry/6
cz = B*(a*A+2*s*a**2*z) + z**2*R/2
cxx = Bxx*P + 2*Bx*(2*A+a*z)*Ax + B*(2*Ax**2+(2*A+a*z)*Axx) + z**3*Rxx/6
cxy = Bxy*P + Bx*(2*A+a*z)*Ay + By*(2*A+a*z)*Ax + B*(2*Ax*Ay+(2*A+a*z)*Axy) + z**3*Rxy/6
cyy = Byy*P + 2*By*(2*A+a*z)*Ay + B*(2*Ay**2+(2*A+a*z)*Ayy) + z**3*Ryy/6
cxz = Bx*(a*A+2*s*a**2*z) + B*a*Ax + z**2*Rx/2
cyz = By*(a*A+2*s*a**2*z) + B*a*Ay + z**2*Ry/2
czz = 2*s*a**2*B + z*R
g3 = sp.Matrix([cx,cy,cz])
H3 = sp.Matrix([[cxx,cxy,cxz],[cxy,cyy,cyz],[cxz,cyz,czz]])
J3 = sp.Poly(sp.expand((g3.T*H3.adjugate()*g3)[0]),z)

# Scalar z^3 tail: R=t, all derivatives zero.
t = sp.symbols("t")
scalar_tail = {R:t,Rx:0,Ry:0,Rxx:0,Rxy:0,Ryy:0}
z8_scalar = sp.factor(J3.coeff_monomial(z**8).subs(scalar_tail))
assert z8_scalar == 9*a**4*s**2*t**2*(Bxx*Byy-Bxy**2)
z6_scalar_zero = sp.factor(J3.coeff_monomial(z**6).subs(scalar_tail).subs(s,0))
Pxx=Axx*B+2*Ax*Bx+A*Bxx
Pxy=Axy*B+Ax*By+Ay*Bx+A*Bxy
Pyy=Ayy*B+2*Ay*By+A*Byy
assert sp.factor(z6_scalar_zero - 9*a**2*t**2*(Pxx*Pyy-Pxy**2)) == 0

# Nonconstant z^3 tail normalized to X^r.
power_tail = {
    R:X**r,Rx:r*X**(r-1),Ry:0,
    Rxx:r*(r-1)*X**(r-2),Rxy:0,Ryy:0,
}
z9 = sp.factor(J3.coeff_monomial(z**9).subs(power_tail))
assert sp.factor(z9 + a**2*s*r**2*(7-3*r)*X**(3*r-2)*Byy/72) == 0
z8_zero = sp.factor(J3.coeff_monomial(z**8).subs(s,0).subs(power_tail))
assert sp.factor(z8_zero + a*r**2*(7-3*r)*X**(3*r-2)*Pyy/72) == 0

# ---------------------------------------------------------------------------
# 4. Boundary d=2o: add scalar z^4 tail t4/24.
# ---------------------------------------------------------------------------
t4 = sp.symbols("t4")
g4 = sp.Matrix([cx,cy,cz+t4*z**3/6])
H4 = sp.Matrix([
    [cxx,cxy,cxz],
    [cxy,cyy,cyz],
    [cxz,cyz,czz+t4*z**2/2],
])
J4 = sp.Poly(sp.expand((g4.T*H4.adjugate()*g4)[0]),z)
assert sp.factor(J4.coeff_monomial(z**12) - t4**2*(Rxx*Ryy-Rxy**2)/1296) == 0
z11 = sp.factor(J4.coeff_monomial(z**11).subs(power_tail))
assert sp.factor(z11 - Byy*X**(r-2)*a**2*r*s*t4**2*(r-1)/216) == 0
z10_zero = sp.factor(J4.coeff_monomial(z**10).subs(s,0).subs(power_tail))
assert sp.factor(z10_zero - X**(r-2)*a*r*t4**2*(r-1)*Pyy/216) == 0

result = {
    "scope": "all-degree h=0 scalar tail descent",
    "highest_tail": "j R (j R det(H_R) - (j+1) grad(R)^T adj(H_R) grad(R))",
    "mixed_tail_after_R=x^r": "-j r(r+j) x^(3r-2) S_yy",
    "closed_region": "h=0 and d>=2o",
}
(ARTIFACT_DIR / "hc4_minimal_excess_tail_descent.json").write_text(
    json.dumps(result,indent=2)+"\n", encoding="utf-8"
)
print(json.dumps(result,indent=2))
