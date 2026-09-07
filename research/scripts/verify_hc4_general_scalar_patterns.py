#!/usr/bin/env python3
"""Verify all-degree identities behind the scalar reverse-Schur compression.

The script deliberately proves symbolic identities in abstract degree/local
parameters rather than enumerating root partitions.
"""
from __future__ import annotations

import json
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT_DIR = ROOT / "artifacts" / "generated-results"
ARTIFACT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# 1. Resonance square identity
# ---------------------------------------------------------------------------
d, m, e, n = sp.symbols("d m e n")
C = (
    d**2 * m + d**2 * n**2 - 2*d*e*m*n - 2*d*e*m
    - d*m**2 - d*m + e**2*m**2 + 2*e*m**2 + m**2
)
square_form = (d*n - e*m)**2 + m*(d-m)*(d-2*e-1)
assert sp.expand(C - square_form) == 0
assert sp.factor(sp.discriminant(C, n)) == -4*d**2*m*(d-m)*(d-2*e-1)

kappa = -C / (d*m*(d-m))
kappa_square = (2*e + 1 - d)/d - (d*n-e*m)**2/(d*m*(d-m))
assert sp.factor(kappa-kappa_square) == 0

# Parity-separated root-value coincidence formulas.
m1, m2 = sp.symbols("m1 m2")
k_even = lambda mm: sp.factor((2*e+1-d)/d - mm*(d-2*e)**2/(4*d*(d-mm)))
k_odd = lambda mm: sp.factor((2*e+1-d)/d - (d+mm*(d-2*e))**2/(4*d*mm*(d-mm)))
assert sp.factor(
    k_even(m1)-k_even(m2)
    + (d-2*e)**2*(m1-m2)/(4*(d-m1)*(d-m2))
) == 0
odd_coincidence = d*(m1+m2-d) + m1*m2*(d-2*e)*(d-2*e+2)
assert sp.factor(
    k_odd(m1)-k_odd(m2)
    + (m1-m2)*odd_coincidence/(4*m1*m2*(d-m1)*(d-m2))
) == 0

# ---------------------------------------------------------------------------
# 2. Transverse-excess bookkeeping
# ---------------------------------------------------------------------------
o, h = sp.symbols("o h")
# W=sum ceil(m_i/2)=(d+o)/2.  With e=W+h:
W = (d+o)/2
e_from_h = W+h
assert sp.expand(e_from_h-W-h) == 0
# q has degree 2e-d; after q=B*C with deg B=o, deg C=2h.
assert sp.expand((2*e_from_h-d)-o-2*h) == 0
# j-th z-tail has binary degree j*e-(j-1)d.
j = sp.symbols("j")
tail_degree = sp.expand(j*e_from_h-(j-1)*d)
assert tail_degree == sp.expand(d + j*(o-d)/2 + j*h)

# ---------------------------------------------------------------------------
# 3. All-even h=0 complete-face obstruction
# ---------------------------------------------------------------------------
# Let A=A(x,y) abstractly.  We only need its first/second derivatives to
# extract the highest z coefficient from
#   c=A^2 + a z A + a^2 z^2/(2d).
# No degree-specific expansion is used.
A, Ax, Ay, Axx, Axy, Ayy, a, z = sp.symbols(
    "A Ax Ay Axx Axy Ayy a z"
)
grad = sp.Matrix([
    (2*A+a*z)*Ax,
    (2*A+a*z)*Ay,
    a*A + a**2*z/d,
])
hess = sp.Matrix([
    [(2*A+a*z)*Axx + 2*Ax**2,
     (2*A+a*z)*Axy + 2*Ax*Ay,
     a*Ax],
    [(2*A+a*z)*Axy + 2*Ax*Ay,
     (2*A+a*z)*Ayy + 2*Ay**2,
     a*Ay],
    [a*Ax, a*Ay, a**2/d],
])
J = sp.Poly(sp.expand((grad.T*hess.adjugate()*grad)[0]), z)
z4 = sp.factor(J.coeff_monomial(z**4))
assert z4 == a**6*(Axx*Ayy-Axy**2)/d**2

result = {
    "scope": "all-degree synchronized scalar reverse-Schur identities",
    "resonance_square_identity": "C=(dn-em)^2+m(d-m)(d-2e-1)",
    "transverse_excess": {
        "W": "(d+o)/2",
        "deg_H": "h=e-W",
        "deg_C": "2h",
    },
    "all_even_h0": {
        "family": "f=A^2, g=aA, q=a^2/d",
        "terminal_coefficient": "[z^4]J=a^6 det(Hess A)/d^2",
        "consequence": "binary Hess(A)=0, hence A is a power of a linear form",
    },
}
output = ARTIFACT_DIR / "hc4_general_scalar_patterns.json"
output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
print(json.dumps(result, indent=2))
