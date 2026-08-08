#!/usr/bin/env python3
"""Exact algebraic checks for HC4RSD71.

This checker verifies the finite matrix identities used in the
 developable-gradient-image reduction.  The classical focal/developable
classification is external geometry and is not represented as a CAS claim.
"""
from __future__ import annotations

import json
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "artifacts" / "generated-results" / "hc4_final_rank_three_developable_image.json"
OUT.parent.mkdir(parents=True, exist_ok=True)

s, t = sp.symbols("s t")

# ---------------------------------------------------------------------------
# 1. Quotient determinant / Krylov identity.
# ---------------------------------------------------------------------------
m = sp.symbols("m0:9")
M = sp.Matrix(3, 3, m)
r = sp.Matrix(sp.symbols("r0:3"))
c = sp.Matrix(sp.symbols("c0:3"))

poly = sp.expand((r.T * (sp.eye(3) + s * M).adjugate() * c)[0])
assert sp.Poly(poly, s).degree() <= 2

c0 = sp.expand(poly).coeff(s, 0)
c1 = sp.expand(poly).coeff(s, 1)
c2 = sp.expand(poly).coeff(s, 2)

assert sp.expand(c0 - (r.T * c)[0]) == 0
assert sp.expand(c1 - (sp.trace(M) * (r.T * c)[0] - (r.T * M * c)[0])) == 0
assert sp.expand(c2 - (r.T * M.adjugate() * c)[0]) == 0

# Cayley-Hamilton form of adj(M) for 3x3 matrices.  Hence after r.c=r.M.c=0,
# the s^2 coefficient is exactly r.M^2.c.
sigma2 = (
    sp.trace(M) ** 2 - sp.trace(M * M)
) / 2
assert sp.simplify(
    M.adjugate() - (M * M - sp.trace(M) * M + sigma2 * sp.eye(3))
) == sp.zeros(3)

# ---------------------------------------------------------------------------
# 2. Cubic orbit formula in the tangent hyperplane.
# ---------------------------------------------------------------------------
# Use a generic strict upper-triangular nilpotent operator on a 3-dimensional
# invariant hyperplane.  This is the universal Jordan-coordinate model for the
# restriction K^T|_{k^perp}; polynomial identities are basis invariant.
a, b, d = sp.symbols("a b d")
L = sp.Matrix([[0, a, b], [0, 0, d], [0, 0, 0]])
p = sp.Matrix(sp.symbols("p0:3"))
I3 = sp.eye(3)
assert L ** 3 == sp.zeros(3)

inverse_on_p = (I3 + t * L).inv() * p
expected_inverse = p - t * L * p + t**2 * L**2 * p
assert sp.simplify(inverse_on_p - expected_inverse) == sp.zeros(3, 1)

H0 = sp.Matrix(sp.symbols("H0:3"))
H_orbit = H0 + t * inverse_on_p
expected_H = H0 + t * p - t**2 * L * p + t**3 * L**2 * p
assert sp.simplify(H_orbit - expected_H) == sp.zeros(3, 1)

# ---------------------------------------------------------------------------
# 3. Riccati solution identity (noncommutative content checked symbolically
#    with exact numerical matrices; the proof in the note is by differentiation).
# ---------------------------------------------------------------------------
T0 = sp.Matrix([[0, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0]])
B = sp.Matrix([[1, 2, 0, 1], [2, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0]])
R = T0 * (sp.eye(4) + t * B * T0).inv()
assert sp.simplify(sp.diff(R, t) + R * B * R) == sp.zeros(4)
assert sp.simplify(R.subs(t, 0) - T0) == sp.zeros(4)

result = {
    "scope": "final rank-three [4] developable-gradient-image reduction",
    "status": "exact orbit and quotient identities verified",
    "identities": {
        "quotient": "r^T adj(I+sM)c has coefficients r.c, tr(M)r.c-r.M.c, r.adj(M).c",
        "krylov": "r.c=r.M.c=0 and top coefficient nonzero imply r.M^2.c nonzero",
        "orbit": "H(phi_t)=H+t p-t^2 K^T p+t^3 (K^T)^2 p on k^perp",
        "riccati": "T'=-T B T"
    },
    "external_inputs": [
        "classification of projective varieties with Gauss rank two",
        "focal scheme of a degenerate Gauss map lies in the singular locus"
    ]
}
OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
print(json.dumps(result, indent=2))
