#!/usr/bin/env python3
"""Verify the exact algebraic bridge in the global smooth-developable HC4 obstruction.

The projective developable-surface classification and classical Bertini
irreducibility theorem are external geometric inputs recorded in the companion
note.  This checker verifies every polynomial-matrix identity used to pass from
a constant-Hessian four-variable scalar packet to the three-variable
bordered-Hessian equation.
"""
from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "hc4_global_smooth_developable_obstruction.json"
)
OUTPUT.parent.mkdir(parents=True, exist_ok=True)

w = sp.symbols("w")
e00, e01, e02, e11, e12, e22 = sp.symbols(
    "e00 e01 e02 e11 e12 e22"
)
g00, g01, g02, g11, g12, g22 = sp.symbols(
    "g00 g01 g02 g11 g12 g22"
)
p0, p1, p2 = sp.symbols("p0 p1 p2")

E = sp.Matrix(
    [
        [e00, e01, e02],
        [e01, e11, e12],
        [e02, e12, e22],
    ]
)
G = sp.Matrix(
    [
        [g00, g01, g02],
        [g01, g11, g12],
        [g02, g12, g22],
    ]
)
p = sp.Matrix([p0, p1, p2])
M = E + w * G
block = sp.Matrix.vstack(
    sp.Matrix.hstack(M, p),
    sp.Matrix([[p0, p1, p2, 0]]),
)

block_det = sp.expand(block.det(method="domain-ge"))
schur_expression = sp.expand(-(p.T * M.adjugate() * p)[0])
assert sp.expand(block_det - schur_expression) == 0

block_polynomial = sp.Poly(block_det, w)
assert block_polynomial.degree() <= 2
assert sp.factor(
    block_polynomial.coeff_monomial(w**2)
    + (p.T * G.adjugate() * p)[0]
) == 0

# If the gradient column vanishes at a point, the last row and column of the
# four-variable Hessian vanish there, so its determinant is zero.
assert sp.expand(block_det.subs({p0: 0, p1: 0, p2: 0})) == 0

# The bordered-Hessian determinant for c is minus the UFE polynomial
# grad(c)^T adj(Hess(c)) grad(c).
bordered_c = sp.Matrix.vstack(
    sp.Matrix.hstack(G, p),
    sp.Matrix([[p0, p1, p2, 0]]),
)
ufe = sp.expand((p.T * G.adjugate() * p)[0])
assert sp.expand(bordered_c.det(method="domain-ge") + ufe) == 0

result = {
    "scope": "global smooth-developable obstruction for scalar HC4 reverse Schur",
    "status": "exact algebraic bridge verified",
    "identities": {
        "block_determinant": "det[[E+wG,p],[p^T,0]]=-p^T adj(E+wG)p",
        "quadratic_coefficient": "[w^2]det=-p^T adj(G)p",
        "gradient_zero": "p=0 forces determinant zero",
        "bordered_hessian": "det[[G,p],[p^T,0]]=-p^T adj(G)p",
    },
    "external_inputs": [
        "classical Bertini irreducibility for non-composite polynomials",
        "classification of irreducible projective developable surfaces in P3",
        "singularity of a non-planar tangent developable along its edge of regression",
    ],
}
OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
print(json.dumps(result, indent=2))
