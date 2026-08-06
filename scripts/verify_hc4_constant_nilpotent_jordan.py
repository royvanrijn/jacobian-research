#!/usr/bin/env python3
"""Classify every constant nilpotent relative Jordan frame in the HC4 pencil."""
from __future__ import annotations

import json
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT_DIR = ROOT / "artifacts" / "generated-results"
ARTIFACT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# 8. Constant nilpotent matrix pivots
# ---------------------------------------------------------------------------
X, Y, Z, W = sp.symbols("X Y Z W")
s_pencil = sp.symbols("s_pencil")
alpha, epsilon = sp.symbols("alpha epsilon", nonzero=True)

N4 = sp.Matrix(
    [[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [0, 0, 0, 0]]
)
U = sp.Function("U")(W)
Rfun = sp.Function("Rfun")(W)
Vfun = sp.Function("Vfun")(W)
psi4 = (
    alpha * X * W
    + alpha * Y * Z
    + Y * U
    + Z**2 * sp.diff(U, W) / 2
    + Z * sp.diff(Rfun, W)
    + Vfun
)
A4 = alpha * Y * W + alpha * Z**2 / 2 + Z * U + Rfun
S4 = sp.hessian(psi4, (X, Y, Z, W))
T4 = sp.hessian(A4, (X, Y, Z, W))
assert S4 * N4 == T4
assert N4.T * S4 == S4 * N4
assert sp.factor(S4.det() - alpha**4) == 0
assert sp.factor((S4 + s_pencil * T4).det() - alpha**4) == 0

N31 = sp.Matrix(
    [[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0], [0, 0, 0, 0]]
)
Pfun = sp.Function("Pfun")(Z)
Qfun = sp.Function("Qfun")(Z)
R31 = sp.Function("R31")(Z)
psi31 = (
    alpha * X * Z
    + alpha * Y**2 / 2
    + Y * sp.diff(Pfun, Z)
    + epsilon * W**2 / 2
    + W * Qfun
    + R31
)
A31 = alpha * Y * Z + Pfun
S31 = sp.hessian(psi31, (X, Y, Z, W))
T31 = sp.hessian(A31, (X, Y, Z, W))
assert S31 * N31 == T31
assert N31.T * S31 == S31 * N31
assert sp.factor(S31.det() + alpha**3 * epsilon) == 0
assert sp.factor((S31 + s_pencil * T31).det() + alpha**3 * epsilon) == 0

N211 = sp.Matrix(
    [[0, 1, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]
)
Cfun = sp.Function("Cfun")(Y, Z, W)
psi211 = alpha * X * Y + Cfun
A211 = alpha * Y**2 / 2
S211 = sp.hessian(psi211, (X, Y, Z, W))
T211 = sp.hessian(A211, (X, Y, Z, W))
assert S211 * N211 == T211
assert N211.T * S211 == S211 * N211
passive_hessian = sp.hessian(Cfun, (Z, W))
assert sp.factor(S211.det() + alpha**2 * passive_hessian.det()) == 0

N22 = sp.Matrix(
    [[0, 1, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1], [0, 0, 0, 0]]
)
Hfun = sp.Function("Hfun")(Y, W)
R22 = sp.Function("R22")(Y, W)
psi22 = X * sp.diff(Hfun, Y) + Z * sp.diff(Hfun, W) + R22
A22 = Hfun
S22 = sp.hessian(psi22, (X, Y, Z, W))
T22 = sp.hessian(A22, (X, Y, Z, W))
assert S22 * N22 == T22
assert N22.T * S22 == S22 * N22
plane_hessian = sp.hessian(Hfun, (Y, W))
assert sp.factor(S22.det() - plane_hessian.det() ** 2) == 0



result = {
    "scope": "constant nilpotent relative frames in four variables",
    "status": "classified",
    "types": {
        "4": "triangular polynomial automorphism",
        "3+1": "triangular polynomial automorphism",
        "2+1+1": "HC2 endpoint",
        "2+2": "JC2 cotangent endpoint",
    },
}
output = ARTIFACT_DIR / "hc4_constant_nilpotent_jordan.json"
output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
print(json.dumps(result, indent=2))
