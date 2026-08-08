#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "artifacts" / "generated-results" / "hc4_final_rank_three_krylov_flag.json"
OUT.parent.mkdir(parents=True, exist_ok=True)

t = sp.symbols("t")
a0,a1,a2 = sp.symbols("a0 a1 a2")
c11,c12,c13,c21,c22,c23,c31,c32,c33 = sp.symbols(
    "c11 c12 c13 c21 c22 c23 c31 c32 c33"
)

M0 = sp.Matrix([
    [0,0,a0],
    [1,0,a1],
    [0,1,a2],
])
C = sp.Matrix([
    [c11,c12,c13],
    [c21,c22,c23],
    [c31,c32,c33],
])
c = sp.Matrix([1,0,0])
r = sp.Matrix([[0,0,1]])
M = M0+t*C

# The s coefficient in the quotient determinant is r^T M c.
linear = sp.expand((r*M*c)[0])
assert linear == c31*t

# The s^2 coefficient is r^T adj(M)c.
quadratic = sp.expand((r*M.adjugate()*c)[0])
expected = 1 + t*(c21+c32) + t**2*(c21*c32-c22*c31)
assert sp.expand(quadratic-expected) == 0

# Once c31=0, t-independence gives c21+c32=0 and c21*c32=0, hence both zero.
u = sp.symbols("u")
assert sp.solve([u + (-u), u*(-u)], [u], dict=True) == [{u: 0}]

result = {
    "scope": "fiberwise Krylov flag in final rank-three [4] HC4 stratum",
    "status": "exact quotient identities verified",
    "identities": {
        "s_coefficient": "r^T(M0+tC)c = t*C31",
        "s2_coefficient": "r^T adj(M0+tC)c = 1+t(C21+C32)+t^2(C21*C32-C22*C31)",
        "conclusion": "C21=C31=C32=0"
    }
}
OUT.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
print(json.dumps(result, indent=2))
