#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "artifacts" / "generated-results" / "hc4_final_frobenius_closure.json"
OUTPUT.parent.mkdir(parents=True, exist_ok=True)

# Canonical S-self-adjoint nilpotent Jordan pair.
N = sp.Matrix([
    [0, 1, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1],
    [0, 0, 0, 0],
])
S = sp.Matrix([
    [0, 0, 0, 1],
    [0, 0, 1, 0],
    [0, 1, 0, 0],
    [1, 0, 0, 0],
])
assert N.T * S == S * N
T = S * N
assert T == sp.Matrix([
    [0, 0, 0, 0],
    [0, 0, 0, 1],
    [0, 0, 1, 0],
    [0, 1, 0, 0],
])

# B = Hess(g)(F), symmetric.
b11,b12,b13,b14,b22,b23,b24,b33,b34,b44 = sp.symbols(
    "b11 b12 b13 b14 b22 b23 b24 b33 b34 b44"
)
B = sp.Matrix([
    [b11,b12,b13,b14],
    [b12,b22,b23,b24],
    [b13,b23,b33,b34],
    [b14,b24,b34,b44],
])
J = B * T

# Exact trace identity.
trace_identity = sp.factor(sp.trace(J))
assert trace_identity == 2*b24 + b33

# lambda = S e1.  Hessian symmetry gives d lambda from the skew part of S J.
SJ = S * J
skew = SJ - SJ.T
# Restrict to E3 = span(e1,e2,e3): entries (12),(13),(23).
frob_entries = [
    sp.factor(skew[0,1]),
    sp.factor(skew[0,2]),
    sp.factor(skew[1,2]),
]
assert frob_entries == [b44, b34, -b24 + b33]

# The Gauss line ell = S e1 = e4* is the radical of II_Y, hence
# b24=b34=b44=0.  On that locus trace(J)=b33 and the sole remaining
# Frobenius coefficient is also b33.
radical_subs = {b24:0, b34:0, b44:0}
assert sp.factor(trace_identity.subs(radical_subs)) == b33
assert [sp.factor(v.subs(radical_subs)) for v in frob_entries] == [0,0,b33]

result = {
    "scope": "final rank-three [4] HC4 Frobenius gap",
    "status": "closed",
    "canonical_pair": "N single nilpotent Jordan block, S anti-diagonal self-adjoint metric",
    "trace_identity": "tr(Hess(g)(F) * T) = 2 b24 + b33",
    "gauss_radical": ["b24=0", "b34=0", "b44=0"],
    "remaining_scalar": "b33 = II_Y(m,m) = Frobenius obstruction for ker N^3",
    "nilpotence_consequence": "tr Jn = 0 => b33 = 0",
}
OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
print(json.dumps(result, indent=2))
