#!/usr/bin/env python3
"""Exact first-order verification for HC4RSD77.

Checks, in an S-adapted regular-nilpotent frame, that Codazzi symmetry of
S and T=SN plus Frobenius of ker N^2, quasi-translation normalization and
unit affine volume force ker N^2 to be autoparallel.  Also verifies that the
three previously visible kernel-motion coefficients contain only two
independent modes.
"""
from __future__ import annotations

import json
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "artifacts" / "generated-results" / "hc4_affine_plane_bridge.json"
OUT.parent.mkdir(parents=True, exist_ok=True)

n = 4
Gamma = {}
variables = []
for i in range(n):
    for j in range(n):
        for k in range(n):
            s = sp.symbols(f"g{i+1}{j+1}{k+1}")
            Gamma[i,j,k] = s
            variables.append(s)

S = sp.zeros(n)
for i in range(n):
    S[i,n-1-i] = 1
N = sp.zeros(n)
for j in range(1,n):
    N[j-1,j] = 1
T = S*N


def cov(M, i, j, k):
    return -sum(
        Gamma[i,j,a]*M[a,k] + Gamma[i,k,a]*M[j,a]
        for a in range(n)
    )


eqs = []
for M in (S,T):
    for i in range(n):
        for j in range(n):
            for k in range(n):
                c = cov(M,i,j,k)
                eqs.append(c-cov(M,j,i,k))
                eqs.append(c-cov(M,i,k,j))

# Frobenius ker N^2=<e1,e2>, ker N^3=<e1,e2,e3>.
for k in (2,3):
    eqs.append(Gamma[0,1,k]-Gamma[1,0,k])
for i in range(3):
    for j in range(i+1,3):
        eqs.append(Gamma[i,j,3]-Gamma[j,i,3])

# Primitive kernel and unit affine volume.
for k in range(n):
    eqs.append(Gamma[0,0,k])
for i in range(n):
    eqs.append(sum(Gamma[i,j,j] for j in range(n)))

eqs = [sp.expand(e) for e in eqs if sp.expand(e) != 0]
A,_ = sp.linear_eq_to_matrix(eqs, variables)
null = A.nullspace()


def forced_zero(symbol):
    idx = variables.index(symbol)
    return all(v[idx] == 0 for v in null)

# Autoparallelity: nabla_{ei} ej has no e3/e4 component for i,j in E2.
for i in (0,1):
    for j in (0,1):
        for k in (2,3):
            assert forced_zero(Gamma[i,j,k]), (i,j,k)

# Previously visible kernel motions.
a = Gamma[2,0,1]  # Gamma^2_{3,1}
b = Gamma[3,0,1]  # Gamma^2_{4,1}
c = Gamma[3,0,2]  # Gamma^3_{4,1}

# Verify c-a is forced zero, while a and b are individually allowed.
idx = {v:i for i,v in enumerate(variables)}
assert all((v[idx[c]]-v[idx[a]]) == 0 for v in null)
assert not forced_zero(a)
assert not forced_zero(b)

result = {
    "scope": "HC4RSD77 affine-plane middle foliation",
    "status": "verified",
    "autoparallel_E2": True,
    "independent_kernel_twists": 2,
    "twist_relation": "Gamma^3_{4,1}=Gamma^2_{3,1}",
}
OUT.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
print(json.dumps(result, indent=2))
