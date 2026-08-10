#!/usr/bin/env python3
"""Exact first-order verification for HC4RSD77 and audit of HC4RSD78.

Checks, in an S-adapted regular-nilpotent frame, that Codazzi symmetry of
S and T=SN plus Frobenius of ker N^2, quasi-translation normalization and
unit affine volume force ker N^2 to be autoparallel.  Also verifies that the
three previously visible kernel-motion coefficients contain only two
independent modes.

The same linear system does not force the Grassmann derivative of ker N^2
into the rank-one tangent cone.  We extract its exact normal form and exhibit
a formal first-order solution with maximal kernel-line motion for which one
Grassmann tangent matrix has rank two.  Thus the former Schubert-dichotomy
formulation needs an additional identity beyond the inputs checked here;
HC4RSD78 is retained only as the exact normal-form statement.
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


def forced_relation(expression):
    coefficients = [sp.expand(expression).coeff(v) for v in variables]
    return all(
        sum(coefficients[i] * vector[i] for i in range(len(variables))) == 0
        for vector in null
    )

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

# The derivative of E2=<e1,e2> in the transverse e3/e4 directions is a pair
# in Hom(E2,V/E2).  The exact first-order normal form is
#
#   A3 = [[0,-(a+q)/2],[0,0]],  A4 = [[a,r],[0,q]].
#
# Here ``a`` controls projective motion of the source kernel line, while q is
# an independent scalar.  The alpha/beta Schubert conclusion requires q=0,
# but the imposed equations do not force it.
q = Gamma[3,1,3]  # Gamma^4_{4,2}
r = Gamma[3,1,2]  # Gamma^3_{4,2}

assert forced_zero(Gamma[2,0,2])
assert forced_zero(Gamma[2,0,3])
assert forced_relation(Gamma[2,1,2] + (a + q) / 2)
assert forced_zero(Gamma[2,1,3])
assert forced_relation(Gamma[3,0,2] - a)
assert forced_zero(Gamma[3,0,3])
assert not forced_zero(q)
assert not forced_zero(r)

# An exact formal witness: take the two nullspace modes supporting a and q.
# It has a=q=1, r=0, so the kernel-line direction map has rank two while A4
# is the identity and hence is not in the rank-one Grassmann tangent cone.
idx = {v:i for i,v in enumerate(variables)}
a_mode = next(vector for vector in null if vector[idx[a]] == 1)
q_mode = next(
    vector
    for vector in null
    if vector[idx[q]] == 1 and vector[idx[a]] == 0
)
witness = a_mode + q_mode
assert witness[idx[a]] == 1
assert witness[idx[q]] == 1
assert witness[idx[r]] == 0
assert witness[idx[Gamma[2,0,1]]] == 1
assert witness[idx[Gamma[3,0,2]]] == 1

A3_witness = sp.Matrix([
    [witness[idx[Gamma[2,0,2]]], witness[idx[Gamma[2,1,2]]]],
    [witness[idx[Gamma[2,0,3]]], witness[idx[Gamma[2,1,3]]]],
])
A4_witness = sp.Matrix([
    [witness[idx[Gamma[3,0,2]]], witness[idx[Gamma[3,1,2]]]],
    [witness[idx[Gamma[3,0,3]]], witness[idx[Gamma[3,1,3]]]],
])
assert A3_witness == sp.Matrix([[0, -1], [0, 0]])
assert A4_witness == sp.eye(2)

result = {
    "scope": "HC4RSD77 affine-plane middle foliation",
    "status": "verified",
    "autoparallel_E2": True,
    "independent_kernel_twists": 2,
    "twist_relation": "Gamma^3_{4,1}=Gamma^2_{3,1}",
    "grassmann_derivative_normal_form": {
        "A3": "[[0,-(a+q)/2],[0,0]]",
        "A4": "[[a,r],[0,q]]",
    },
    "schubert_gate": "q=0",
    "schubert_gate_forced": False,
    "formal_non_schubert_witness": {
        "parameters": "a=1, q=1, r=0",
        "A3": "[[0,-1],[0,0]]",
        "A4": "[[1,0],[0,1]]",
        "kernel_line_motion_rank": 2,
    },
    "conclusion": (
        "HC4RSD77 and the exact first-order normal form HC4RSD78 are "
        "verified; the Schubert gate q=0 is not forced"
    ),
}
OUT.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
print(json.dumps(result, indent=2))
