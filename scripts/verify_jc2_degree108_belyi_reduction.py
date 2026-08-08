#!/usr/bin/env python3
from __future__ import annotations

import itertools
import json
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "artifacts" / "generated-results" / "jc2_degree108_belyi_reduction.json"
OUT.parent.mkdir(parents=True, exist_ok=True)

X, s = sp.symbols("X s")
A, B, C, D, E, F, G = [sp.Function(name)(X) for name in "ABCDEFG"]
P = A*s**2 + B*s + C
Q = D*s**3 + E*s**2 + F*s + G
J = sp.expand(sp.diff(P, X)*sp.diff(Q, s) - sp.diff(P, s)*sp.diff(Q, X))
coeffs = [sp.factor(J.coeff(s, k)) for k in range(5)]
expected = [
    -B*sp.diff(G,X) + F*sp.diff(C,X),
    -2*A*sp.diff(G,X) - B*sp.diff(F,X) + 2*E*sp.diff(C,X) + F*sp.diff(B,X),
    -2*A*sp.diff(F,X) - B*sp.diff(E,X) + 3*D*sp.diff(C,X) + 2*E*sp.diff(B,X) + F*sp.diff(A,X),
    -2*A*sp.diff(E,X) - B*sp.diff(D,X) + 3*D*sp.diff(B,X) + 2*E*sp.diff(A,X),
    -2*A*sp.diff(D,X) + 3*D*sp.diff(A,X),
]
assert all(sp.expand(a-b) == 0 for a,b in zip(coeffs, expected))

# The Laurent coordinate change X=x*y^2, s=1/y has determinant -1.
x, y = sp.symbols("x y")
coord_jac = sp.det(sp.Matrix([
    [sp.diff(x*y**2, x), sp.diff(x*y**2, y)],
    [sp.diff(1/y, x), sp.diff(1/y, y)],
]))
assert sp.simplify(coord_jac + 1) == 0

# Exact dessin enumeration for passport (2^10,1), (3^7), (17,1^4).
def compose(p, q):
    return tuple(p[q[i]] for i in range(len(p)))

def cycle_type(p):
    seen = [False]*len(p)
    out = []
    for i in range(len(p)):
        if not seen[i]:
            j=i; length=0
            while not seen[j]:
                seen[j]=True; length+=1; j=p[j]
            out.append(length)
    return sorted(out, reverse=True)

n = 21
sigma_inf = list(range(n))
for i in range(17):
    sigma_inf[i] = (i+1) % 17
sigma_inf = tuple(sigma_inf)
fixed_inf = list(range(17,21))
solutions = []
for centers in itertools.combinations(range(17), 4):
    forced = {}
    good = True
    for f, r in zip(fixed_inf, centers):
        for u, v in [(f, r), ((r+1)%17, (r-1)%17)]:
            if u == v or u in forced or v in forced:
                good = False
                break
            forced[u] = v
            forced[v] = u
        if not good:
            break
    if not good:
        continue
    rem = [i for i in range(17) if i not in forced]
    if len(rem) != 5:
        continue
    for fixed0 in rem:
        rest = [i for i in rem if i != fixed0]
        a,b,c,d = rest
        matchings = [
            [(a,b),(c,d)],
            [(a,c),(b,d)],
            [(a,d),(b,c)],
        ]
        for pairs in matchings:
            sigma0 = list(range(n))
            for u,v in forced.items(): sigma0[u]=v
            sigma0[fixed0]=fixed0
            for u,v in pairs:
                sigma0[u]=v; sigma0[v]=u
            sigma0=tuple(sigma0)
            sigma1=compose(sigma0, sigma_inf)
            if cycle_type(sigma0)==[2]*10+[1] and cycle_type(sigma1)==[3]*7:
                solutions.append((centers, fixed0, pairs))

assert len(solutions) == 85
center_sets = {tuple(sol[0]) for sol in solutions}
assert len(center_sets) == 85

def canonical_rotation(cs):
    return min(tuple(sorted((x+r)%17 for x in cs)) for r in range(17))

orbits = sorted({canonical_rotation(cs) for cs in center_sets})
expected_orbits = sorted([
    (0,3,7,11),
    (0,3,7,12),
    (0,3,8,11),
    (0,3,8,13),
    (0,3,9,13),
])
assert orbits == expected_orbits

result = {
    "scope": "open (72,108) / (8,28) no-vertical-edge JC2 residue",
    "status": "reduced to finite Belyi deformation problem",
    "wronskian_equations": [str(v) for v in coeffs[::-1]],
    "belyi_passport": ["2^10,1", "3^7", "17,1,1,1,1"],
    "labelled_center_sets": len(center_sets),
    "rotation_orbits": [list(v) for v in orbits],
}
OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
print(json.dumps(result, indent=2))
