#!/usr/bin/env python3
"""Exact algebraic checks for HC4RSD66--67."""
from __future__ import annotations

import json
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "artifacts" / "generated-results" / "hc4_rank_three_null_fiber.json"
OUT.parent.mkdir(parents=True, exist_ok=True)

s = sp.symbols("s")
# Use the canonical self-adjoint Jordan pair.  This is enough to verify the
# polynomial inverse/adjugate coefficient formulas and rank pattern.
a,b,c,d = sp.symbols("a b c d", nonzero=True)
N = sp.Matrix([[0,1,0,0],[0,0,1,0],[0,0,0,1],[0,0,0,0]])
S = sp.Matrix([[0,0,0,a],[0,0,a,b],[0,a,b,c],[a,b,c,d]])
assert N.T * S == S * N
T = S * N
assert T == T.T
assert N**4 == sp.zeros(4)
delta = sp.factor(S.det())
assert delta == a**4
M = S + s*T
assert sp.factor(M.det() - delta) == 0
G = sp.simplify(M.inv())
G_expected = (sp.eye(4)-s*N+s**2*N**2-s**3*N**3)*S.inv()
assert sp.simplify(G-G_expected) == sp.zeros(4)
Adj = sp.simplify(M.adjugate())
assert sp.simplify(Adj-delta*G_expected) == sp.zeros(4)
coeffs=[]
for j in range(4):
    Cj=Adj.applyfunc(lambda entry: sp.Poly(sp.expand(entry),s).coeff_monomial(s**j))
    expected=sp.simplify(delta*((-1)**j)*(N**j)*S.inv())
    assert sp.simplify(Cj-expected)==sp.zeros(4)
    coeffs.append(Cj)
assert coeffs[3] == T.adjugate()
assert coeffs[3].rank() == 1
assert coeffs[2].rank() == 2
assert coeffs[1].rank() == 3

# Generic determinant identity: for a 4x4 matrix pencil the s^3 coefficient
# is tr(adj(T) S).  Verify on a dense symbolic symmetric rank-three chart
# T=diag(h1,h2,h3,0), with a general symmetric S.
h1,h2,h3=sp.symbols("h1 h2 h3", nonzero=True)
t00,t01,t02,t03,t11,t12,t13,t22,t23,t33=sp.symbols("t00 t01 t02 t03 t11 t12 t13 t22 t23 t33")
Sg=sp.Matrix([[t00,t01,t02,t03],[t01,t11,t12,t13],[t02,t12,t22,t23],[t03,t13,t23,t33]])
Tg=sp.diag(h1,h2,h3,0)
poly=sp.Poly(sp.expand((Sg+s*Tg).det()),s)
coef3=sp.factor(poly.coeff_monomial(s**3))
trace_expr=sp.factor(sp.trace(Tg.adjugate()*Sg))
assert sp.factor(coef3-trace_expr)==0
assert coef3 == h1*h2*h3*t33

result={
  "scope":"final rank-three [4] null-fiber/cofactor reduction",
  "status":"exact algebraic identities verified",
  "identities":{
    "inverse_pencil":"(S+sT)^-1=(I-sN+s^2N^2-s^3N^3)S^-1",
    "cofactor_flag":"C_j=delta(-1)^j N^j S^-1",
    "top_cofactor":"C_3=adj(T), rank 1",
    "top_determinant":"[s^3] det(S+sT)=tr(adj(T)S)"
  }
}
OUT.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
print(json.dumps(result,indent=2))
