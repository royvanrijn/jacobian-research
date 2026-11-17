#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "artifacts" / "generated-results" / "hc4_final_rank_three_smooth_chart.json"
OUT.parent.mkdir(parents=True, exist_ok=True)

# Canonical generic point of the Gauss-rank-two graph chart.
a,b,c,delta,z,s = sp.symbols("a b c delta z s")
A11,A12,A13,A22,A23 = sp.symbols("A11 A12 A13 A22 A23")

A = sp.Matrix([
    [A11,A12,A13],
    [A12,A22,A23],
    [A13,A23,0],
])
B = sp.Matrix([
    [0,1,0],
    [1,0,0],
    [0,0,0],
])
U = sp.Matrix([
    [a,b,0],
    [c,a,0],
    [0,0,0],
])

# Hessian-integrability of the moving kernel: B*U must be symmetric.
assert B*U == (B*U).T

v11,v12,v13,v21,v22,v23,v31,v32,v33 = sp.symbols(
    "v11 v12 v13 v21 v22 v23 v31 v32 v33"
)
V0 = sp.Matrix([
    [v11,v12,v13],
    [v21,v22,v23],
    [v31,v32,v33],
])

# Symmetry of M*(V0+zU), with M=A-zB.
def antisym_entries(M):
    return [sp.expand(M[i,j]-M[j,i]) for i in range(3) for j in range(i+1,3)]

sym0 = antisym_entries(A*V0)
sym1 = antisym_entries(A*U-B*V0)
assert antisym_entries(B*U) == [0,0,0]

# The second set immediately gives the two entries that enter the bordered
# determinant coefficients.
sol_simple = sp.solve(sym1, [v11, v13, v23], dict=True, simplify=False)
assert sol_simple
simple = sol_simple[0]
assert sp.factor(simple[v13] + A13*b + A23*a) == 0
assert sp.factor(simple[v23] + A13*a + A23*c) == 0

M = A-z*B
V = V0+z*U
e3 = sp.Matrix([0,0,1])
identity = sp.expand(
    -(e3.T*M*(V+s*sp.eye(3)).adjugate()*e3)[0]
    - delta*M.det()
)
poly = sp.Poly(identity, s, z)

# Substitute the simple symmetry relations.
identity_sub = sp.expand(identity.subs(simple))
poly_sub = sp.Poly(identity_sub, s, z)
coeff_s = sp.factor(poly_sub.coeff_monomial(s))
coeff_z = sp.factor(poly_sub.coeff_monomial(z))

expected_s = -(A13**2*b + 2*A13*A23*a + A23**2*c)
assert sp.factor(coeff_s-expected_s) == 0

expected_z = -2*A13*A23*(a**2-b*c-delta)
assert sp.factor(coeff_z-expected_z) == 0

# Generic case A13*A23 != 0: coefficient z gives det(projective derivative)=delta.
assert sp.factor(expected_z / (-2*A13*A23) - (a**2-b*c-delta)) == 0

# Special case A23=0.  Nonsingularity then implies A13*A22 != 0.
# coeff_s forces b=0.  Use all remaining symmetry equations and verify that
# the constant term becomes -A13^2*A22*(a^2-delta).
def special_constant(substitutions, solve_vars):
    equations = [sp.expand(eq.subs(simple).subs(substitutions)) for eq in sym0]
    solution = sp.solve(equations, solve_vars, dict=True, simplify=False)
    assert solution
    expr = sp.factor(
        poly_sub.coeff_monomial(1)
        .subs(substitutions)
        .subs(solution[0])
    )
    return expr

case_beta0 = special_constant({A23:0,b:0}, [v12,v21,v31])
assert sp.factor(case_beta0 + A13**2*A22*(a**2-delta)) == 0

case_alpha0 = special_constant({A13:0,c:0}, [v12,v21,v31])
assert sp.factor(case_alpha0 + A11*A23**2*(a**2-delta)) == 0

result = {
    "scope": "final rank-three [4] smooth Gauss-rank-two chart",
    "status": "canonical local motion identities verified",
    "identities": {
        "projective_kernel_jet": "U=[[a,b,0],[c,a,0],[0,0,0]]",
        "s_coefficient": "A13^2 b + 2 A13 A23 a + A23^2 c = 0",
        "z_coefficient": "2 A13 A23 (delta-(a^2-bc)) = 0",
        "special_A23_zero": "b=0 and A13^2 A22 (delta-a^2)=0",
        "special_A13_zero": "c=0 and A11 A23^2 (delta-a^2)=0",
        "conclusion": "det(projective kernel differential)=a^2-bc=delta != 0"
    }
}
OUT.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
print(json.dumps(result, indent=2))
