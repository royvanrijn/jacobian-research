#!/usr/bin/env python3
"""Exact first-order audit of the final HC4 affine-straightening bridge.

We work in an S-adapted regular-nilpotent frame e1,...,e4 with
N e2=e1, N e3=e2, N e4=e3 and anti-diagonal S.  The frame is normalized so
S and T=SN have constant component matrices.  Unknowns Gamma[i,j,k] are the
coefficients of the ambient flat affine connection in this moving frame.

Imposed linear identities:
  * Hessian/Codazzi symmetry for S;
  * Hessian/Codazzi symmetry for T;
  * Frobenius integrability of ker N^2 and ker N^3;
  * quasi-translation normalization nabla_{e1} e1=0;
  * constant affine volume of the normalized frame (trace Gamma_i=0).

The audit verifies that these identities do NOT force the kernel line to be
parallel for the affine connection.  Exactly the upper-triangular transverse
motions Gamma^2_{31}, Gamma^2_{41}, Gamma^3_{41} can remain nonzero.
Therefore the passage from the local Frobenius theorem HC4RSD75 to a constant
affine flag requires a genuinely global/polynomial argument.
"""
from __future__ import annotations

import json
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "artifacts" / "generated-results" / "hc4_affine_bridge_first_order.json"
OUT.parent.mkdir(parents=True, exist_ok=True)

n = 4
Gamma: dict[tuple[int, int, int], sp.Symbol] = {}
variables: list[sp.Symbol] = []
for i in range(n):
    for j in range(n):
        for k in range(n):
            symbol = sp.symbols(f"g{i+1}{j+1}{k+1}")
            Gamma[i, j, k] = symbol
            variables.append(symbol)

S = sp.zeros(n)
for i in range(n):
    S[i, n - 1 - i] = 1
N = sp.zeros(n)
for j in range(1, n):
    N[j - 1, j] = 1
T = S * N


def covariant_derivative(matrix: sp.Matrix, i: int, j: int, k: int) -> sp.Expr:
    # Components of (nabla_{e_i} matrix)(e_j,e_k) when the matrix has
    # constant components in the chosen moving frame.
    return -sum(
        Gamma[i, j, a] * matrix[a, k] + Gamma[i, k, a] * matrix[j, a]
        for a in range(n)
    )


equations: list[sp.Expr] = []
for matrix in (S, T):
    for i in range(n):
        for j in range(n):
            for k in range(n):
                c = covariant_derivative(matrix, i, j, k)
                equations.append(c - covariant_derivative(matrix, j, i, k))
                equations.append(c - covariant_derivative(matrix, i, k, j))

# ker N^2 = <e1,e2> is Frobenius.
for k in (2, 3):
    equations.append(Gamma[0, 1, k] - Gamma[1, 0, k])

# ker N^3 = <e1,e2,e3> is Frobenius.
for i in range(3):
    for j in range(i + 1, 3):
        equations.append(Gamma[i, j, 3] - Gamma[j, i, 3])

# Primitive quasi-translation kernel.
for k in range(n):
    equations.append(Gamma[0, 0, k])

# The S-normalized frame has constant affine volume because det S is a unit.
for i in range(n):
    equations.append(sum(Gamma[i, j, j] for j in range(n)))

equations = [sp.expand(eq) for eq in equations if sp.expand(eq) != 0]
linear_matrix, _ = sp.linear_eq_to_matrix(equations, variables)
nullspace = linear_matrix.nullspace()


def forced_zero(symbol: sp.Symbol) -> bool:
    index = variables.index(symbol)
    return all(vector[index] == 0 for vector in nullspace)

kernel_motion = {
    f"Gamma^{k+1}_{i+1},1": forced_zero(Gamma[i, 0, k])
    for i in range(n)
    for k in range(1, n)
}
allowed = [name for name, is_zero in kernel_motion.items() if not is_zero]
expected_allowed = ["Gamma^2_3,1", "Gamma^2_4,1", "Gamma^3_4,1"]
assert allowed == expected_allowed, (allowed, expected_allowed)

result = {
    "scope": "first-order affine bridge after HC4RSD75",
    "linear_unknowns": len(variables),
    "linear_rank": int(linear_matrix.rank()),
    "solution_dimension": len(variables) - int(linear_matrix.rank()),
    "allowed_nonparallel_kernel_motion": allowed,
    "conclusion": (
        "Hessian/Codazzi + full Jordan-flag Frobenius + quasi-translation + "
        "unit Hessian volume do not by themselves force a constant affine kernel line"
    ),
}
OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
print(json.dumps(result, indent=2))
