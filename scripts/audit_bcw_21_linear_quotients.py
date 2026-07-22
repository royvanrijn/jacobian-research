#!/usr/bin/env python3
"""Classify invariant row modules of the essential 21D BCW map."""

import json
from pathlib import Path

import sympy as sp

from audit_bcw_22_linear_quotients import (
    decode_h,
    generated_algebra_dimension,
    jacobian_coefficient_matrices,
    modular_flat,
)


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "artifacts" / "generated-results" / "essential_bcw_21_counterexample.json"
stored = json.loads(SOURCE.read_text())
variables, h = decode_h(stored)
assert len(variables) == 21
matrices = jacobian_coefficient_matrices(h, variables)
assert len(matrices) == 65

# The homogenizing covector L is invariant and constant on the collision.
last_row = sp.zeros(1, 21)
last_row[0, 20] = 1
assert all(last_row * matrix == sp.zeros(1, 21) for matrix in matrices)
points = [
    sp.Matrix([sp.Rational(value) for value in point])
    for point in stored["collision_points"]
]
assert len({point[20] for point in points}) == 1

# The generated coefficient algebra on V/L is the full M_20.  Full rank at
# one good prime certifies full rank over QQ, so V/L is a simple module.
quotient_generators = [modular_flat(matrix[:20, :20]) for matrix in matrices]
assert generated_algebra_dimension(quotient_generators, 20) == 20 * 20

# A proper invariant subspace surjecting onto V/L would be an invariant
# hyperplane.  Its annihilator would be a common invariant column line.  The
# 64 nilpotent coefficient matrices act by zero on such a line, but their
# common kernel is zero, excluding it.
nilpotent = [matrix for matrix in matrices if matrix**21 == sp.zeros(21)]
exceptional = [matrix for matrix in matrices if matrix**21 != sp.zeros(21)]
assert len(nilpotent) == 64 and len(exceptional) == 1
assert sp.Matrix.vstack(*nilpotent).nullspace() == []
lam = sp.Symbol("lambda")
assert sp.factor(exceptional[0].charpoly().as_expr()) == lam**19 * (lam - 6) * (lam + 6)

print("PASS 21D module audit: 65 exact Jacobian coefficient matrices")
print("PASS 21D module audit: quotient algebra mod <e_20^T> is full M_20(Q)")
print("PASS 21D module audit: no common invariant column line")
print("PASS 21D module audit: only proper row module is constant on the collision")
print("PASS 21D module audit: no further collision-preserving linear quotient")
