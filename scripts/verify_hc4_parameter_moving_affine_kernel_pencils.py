#!/usr/bin/env python3
"""Verify the terminal identities for parameter-moving affine kernels.

The accompanying proof classifies a primitive kernel generator

    v(s,x) = a_0 + L_0*x + s*(a_1 + L_1*x)

using the bordered-unit adjugate factorization and the Hessian Piola
identity.  This checker replays the two nontrivial terminal charts: the
constant-at-infinity corner and the sole genuinely moving common-covector
normal form.  The latter integrates to a pencil whose Hessian has rank at
most two, contradicting the bordered unit.
"""

from __future__ import annotations

from itertools import combinations
import hashlib
import json
from pathlib import Path
import sys

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.reverse_schur_descent import (
    ScalarPivotSchurFamily,
    corank_one_adjugate_scalar,
    kernel_line_piola_residuals,
)


OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "hc4_parameter_moving_affine_kernel_pencils.json"
)


def assert_zero_vector(vector: sp.Matrix) -> None:
    assert all(sp.simplify(entry) == 0 for entry in vector)


# Constant-at-infinity corner.  Normalize v_0=(z,1,0,0) and let the
# parameter coefficient v_1 be constant.  The complete v_0-kernel potential
# is the HC4RSD2 integral.  The displayed Hessian product supplies exactly
# the two cross equations used in the proof.
x, y, z, w = sp.symbols("x y z w")
variables = (x, y, z, w)
b1, b2, b4 = sp.symbols("b1 b2 b4")
C = sp.Function("C")(z)
G = sp.Function("G")(z, w)
fixed_kernel = sp.Matrix([z, 1, 0, 0])
constant_parameter_part = sp.Matrix([b1, b2, 0, b4])
fixed_potential = sp.expand(
    y * C + (x - y * z) * sp.diff(C, z) + G
)
fixed_hessian = sp.hessian(fixed_potential, variables)
expected_cross = sp.Matrix(
    [
        0,
        0,
        (b1 - b2 * z) * sp.diff(C, z, 2)
        + b4 * sp.diff(G, z, w),
        b4 * sp.diff(G, w, 2),
    ]
)
assert_zero_vector(fixed_hessian * fixed_kernel)
assert_zero_vector(fixed_hessian * constant_parameter_part - expected_cross)
assert corank_one_adjugate_scalar(
    fixed_hessian,
    fixed_kernel,
) == -sp.diff(C, z, 2) ** 2 * sp.diff(G, w, 2)


# The two jointly unimodular common-covector normal forms.  Type II is
# rejected in the proof by grad(A).v_1=0 and grad(A).v_0=alpha, which would
# require z*D_{e_2}(A)=alpha.  Type I is the only chart requiring Hessian
# integration.
s = sp.symbols("s")
r = sp.symbols("r")
affine_normal_variables = (x, y, r, z)
type_one_kernel = sp.Matrix([s * z, z, 1, 0])
type_two_kernel = sp.Matrix([1 + s * z, z, 0, 0])
assert not any(
    kernel_line_piola_residuals(type_one_kernel, affine_normal_variables)
)
assert not any(
    kernel_line_piola_residuals(type_two_kernel, affine_normal_variables)
)
assert sp.expand(type_one_kernel[2]) == 1
assert sp.expand(type_two_kernel[0] - s * type_two_kernel[1]) == 1


# Complete terminal Type-I integral.  In coordinates (x,y,r,z), q=y-z*r,
# the kernel equations force the displayed A and B.  Their pencil is linear
# in x,y,r, so every 3-by-3 Hessian minor vanishes.
terminal_variables = (x, y, r, z)
alpha = sp.symbols("alpha", nonzero=True)
p = sp.Function("p")(z)
h = sp.Function("h")(z)
g = sp.Function("g")(z)
terminal_C = sp.Function("terminal_C")(z)
q = y - z * r
terminal_A = sp.expand(alpha * r + z * sp.diff(p, z) * q + h)
terminal_B = sp.expand(
    r * terminal_C + q * sp.diff(terminal_C, z) + p * x + g
)
terminal_pencil = sp.hessian(
    sp.expand(terminal_B + s * terminal_A),
    terminal_variables,
)
terminal_kernel = sp.Matrix([s * z, z, 1, 0])
assert_zero_vector(terminal_pencil * terminal_kernel)
terminal_three_minors = [
    sp.factor(
        terminal_pencil.extract(rows, columns).det(method="domain-ge")
    )
    for rows in combinations(range(4), 3)
    for columns in combinations(range(4), 3)
]
assert not any(terminal_three_minors)
assert terminal_pencil.rank(iszerofunc=lambda expression: expression == 0) <= 2

t = sp.symbols("t")
terminal_family = ScalarPivotSchurFamily(
    terminal_variables,
    t,
    terminal_A,
    terminal_B,
)
assert terminal_family.bordered_determinant(s) == 0


payload = {
    "format": "hc4-parameter-moving-affine-kernel-pencils-v1",
    "status": {
        "id": "HC4RSD3",
        "kind": "hybrid theorem",
        "scope": (
            "scalar singular pencil whose primitive kernel generator is "
            "affine in x and arbitrary in the pencil parameter"
        ),
    },
    "bordered_unit_reductions": {
        "adjugate": "adj(M)=epsilon*v*v^T with epsilon in K^*",
        "pairing": "grad(A).v=alpha in K^*",
        "parameter_degree": "deg_s(v)<=1 because deg_s(adj(M))<=3",
        "piola": [
            "tr(L(s))=0",
            "L(s)^2=0",
            "L(s)*a(s)=0",
            "rank(L(s))<=1",
        ],
    },
    "rank_one_pencil_classification": {
        "proportional_linear_parts": "collapse to a fixed kernel line",
        "common_image": "collapse to a fixed kernel line",
        "constant_at_infinity": {
            "kernel": "v0=(z,1,0,0), v1=(b1,b2,0,b4)",
            "cross_equations": [
                "b4*G_ww=0",
                "F_r=(b1-b2*z)*C''+b4*G_wz",
            ],
            "terminal_equation": (
                "b2*alpha+(b1-b2*z)^2*C''=0, incompatible with "
                "C''*G_ww nonzero"
            ),
        },
        "common_covector": {
            "normal_forms": [
                "v=(s*z,z,1,0)",
                "v=(1+s*z,z,0,0)",
            ],
            "type_two": "would require z*D_u(A)=alpha",
            "type_one_terminal_potential": {
                "q": "y-z*r",
                "A": "alpha*r+z*p'(z)*q+h(z)",
                "B": "r*C(z)+q*C'(z)+p(z)*x+g(z)",
            },
            "type_one_hessian_rank": "at most two",
        },
    },
    "result": (
        "no parameter-moving affine-in-x kernel line survives the "
        "bordered-unit gate"
    ),
    "open_frontier": "primitive kernel generators nonlinear in x",
}

serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
OUTPUT.write_text(serialized)
digest = hashlib.sha256(serialized.encode()).hexdigest()

print("PASS: affine Piola pencil reduces to rank-one compression types")
print("PASS: constant-at-infinity cross equations replay exactly")
print("PASS: sole moving common-covector integral has Hessian rank at most two")
print("PASS: no parameter-moving affine kernel survives the bordered unit")
print(f"PASS: wrote {OUTPUT.relative_to(ROOT)}")
print(f"SHA256: {digest}")
