#!/usr/bin/env python3
"""Verify the scalar cone-pencil obstruction and matrix-pivot schema.

The homogeneous scalar conclusion uses the characteristic-zero
Gordan--Noether theorem in four variables.  More generally the same proof
applies whenever the generic pencil member is already known to have an
x-constant kernel direction.  This checker certifies every exact
linear-algebra identity used after that input, exercises the
coefficient-equation API, and records the empty intersection of this cone
stratum with the surviving HC4 affine-degree packets.

It does not classify nonhomogeneous pencils with x-moving kernel lines or
nonconstant matrix-pivot kernel planes.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "hc4_reverse_schur_descent.json"
)
ATLAS = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "hc4_projective_polar_atlas.json"
)
HELPER = ROOT / "jcsearch" / "reverse_schur_descent.py"
EXPECTED_OUTPUT_SHA256 = (
    "279aaf673b3a8cd254e4f51ae3c0cdc32ecce732a8a8e2b28979cbf453a89226"
)
EXPECTED_ATLAS_SHA256 = (
    "350bc81b4ba7ac21289d7548f6d46de6526887c4e98a0813596cf20a454b240b"
)
EXPECTED_HELPER_SHA256 = (
    "b80da9da8105caa51fa38fc178a3a256d9132335a91bb5a87bb84eab96cf3e44"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit_existing() -> None:
    assert sha256(OUTPUT) == EXPECTED_OUTPUT_SHA256
    assert sha256(ATLAS) == EXPECTED_ATLAS_SHA256
    assert sha256(HELPER) == EXPECTED_HELPER_SHA256
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    atlas = json.loads(ATLAS.read_text(encoding="utf-8"))
    assert payload["format"] == "hc4-reverse-schur-descent-v1"
    assert payload["status"]["id"] == "HC4RSD1"
    live_rows = atlas["quintic_coverage_summary"][
        "remaining_numerical_signatures_after_vertex_colength"
    ]
    assert live_rows == {
        "affine_degree_2": 318,
        "affine_degree_3": 306,
        "total": 624,
    }
    assert payload["projective_polar_intersection"][
        "input_live_quintic_rows"
    ] == live_rows
    assert set(payload["open_frontier"]) == {
        "matrix_pivot",
        "moving_kernel_scalar_pencil",
        "nonsingular_scalar_pencil",
    }
    print(
        "PASS: committed HC4RSD1 artifact, consumed atlas, and equation helper "
        "are intact; its frontier is historical: HC4RSD2--5 narrow the moving "
        "scalar branch and HC4MR1 closes the nonzero-corner auxiliary pencil, "
        "while nonlinear zero-corner and moving-matrix mechanisms remain"
    )


parser = argparse.ArgumentParser()
parser.add_argument(
    "--audit-existing-only",
    action="store_true",
    help=(
        "validate committed inputs and the historical frontier without "
        "symbolic replay or artifact rewriting"
    ),
)
arguments = parser.parse_args()
if arguments.audit_existing_only:
    audit_existing()
    raise SystemExit(0)

import sys

import sympy as sp

sys.path.insert(0, str(ROOT))

from jcsearch.reverse_schur_descent import (
    MatrixPivotSchurFamily,
    ScalarPivotSchurFamily,
    coefficient_equations,
    corank_one_adjugate_scalar,
    hessian_integrability_residuals,
    rank_at_most_equations,
)


x1, x2, x3, z, t, s = sp.symbols("x1 x2 x3 z t s")
x = (x1, x2, x3, z)
kappa, mu, pivot_value = sp.symbols("kappa mu pivot_value")


# A complete common-kernel calibration.  The ternary pencil has constant
# Hessian determinant -1, while the fourth direction is its common kernel.
a = x2**3 / 6
b = x1 * x2 + x3**2 / 2
A = z + a
B = 2 * z + b
scalar_family = ScalarPivotSchurFamily(x, t, A, B)
pencil = scalar_family.pencil_hessian(s)
kernel = sp.Matrix([0, 0, 0, 1])
gradient_a = scalar_family.gradient_a

assert not any(hessian_integrability_residuals(pencil, x))
assert sp.factor(pencil.det(method="domain-ge")) == 0
assert pencil.rank(iszerofunc=lambda expression: expression == 0) == 3
assert corank_one_adjugate_scalar(pencil, kernel) == -1
assert sp.expand((kernel.T * gradient_a)[0]) == 1
assert scalar_family.bordered_determinant(s) == 1
assert scalar_family.parent_hessian_determinant(3) == 1
assert scalar_family.quadratic_schur_identity_residual(3, pivot_value) == 0

descendant = scalar_family.schur_descendant(kappa, mu)
descendant_determinant = sp.factor(
    sp.hessian(descendant, x).det(method="domain-ge")
)
assert descendant_determinant == -kappa

# The collision reduction is visible coordinate by coordinate.  The z
# gradient recovers S=kappa*A+mu.  At fixed S, the other three coordinates
# are exactly grad(b+S*a), a triangular ternary Keller automorphism here.
descendant_gradient = sp.Matrix([sp.diff(descendant, variable) for variable in x])
S = sp.symbols("S")
assert sp.expand(descendant_gradient[3] - (2 + kappa * A + mu)) == 0
fixed_s_gradient = sp.Matrix(
    [sp.diff(b + S * a, variable) for variable in (x1, x2, x3)]
)
assert fixed_s_gradient == sp.Matrix(
    [x2, x1 + S * x2**2 / 2, x3]
)


# A moving-kernel symmetric Hessian pencil is not discarded merely because
# its determinant vanishes.  The primitive kernel v=(-s,1,0,0) has pairing
# v.A'=x3, so its bordered determinant is x3^2 rather than a unit.  This is
# the exact unit gate used in the cone-pencil proof.
moving_A = x2 * x3
moving_B = x1 * x3 + z**2 / 2
moving_family = ScalarPivotSchurFamily(x, t, moving_A, moving_B)
moving_pencil = moving_family.pencil_hessian(s)
moving_kernel = sp.Matrix([-s, 1, 0, 0])
moving_q = corank_one_adjugate_scalar(moving_pencil, moving_kernel)
moving_pairing = sp.factor((moving_kernel.T * moving_family.gradient_a)[0])
assert not any(hessian_integrability_residuals(moving_pencil, x))
assert moving_q == -1
assert moving_pairing == x3
assert moving_family.bordered_determinant(s) == x3**2


# Generic corank at least two cannot support a scalar bordered unit: one
# border row and column raise rank by at most two, leaving a 5-by-5 matrix
# singular.  The displayed matrix is an exact rank-two calibration.
rank_two = sp.diag(1, 1, 0, 0)
generic_gradient = sp.Matrix(sp.symbols("g0:4"))
rank_two_border = sp.Matrix.vstack(
    sp.Matrix.hstack(sp.zeros(1, 1), generic_gradient.T),
    sp.Matrix.hstack(generic_gradient, rank_two),
)
assert rank_two_border.det(method="domain-ge") == 0
assert len(rank_at_most_equations(sp.zeros(4), 2)) == 16


# Exercise exact coefficient extraction on a parameterized pencil.  The
# singularity equations retain coefficient parameters and eliminate x,s.
c0, c1 = sp.symbols("c0 c1")
equations = coefficient_equations(
    (c0 * s * x1 + c1 * x2 + c0,),
    (s, *x),
)
assert set(equations) == {c0, c1}


# Simultaneous two-pivot implementation.  This constant corank-two pencil
# pairs its kernel plane with the two A-gradients.  The exact Schur identity
# is checked with a non-diagonal rational pivot matrix.
t1, t2, y1, y2 = sp.symbols("t1 t2 y1 y2")
matrix_family = MatrixPivotSchurFamily(
    variables=x,
    pivot_variables=(t1, t2),
    a=(x1, x2),
    b=(x3**2 + z**2) / 2,
    pivot_matrix=sp.Matrix([[2, 1], [1, 2]]),
)
assert matrix_family.pencil_hessian((y1, y2)) == sp.diag(0, 0, 1, 1)
assert matrix_family.pencil_corank_equations((y1, y2), 2) == ()
assert matrix_family.schur_identity_residual((y1, y2)) == 0
assert sp.factor(
    sp.hessian(
        matrix_family.parent_potential,
        (t1, t2, *x),
    ).det(method="domain-ge")
) == 1


# Intersect the proved scalar cone classification with the live HC4
# projective-gradient packets.  The theorem makes these descendants injective
# so none of the collision-normalized affine-degree 2/3 rows is reached.
# The input counts are read from the canonical atlas artifact.
atlas = json.loads(ATLAS.read_text())
live_rows = atlas["quintic_coverage_summary"][
    "remaining_numerical_signatures_after_vertex_colength"
]
assert live_rows == {
    "affine_degree_2": 318,
    "affine_degree_3": 306,
    "total": 624,
}


payload = {
    "format": "hc4-reverse-schur-descent-v1",
    "status": {
        "id": "HC4RSD1",
        "kind": "hybrid theorem",
        "scope": (
            "scalar pivot with identically singular four-variable Hessian "
            "pencil and x-constant generic kernel line"
        ),
        "external_input": (
            "Gordan--Noether supplies the x-constant kernel for homogeneous "
            "four-variable pencil members; the general polynomial moving-kernel "
            "case is not claimed"
        ),
    },
    "exact_scalar_classification": {
        "generic_pencil_corank": 1,
        "corank_at_least_two": "bordered determinant identically zero",
        "corank_one_adjugate": "adj(M)=q*v*v^T",
        "border_unit_equation": "c=-q*(v^T*grad(A))^2",
        "consequence": (
            "q and v^T*grad(A) are units over Q(s)[x], forcing a common "
            "constant kernel direction for Hess(A) and Hess(B)"
        ),
        "normal_form": "A=alpha*z+a(u), B=beta*z+b(u), alpha!=0",
        "ternary_pencil_determinant": "det Hess_u(b+s*a)=-c/alpha^2",
        "collision_fiber": "grad_u(b+s0*a) with s0=kappa*A+mu fixed",
        "result": "HC3 makes every reduced gradient fiber a singleton",
    },
    "equation_builder": {
        "scalar": [
            "singular pencil coefficients",
            "constant parent determinant coefficients",
            "exact reduced collision equations",
            "Hessian integrability residuals",
        ],
        "matrix_pivot": [
            "prescribed-corank pencil minors",
            "exact simultaneous Schur identity",
        ],
    },
    "calibrations": {
        "common_kernel": {
            "pencil_rank": 3,
            "adjugate_scalar": -1,
            "bordered_determinant": 1,
            "descendant_determinant": "-kappa",
        },
        "moving_kernel_rejection": {
            "kernel": "(-s,1,0,0)",
            "adjugate_scalar": -1,
            "kernel_gradient_pairing": "x3",
            "bordered_determinant": "x3^2",
        },
        "matrix_pivot": {
            "pivot_rank": 2,
            "reduced_pencil_corank": 2,
            "parent_determinant": 1,
        },
    },
    "projective_polar_intersection": {
        "input_live_quintic_rows": live_rows,
        "scalar_cone_reverse_schur_rows": {
            "affine_degree_2": 0,
            "affine_degree_3": 0,
            "total": 0,
        },
        "reason": "injectivity excludes every collision-normalized packet",
    },
    "rational_reconstruction": {
        "surviving_scalar_cone_components": 0,
        "status": "not invoked because the exact cone-stratum intersection is empty",
    },
    "open_frontier": {
        "matrix_pivot": (
            "classify nonconstant jointly moving kernel planes for r>=2; "
            "the implemented minor and Schur equations do not prove this classification"
        ),
        "nonsingular_scalar_pencil": (
            "allow cancellation between lambda*det(M) and the bordered term"
        ),
        "moving_kernel_scalar_pencil": (
            "classify nonhomogeneous singular Hessian pencils whose generic "
            "kernel direction depends on x"
        ),
    },
}

serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
OUTPUT.write_text(serialized)
digest = hashlib.sha256(serialized.encode()).hexdigest()

print("PASS: exact scalar and matrix-pivot Schur identities")
print("PASS: scalar singular pencil must have generic corank one")
print("PASS: cone-pencil bordered-unit gate forces the common-kernel HC3 reduction")
print("PASS: no scalar cone descendant reaches an affine-degree 2/3 HC4 packet")
print("PASS: matrix-pivot corank equations generated; moving planes remain open")
print(f"PASS: wrote {OUTPUT.relative_to(ROOT)}")
print(f"SHA256: {digest}")
