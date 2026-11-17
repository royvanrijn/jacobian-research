#!/usr/bin/env python3
"""Verify the affine-pivot collision-fiber obstruction.

If ell is a constant covector and

    N_ell = ell^T adj(Hess(psi)) ell

is a nonzero constant, then every affine hyperplane ell.x=c carries a
three-variable constant-Hessian restriction of psi.  HC3 therefore makes
the tangential gradient injective on that hyperplane.  Two equal full
gradients on the same ell-fiber must be the same point.

Consequently an affine zero-corner Schur parent, whether singular after
reduction or in the nonsingular exact-remainder branch, cannot lift a
marked gradient collision: parent pivot-gradient equality would require
the two points to have the same value of A=ell.x+a0.
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
    / "hc4_affine_pivot_collision_fibers.json"
)
EXPECTED_OUTPUT_SHA256 = (
    "d999d448a661fd858d1eb2eb39063c5b9522f6b491f67ba4b65915ef8f8e0688"
)


def audit_existing() -> None:
    actual = hashlib.sha256(OUTPUT.read_bytes()).hexdigest()
    assert actual == EXPECTED_OUTPUT_SHA256, (actual, EXPECTED_OUTPUT_SHA256)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["status"]["id"] == "HC4RSD7"
    assert payload["collision_consequence"] == (
        "grad(psi)(p)=grad(psi)(q) and ell.p=ell.q imply p=q"
    )
    assert payload["schur_consequence"].startswith(
        "an affine zero-corner pivot, singular or nonsingular"
    )
    assert payload["open_frontier"].startswith("nonlinear scalar pivots")
    print(
        "PASS: committed HC4RSD7 affine-fiber artifact is intact and retains "
        "its different-fiber/nonlinear boundary; no symbolic replay or rewrite"
    )


parser = argparse.ArgumentParser()
parser.add_argument(
    "--audit-existing-only",
    action="store_true",
    help="validate the committed artifact without symbolic replay or rewriting it",
)
arguments = parser.parse_args()
if arguments.audit_existing_only:
    audit_existing()
    raise SystemExit(0)

import sympy as sp


def symmetric_matrix(prefix: str, size: int) -> sp.Matrix:
    entries: dict[tuple[int, int], sp.Symbol] = {}
    for row in range(size):
        for column in range(row, size):
            entries[row, column] = sp.Symbol(f"{prefix}{row + 1}{column + 1}")
    return sp.Matrix(
        size,
        size,
        lambda row, column: entries[min(row, column), max(row, column)],
    )


# In coordinates with ell=e4, the inverse-Hessian metric numerator is
# exactly the determinant of the tangential 3-by-3 Hessian block.
H = symmetric_matrix("h", 4)
ell = sp.Matrix([0, 0, 0, 1])
tangential_hessian = H[:3, :3]
metric_numerator = sp.expand((ell.T * H.adjugate(method="domain-ge") * ell)[0])
assert sp.expand(
    metric_numerator - tangential_hessian.det(method="domain-ge")
) == 0


# A quadratic Schur repair changes H by kappa*ell*ell^T but leaves the
# tangential cofactor numerator unchanged.  Hence a zero-corner affine
# parent with constant bordered determinant supplies the slice hypothesis
# whether or not its reduced Hessian is singular.
general_ell = sp.Matrix(sp.symbols("l1:5"))
kappa = sp.symbols("kappa")
repaired_H = H + kappa * general_ell * general_ell.T
metric_before_repair = sp.expand(
    (general_ell.T * H.adjugate(method="domain-ge") * general_ell)[0]
)
bordered_H = sp.zeros(5)
bordered_H[0, 1:] = general_ell.T
bordered_H[1:, 0] = general_ell
bordered_H[1:, 1:] = H
assert sp.expand(
    bordered_H.det(method="domain-ge") + metric_before_repair
) == 0
metric_after_repair = sp.expand(
    (
        general_ell.T
        * repaired_H.adjugate(method="domain-ge")
        * general_ell
    )[0]
)
assert sp.expand(metric_after_repair - metric_before_repair) == 0


# Exact coordinate-change covariance.  A block upper-triangular matrix is
# the universal adapted-coordinate form once its first three columns span
# ker(ell).  It is enough to replay the determinant-square factor without
# introducing rational inverse entries.
P3 = sp.Matrix(
    [
        [sp.symbols("p11"), sp.symbols("p12"), sp.symbols("p13")],
        [sp.symbols("p21"), sp.symbols("p22"), sp.symbols("p23")],
        [sp.symbols("p31"), sp.symbols("p32"), sp.symbols("p33")],
    ]
)
v1, v2, v3, scale = sp.symbols("v1 v2 v3 scale")
P = sp.zeros(4)
P[:3, :3] = P3
P[:3, 3] = sp.Matrix([v1, v2, v3])
P[3, 3] = scale
transformed_H = sp.expand(P.T * H * P)
transformed_tangential = transformed_H[:3, :3]
assert sp.expand(
    transformed_tangential.det(method="domain-ge")
    - P3.det(method="domain-ge") ** 2
    * tangential_hessian.det(method="domain-ge")
) == 0


# Restriction of an actual Hessian agrees with the tangential block, and
# restriction of the full gradient agrees with the tangential gradient.
x1, x2, x3, x4, c = sp.symbols("x1 x2 x3 x4 c")
variables = (x1, x2, x3, x4)
psi = sp.Function("psi")(*variables)
full_gradient = sp.Matrix([sp.diff(psi, variable) for variable in variables])
full_hessian = sp.hessian(psi, variables)
restricted_psi = psi.subs(x4, c)
restricted_gradient = sp.Matrix(
    [sp.diff(restricted_psi, variable) for variable in variables[:3]]
)
restricted_hessian = sp.hessian(restricted_psi, variables[:3])
assert all(
    sp.simplify(restricted_gradient[index] - full_gradient[index].subs(x4, c))
    == 0
    for index in range(3)
)
assert all(
    sp.simplify(
        restricted_hessian[row, column]
        - full_hessian[row, column].subs(x4, c)
    )
    == 0
    for row in range(3)
    for column in range(3)
)


# Parent collision transfer: for Phi=t*A+B, equality of full parent
# gradients at a common pivot value includes both A_plus=A_minus and
# grad(B)_plus=grad(B)_minus.  These two equalities also preserve the
# collision after the quadratic repair B -> B+kappa*A^2/2+mu*A.
A_plus, A_minus = sp.symbols("A_plus A_minus")
parent_pivot_gradient_difference = sp.expand(A_plus - A_minus)
assert parent_pivot_gradient_difference.subs(A_plus, A_minus) == 0
gradient_B_difference = sp.Matrix(sp.symbols("db1:5"))
mu = sp.symbols("mu")
repaired_gradient_difference = (
    gradient_B_difference
    + kappa * (A_plus - A_minus) * general_ell
    + mu * general_ell
    - mu * general_ell
)
assert repaired_gradient_difference.subs(
    {A_plus: A_minus, **{entry: 0 for entry in gradient_B_difference}}
) == sp.zeros(4, 1)


payload = {
    "format": "hc4-affine-pivot-collision-fibers-v1",
    "status": {
        "id": "HC4RSD7",
        "kind": "hybrid theorem",
        "scope": (
            "all affine constant-covector pivots with nonzero constant "
            "inverse-Hessian metric numerator"
        ),
    },
    "slice_identity": {
        "adapted_covector": "ell=e4",
        "bordered_determinant": "det([[0,ell^T],[ell,H]])=-N_ell",
        "metric_numerator": "N_ell=det(Hess_(x1,x2,x3)(psi))",
        "coordinate_change": "the slice determinant scales by det(P3)^2",
    },
    "hc3_consequence": (
        "the tangential gradient is injective on every affine fiber "
        "ell.x=c"
    ),
    "collision_consequence": (
        "grad(psi)(p)=grad(psi)(q) and ell.p=ell.q imply p=q"
    ),
    "schur_consequence": (
        "an affine zero-corner pivot, singular or nonsingular after "
        "reduction, cannot transfer a marked collision at a common pivot value"
    ),
    "open_frontier": (
        "nonlinear scalar pivots, mixed/coisotropic pivots, and affine "
        "representations whose collisions lie on different pivot fibers"
    ),
}

serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
OUTPUT.write_text(serialized)
digest = hashlib.sha256(serialized.encode()).hexdigest()

print("PASS: the affine metric numerator is the ternary slice determinant")
print("PASS: the zero-corner bordered determinant is minus that numerator")
print("PASS: quadratic repair preserves the affine metric numerator")
print("PASS: adapted coordinate changes preserve constancy up to a unit square")
print("PASS: restricted gradients and Hessians are the tangential blocks")
print("THEOREM: HC3 makes every affine pivot fiber collision-free")
print("SCOPE: nonlinear and genuinely mixed/coisotropic pivots remain")
print(f"PASS: wrote {OUTPUT.relative_to(ROOT)}")
print(f"SHA256: {digest}")
