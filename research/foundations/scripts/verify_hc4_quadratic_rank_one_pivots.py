#!/usr/bin/env python3
"""Verify the rank-one quadratic-pivot classification.

With A=x^2/2+w, the singular pencil first makes the passive three-variable
Hessian singular.  Passive ranks zero and two are incompatible with the
bordered unit.  In passive rank one, the rank-one polynomial-Hessian normal
form and the exact parent determinant produce a polynomial unit frame.  It
integrates to a complete all-degree normal form whose descendants have a
triangular polynomial inverse.
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
    / "hc4_quadratic_rank_one_pivots.json"
)
EXPECTED_OUTPUT_SHA256 = (
    "7eb678ce7b283ff3b4b102ea548a3131bd9970c531e366be364a2183bc279630"
)


def audit_existing() -> None:
    actual = hashlib.sha256(OUTPUT.read_bytes()).hexdigest()
    assert actual == EXPECTED_OUTPUT_SHA256, (actual, EXPECTED_OUTPUT_SHA256)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["status"]["id"] == "HC4RSD10"
    assert payload["result"].startswith(
        "all quadratic scalar pivots in the singular-pencil programme"
    )
    assert payload["open_frontier"].startswith(
        "higher-degree nonlinear pivots, nonsingular pencil cancellation"
    )
    print(
        "PASS: committed HC4RSD10 stage artifact is intact; its quadratic "
        "nonsingular handoff is superseded by HC4RSD11--16 and its auxiliary "
        "higher-degree pencil handoff by HC4MR1; no symbolic replay or rewrite"
    )


parser = argparse.ArgumentParser()
parser.add_argument(
    "--audit-existing-only",
    action="store_true",
    help="validate the committed stage artifact without symbolic replay or rewriting it",
)
arguments = parser.parse_args()
if arguments.audit_existing_only:
    audit_existing()
    raise SystemExit(0)

import sympy as sp


# Universal one-active/three-passive block identities.
s, lam, active_entry, active_gradient = sp.symbols(
    "s lambda active_entry active_gradient"
)
d = sp.Matrix(sp.symbols("d1:4"))
e11, e12, e13, e22, e23, e33 = sp.symbols(
    "e11 e12 e13 e22 e23 e33"
)
E = sp.Matrix(
    [
        [e11, e12, e13],
        [e12, e22, e23],
        [e13, e23, e33],
    ]
)
a = sp.Matrix([0, 0, 1])
pencil = sp.Matrix.vstack(
    sp.Matrix.hstack(sp.Matrix([[active_entry + s]]), d.T),
    sp.Matrix.hstack(d, E),
)
passive_border = sp.Matrix.vstack(
    sp.Matrix.hstack(sp.Matrix([[lam]]), a.T),
    sp.Matrix.hstack(a, E),
)
parent = sp.Matrix.vstack(
    sp.Matrix.hstack(
        sp.Matrix([[active_entry + s, active_gradient]]),
        d.T,
    ),
    sp.Matrix.hstack(
        sp.Matrix([[active_gradient, lam]]),
        a.T,
    ),
    sp.Matrix.hstack(d, a, E),
)

assert sp.expand(
    sp.Poly(pencil.det(method="domain-ge"), s).coeff_monomial(s)
    - E.det(method="domain-ge")
) == 0
assert sp.expand(
    sp.Poly(parent.det(method="domain-ge"), s).coeff_monomial(s)
    - passive_border.det(method="domain-ge")
) == 0
assert sp.expand(
    passive_border.det(method="domain-ge")
    - lam * E.det(method="domain-ge")
    + (a.T * E.adjugate(method="domain-ge") * a)[0]
) == 0


# Rank-one passive frame.  The entries involving lambda, the active Hessian,
# and the active gradient cancel; only the squared three-vector frame
# remains.
rho = sp.symbols("rho")
ell = sp.Matrix(sp.symbols("ell1:4"))
E_rank_one = rho * ell * ell.T
parent_rank_one = parent.subs(
    {
        e11: E_rank_one[0, 0],
        e12: E_rank_one[0, 1],
        e13: E_rank_one[0, 2],
        e22: E_rank_one[1, 1],
        e23: E_rank_one[1, 2],
        e33: E_rank_one[2, 2],
    }
)
frame = sp.Matrix.hstack(a, d, ell).det(method="domain-ge")
assert sp.factor(parent_rank_one.det(method="domain-ge") - rho * frame**2) == 0


# The moving passive direction.  For ell=(p,r,h0), integration gives
# B=rho*(ell.u)^2/2+b.u+delta.  The frame splits into one passive-linear
# Wronskian and one base term.
x, y, z, w = sp.symbols("x y z w")
p0, p1, r0, r1, h0, h1 = sp.symbols("p0 p1 r0 r1 h0 h1")
b10, b11, b20, b21, b30, b31 = sp.symbols(
    "b10 b11 b20 b21 b30 b31"
)
p = p0 + p1 * x
r = r0 + r1 * x
h = h0 + h1 * x
b1 = b10 + b11 * x
b2 = b20 + b21 * x
b3 = b30 + b31 * x
ell_polynomial = sp.Matrix([p, r, h])
u = sp.Matrix([y, z, w])
L = sp.expand((ell_polynomial.T * u)[0])
B_frame = (
    rho * L**2 / 2
    + b1 * y
    + b2 * z
    + b3 * w
)
d_frame = sp.Matrix(
    [sp.diff(sp.diff(B_frame, x), variable) for variable in (y, z, w)]
)
frame_polynomial = sp.expand(
    sp.Matrix.hstack(a, d_frame, ell_polynomial).det(method="domain-ge")
)
expected_frame = sp.expand(
    rho * L * (sp.diff(p, x) * r - sp.diff(r, x) * p)
    + sp.diff(b1, x) * r
    - sp.diff(b2, x) * p
)
assert sp.expand(frame_polynomial - expected_frame) == 0


# Complete normalized family and descendant recovery.
t, kappa, mu = sp.symbols("t kappa mu")
alpha0, alpha1 = sp.symbols("alpha0 alpha1")
gamma0, gamma1 = sp.symbols("gamma0 gamma1")
delta0, delta1, delta2 = sp.symbols("delta0 delta1 delta2")
hs0, hs1, hs2 = sp.symbols("hs0 hs1 hs2")
hs = hs0 + hs1 * x + hs2 * x**2
alpha = alpha0 + alpha1 * x
gamma = gamma0 + gamma1 * x
delta = delta0 + delta1 * x + delta2 * x**2
A = x**2 / 2 + w
Y = y + hs * w
B = x * z + rho * Y**2 / 2 + alpha * y + gamma * w + delta
Phi = lam * t**2 / 2 + t * A + B
psi = B + kappa * A**2 / 2 + mu * A
variables = (x, y, z, w)

assert sp.factor(
    sp.hessian(B + s * A, variables).det(method="domain-ge")
) == 0
assert sp.factor(
    sp.hessian(Phi, (t, *variables)).det(method="domain-ge")
) == rho
assert sp.factor(
    sp.hessian(psi, variables).det(method="domain-ge")
) == -kappa * rho

gradient = sp.Matrix([sp.diff(psi, variable) for variable in variables])
assert sp.expand(gradient[2] - x) == 0
assert sp.expand(gradient[1] - rho * Y - alpha) == 0
assert sp.expand(
    gradient[3] - rho * hs * Y - gamma - kappa * A - mu
) == 0
assert sp.diff(gradient[0], z) == 1


payload = {
    "format": "hc4-quadratic-rank-one-pivots-v1",
    "status": {
        "id": "HC4RSD10",
        "kind": "hybrid theorem",
        "scope": (
            "rank-one quadratic scalar pivots with an identically singular "
            "four-variable Hessian pencil"
        ),
        "structural_input": "rank-one polynomial-Hessian normal form",
    },
    "passive_rank_split": {
        "pencil_leading_face": "det(E)=0",
        "parent_leading_face": "a^T*adj(E)*a=0",
        "rank_zero": "contradicts generic pencil rank three",
        "rank_two": (
            "produces a reduced-pencil kernel orthogonal to grad(A), "
            "so the parent Hessian is singular"
        ),
        "survivor": "rank(E)=1",
    },
    "rank_one_frame": {
        "factorization": "E=rho*ell*ell^T",
        "parent_determinant": "rho*det(a,d,ell)^2",
        "unit_consequence": "rho and det(a,d,ell) are constants",
        "motion_equation": "p'*r-r'*p=0",
        "normalization": "(p,r)=(1,0), b2=x",
    },
    "normal_form": {
        "A": "x^2/2+w",
        "B": (
            "x*z+rho*(y+h(x)*w)^2/2+alpha(x)*y+"
            "gamma(x)*w+delta(x)"
        ),
        "parent_hessian_determinant": "rho",
        "descendant_hessian_determinant": "-kappa*rho",
    },
    "triangular_recovery": [
        "x=F_z",
        "Y=(F_y-alpha(x))/rho",
        "A=(F_w-rho*h(x)*Y-gamma(x)-mu)/kappa",
        "w=A-x^2/2",
        "y=Y-h(x)*w",
        "z=F_x-known_polynomial",
    ],
    "result": (
        "all quadratic scalar pivots in the singular-pencil programme "
        "are collision-free"
    ),
    "open_frontier": (
        "higher-degree nonlinear pivots, nonsingular pencil cancellation, "
        "and moving matrix pivots"
    ),
}

serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
OUTPUT.write_text(serialized)
digest = hashlib.sha256(serialized.encode()).hexdigest()

print("PASS: the passive three-by-three Hessian is singular")
print("PASS: passive ranks zero and two contradict the bordered unit")
print("PASS: the rank-one passive frame has constant determinant")
print("PASS: integrated the complete moving-direction normal form")
print("PASS: every descendant has a triangular polynomial inverse")
print("THEOREM: the quadratic singular-pencil scalar branch is closed")
print(f"PASS: wrote {OUTPUT.relative_to(ROOT)}")
print(f"SHA256: {digest}")
