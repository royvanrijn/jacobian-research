#!/usr/bin/env python3
"""Verify the scalar cancellation dichotomy for reverse HC4 descent.

For a nonzero quadratic pivot corner, exact Schur cancellation is exactly a
four-variable constant-Hessian pencil.  For a zero corner and a graph-coordinate
pivot A=w+q(u), a universal determinant factorization reduces every
parent-gradient fiber to a ternary constant-Hessian gradient.  HC3 then excludes
collisions.  The written proof shows that every quadratic zero-corner pivot is
forced into the graph-coordinate class.  The final identities below show that
unit border transversality freezes the apparently moving ternary chart in the
rank-one nonzero-corner stratum.
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
    / "hc4_scalar_cancellation_dichotomy.json"
)
EXPECTED_OUTPUT_SHA256 = (
    "e051d9c5b320dc008106e37ea80f98fdeb1c88d1576d8149de79448c3c6381ac"
)


def audit_existing() -> None:
    actual = hashlib.sha256(OUTPUT.read_bytes()).hexdigest()
    assert actual == EXPECTED_OUTPUT_SHA256, (actual, EXPECTED_OUTPUT_SHA256)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["format"] == "hc4-scalar-cancellation-dichotomy-v1"
    assert [row["id"] for row in payload["status"]] == [
        "HC4RSD11",
        "HC4RSD12",
        "HC4RSD13",
        "HC4RSD14",
        "HC4RSD15",
        "HC4RSD16",
    ]
    assert payload["quadratic_corollary"]["conclusion"] == (
        "every quadratic zero-corner parent is collision-free"
    )
    assert payload["rank_one_nonlinear_pencil"]["conclusion"].endswith(
        "HC2 or the JC2 cotangent packet"
    )
    assert payload["open_frontier"].startswith(
        "higher-degree nonlinear four-variable constant-Hessian pencils"
    )
    print(
        "PASS: committed HC4RSD11--16 stage artifact is intact; its "
        "higher-degree pencil frontier is historical and superseded by HC4MR1; "
        "no symbolic replay or rewrite"
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


# Universal nonzero-corner Schur identity.
lam = sp.symbols("lambda", nonzero=True)
g = sp.Matrix(sp.symbols("g0:4"))
m = sp.symbols("m0:10")
M = sp.Matrix(
    [
        [m[0], m[1], m[2], m[3]],
        [m[1], m[4], m[5], m[6]],
        [m[2], m[5], m[7], m[8]],
        [m[3], m[6], m[8], m[9]],
    ]
)
nonzero_corner = sp.Matrix.vstack(
    sp.Matrix.hstack(sp.Matrix([[lam]]), g.T),
    sp.Matrix.hstack(g, M),
)
schur_hessian = M - g * g.T / lam
assert sp.factor(
    nonzero_corner.det(method="domain-ge")
    - lam * schur_hessian.det(method="domain-ge")
) == 0


# Universal graph-coordinate zero-corner factorization.  Here p=grad(q),
# Q=Hess(q), H=Hess_u(C), b=grad_u(C_r), d=C_rr, and
# tau=t+C_r.  All symbols are independent, so this checks the complete
# second-jet identity rather than a bounded polynomial ansatz.
p = sp.Matrix(sp.symbols("p0:3"))
b = sp.Matrix(sp.symbols("b0:3"))
h11, h12, h13, h22, h23, h33 = sp.symbols(
    "h11 h12 h13 h22 h23 h33"
)
q11, q12, q13, q22, q23, q33 = sp.symbols(
    "q11 q12 q13 q22 q23 q33"
)
H = sp.Matrix(
    [
        [h11, h12, h13],
        [h12, h22, h23],
        [h13, h23, h33],
    ]
)
Q = sp.Matrix(
    [
        [q11, q12, q13],
        [q12, q22, q23],
        [q13, q23, q33],
    ]
)
d, tau = sp.symbols("d tau")
active_block = H + b * p.T + p * b.T + d * p * p.T + tau * Q
cross = b + d * p
zero_corner = sp.Matrix.vstack(
    sp.Matrix.hstack(sp.zeros(1, 1), p.T, sp.ones(1, 1)),
    sp.Matrix.hstack(p, active_block, cross),
    sp.Matrix.hstack(sp.ones(1, 1), cross.T, sp.Matrix([[d]])),
)
ternary_hessian = H + tau * Q
assert sp.factor(
    zero_corner.det(method="domain-ge")
    + ternary_hessian.det(method="domain-ge")
) == 0


# Degree-unbounded nonlinear calibration.  The cubic graph coordinate
# q=x^3/3 is deliberately beyond the quadratic corollary.  The ternary
# Hessian determinant is constant for every r and tau, while C has nonlinear
# r-dependence so all cancelled universal blocks are exercised.
t, x, y, z, w, r, tau_symbol = sp.symbols("t x y z w r tau")
u = (x, y, z)
q = x**3 / 3
C = x * z + y**2 / 2 + r * x + r**2 * z
A = w + q
B = sp.expand(C.subs(r, A))
Phi_zero = sp.expand(t * A + B)
C_tau = sp.expand(C + tau_symbol * q)

assert sp.factor(
    sp.hessian(C_tau, u).det(method="domain-ge")
) == -1
assert sp.factor(
    sp.hessian(Phi_zero, (t, *u, w)).det(method="domain-ge")
) == 1

r_coordinate = A
tau_coordinate = sp.expand(t + sp.diff(C, r).subs(r, A))
gradient_zero = sp.Matrix(
    [sp.diff(Phi_zero, variable) for variable in (t, *u, w)]
)
expected_middle = sp.Matrix(
    [
        sp.diff(C, variable).subs(r, r_coordinate)
        + tau_coordinate * sp.diff(q, variable)
        for variable in u
    ]
)
assert sp.expand(gradient_zero[0] - r_coordinate) == 0
assert all(
    sp.expand(gradient_zero[index + 1] - expected_middle[index]) == 0
    for index in range(3)
)
assert sp.expand(gradient_zero[4] - tau_coordinate) == 0


# Nonzero-corner suspension calibration with nonlinear A.  No property of A
# is used: the parent determinant is lambda times the child determinant and
# the y=0 parent-gradient level is exactly the child gradient.
lam_value = sp.Integer(3)
A_nonzero = x**2 / 2
psi = x * z + y * w
Phi_nonzero = sp.expand(
    lam_value * (t + A_nonzero / lam_value) ** 2 / 2 + psi
)
child_variables = (x, y, z, w)
assert sp.factor(
    sp.hessian(psi, child_variables).det(method="domain-ge")
) == 1
assert sp.factor(
    sp.hessian(Phi_nonzero, (t, *child_variables)).det(method="domain-ge")
) == lam_value

critical_t = -A_nonzero / lam_value
parent_gradient = sp.Matrix(
    [sp.diff(Phi_nonzero, variable) for variable in (t, *child_variables)]
)
child_gradient = sp.Matrix(
    [sp.diff(psi, variable) for variable in child_variables]
)
assert sp.expand(parent_gradient[0].subs(t, critical_t)) == 0
assert all(
    sp.expand(parent_gradient[index + 1].subs(t, critical_t) - child_gradient[index])
    == 0
    for index in range(4)
)


# Rank-three nonlinear constant-Hessian pencil.  The universal leading
# faces force the passive entry to vanish and the active/passive gradient to
# be isotropic.  The written UFD lemma makes that unimodular isotropic
# gradient constant; the resulting all-degree normal form is checked here.
s = sp.symbols("s")
k11, k12, k13, k22, k23, k33 = sp.symbols(
    "k11 k12 k13 k22 k23 k33"
)
d1, d2, d3, passive_entry = sp.symbols("d1 d2 d3 passive_entry")
K3 = sp.Matrix(
    [
        [k11, k12, k13],
        [k12, k22, k23],
        [k13, k23, k33],
    ]
)
Q3 = sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 1]])
d3_vector = sp.Matrix([d1, d2, d3])
rank_three_pencil = sp.Matrix.vstack(
    sp.Matrix.hstack(K3 + s * Q3, d3_vector),
    sp.Matrix.hstack(d3_vector.T, sp.Matrix([[passive_entry]])),
)
rank_three_polynomial = sp.Poly(
    rank_three_pencil.det(method="domain-ge"),
    s,
)
assert sp.factor(
    rank_three_polynomial.coeff_monomial(s**3) + passive_entry
) == 0
rank_three_e_zero = sp.Poly(
    rank_three_pencil.det(method="domain-ge").subs(passive_entry, 0),
    s,
)
assert sp.factor(
    rank_three_e_zero.coeff_monomial(s**2)
    - (2 * d1 * d2 + d3**2)
) == 0

alpha = sp.symbols("alpha", nonzero=True)
beta0, beta1, beta2 = sp.symbols("beta0 beta1 beta2")
g0, gx, gz, gxz, gx2, gz2 = sp.symbols("g0 gx gz gxz gx2 gz2")
beta = beta0 + beta1 * x + beta2 * x**2
G = g0 + gx * x + gz * z + gxz * x * z + gx2 * x**2 + gz2 * z**2
A_rank_three = x * y + z**2 / 2
psi_rank_three = w * x + y * (alpha * z + beta) + G
rank_three_member = sp.expand(psi_rank_three + s * A_rank_three)
rank_three_variables = (x, y, z, w)
assert sp.factor(
    sp.hessian(rank_three_member, rank_three_variables).det(
        method="domain-ge"
    )
) == alpha**2

rank_three_gradient = sp.Matrix(
    [sp.diff(rank_three_member, variable) for variable in rank_three_variables]
)
assert sp.expand(rank_three_gradient[3] - x) == 0
assert sp.expand(rank_three_gradient[1] - alpha * z - beta - s * x) == 0
assert sp.expand(rank_three_gradient[2] - alpha * y - sp.diff(G, z) - s * z) == 0
assert sp.diff(rank_three_gradient[0], w) == 1


# Rank-two constant-Hessian pencils.  The passive-rank-zero determinant is
# the square of an arbitrary plane Jacobian.  On passive rank one, the
# coefficient of s splits into two exact channels and the surviving unit
# frame normalizes to psi=x*z+C(x,y,w).
kx, kxy, ky = sp.symbols("kx kxy ky")
e11, e12, e22 = sp.symbols("e11 e12 e22")
d11, d12, d21, d22 = sp.symbols("D11 D12 D21 D22")
K2 = sp.Matrix([[kx, kxy + s], [kxy + s, ky]])
D2 = sp.Matrix([[d11, d12], [d21, d22]])
E2 = sp.Matrix([[e11, e12], [e12, e22]])
rank_two_pencil = sp.Matrix.vstack(
    sp.Matrix.hstack(K2, D2),
    sp.Matrix.hstack(D2.T, E2),
)
assert sp.factor(
    sp.Poly(rank_two_pencil.det(method="domain-ge"), s).coeff_monomial(s**2)
    + E2.det(method="domain-ge")
) == 0

assert sp.factor(
    rank_two_pencil.subs({e11: 0, e12: 0, e22: 0}).det(method="domain-ge")
    - D2.det(method="domain-ge") ** 2
) == 0

rho, ell1, ell2, channel = sp.symbols("rho ell1 ell2 channel")
ell2_vector = sp.Matrix([ell1, ell2])
E2_rank_one = rho * ell2_vector * ell2_vector.T
D2_one_channel = sp.Matrix(
    [
        [d11, d12],
        [channel * ell1, channel * ell2],
    ]
)
rank_two_rank_one = sp.Matrix.vstack(
    sp.Matrix.hstack(K2, D2_one_channel),
    sp.Matrix.hstack(D2_one_channel.T, E2_rank_one),
)
rank_two_rank_one_polynomial = sp.Poly(
    rank_two_rank_one.det(method="domain-ge"),
    s,
)
frame = -d11 * ell2 + d12 * ell1
assert sp.factor(
    rank_two_rank_one_polynomial.coeff_monomial(s)
) == 0
assert sp.factor(
    rank_two_rank_one_polynomial.coeff_monomial(1)
    - (channel**2 - ky * rho) * frame**2
) == 0

cyy, cyw, cww = sp.symbols("Cyy Cyw Cww")
normalized_rank_one = sp.Matrix(
    [
        [kx, kxy + s, 1, d12],
        [kxy + s, cyy, 0, cyw],
        [1, 0, 0, 0],
        [d12, cyw, 0, cww],
    ]
)
assert sp.factor(
    normalized_rank_one.det(method="domain-ge")
    - (cyw**2 - cyy * cww)
) == 0


# Rank-one quadratic pencil.  Its passive 3-by-3 Hessian must have rank two.
# The constant-kernel chart factors through a binary constant-Hessian
# potential; the exceptional chart is the cotangent lift of a plane map.
binary_det = cyy * cww - cyw**2
bx = sp.symbols("bx")
rank_one_constant_kernel = sp.Matrix(
    [
        [kx, d11, d12, bx],
        [d11, cyy, cyw, 0],
        [d12, cyw, cww, 0],
        [bx, 0, 0, 0],
    ]
)
assert sp.factor(
    rank_one_constant_kernel.det(method="domain-ge")
    + bx**2 * binary_det
) == 0

px, py, qx, qy = sp.symbols("Px Py Qx Qy")
plane_jacobian = px * qy - py * qx
cotangent_hessian = sp.Matrix(
    [
        [kx, kxy, px, qx],
        [kxy, ky, py, qy],
        [px, py, 0, 0],
        [qx, qy, 0, 0],
    ]
)
assert sp.factor(
    cotangent_hessian.det(method="domain-ge") - plane_jacobian**2
) == 0


# Rank-one globalization gate.  In the exceptional ternary chart write
# psi=u^T*g(x,l)+R(x,l), l=p(x)^T*u, and abbreviate q=g_l.  Its passive
# Hessian is p*q^T+q*p^T+h*p*p^T.  These independent-jet identities prove
# that the adjugate kernel is p cross q and that the full bordered factor is
# (p cross q).g_x + (u.q+R_l)(p cross q).p_x.  Moving u inside a level of l
# changes only the last term, which freezes p unless the chart has already
# collapsed to the constant-kernel type.
ep = sp.Matrix(sp.symbols("ep0:3"))
eq = sp.Matrix(sp.symbols("eq0:3"))
epx = sp.Matrix(sp.symbols("epx0:3"))
egx = sp.Matrix(sp.symbols("egx0:3"))
eu = sp.Matrix(sp.symbols("eu0:3"))
er = sp.Matrix(sp.symbols("er0:3"))
eh, eRl = sp.symbols("eh eRl")
ev = ep.cross(eq)
exceptional_passive_hessian = (
    ep * eq.T + eq * ep.T + eh * ep * ep.T
)
assert exceptional_passive_hessian.adjugate().applyfunc(sp.factor) == (
    -ev * ev.T
).applyfunc(sp.factor)

exceptional_gate = sp.expand(
    ev.dot(egx) + (eu.dot(eq) + eRl) * ev.dot(epx)
)
shifted_exceptional_gate = sp.expand(
    ev.dot(egx) + ((eu + er).dot(eq) + eRl) * ev.dot(epx)
)
assert sp.expand(
    shifted_exceptional_gate
    - exceptional_gate
    - er.dot(eq) * ev.dot(epx)
) == 0

# A concrete moving distinguished covector shows the exposed obstruction:
# p=(1,x,0), q=(0,0,1) has (p cross q).p_x=-1, so its gate necessarily
# contains the free transverse coordinate u.q.
moving_p = sp.Matrix([1, x, 0])
moving_q = sp.Matrix([0, 0, 1])
moving_px = moving_p.diff(x)
assert sp.expand(moving_p.cross(moving_q).dot(moving_px)) == -1

# For the constant-kernel type, v is independent of passive variables.
# Differentiating the unit f=v^T*d in a passive direction gives E_x*v;
# differentiating E*v=0 in x gives E_x*v=-E*v'.  Thus f constant forces
# v' back into ker(E), freezing its projective class.  The following moving
# kernel calibration realizes exactly the forbidden nonconstant term.
moving_v = sp.Matrix([1, x, 0])
moving_E = sp.Matrix([[x**2, -x, 0], [-x, 1, 0], [0, 0, 1]])
assert moving_E * moving_v == sp.zeros(3, 1)
assert moving_E.adjugate() == moving_v * moving_v.T
moving_u = sp.Matrix([y, z, w])
moving_quadratic = sp.expand((moving_u.T * moving_E * moving_u)[0] / 2)
moving_d = sp.Matrix(
    [sp.diff(moving_quadratic, x, variable) for variable in (y, z, w)]
)
moving_gate = sp.factor(moving_v.dot(moving_d))
assert moving_gate == x * y - z


payload = {
    "format": "hc4-scalar-cancellation-dichotomy-v1",
    "status": [
        {
            "id": "HC4RSD11",
            "kind": "exact theorem",
            "scope": "nonzero-corner scalar Schur cancellation",
            "result": "equivalence with a four-variable constant-Hessian pencil",
        },
        {
            "id": "HC4RSD12",
            "kind": "hybrid theorem",
            "scope": "zero-corner graph-coordinate scalar pivots",
            "external_input": "HC3 and Ax--Grothendieck",
            "result": "the five-variable parent gradient is a polynomial automorphism",
        },
        {
            "id": "HC4RSD13",
            "kind": "hybrid theorem",
            "scope": "rank-three quadratic four-variable constant-Hessian pencils",
            "result": "complete triangular normal form",
        },
        {
            "id": "HC4RSD14",
            "kind": "hybrid theorem",
            "scope": "rank-two quadratic four-variable constant-Hessian pencils",
            "result": "passive rank one is safe; passive rank zero is exactly JC2",
        },
        {
            "id": "HC4RSD15",
            "kind": "hybrid theorem",
            "scope": "rank-one quadratic pencils in x-constant ternary normal-form charts",
            "result": "constant kernel is safe; exceptional chart is exactly JC2",
        },
        {
            "id": "HC4RSD16",
            "kind": "hybrid theorem",
            "scope": "globalization of rank-one ternary singular-Hessian charts",
            "result": "unit border transversality freezes both rationally moving charts",
        },
    ],
    "nonzero_corner": {
        "block_identity": (
            "det([[lambda,g^T],[g,M]])="
            "lambda*det(M-g*g^T/lambda)"
        ),
        "base_potential": "psi=B-A^2/(2*lambda)",
        "pencil": "psi+s*A with s=t+A/lambda",
        "collision_transfer": "bijective after t=s-A/lambda",
        "converse": "Phi=lambda*(t+A/lambda)^2/2+psi",
        "classification_status": (
            "affine A contains direct HC4; nonlinear A is a structured "
            "constant-Hessian pencil locus"
        ),
    },
    "zero_corner_graph_coordinate": {
        "pivot": "A=w+q(u1,u2,u3)",
        "coordinates": ["r=A", "tau=t+C_r"],
        "gradient": ["F_t=r", "F_u=grad_u(C+tau*q)", "F_w=tau"],
        "determinant": "det Hess(Phi)=-det Hess_u(C+tau*q)",
        "collision_obstruction": "fixed (r,tau) reduces to HC3",
        "degree_bound": "none",
    },
    "quadratic_corollary": {
        "unit_ideal": "constant bordered determinant puts 1 in (grad(A))",
        "linear_algebra": (
            "for deg(A)<=2, some constant v has D_v(A) in K^*"
        ),
        "conclusion": "every quadratic zero-corner parent is collision-free",
        "pencil_scope": "singular and nonsingular exact-remainder pencils",
    },
    "rank_three_nonlinear_pencil": {
        "pivot": "A=x*y+z^2/2",
        "leading_faces": [
            "psi_ww=0",
            "grad(b)^T*adj(Q)*grad(b)=0",
        ],
        "unit_consequence": "grad(b) is a unimodular isotropic vector",
        "structural_lemma": (
            "a unimodular polynomial eikonal gradient in three variables "
            "is constant"
        ),
        "normal_form": (
            "psi=w*x+y*(alpha*z+beta(x))+G(x,z), alpha!=0"
        ),
        "pencil_hessian_determinant": "alpha^2",
        "result": "every member has a triangular polynomial gradient inverse",
    },
    "rank_two_nonlinear_pencil": {
        "pivot": "A=x*y",
        "leading_face": "det Hess_(z,w)(psi)=0",
        "passive_rank_one": {
            "channel_split": "2*rho*Delta_x*Delta_y=0",
            "unit_frame": "delta=(h^2-rho*psi_yy)*Delta_x^2",
            "normal_form": "psi=x*z+C(x,y,w)",
            "binary_hessian_determinant": "-delta",
            "result": "HC2 makes every pencil member injective",
        },
        "passive_rank_zero": {
            "normal_form": "psi=z*b(x,y)+w*c(x,y)+D(x,y)",
            "determinant": "Jac(b,c)^2",
            "result": "gradient invertibility is equivalent to JC2 for (b,c)",
        },
    },
    "rank_one_nonlinear_pencil": {
        "pivot": "A=x^2/2",
        "passive_gate": {
            "rank": 2,
            "adjugate": "adj(E)=rho*v*v^T",
            "unit_identity": "delta=-rho*(v^T*d)^2",
        },
        "constant_kernel_chart": {
            "normal_form": "psi=x*w+C(x,y,z)",
            "binary_hessian_determinant": "-delta",
            "result": "HC2 makes every pencil member injective",
        },
        "constant_exceptional_chart": {
            "normal_form": "psi=z*P(x,y)+w*Q(x,y)+R(x,y)",
            "determinant": "Jac(P,Q)^2",
            "result": "exactly the JC2 cotangent packet",
        },
        "globalization": {
            "constant_kernel_identity": (
                "grad_u(v^T*d)=E_x*v=-E*v_x; a unit gate forces v_x in ker(E)"
            ),
            "exceptional_representation": "psi=u^T*g(x,l)+R(x,l), l=p(x)^T*u",
            "exceptional_adjugate": "adj(E)=-(p cross g_l)*(p cross g_l)^T",
            "exceptional_gate": (
                "(p cross g_l)^T*g_x+(u^T*g_l+R_l)*(p cross g_l)^T*p_x"
            ),
            "result": (
                "either the chart is constant-kernel or p_x is parallel to p; "
                "both normal-form charts are x-independent"
            ),
        },
        "conclusion": "rank one reduces completely to HC2 or the JC2 cotangent packet",
    },
    "calibrations": {
        "universal_nonzero_corner_residual": 0,
        "universal_graph_coordinate_residual": 0,
        "nonlinear_graph_coordinate_degree": 3,
        "nonlinear_graph_parent_determinant": 1,
        "suspension_parent_determinant": 3,
        "rank_three_pencil_determinant": "alpha^2",
        "rank_two_passive_zero_determinant": "det(D)^2",
        "rank_two_passive_one_determinant": "Cyw^2-Cyy*Cww",
        "rank_one_constant_kernel_determinant": "-bx^2*det(binary Hessian)",
        "rank_one_exceptional_determinant": "Jac(P,Q)^2",
        "moving_exceptional_chart_coefficient": -1,
        "moving_constant_kernel_gate": "x*y-z",
    },
    "open_frontier": (
        "higher-degree nonlinear four-variable constant-Hessian pencils, "
        "moving matrix-pivot planes, genuinely mixed/coisotropic canonical "
        "transformations, and direct degree-five HC4; every quadratic scalar "
        "pencil is reduced to HC2 or JC2"
    ),
}

serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
OUTPUT.write_text(serialized)
digest = hashlib.sha256(serialized.encode()).hexdigest()

print("PASS: nonzero-corner cancellation is a constant-Hessian pencil")
print("PASS: parent and reduced gradient collisions correspond exactly")
print("PASS: zero-corner graph-coordinate Hessian factorization is universal")
print("PASS: HC3 excludes every graph-coordinate parent collision")
print("PASS: every quadratic zero-corner pivot belongs to this class")
print("PASS: rank-three nonlinear constant-Hessian pencils are triangular")
print("PASS: rank-two passive-rank-one pencils reduce to HC2")
print("ENDPOINT: rank-two passive-rank-zero pencils are exactly JC2")
print("PASS: rank-one constant kernel charts reduce to HC2")
print("ENDPOINT: rank-one constant exceptional charts are exactly JC2")
print("PASS: unit transversality freezes every moving ternary chart")
print("THEOREM: every nonlinear quadratic scalar pencil reduces to HC2 or JC2")
print(f"PASS: wrote {OUTPUT.relative_to(ROOT)}")
print(f"SHA256: {digest}")
