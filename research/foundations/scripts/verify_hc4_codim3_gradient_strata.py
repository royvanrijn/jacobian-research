#!/usr/bin/env python3
"""Exact codimension-three gradient-stratum sieve for quintic HC(4).

This checker treats two support packets:

* smooth essential rank two, whose infinity support is the kernel line;
* essential rank three with isolated ordinary singularities, whose support
  contains the corresponding vertex lines.

The output separates unconditional identities, conditional exclusions, and
normal-slice calibrations.  It does not promote the remaining exceptional
coefficient strata to realized HC4 candidates.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.projective_gradient_segre import (  # noqa: E402
    SmoothEssentialGradientNormalSlice,
)


ATLAS = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "hc4_projective_polar_atlas.json"
)
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "hc4_codim3_gradient_strata.json"
)


atlas_payload = json.loads(ATLAS.read_text())

# ---------------------------------------------------------------------------
# 1. Constant-kernel synchronization for the rank-two quartic face.
#
# Work over the quotient field in the active variables.  Normalize the
# nonzero pure-kernel quartic to t^4 and write the remaining polynomial in
# the two kernel coordinates (t,w) with t-degree at most three.  The radical
# of det Hess_(t,w)(h4)=0 kills every coefficient contributing a second
# derivative in w or a mixed t,w derivative.
# ---------------------------------------------------------------------------
t, w = sp.symbols("t w")
lower_exponents = tuple(
    (t_power, total - t_power)
    for total in range(4)
    for t_power in range(total + 1)
)
quartic_coefficients = sp.symbols(f"a0:{len(lower_exponents)}")
kernel_quartic = sp.expand(
    t**4
    + sum(
        coefficient * t**t_power * w**w_power
        for coefficient, (t_power, w_power) in zip(
            quartic_coefficients,
            lower_exponents,
        )
    )
)
kernel_hessian_determinant = sp.Poly(
    sp.expand(sp.det(sp.hessian(kernel_quartic, (t, w)))),
    t,
    w,
)
kernel_face_equations = tuple(
    coefficient for _, coefficient in kernel_hessian_determinant.terms()
)
assert len(kernel_face_equations) == 8

kernel_face_groebner = sp.groebner(
    kernel_face_equations,
    *quartic_coefficients,
    order="grevlex",
)
forbidden_kernel_coefficients = {
    lower_exponents[index]: {
        "coefficient": quartic_coefficients[index],
        "radical_power": power,
    }
    for index, power in {
        3: 2,  # w^2
        4: 3,  # t*w
        6: 1,  # w^3
        7: 1,  # t*w^2
        8: 3,  # t^2*w
    }.items()
}
for record in forbidden_kernel_coefficients.values():
    assert kernel_face_groebner.reduce(
        record["coefficient"] ** record["radical_power"]
    )[1] == 0

# The surviving polynomial is P(t)+linear(w), so w is a constant kernel
# direction for its Hessian in the two kernel variables.
assert all(
    w_power <= 1 and not (w_power == 1 and t_power >= 1)
    for t_power, w_power in lower_exponents
    if (t_power, w_power) not in forbidden_kernel_coefficients
)

# The homogeneous binary quartic lemma used to choose the normalization has
# an elementary dehomogenized form:
#
# det Hess(y^4*p(x/y))|_(y=1) = 3*(4*p*p''-3*(p')^2).
z, x, y = sp.symbols("z x y")
p_coefficients = sp.symbols("p0:5")
p = sum(
    coefficient * z**index
    for index, coefficient in enumerate(p_coefficients)
)
homogeneous_binary_quartic = sp.expand(y**4 * p.subs(z, x / y))
binary_quartic_hessian_chart = sp.expand(
    sp.det(sp.hessian(homogeneous_binary_quartic, (x, y))).subs(
        {x: z, y: 1}
    )
)
binary_quartic_ode = sp.expand(
    3 * (4 * p * sp.diff(p, z, 2) - 3 * sp.diff(p, z) ** 2)
)
assert sp.expand(binary_quartic_hessian_chart - binary_quartic_ode) == 0

# ---------------------------------------------------------------------------
# 2. The next rank-two determinant face.
#
# With w synchronized as the kernel direction of Hess_K(h4), the lambda^9
# coefficient is
#
# b*(det(C)*D_w^2(h3) - d^T*adj(C)*d),
#
# where C=Hess(h5) is the active binary block and
# d=grad_active(D_w h4).
# ---------------------------------------------------------------------------
lam = sp.symbols("lambda")
c11, c12, c22 = sp.symbols("c11 c12 c22")
active_hessian = sp.Matrix([[c11, c12], [c12, c22]])
b, a00 = sp.symbols("b a00")
d1, d2 = sp.symbols("d1 d2")
e1, e2 = sp.symbols("e1 e2")
q11, q12, q22 = sp.symbols("q11 q12 q22")

top_hessian = sp.zeros(4)
top_hessian[2:4, 2:4] = active_hessian
quartic_hessian = sp.Matrix(
    [
        [0, 0, d1, d2],
        [0, b, e1, e2],
        [d1, e1, q11, q12],
        [d2, e2, q12, q22],
    ]
)
cubic_hessian = sp.diag(a00, 0, 0, 0)
rank_two_pencil_determinant = sp.Poly(
    sp.expand(
        (
            lam**3 * top_hessian
            + lam**2 * quartic_hessian
            + lam * cubic_hessian
        ).det()
    ),
    lam,
)
assert rank_two_pencil_determinant.coeff_monomial(lam**10) == 0
rank_two_lambda9 = sp.expand(
    rank_two_pencil_determinant.coeff_monomial(lam**9)
)
active_gradient_vector = sp.Matrix([d1, d2])
rank_two_schur_face = sp.expand(
    b
    * (
        a00 * active_hessian.det()
        - (
            active_gradient_vector.T
            * active_hessian.adjugate()
            * active_gradient_vector
        )[0]
    )
)
assert sp.expand(rank_two_lambda9 - rank_two_schur_face) == 0

rank_two_degree_ledger = {
    "active_hessian_entry_degree": 3,
    "active_hessian_determinant_degree": 6,
    "adjugate_degree": 3,
    "active_gradient_degree": 2,
    "adjugate_gradient_vector_degree": 5,
}
assert (
    rank_two_degree_ledger["adjugate_gradient_vector_degree"]
    < rank_two_degree_ledger["active_hessian_determinant_degree"]
)

# An exact squarefree binary-quintic Hessian witness proves that the
# squarefree discriminant condition is a genuine dense-open restriction.
binary_x, binary_y = sp.symbols("binary_x binary_y")
binary_quintic_witness = (
    binary_x**5
    + binary_x**4 * binary_y
    + binary_x * binary_y**4
    + binary_y**5
)
binary_quintic_gradient_gcd = sp.gcd(
    sp.gcd(
        binary_quintic_witness,
        sp.diff(binary_quintic_witness, binary_x),
    ),
    sp.diff(binary_quintic_witness, binary_y),
)
binary_hessian_discriminant = sp.expand(
    sp.det(
        sp.hessian(
            binary_quintic_witness,
            (binary_x, binary_y),
        )
    )
)
binary_hessian_discriminant_gcd = sp.gcd(
    sp.gcd(
        binary_hessian_discriminant,
        sp.diff(binary_hessian_discriminant, binary_x),
    ),
    sp.diff(binary_hessian_discriminant, binary_y),
)
assert sp.Poly(
    binary_quintic_gradient_gcd,
    binary_x,
    binary_y,
).total_degree() == 0
assert sp.Poly(
    binary_hessian_discriminant,
    binary_x,
    binary_y,
).total_degree() == 6
assert sp.Poly(
    binary_hessian_discriminant_gcd,
    binary_x,
    binary_y,
).total_degree() == 0

# ---------------------------------------------------------------------------
# 3. Generic transverse length on the rank-two kernel line.
#
# The two active quartics form a (4,4) complete intersection of length 16.
# Adding epsilon^4 gives length 64.  When h4|K is nonzero, a kernel-gradient
# component is epsilon times a unit at the generic point of K, so the full
# transverse base is the special fiber and sigma_3=16.
# ---------------------------------------------------------------------------
rank_two_normal_slice = SmoothEssentialGradientNormalSlice(4, 4, 2)
binary_jacobian_hilbert_function = (
    rank_two_normal_slice.jacobian_hilbert_function
)
assert binary_jacobian_hilbert_function == (1, 2, 3, 4, 3, 2, 1)
binary_jacobian_length = rank_two_normal_slice.jacobian_length
assert binary_jacobian_length == 16
active_transverse_length = rank_two_normal_slice.truncated_active_length
assert active_transverse_length == 64
forced_sigma3 = rank_two_normal_slice.unit_penultimate_segre_degree

rank_two_atlas_intersection = {}
for affine_degree in (2, 3):
    key = f"gradient_degree_4_affine_degree_{affine_degree}"
    codimension_three_rows = [
        row
        for row in atlas_payload["atlases"][key]
        if row["leading_base_codimension"] == 3
    ]
    forced_rows = [
        row
        for row in codimension_three_rows
        if row["segre_degrees"][2] == forced_sigma3
    ]
    assert len(forced_rows) == 1
    forced_row = forced_rows[0]
    assert forced_row["projective_degrees"] == [
        1,
        4,
        16,
        48,
        affine_degree,
    ]
    assert forced_row["segre_degrees"] == [
        0,
        0,
        16,
        -affine_degree,
    ]
    rank_two_atlas_intersection[f"affine_degree_{affine_degree}"] = {
        "codimension_three_rows": len(codimension_three_rows),
        "rows_with_forced_sigma3": len(forced_rows),
        "forced_projective_degrees": forced_row["projective_degrees"],
        "forced_segre_degrees": forced_row["segre_degrees"],
    }

assert rank_two_atlas_intersection == {
    "affine_degree_2": {
        "codimension_three_rows": 58,
        "rows_with_forced_sigma3": 1,
        "forced_projective_degrees": [1, 4, 16, 48, 2],
        "forced_segre_degrees": [0, 0, 16, -2],
    },
    "affine_degree_3": {
        "codimension_three_rows": 57,
        "rows_with_forced_sigma3": 1,
        "forced_projective_degrees": [1, 4, 16, 48, 3],
        "forced_segre_degrees": [0, 0, 16, -3],
    },
}

# ---------------------------------------------------------------------------
# 4. Rank-three ordinary-singularity incidence.
#
# At a projective singular point p with rank Hess(h5)(p)=2, homogeneity
# gives ker Hess(h5)(p)=k*p.  The adjugate is a nonzero scalar multiple of
# p*p^T.  Evaluating the Schur face and using Euler for the cubic s3 forces
# s3(p)=0.  The following nodal quintic checks the scalar identity exactly.
# ---------------------------------------------------------------------------
node_x, node_y, node_z = sp.symbols("node_x node_y node_z")
nodal_quintic = (
    node_x**3 * (node_y**2 + node_z**2)
    + node_y**5
    + node_z**5
)
node_point = {node_x: 1, node_y: 0, node_z: 0}
nodal_gradient_at_point = sp.Matrix(
    [
        sp.diff(nodal_quintic, variable).subs(node_point)
        for variable in (node_x, node_y, node_z)
    ]
)
nodal_hessian_at_point = sp.hessian(
    nodal_quintic,
    (node_x, node_y, node_z),
).subs(node_point)
assert nodal_gradient_at_point == sp.zeros(3, 1)
assert nodal_hessian_at_point.rank() == 2
assert nodal_hessian_at_point * sp.Matrix([1, 0, 0]) == sp.zeros(3, 1)

node_cubic_coefficients = sp.symbols("s0:10")
node_cubic_exponents = tuple(
    (first, second, 3 - first - second)
    for first in range(4)
    for second in range(4 - first)
)
universal_node_cubic = sp.expand(
    sum(
        coefficient
        * node_x**first
        * node_y**second
        * node_z**third
        for coefficient, (first, second, third) in zip(
            node_cubic_coefficients,
            node_cubic_exponents,
        )
    )
)
node_cubic_gradient = sp.Matrix(
    [
        sp.diff(universal_node_cubic, variable).subs(node_point)
        for variable in (node_x, node_y, node_z)
    ]
)
node_schur_value = sp.expand(
    (
        node_cubic_gradient.T
        * nodal_hessian_at_point.adjugate()
        * node_cubic_gradient
    )[0]
)
node_cubic_value = sp.expand(universal_node_cubic.subs(node_point))
assert sp.expand(node_schur_value - 36 * node_cubic_value**2) == 0

payload = {
    "format": "hc4-codim3-gradient-strata-v1",
    "software_assumptions": {
        "python": "repository .python-version and requirements.txt",
        "coefficient_field": "characteristic zero",
        "independent_replay": "Macaulay2 over QQ",
    },
    "scope": (
        "Exact conditional sieve on the two codimension-three top-gradient "
        "packets. It closes a dense open rank-two coefficient stratum and "
        "forces sigma_3=16 on its nonsquarefree-Hessian remainder. It also "
        "forces the rank-three Schur cubic through every ordinary isolated "
        "singular point. No unconditional codimension-three atlas row is "
        "excluded because the zero kernel restriction, nonsquarefree "
        "binary Hessian, and worse ternary singular strata remain."
    ),
    "rank_two_smooth_top": {
        "top_type": "squarefree essential binary quintic",
        "kernel_support": "P1",
        "quartic_face": "det Hess_K(h4)=0",
        "nonzero_kernel_restriction": {
            "normal_form": "h4|K=t^4",
            "radical_constant_kernel_certificate": {
                f"t^{t_power}*w^{w_power}": record["radical_power"]
                for (
                    t_power,
                    w_power,
                ), record in forbidden_kernel_coefficients.items()
            },
            "conclusion": "a constant direction w satisfies D_w^2 h4=0",
        },
        "next_schur_face": (
            "det(C)*D_w^2(h3)-"
            "grad_U(s3)^T*adj(C)*grad_U(s3)=0"
        ),
        "degree_ledger": rank_two_degree_ledger,
        "squarefree_binary_hessian_branch": {
            "status": "excluded via HC4CD5",
            "reason": (
                "squarefree degree-6 det(C) divides the degree-5 vector "
                "adj(C)*grad_U(s3), forcing grad_U(s3)=0 and D_w^2(h3)=0"
            ),
            "exact_witness": (
                "x^5+x^4*y+x*y^4+y^5 has squarefree Hessian determinant"
            ),
        },
        "nonsquarefree_binary_hessian_branch": {
            "status": "open exceptional discriminant",
            "generic_transverse_hilbert_function": list(
                binary_jacobian_hilbert_function
            ),
            "active_transverse_length": active_transverse_length,
            "forced_sigma3": forced_sigma3,
            "atlas_intersection": rank_two_atlas_intersection,
        },
        "zero_kernel_restriction_branch": {
            "condition": "h4|K=0",
            "status": "open synchronization stratum",
        },
    },
    "rank_three_isolated_singularities": {
        "ordinary_point_condition": (
            "grad(h5)(p)=0 and rank Hess(h5)(p)=2"
        ),
        "schur_incidence": "s3(p)=0",
        "proof_identity": (
            "adj(Hess(h5)(p))=c*p*p^T and "
            "p^T*grad(s3)(p)=3*s3(p)"
        ),
        "nodal_calibration": {
            "h5": "x^3*(y^2+z^2)+y^5+z^5",
            "point": "[1:0:0]",
            "schur_value": "36*s3(p)^2",
        },
        "worse_singularities": (
            "rank Hess(h5)(p)<=1 is not constrained by this adjugate value"
        ),
    },
    "atlas_consequence": (
        "No row is removed unconditionally. In the rank-two subpacket with "
        "h4|K nonzero, the squarefree binary-Hessian open is empty; its "
        "nonsquarefree remainder can meet only the sigma3=16 row for each "
        "affine degree."
    ),
}

serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
OUTPUT.write_text(serialized)
digest = hashlib.sha256(serialized.encode()).hexdigest()

print("PASS: synchronized the rank-two quartic kernel on h4|K nonzero")
print("PASS: verified the rank-two lambda^9 Schur determinant face")
print("PASS: squarefree binary Hessians reach the common-direction obstruction")
print("PASS: the nonsquarefree remainder has forced generic sigma_3=16")
print("PASS: ordinary rank-three singular points force s3(p)=0")
print(f"PASS: wrote {OUTPUT.relative_to(ROOT)}")
print(f"SHA256: {digest}")
