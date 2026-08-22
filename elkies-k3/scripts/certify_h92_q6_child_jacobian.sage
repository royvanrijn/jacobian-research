#!/usr/bin/env sage -python
"""Compute the Jacobian fibre signature of the first exact H3 q=6 pencil."""

import argparse
import hashlib
import json
from importlib.machinery import SourceFileLoader
from pathlib import Path

from sage.all import PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[2]
RR = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-global-rr.json"
RR_COVER = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-actual-resolved-rr-cover.json"
ELIMINATION = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-global-pencil-elimination.json"
SECTION = ROOT / "artifacts/generated-results/elkies-k3-h92-p1-lift.json"
CORE = ROOT / "elkies-k3/scripts/elliptic_neighbor_compiler.sage"
SECTION_SHA256 = "c323bf6346bb239934a5a2d8b1a3f4067e70e993d2e4eb32aaa30f469fca6397"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-jacobian.json"

exec(compile(CORE.read_text(), str(CORE), "exec"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial(ring, coefficients):
    return ring([QQ(value) for value in coefficients])


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--rr", type=Path, default=RR)
parser.add_argument("--rr-cover", type=Path, default=RR_COVER)
parser.add_argument("--elimination", type=Path, default=ELIMINATION)
parser.add_argument("--section", type=Path, default=SECTION)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

rr = json.loads(args.rr.read_text())
rr_cover = json.loads(args.rr_cover.read_text())
elimination = json.loads(args.elimination.read_text())
assert rr["status"] == "PASS_EXACT_GLOBAL_RR_KERNEL"
assert rr_cover["status"] == "PASS_EXACT_Q6_ACTUAL_RESOLVED_RR_COVER"
assert rr_cover["inputs"]["global_rr"]["sha256"] == digest(args.rr)
assert rr_cover["vertical_condition_matrix"]["h0_D"] == 2
assert elimination["status"] == "PASS_EXACT_GENUS_ONE_GATE"
assert digest(args.section) == SECTION_SHA256
section = json.loads(args.section.read_text())

u_ring = PolynomialRing(QQ, "u")
u = u_ring.gen()
u_field = u_ring.fraction_field()
h = polynomial(u_ring, section["structured_denominator"]["Z4_coefficients"])
x_p = u_field(polynomial(
    u_ring, section["x_entrance_base"]["numerator_coefficients"]
)) / u_field(polynomial(
    u_ring, section["x_entrance_base"]["denominator_coefficients"]
))
y_p = u_field(polynomial(
    u_ring, section["y_entrance_base"]["numerator_coefficients"]
)) / u_field(polynomial(
    u_ring, section["y_entrance_base"]["denominator_coefficients"]
))

ANCHOR = ROOT / "elkies-k3/scripts/verify_h3_noncm_q6_source_anchor.sage"
H92 = ROOT / "artifacts/local/humbert-inputs/92/igusa92.txt"
anchor = SourceFileLoader("h92_q6_child_anchor", str(ANCHOR)).load_module()
h92_ring, h92_formulas = anchor.parse_h92(H92)
r92, s92 = anchor.EXPECTED_H92
A1, A, B1, B, B2 = tuple(QQ(value(r92, s92)) for value in h92_formulas)
old_a = A1 / u**3 + A / u**4
old_b = B1 / u**5 + B / u**6 + B2 / u**7

T_ring = PolynomialRing(QQ, "T")
T = T_ring.gen()
T_field = T_ring.fraction_field()
u_over_T = PolynomialRing(T_field, "u")
uu = u_over_T.gen()
uT_field = u_over_T.fraction_field()

def transport(value):
    return uT_field(u_over_T([T_field(entry) for entry in value.numerator().list()])) / uT_field(
        u_over_T([T_field(entry) for entry in value.denominator().list()])
    )

# Reconstruct the same complete q6 RR matrix used by the exact resolved-cover
# certificate.  The child equation is compiled from its certified kernel, not
# from an independently substituted pair of global sections.
ambient_basis = tuple(tuple(entry) for entry in rr["ambient"]["basis"])
collision_matrix = matrix(
    QQ, [[QQ(value) for value in row] for row in rr["collision_condition"]["matrix"]]
)
collision_block = quotient_condition(
    "smooth_P1_O_collision",
    ambient_basis,
    lambda column: tuple(
        collision_matrix[row, ambient_basis.index(column)]
        for row in range(collision_matrix.nrows())
    ),
    tuple("collision_{}".format(row) for row in range(collision_matrix.nrows())),
    "pinned q6 smooth collision quotient",
)
rr_compilation = compile_resolved_conditions(
    ambient_basis, (collision_block,), complete=True, compute_kernel=False
)
assert rr_compilation["condition_matrix"] == collision_matrix
assert rr_compilation["h0_certified"]
kernel_basis = matrix(QQ, [
    [QQ(value) for value in row] for row in rr["kernel"]["basis_matrix"]
])
chord_expansions = tuple(
    (u_field(u**power) / u_field(h**2), u_field(0))
    if kind == "A" else
    (u_field(0), u_field(u**power) / u_field(h))
    for kind, power in ambient_basis
)
child_compilation = compile_resolved_degree_two_child_jacobian(
    rr_compilation,
    kernel_basis,
    tuple((transport(a), transport(b)) for a, b in chord_expansions),
    u_over_T,
    T_ring,
    T,
    transport(x_p),
    transport(y_p),
    transport(old_a),
    old_b=transport(old_b),
)
conversion = child_compilation["resolved_hop"]["conversion"]
quartic = conversion["binary_quartic"]
assert quartic.degree() == 4
jacobian_A = child_compilation["jacobian_a"]
jacobian_B = child_compilation["jacobian_b"]
classification = child_compilation["finite_classification"]
finite_data = [{
    "factor": str(item["factor"]), "degree": item["degree"],
    "raw_orders": list(item["raw_orders"]), "scaling": item["scaling"],
    "minimal_orders": list(item["minimal_orders"]), "kodaira": item["kodaira"],
} for item in classification["finite_fibres"]]
root_rank = classification["finite_root_rank"]
root_euler = classification["finite_euler_number"]
root_determinant = classification["finite_root_determinant"]
infinity = classification["infinity_boundary"]
raw_infinity = infinity["raw_orders"]
infinity_scaling = infinity["scaling"]
infinity_orders = infinity["normalized_orders"]
infinity_kind = "smooth"
if infinity_orders[2] > 0:
    rank, euler, determinant, infinity_kind = kodaira_data_from_short_orders(
        *infinity_orders
    )
    root_rank += rank
    root_euler += euler
    root_determinant *= determinant
print(
    "H92Q6CHILD|stage=classification|finite={}|infinity={}|root_rank={}|root_euler={}|root_det={}".format(
        [(item["degree"], item["minimal_orders"], item["kodaira"]) for item in finite_data],
        (infinity_orders, infinity_kind), root_rank, root_euler, root_determinant
    ),
    flush=True,
)
assert root_euler == 24
assert (root_rank, root_determinant) == (14, 3)
minimal_A = classification["finite_minimization"]["minimal_a"]
minimal_B = classification["finite_minimization"]["minimal_b"]
minimal_delta = classification["finite_minimization"]["minimal_discriminant"]
assert minimal_A.degree() <= 8 and minimal_B.degree() <= 12 and minimal_delta.degree() <= 24

payload = {
    "schema": "elkies-k3.h92-q6-child-jacobian.v1",
    "status": "PASS_EXACT_E8_E6_CHILD_JACOBIAN",
    "inputs": {
        "global_rr": {"path": str(args.rr.relative_to(ROOT)), "sha256": digest(args.rr)},
        "actual_rr_cover": {"path": str(args.rr_cover.relative_to(ROOT)), "sha256": digest(args.rr_cover)},
        "elimination": {"path": str(args.elimination.relative_to(ROOT)), "sha256": digest(args.elimination)},
        "marked_section": {"path": str(args.section.relative_to(ROOT)), "sha256": SECTION_SHA256},
    },
    "binary_quartic_degree": int(quartic.degree()),
    "minimal_short_weierstrass": {
        "equation": "Y^2=X^3+A(T)*X+B(T)",
        "A_coefficients_low_to_high": [str(value) for value in minimal_A.list()],
        "B_coefficients_low_to_high": [str(value) for value in minimal_B.list()],
        "Delta_coefficients_low_to_high": [str(value) for value in minimal_delta.list()],
    },
    "finite_fibres": finite_data,
    "infinity": {"raw_orders": [int(value) for value in raw_infinity], "scaling": int(infinity_scaling), "minimal_orders": [int(value) for value in infinity_orders], "kodaira": infinity_kind},
    "root_data": {"rank": int(root_rank), "determinant": int(root_determinant), "type": "E8+E6"},
    "mordell_weil_rank_if_rho_19": 3,
    "boundary": "This is a finite-place Jacobian root-signature certificate. It does not assign section components; the transported-section and Shioda-Gram certificates are separate resolved-chart evidence.",
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92Q6CHILD|roots=E8+E6|root_rank=14|root_det=3|MW_rank=3|"
    "status=PASS_EXACT_E8_E6_CHILD_JACOBIAN",
    flush=True,
)
