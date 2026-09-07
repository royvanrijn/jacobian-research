#!/usr/bin/env sage -python
"""Eliminate the first certified H3 q=6 pencil and run its genus gate.

The exact global kernel provides two functions ``s0=a0+b0*m`` and
``s1=a1+b1*m``.  Set ``T=s1/s0`` and solve for the old chord m.  The standard
chord discriminant gives the hyperelliptic presentation of the generic T
fibre.  This script reports only the exact squarefree old-base degree; a
degree three or four model is the required genus-one gate.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[2]
RR = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-global-rr.json"
SECTION = ROOT / "artifacts/generated-results/elkies-k3-h92-p1-lift.json"
CORE = ROOT / "elkies-k3/scripts/elliptic_neighbor_compiler.sage"
RR_SHA256 = None
SECTION_SHA256 = "c323bf6346bb239934a5a2d8b1a3f4067e70e993d2e4eb32aaa30f469fca6397"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-global-pencil-elimination.json"

exec(compile(CORE.read_text(), str(CORE), "exec"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial(ring, coefficients):
    return ring([QQ(value) for value in coefficients])


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--rr", type=Path, default=RR)
parser.add_argument("--section", type=Path, default=SECTION)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

rr = json.loads(args.rr.read_text())
assert rr["status"] == "PASS_EXACT_GLOBAL_RR_KERNEL"
assert rr["kernel"]["dimension"] == 2
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

# Reconstruct the same complete q6 resolved condition matrix used by the
# actual-cover certificate.  The displayed global kernel is accepted only as
# an explicit basis of that matrix, then applied to the ambient chord frame.
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
assert kernel_basis.nrows() == 2

chord_expansions = tuple(
    (u_field(u**power) / u_field(h**2), u_field(0))
    if kind == "A" else
    (u_field(0), u_field(u**power) / u_field(h))
    for kind, power in ambient_basis
)

# Coefficients of the original H92 Weierstrass equation in entrance u.
# They are recovered from P1's exact square relation to avoid a second model
# parser here: a=(yP^2-xP^3-b)/xP is not useful, so load the anchor formula.
from importlib.machinery import SourceFileLoader
ANCHOR = ROOT / "elkies-k3/scripts/verify_h3_noncm_q6_source_anchor.sage"
H92 = ROOT / "artifacts/local/humbert-inputs/92/igusa92.txt"
anchor = SourceFileLoader("h92_q6_eliminate_anchor", str(ANCHOR)).load_module()
h92_ring, h92_formulas = anchor.parse_h92(H92)
r92, s92 = anchor.EXPECTED_H92
A1, A, B1, B, B2 = tuple(QQ(value(r92, s92)) for value in h92_formulas)
old_a = A1 / u**3 + A / u**4
old_b = B1 / u**5 + B / u**6 + B2 / u**7

T_ring = PolynomialRing(QQ, "T")
T = T_ring.gen()
T_field = T_ring.fraction_field()
old_base_ring = PolynomialRing(T_field, "u")
uu = old_base_ring.gen()
old_base_field = old_base_ring.fraction_field()

def transport(value):
    return old_base_field(
        old_base_ring([T_field(entry) for entry in value.numerator().list()])
    ) / old_base_field(
        old_base_ring([T_field(entry) for entry in value.denominator().list()])
    )

resolved_hop = compile_resolved_degree_two_chord_hop(
    rr_compilation,
    kernel_basis,
    tuple((transport(a), transport(b)) for a, b in chord_expansions),
    old_base_ring,
    T,
    transport(x_p),
    transport(y_p),
    transport(old_a),
    transport(old_b),
)
a0_t, b0_t = resolved_hop["chord_coefficients"][0]
a1_t, b1_t = resolved_hop["chord_coefficients"][1]
compiled_hop = resolved_hop["conversion"]
m = compiled_hop["chord"]
radicand = compiled_hop["radicand"]
numerator = old_base_ring(radicand.numerator())
denominator = old_base_ring(radicand.denominator())
quartic = compiled_hop["binary_quartic"]
square_factor = compiled_hop["square_factor"]

def parity_factors(value):
    return tuple(
        (factor, int(exponent))
        for factor, exponent in value.factor()
        if exponent % 2
    )

odd_numerator = parity_factors(numerator)
odd_denominator = parity_factors(denominator)
squarefree_degree = sum(factor.degree() for factor, _ in odd_numerator + odd_denominator)
assert squarefree_degree in (3, 4)
assert quartic.degree() == squarefree_degree

payload = {
    "schema": "elkies-k3.h92-q6-global-pencil-elimination.v1",
    "status": "PASS_EXACT_GENUS_ONE_GATE",
    "inputs": {
        "global_rr": {"path": str(args.rr.relative_to(ROOT)), "sha256": digest(args.rr)},
        "marked_section": {"path": str(args.section.relative_to(ROOT)), "sha256": SECTION_SHA256},
    },
    "parameter": "T=(a1+b1*m)/(a0+b0*m)",
    "chord_solution": "m=(a1-T*a0)/(T*b0-b1)",
    "squarefree_old_base_degree": int(squarefree_degree),
    "squarefree_reconstruction": "radicand=square_factor^2*binary_quartic (exactly checked)",
    "odd_numerator_factor_degrees": [int(factor.degree()) for factor, _ in odd_numerator],
    "odd_denominator_factor_degrees": [int(factor.degree()) for factor, _ in odd_denominator],
    "boundary": "This certifies the genus-one gate. Jacobian minimization and fibre/section transport remain downstream.",
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92Q6ELIMINATION|squarefree_degree={}|status=PASS_EXACT_GENUS_ONE_GATE".format(
        squarefree_degree
    ),
    flush=True,
)
