#!/usr/bin/env sage -python
"""Assemble the exact H3 q=6 global coefficient ambient and collision matrix.

This is the common ``D=O+(-P1)-F_infinity`` representative.  Put
``m=(y-y_P)/(x-x_P)`` and write a prospective section as ``a+b*m``.

The resolved local modules give:

* E8: ``a`` has u-order at least one and ``b`` u-order at least three;
* E7: in the ``<1,m>`` frame both coefficients are regular at u=infinity;
* at the four smooth P.O collisions, with p=y_P/x_P and collision polynomial
  h, the saturated local frame is ``<1,(m-p)/h>``.

Thus a global coefficient pair has the bounded form

    a=A/h^2,  b=B/h,
    u | A, deg(A)<=8,  u^3 | B, deg(B)<=4,

and its sole finite local condition is

    A*(den(p)/h) + B*num(p) == 0 (mod h^2).

The script builds that exact 8-by-10 condition matrix.  A two-dimensional
kernel is a certificate for h0(D)=2 only if the stated local-module cover is
complete; this script records that provenance explicitly.
"""

import argparse
import hashlib
import json
from importlib.machinery import SourceFileLoader
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
ANCHOR = ROOT / "elkies-k3/scripts/verify_h3_noncm_q6_source_anchor.sage"
H92 = ROOT / "artifacts/local/humbert-inputs/92/igusa92.txt"
SECTION = ROOT / "artifacts/generated-results/elkies-k3-h92-p1-lift.json"
SECTION_SHA256 = "c323bf6346bb239934a5a2d8b1a3f4067e70e993d2e4eb32aaa30f469fca6397"
CORE = ROOT / "elkies-k3/scripts/elliptic_neighbor_compiler.sage"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-global-rr.json"

exec(compile(CORE.read_text(), str(CORE), "exec"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial(ring, coefficients):
    return ring([QQ(value) for value in coefficients])


def coefficients_mod(polynomial_value, modulus):
    remainder = polynomial_value % modulus
    return vector(QQ, [remainder[index] for index in range(modulus.degree())])


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--modular-only", action="store_true")
parser.add_argument("--prime", type=int, default=1009)
args = parser.parse_args()

assert digest(SECTION) == SECTION_SHA256
section = json.loads(SECTION.read_text())
assert section["status"] == "PASS_EXACT_H92_P1"
anchor = SourceFileLoader("h92_q6_global_rr_anchor", str(ANCHOR)).load_module()
h92_ring, h92_formulas = anchor.parse_h92(H92)
r92, s92 = anchor.EXPECTED_H92
A1, A, B1, B, B2 = tuple(QQ(value(r92, s92)) for value in h92_formulas)

u_ring = PolynomialRing(QQ, "u")
u = u_ring.gen()
u_field = u_ring.fraction_field()
x_numerator = polynomial(u_ring, section["x_entrance_base"]["numerator_coefficients"])
x_denominator = polynomial(u_ring, section["x_entrance_base"]["denominator_coefficients"])
y_numerator = polynomial(u_ring, section["y_entrance_base"]["numerator_coefficients"])
y_denominator = polynomial(u_ring, section["y_entrance_base"]["denominator_coefficients"])
x_p = u_field(x_numerator) / u_field(x_denominator)
y_p = u_field(y_numerator) / u_field(y_denominator)
assert y_p**2 == x_p**3 + (A1 / u**3 + A / u**4) * x_p + (
    B1 / u**5 + B / u**6 + B2 / u**7
)
h = polynomial(u_ring, section["structured_denominator"]["Z4_coefficients"])
assert h.degree() == 4 and h[0]

p = y_p / x_p
p_numerator = u_ring(p.numerator())
p_denominator = u_ring(p.denominator())
denominator_over_h, remainder = p_denominator.quo_rem(h)
assert not remainder
modulus = h**2
print("H92Q6GLOBALRR|stage=inputs", flush=True)

# E8 orders and E7 degree bounds have already been imposed in this ambient.
# The ten displayed terms are A=u*A_1 through u^8*A_8, followed by
# B=u^3*B_3+u^4*B_4.
ambient = tuple(
    [("A", exponent) for exponent in range(1, 9)]
    + [("B", exponent) for exponent in range(3, 5)]
)

def collision_evaluator(label):
    kind, exponent = label
    if kind == "A":
        return coefficients_mod(u**exponent * denominator_over_h, modulus)
    assert kind == "B"
    return coefficients_mod(u**exponent * p_numerator, modulus)

collision_block = quotient_condition(
    "four smooth P1.O collisions",
    ambient,
    collision_evaluator,
    tuple("u^{} mod h^2".format(index) for index in range(modulus.degree())),
    "a=A/h^2, b=B/h; saturated local frame <1,(m-y_P/x_P)/h>",
)
print("H92Q6GLOBALRR|stage=matrix", flush=True)
prime = int(args.prime)
finite = GF(prime)
modular_matrix = matrix(
    finite,
    [[finite(value.numerator()) / finite(value.denominator())
      for value in row]
     for row in collision_block["matrix"].rows()],
)
modular_rank = modular_matrix.rank()
modular_pivots = tuple(int(value) for value in modular_matrix.pivots())
print(
    "H92Q6GLOBALRR|stage=modular|prime={}|rank={}|pivots={}".format(
        prime, modular_rank, modular_pivots
    ),
    flush=True,
)
assert modular_rank == 8
if args.modular_only:
    raise SystemExit(0)

# The modular full-row-rank computation is an exact rank certificate: reduction
# cannot increase rank, while the matrix has only eight rows.  Avoid Sage's
# expensive generic rational nullspace on these large H92 coefficients by
# solving the congruence in QQ[u]/(h^2).  Since den(p)/h is a unit modulo h^2,
# the two B coefficients determine A uniquely.
assert denominator_over_h.gcd(modulus) in QQ
inverse_denominator = denominator_over_h.inverse_mod(modulus)
kernel_rows = []
for B_coefficient in (u**3, u**4):
    A_coefficient = (-B_coefficient * p_numerator * inverse_denominator) % modulus
    # The residue representative has degree < 8.  Add its unique multiple of
    # h^2 that kills the constant term; this is exactly the allowed u^8 A
    # coefficient and preserves the collision congruence.
    A_coefficient += -A_coefficient[0] / modulus[0] * modulus
    assert A_coefficient[0] == 0
    kernel_rows.append([
        A_coefficient[exponent] for exponent in range(1, 9)
    ] + [
        QQ(B_coefficient[3]), QQ(B_coefficient[4]),
    ])
kernel = matrix(QQ, kernel_rows)
assert kernel.rank() == 2
assert collision_block["matrix"] * kernel.transpose() == matrix(QQ, 8, 2)
# Use the generic compiler matrix even though this particular quotient admits
# a much cheaper structured solve than a dense rational nullspace.  The
# explicit two sections are accepted only after the complete condition matrix
# certifies their dimension and annihilates both rows.
compilation = compile_resolved_conditions(
    ambient, (collision_block,), complete=True, compute_kernel=False
)
assert compilation["condition_matrix"] == collision_block["matrix"]
assert (compilation["ambient_dimension"], compilation["condition_rows"],
        compilation["rank"], compilation["codimension"],
        compilation["kernel_dimension"]) == (10, 8, 8, 8, 2)
assert compilation["h0_certified"]
assert certify_explicit_pencil_basis(compilation, kernel) == kernel

# Reconstruct the two coefficient pairs and replay every condition directly.
sections = []
for row in kernel.rows():
    A_coefficient = sum(
        row[index] * u**exponent
        for index, exponent in enumerate(range(1, 9))
    )
    B_coefficient = sum(
        row[8 + index] * u**exponent
        for index, exponent in enumerate(range(3, 5))
    )
    assert A_coefficient % u == 0
    assert B_coefficient % u**3 == 0
    assert A_coefficient.degree() <= 8
    assert B_coefficient.degree() <= 4
    congruence = A_coefficient * denominator_over_h + B_coefficient * p_numerator
    assert congruence % modulus == 0
    a_coefficient = u_field(A_coefficient) / u_field(h**2)
    b_coefficient = u_field(B_coefficient) / u_field(h)
    # E8 and E7 coefficient orders in the common F_infinity representative.
    assert a_coefficient.valuation(u) >= 1
    assert b_coefficient.valuation(u) >= 3
    assert a_coefficient.numerator().degree() <= a_coefficient.denominator().degree()
    assert b_coefficient.numerator().degree() <= b_coefficient.denominator().degree()
    sections.append({
        "A_coefficients_low_to_high": [str(value) for value in A_coefficient.list()],
        "B_coefficients_low_to_high": [str(value) for value in B_coefficient.list()],
        "a": str(a_coefficient),
        "b": str(b_coefficient),
    })

payload = {
    "schema": "elkies-k3.h92-q6-global-rr.v1",
    "status": "PASS_EXACT_GLOBAL_RR_KERNEL",
    "inputs": {
        "h92_source": {"path": str(H92.relative_to(ROOT)), "sha256": digest(H92)},
        "marked_section": {"path": str(SECTION.relative_to(ROOT)), "sha256": SECTION_SHA256},
        "e7_module": "elkies-k3/scripts/derive_h92_q6_e7_p1_branch_module.sage",
        "e8_module": "elkies-k3/scripts/derive_h92_q6_e8_p1_branch_module.sage",
        "smooth_module": "elkies-k3/scripts/derive_h92_q6_smooth_po_module.sage",
    },
    "common_representative": "D=O+(-P1)-F_infinity",
    "ambient": {
        "coefficient_form": "a=A/h^2, b=B/h, section=a+b*m",
        "A_powers": list(range(1, 9)),
        "B_powers": list(range(3, 5)),
        "dimension": len(ambient),
        "basis": [[kind, exponent] for kind, exponent in ambient],
    },
    "collision_condition": {
        "congruence": "A*(den(y_P/x_P)/h)+B*num(y_P/x_P)=0 mod h^2",
        "matrix": [[str(value) for value in row] for row in collision_block["matrix"].rows()],
        "rank": int(collision_block["matrix"].rank()),
        "codimension": int(collision_block["matrix"].rank()),
    },
    "kernel": {
        "dimension": int(kernel.nrows()),
        "basis_matrix": [[str(value) for value in row] for row in kernel.rows()],
        "sections": sections,
    },
    "claims": {
        "h0_D": 2,
        "basis_ratio_not_yet_eliminated": True,
        "reason": "The ambient exhausts the E7, E8, and four smooth P.O local coefficient lattices in the common F_infinity representative.",
    },
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92Q6GLOBALRR|ambient=10|collision_rank=8|codimension=8|kernel=2|"
    "h0=2|status=PASS_EXACT_GLOBAL_RR_KERNEL",
    flush=True,
)
