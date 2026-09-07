#!/usr/bin/env sage -python
"""Compile the physical q4/orbit164 presentation over QQ.

The selected zero is second_I6_affine and the horizontal is the exact
old_A11_component_0 section on the pointed q4/orbit1584 model.  The vertical
correction uses the two finite I2 fibres where that section meets the identity
component.  At those two supports the chord quartic has a unique rational
double branch.  Linear interpolation through the two branch values gives

    u = (m-ell(T))/((T-r0)(T-r1)).

Substitution removes the exact square of the denominator and leaves a quartic
in T.  This realizes the 4 -> 2 -> 2 RR calculation without a Groebner basis.
"""

import hashlib
import json
import time
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"
PARENT = LOCAL / "q4o208-physical-q4o1584-rr-qq.json"
MARKING = LOCAL / "q4o1584-second-affine-equation-marking-qq.json"
SOURCE_MARKING = GENERATED / "elkies-k3-h3-q4o208-physical-q4o1584-second_i6_affine_component-marking.json"
SCORE = GENERATED / "elkies-k3-h3-q4o208-q4o1584-second_i6_affine_component-q4d2-equation-cost.json"
EDGE = GENERATED / "elkies-k3-h3-q4o208-q4o1584-q4o164-a1a1a3a3-certificate.json"
OUTPUT = LOCAL / "q4o1584-physical-q4o164-rr-qq.json"
INPUTS = (PARENT, MARKING, SOURCE_MARKING, SCORE, EDGE)

started = time.monotonic()
stage_times = {}


def log(stage, **fields):
    elapsed = time.monotonic() - started
    stage_times[stage] = elapsed
    tail = "|".join(f"{key}={value}" for key, value in fields.items())
    print(
        f"Q4O164RRQQ|stage={stage}|elapsed={elapsed:.3f}"
        + (f"|{tail}" if tail else ""),
        flush=True,
    )


def rational_bits(values):
    entries = [QQ(value) for value in values]
    return max(
        max(abs(value.numerator()).nbits(), value.denominator().nbits())
        for value in entries
    )


def polynomial_payload(value):
    return [str(entry) for entry in value.list()]


parent = json.loads(PARENT.read_text())
marking = json.loads(MARKING.read_text())
source_marking = json.loads(SOURCE_MARKING.read_text())
score = json.loads(SCORE.read_text())
edge = json.loads(EDGE.read_text())
assert parent["status"] == "PASS_EXACT_QQ_PHYSICAL_Q4O1584_D4_A3_3A1_RR_AND_JACOBIAN"
assert marking["status"] == "PASS_EXACT_QQ_Q4O1584_SECOND_AFFINE_POINTING_AND_C0_MARKING"
assert source_marking["status"] == "PASS_EXACT_Q4O1584_PHYSICAL_EFFECTIVE_ZERO_MARKING"
assert score["status"] == "PASS_EXACT_MARKED_FRONTIER_EQUATION_COST_SCORING"
assert edge["status"] == "PASS_EXACT_MARKED_DEGREE_TWO_CANDIDATE_CERTIFICATE"

candidate = score["best_candidate"]
assert candidate["candidate_id"] == {"q": 4, "old_fibre_degree": 2, "orbit_index": 164}
assert candidate["expected_RR_ambient"] == 4
assert candidate["horizontal"]["P_dot_O"] == 0
assert candidate["horizontal"]["vertical"] == [0, 0, 0, 0, 0, 0, 0, 1, 0, 1]
assert candidate["horizontal"]["fibre_twist"] == 2
assert edge["candidate_id"] == {
    "label": "q4o164_after_q4o1584", "q": 4, "old_fibre_degree": 2,
}

curves = {
    name: vector(ZZ, coordinates)
    for name, coordinates in source_marking["equation_explicit_curves_in_child"].items()
}
horizontal = vector(ZZ, candidate["horizontal"]["section"])
assert horizontal == curves["old_A11_component_0"]
old_fibre = vector(ZZ, [1, 0] + [0] * 17)
old_zero = vector(ZZ, [-1, 1] + [0] * 17)
vertical = vector(ZZ, [0, 0] + candidate["horizontal"]["vertical"] + [0] * 7)
fibre = vector(ZZ, candidate["fibre"])
assert fibre == old_zero + horizontal + 2 * old_fibre + vertical
log("LITERAL_DIVISOR", zero="second_I6_affine", horizontal="old_A11_component_0", vertical="I2_0+I2_2")

RT = PolynomialRing(QQ, "T")
T = RT.gen()
RU = PolynomialRing(QQ, "U")
U = RU.gen()
SU = PolynomialRing(RU, "T")
TU = SU.gen()
A = RT([QQ(value) for value in parent["child"]["minimal_A_coefficients_low_to_high"]])
B = RT([QQ(value) for value in parent["child"]["minimal_B_coefficients_low_to_high"]])


def read_polynomial(record):
    numerator = RT([QQ(value) for value in record["numerator_coefficients_low_to_high"]])
    denominator = RT([QQ(value) for value in record["denominator_coefficients_low_to_high"]])
    quotient = numerator / denominator
    if quotient.denominator() != 1:
        raise ArithmeticError("q4/o164 horizontal is not polynomial on the stored short model")
    return RT(quotient)


point = marking["old_A11_component_0_on_second_affine_pointed_child"]
PX = read_polynomial(point["x"])
PY = read_polynomial(point["y"])
if PY**2 != PX**3 + A * PX + B:
    raise ArithmeticError("old_A11_component_0 misses q4/o1584 short model")
assert (A.degree(), B.degree(), PX.degree(), PY.degree()) == (6, 9, 3, 4)

i2_records = [
    record for record in parent["child"]["finite_reducible_fibres"]
    if record["kodaira"] == "I2"
]
i2_profile = point["finite_I2_nonidentity_profile"]
if len(i2_records) != 3 or len(i2_profile) != 3:
    raise ArithmeticError("q4/o164 input does not expose three marked finite I2 fibres")
supports = [QQ(record["support"]) for record, hit in zip(i2_records, i2_profile) if not hit]
if len(supports) != 2:
    raise ArithmeticError("C0 does not select the two identity-component I2 supports")

M = PolynomialRing(QQ, "m")
m = M.gen()
branch_values = []
for support in supports:
    special = M(
        m**4 - 6 * PX(support) * m**2 - 8 * PY(support) * m
        - 3 * PX(support)**2 - 4 * A(support)
    )
    factors = list(special.factor())
    double_linear = [
        factor for factor, exponent in factors
        if factor.degree() == 1 and exponent == 2
    ]
    if len(double_linear) != 1:
        raise ArithmeticError("selected I2 support lacks a unique rational double branch")
    branch_values.append(QQ(-double_linear[0][0] / double_linear[0][1]))

r0, r1 = supports
v0, v1 = branch_values
linear = RT(v0 + (v1 - v0) * (T - r0) / (r1 - r0))
denominator = RT((T - r0) * (T - r1))
assert linear(r0) == v0 and linear(r1) == v1
log(
    "RR_CONDITIONS", supports_mod131=",".join(
        str(int(GF(131)(value.numerator()) / GF(131)(value.denominator())))
        for value in supports
    ), branch_bits=rational_bits(branch_values),
)


def lift(poly):
    return SU([RU(value) for value in RT(poly).list()])


m_u = lift(linear) + U * lift(denominator)
px_u, py_u, a_u = lift(PX), lift(PY), lift(A)
radicand = m_u**4 - 6 * px_u * m_u**2 - 8 * py_u * m_u - 3 * px_u**2 - 4 * a_u
quartic, remainder = radicand.quo_rem(lift(denominator)**2)
if remainder or quartic.degree() != 4:
    raise ArithmeticError("two-I2 RR quotient does not produce an exact quartic")
if radicand != lift(denominator)**2 * quartic:
    raise ArithmeticError("two-I2 divisor removal failed")
quartic_coefficients = list(quartic.list()) + [RU.zero()] * 5
coefficient_degrees = [int(value.degree()) for value in quartic_coefficients[:5]]
log("QUARTIC", degree=4, coefficient_degrees=",".join(map(str, coefficient_degrees)))

e, d, c, b, a = quartic_coefficients[:5]
I = 12 * a * e - 3 * b * d + c**2
J = 72 * a * c * e + 9 * b * c * d - 27 * a * d**2 - 27 * b**2 * e - 2 * c**3
A_child = RU(-27 * I)
B_child = RU(-27 * J)
Delta_child = RU(-16 * (4 * A_child**3 + 27 * B_child**2))
log(
    "JACOBIAN_DEGREES", A=A_child.degree(), B=B_child.degree(),
    Delta=Delta_child.degree(),
)
if (A_child.degree(), B_child.degree(), Delta_child.degree()) != (8, 12, 20):
    raise ArithmeticError("q4/o164 child is not a minimal elliptic K3 model")
# In the degree-24 homogeneous discriminant the degree deficit is four,
# while c4 has order zero.  Hence infinity is multiplicative I4.
if not A_child.leading_coefficient() or 24 - Delta_child.degree() != 4:
    raise ArithmeticError("q4/o164 infinity is not multiplicative I4")

gcd_started = time.monotonic()
repeated_gcd = Delta_child.gcd(Delta_child.derivative()).monic()
gcd_seconds = time.monotonic() - gcd_started
if repeated_gcd.degree() != 5:
    raise ArithmeticError("q4/o164 repeated discriminant degree changed")
repeated_support = repeated_gcd.squarefree_part().monic()
if repeated_support.degree() != 3:
    raise ArithmeticError("q4/o164 reducible support degree changed")
support_factors = list(repeated_support.factor())
if len(support_factors) != 3 or any(factor.degree() != exponent for factor, exponent in support_factors):
    raise ArithmeticError("q4/o164 reducible supports are not three rational simple roots")

finite_fibres = []
repeated_product = RU.one()
for factor, unused in support_factors:
    factor = factor.monic()
    root = QQ(-factor[0])
    multiplicity = 0
    quotient = Delta_child
    while not quotient(root):
        quotient, rem = quotient.quo_rem(factor)
        if rem:
            raise ArithmeticError("discriminant multiplicity division failed")
        multiplicity += 1
    if multiplicity not in (2, 4) or not A_child(root) or not B_child(root):
        raise ArithmeticError("q4/o164 reducible fibre is not multiplicative I2/I4")
    repeated_product *= factor**multiplicity
    finite_fibres.append({
        "support": str(root),
        "kodaira": "I4" if multiplicity == 4 else "I2",
        "delta_order": multiplicity,
    })
if sorted(item["delta_order"] for item in finite_fibres) != [2, 2, 4]:
    raise ArithmeticError("q4/o164 finite ADE profile is not 2A1+A3")
nodal, remainder = Delta_child.quo_rem(repeated_product)
if remainder or nodal.degree() != 12 or not nodal.is_squarefree():
    raise ArithmeticError("q4/o164 residual nodal discriminant changed")
if nodal.gcd(repeated_support) != 1 or A_child.gcd(nodal) != 1 or B_child.gcd(nodal) != 1:
    raise ArithmeticError("q4/o164 nodal and reducible supports collide")
log("FIBRES", profile="I4_inf+I4+2I2+12I1", gcd_seconds=f"{gcd_seconds:.3f}")

payload = {
    "schema": "elkies-k3.q4o1584-physical-q4o164-rr-qq.v1",
    "status": "PASS_EXACT_QQ_PHYSICAL_Q4O164_2A3_2A1_RR_AND_JACOBIAN",
    "divisor": {
        "identity": (
            "F_q4o164=O(second_I6_affine)+old_A11_component_0+2*old_fibre+"
            "(F-r7)+(F-r9)"
        ),
        "zero": "second_I6_affine_component",
        "horizontal": "old_A11_component_0",
        "P_dot_O": 0,
        "old_fibre_degree": 2,
        "q": 4,
        "orbit_index": 164,
        "vertical_root_indices": [7, 9],
    },
    "resolved_RR": {
        "ambient_basis": ["1", "T", "T^2", "m"],
        "condition_supports": [str(value) for value in supports],
        "double_branch_values": [str(value) for value in branch_values],
        "interpolating_linear_polynomial": polynomial_payload(linear),
        "common_denominator": polynomial_payload(denominator),
        "kernel_basis_in_[1,T,T^2,m]": [
            [str(denominator[0]), str(denominator[1]), "1", "0"],
            [str(-linear[0]), str(-linear[1]), "0", "1"],
        ],
        "dimensions": {"ambient": 4, "condition_rank": 2, "h0": 2},
        "maximum_condition_rational_bits": rational_bits(supports + branch_values + list(linear)),
    },
    "quartic": {
        "coefficients_in_T_low_to_high": [polynomial_payload(value) for value in quartic.list()],
        "coefficient_degrees_in_U": coefficient_degrees,
        "degree_in_T": 4,
        "exact_divisor_removed": "((T-r0)(T-r1))^2",
        "maximum_rational_bits": rational_bits(
            coefficient for value in quartic.list() for coefficient in value.list()
        ),
    },
    "child": {
        "minimal_A_coefficients_low_to_high": polynomial_payload(A_child),
        "minimal_B_coefficients_low_to_high": polynomial_payload(B_child),
        "degrees_A_B_Delta": [8, 12, 20],
        "finite_reducible_fibres": sorted(
            finite_fibres, key=lambda item: (item["delta_order"], item["support"])
        ),
        "finite_nodal_factor_degree": 12,
        "infinity": {"kodaira": "I4", "orders_A_B_Delta": [0, 0, 4]},
        "ADE": "2A3+2A1",
        "root_rank": 8,
        "MW_rank_if_rho19": 9,
        "euler_number": 24,
        "maximum_A_B_rational_bits": max(
            rational_bits(A_child.list()), rational_bits(B_child.list())
        ),
        "maximum_delta_rational_bits": rational_bits(Delta_child.list()),
    },
    "method": {
        "large_Groebner_required": False,
        "full_discriminant_factorization_required": False,
        "finite_discriminant_gcd_seconds": gcd_seconds,
        "stage_elapsed_seconds": stage_times,
        "runtime_seconds": time.monotonic() - started,
    },
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in INPUTS
        },
    },
    "proof_boundary": (
        "Exact QQ 4-to-2-to-2 RR plane from two rational I2 double-branch values, "
        "exact quartic and binary-quartic Jacobian, finite I4+2I2+12I1 fibres, "
        "I4 infinity, Euler 24, and 2A3+2A1/MW9 conditional on rho=19. "
        "The separate marked edge certificate supplies the bidirectional unimodular "
        "NS transport. Later promoted-route equation edges remain separate gates."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
log(
    "DONE", ambient=4, rank=2, h0=2, quartic=4, ADE="2A3+2A1", MW=9,
    status=payload["status"], output=OUTPUT,
)
