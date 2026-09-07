#!/usr/bin/env sage -python
"""Compile the physical q4/orbit1584 presentation loop over QQ.

status: ACTIVE_PROOF
claim: exact QQ H0 plane, quartic, Jacobian, and fibre profile
outputs: artifacts/local/elkies-k3/q4o208-physical-q4o1584-rr-qq.json

The selected horizontal class is literally the already-attached
``first_I6_affine_component``.  If ``P=(x_P,y_P)`` and

    m = (y+y_P)/(x-x_P),

then the generic degree-two branch equation is

    w^2 = m^4 - 6*x_P*m^2 - 8*y_P*m - 3*x_P^2 - 4*A.

At the exact second inherited I4 support ``s`` this quartic has a unique
rational double branch.  Its unique first jet ``c0+c1*(T-s)`` gives

    u = (m-c0-c1*(T-s))/(T-s)^2.

Exact divisibility by ``(T-s)^4`` leaves a quartic in ``T``.  Thus this is a
resolved one-place jet calculation, not a Groebner or multivariable solve.
"""

import hashlib
import json
import time
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"
PARENT = LOCAL / "q24-2a5-physical-q4o208-rr-qq.json"
MARKING = LOCAL / "q24-2a5-physical-q4o208-equation-marking-qq.json"
PHYSICAL = GENERATED / "elkies-k3-h3-q4o208-physical-3a3-marking.json"
SCORE = GENERATED / "elkies-k3-h3-q4o208-physical-3a3-q4d2-equation-cost.json"
EDGE = GENERATED / "elkies-k3-h3-q4o208-physical-q4o1584-certificate.json"
OUTPUT = LOCAL / "q4o208-physical-q4o1584-rr-qq.json"
INPUTS = (PARENT, MARKING, PHYSICAL, SCORE, EDGE)

started = time.monotonic()
stage_times = {}


def log(stage, **fields):
    elapsed = time.monotonic() - started
    stage_times[stage] = elapsed
    suffix = "|".join(f"{key}={value}" for key, value in fields.items())
    print(
        f"Q4O1584RRQQ|stage={stage}|elapsed={elapsed:.3f}"
        + (f"|{suffix}" if suffix else ""),
        flush=True,
    )


def rational_bits(values):
    values = [QQ(value) for value in values]
    return int(max(
        max(abs(value.numerator()).nbits(), value.denominator().nbits())
        for value in values
    ))


def polynomial_payload(poly):
    poly = poly.parent()(poly)
    return [str(value) for value in poly.list()]


parent = json.loads(PARENT.read_text())
marking = json.loads(MARKING.read_text())
physical = json.loads(PHYSICAL.read_text())
score = json.loads(SCORE.read_text())
edge = json.loads(EDGE.read_text())
assert parent["status"] == "PASS_EXACT_QQ_PHYSICAL_Q4O208_3A3_RR_AND_JACOBIAN"
assert marking["status"] == "PASS_EXACT_QQ_PHYSICAL_Q4O208_C5_EQUATION_MARKING"
assert physical["status"] == "PASS_EXACT_Q4O208_PHYSICAL_3A3_MARKING"
assert score["status"] == "PASS_EXACT_MARKED_FRONTIER_EQUATION_COST_SCORING"
assert edge["status"] == "PASS_EXACT_MARKED_DEGREE_TWO_CANDIDATE_CERTIFICATE"

candidate = score["best_candidate"]
assert candidate["candidate_id"] == {"q": 4, "old_fibre_degree": 2, "orbit_index": 1584}
assert edge["candidate_id"] == {
    "label": "q4o1584-physical-presentation-loop",
    "q": 4,
    "old_fibre_degree": 2,
}
assert candidate["expected_RR_ambient"] == 4
assert candidate["horizontal"]["P_dot_O"] == 0
assert candidate["horizontal"]["vertical"] == [0, 0, 0, 1, 2, 1, 0, 0, 0]
assert candidate["horizontal"]["fibre_twist"] == 2

# Literal divisor identity in the physical 3A3 frame.
fibre = vector(ZZ, candidate["fibre"])
horizontal = vector(ZZ, candidate["horizontal"]["section"])
known_horizontal = vector(
    ZZ, physical["equation_explicit_curves_in_child"]["first_I6_affine_component"]
)
assert horizontal == known_horizontal
old_fibre = vector(ZZ, [1, 0] + [0] * 17)
old_zero = vector(ZZ, [-1, 1] + [0] * 17)
vertical = vector(ZZ, [0, 0] + candidate["horizontal"]["vertical"] + [0] * 8)
assert fibre == old_zero + horizontal + 2 * old_fibre + vertical
log("LITERAL_DIVISOR", horizontal="first_I6_affine_component", vertical="A3_highest_root")


R = PolynomialRing(QQ, "T")
T = R.gen()
A = R([QQ(value) for value in parent["child"]["minimal_A_coefficients_low_to_high"]])
B = R([QQ(value) for value in parent["child"]["minimal_B_coefficients_low_to_high"]])
point = marking["first_I6_affine_component_on_C5_pointed_child"]


def read_rational_polynomial(record):
    numerator = R([QQ(value) for value in record["numerator_coefficients_low_to_high"]])
    denominator = R([QQ(value) for value in record["denominator_coefficients_low_to_high"]])
    quotient = numerator / denominator
    assert quotient.denominator() == 1
    return R(quotient)


PX = read_rational_polynomial(point["x"])
PY = read_rational_polynomial(point["y"])
assert (A.degree(), B.degree(), PX.degree(), PY.degree()) == (8, 12, 4, 6)
assert PY**2 == PX**3 + A * PX + B

support_record = marking["physical_fibres"]["second_old_I6_I4"]
s = QQ(support_record["support"])
assert GF(103)(s.numerator()) / GF(103)(s.denominator()) == GF(103)(38)

# Recover the unique rational double branch c0 at T=s.
M = PolynomialRing(QQ, "m")
m = M.gen()
special = M(m**4 - 6 * PX(s) * m**2 - 8 * PY(s) * m - 3 * PX(s)**2 - 4 * A(s))
special_factors = list(special.factor())
double_linear = [factor for factor, exponent in special_factors if factor.degree() == 1 and exponent == 2]
assert len(double_linear) == 1
assert sorted((factor.degree(), int(exponent)) for factor, exponent in special_factors) == [(1, 2), (2, 1)]
c0 = QQ(-double_linear[0][0] / double_linear[0][1])

# Recover the unique first jet.  The order-two equation is a square of one
# rational linear polynomial; the resulting c1 automatically kills order 3.
C = PolynomialRing(QQ, "c")
c = C.gen()
SL = PolynomialRing(C, "L")
L = SL.gen()


def translate(poly):
    return SL([C(value) for value in R(poly)(L + s).list()])


px_local = translate(PX)
py_local = translate(PY)
a_local = translate(A)
m_jet = C(c0) + c * L
jet_branch = (
    m_jet**4 - 6 * px_local * m_jet**2 - 8 * py_local * m_jet
    - 3 * px_local**2 - 4 * a_local
)
assert jet_branch[0] == jet_branch[1] == 0
order_two_factors = list(jet_branch[2].factor())
jet_linear = [factor for factor, unused in order_two_factors if factor.degree() == 1]
assert len(jet_linear) == 1
assert order_two_factors[0][1] == 2
c1 = QQ(-jet_linear[0][0] / jet_linear[0][1])
assert jet_branch[2](c1) == jet_branch[3](c1) == 0
assert jet_branch[4](c1) != 0
log("RESOLVED_JET", c0_bits=rational_bits([c0]), c1_bits=rational_bits([c1]))


# Compile the exact quartic after m=c0+c1*L+U*L^2.
UR = PolynomialRing(QQ, "U")
U = UR.gen()
SU = PolynomialRing(UR, "L")
LU = SU.gen()


def translate_u(poly):
    return SU([UR(value) for value in R(poly)(LU + s).list()])


px_u = translate_u(PX)
py_u = translate_u(PY)
a_u = translate_u(A)
m_u = UR(c0) + UR(c1) * LU + U * LU**2
radicand = m_u**4 - 6 * px_u * m_u**2 - 8 * py_u * m_u - 3 * px_u**2 - 4 * a_u
quartic, remainder = radicand.quo_rem(LU**4)
assert not remainder and quartic.degree() == 4
assert radicand == LU**4 * quartic
quartic_coefficients = list(quartic.list()) + [UR.zero()] * 5
assert [value.degree() for value in quartic_coefficients[:5]] == [2, 2, 3, 3, 4]
log("QUARTIC", degree=4, coefficient_degrees="2,2,3,3,4")


# Standard binary-quartic invariants.
e, d, q2, b, a = quartic_coefficients[:5]
I = 12 * a * e - 3 * b * d + q2**2
J = 72 * a * q2 * e + 9 * b * q2 * d - 27 * a * d**2 - 27 * b**2 * e - 2 * q2**3
A_child = UR(-27 * I)
B_child = UR(-27 * J)
Delta_child = UR(-16 * (4 * A_child**3 + 27 * B_child**2))
assert (A_child.degree(), B_child.degree(), Delta_child.degree()) == (6, 9, 18)
assert (8 - A_child.degree(), 12 - B_child.degree(), 24 - Delta_child.degree()) == (2, 3, 6)

# The leading residual cubic at infinity is separable, hence the exact local
# type is I0* (D4), not a more degenerate additive fibre.
a2 = A_child.leading_coefficient()
b3 = B_child.leading_coefficient()
assert 4 * a2**3 + 27 * b3**2 != 0
log("INFINITY", orders="2,3,6", kodaira="I0star", root="D4")


# Exact finite repeated support without factoring the full discriminant.
gcd_started = time.monotonic()
repeated_gcd = Delta_child.gcd(Delta_child.derivative()).monic()
gcd_seconds = time.monotonic() - gcd_started
assert repeated_gcd.degree() == 6
repeated_support = repeated_gcd.squarefree_part().monic()
assert repeated_support.degree() == 4
support_factors = list(repeated_support.factor())
assert len(support_factors) == 4
assert all(factor.degree() == exponent == 1 for factor, exponent in support_factors)

finite_fibres = []
repeated_product = UR.one()
for factor, unused in support_factors:
    factor = factor.monic()
    root = QQ(-factor[0])
    multiplicity = 0
    quotient = Delta_child
    while not quotient(root):
        quotient, remainder = quotient.quo_rem(factor)
        assert not remainder
        multiplicity += 1
    assert multiplicity in (2, 4)
    assert A_child(root) and B_child(root)
    repeated_product *= factor**multiplicity
    finite_fibres.append({
        "support": str(root),
        "kodaira": "I4" if multiplicity == 4 else "I2",
        "delta_order": multiplicity,
    })

assert sorted(record["delta_order"] for record in finite_fibres) == [2, 2, 2, 4]
nodal, remainder = Delta_child.quo_rem(repeated_product)
assert not remainder and nodal.degree() == 8 and nodal.is_squarefree()
assert nodal.gcd(repeated_support) == 1
assert A_child.gcd(nodal) == B_child.gcd(nodal) == 1
log("FINITE_FIBRES", profile="I4+3I2+8I1", gcd_seconds=f"{gcd_seconds:.3f}")


payload = {
    "schema": "elkies-k3.q4o208-physical-q4o1584-rr-qq.v1",
    "status": "PASS_EXACT_QQ_PHYSICAL_Q4O1584_D4_A3_3A1_RR_AND_JACOBIAN",
    "divisor": {
        "identity": (
            "F_q4o1584=old_zero+first_I6_affine_component+2*old_fibre+"
            "second_I4_chain_[1,2,1]"
        ),
        "horizontal": "first_I6_affine_component",
        "P_dot_O": 0,
        "old_fibre_degree": 2,
        "q": 4,
        "orbit_index": 1584,
    },
    "resolved_RR": {
        "ambient_basis": ["1", "T", "T^2", "m"],
        "common_denominator": "(T-s)^2",
        "chord": "m=(y+y_P)/(x-x_P)",
        "support": str(s),
        "support_mod_103": 38,
        "double_branch_c0": str(c0),
        "first_jet_c1": str(c1),
        "kernel_basis_in_[1,T,T^2,m]": [
            [str(s**2), str(-2 * s), "1", "0"],
            [str(c1 * s - c0), str(-c1), "0", "1"],
        ],
        "dimensions": {"ambient": 4, "condition_rank": 2, "h0": 2},
        "maximum_jet_rational_bits": rational_bits([s, c0, c1]),
    },
    "quartic": {
        "coefficients_in_L_low_to_high": [polynomial_payload(value) for value in quartic.list()],
        "coefficient_degrees_in_U": [int(value.degree()) for value in quartic.list()],
        "degree_in_L": 4,
        "exact_divisor_removed": "(T-s)^4",
        "maximum_rational_bits": rational_bits(
            coefficient for value in quartic.list() for coefficient in value.list()
        ),
    },
    "child": {
        "minimal_A_coefficients_low_to_high": polynomial_payload(A_child),
        "minimal_B_coefficients_low_to_high": polynomial_payload(B_child),
        "degrees_A_B_Delta": [6, 9, 18],
        "finite_reducible_fibres": sorted(finite_fibres, key=lambda item: (item["delta_order"], item["support"])),
        "finite_nodal_factor_degree": 8,
        "infinity": {"kodaira": "I0*", "orders_A_B_Delta": [2, 3, 6]},
        "ADE": "D4+A3+3A1",
        "root_rank": 10,
        "MW_rank_if_rho19": 7,
        "euler_number": 24,
        "maximum_A_B_rational_bits": max(rational_bits(A_child.list()), rational_bits(B_child.list())),
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
        "Exact QQ H0 plane from the unique resolved double-branch jet, exact quartic, "
        "binary-quartic Jacobian, finite I4+3I2+8I1 fibres, I0* infinity, Euler 24, "
        "and D4+A3+3A1/MW7 conditional on rho=19. The full physical marked-NS edge "
        "and bidirectional unimodular transports are supplied by the separate orbit1584 "
        "candidate certificate. Later route edges remain separate gates."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "Q4O1584RRQQ|ambient=4|rank=2|h0=2|quartic=4|"
    "fibres=I0star+I4+3I2+8I1|ADE=D4+A3+3A1|MW=7|status={}|output={}".format(
        payload["status"], OUTPUT
    ),
    flush=True,
)
