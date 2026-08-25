#!/usr/bin/env sage -python
"""Lift the q4/orbit164 zero-node rank-eight branch to QQ.

The compact polynomial-section chart has twelve coefficients.  Its ordinary
thirteen-equation Weierstrass Jacobian has full column rank at the recovered
mod-131 seed, so direct Newton lifting and rational reconstruction apply with
no node-incidence rows.  The exact inverse C8-pointed quartic map then proves
old-parent degree two.  These coarse equation fingerprints do not determine a
marked component orientation or NS class.  No Groebner basis is used.
"""

import argparse
import hashlib
import json
import time
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, Qp, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
MODEL = LOCAL / "q4o164-compact-weierstrass-qq.json"
QUARTIC = LOCAL / "q4o164-compact-binary-quartic-qq.json"
POINTING = LOCAL / "q4o164-c8-equation-marking-qq.json"
MODULAR = LOCAL / "q4o164-zero-node-section-mod131.json"
AUDIT = LOCAL / "q4o164-zero-one-node-parent-degree-audit.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--precision", type=int, default=160)
parser.add_argument(
    "--output", default="artifacts/local/elkies-k3/q4o164-zero-node-section-qq.json",
)
args = parser.parse_args()
OUTPUT = Path(args.output)
if not OUTPUT.is_absolute():
    OUTPUT = ROOT / OUTPUT
INPUTS = (MODEL, QUARTIC, POINTING, MODULAR, AUDIT)

started = time.monotonic()
prime = ZZ(131)
precision = int(args.precision)
assert precision >= 80


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coefficient_bits(value):
    value = QQ(value)
    return max(abs(ZZ(value.numerator())).nbits(), ZZ(value.denominator()).nbits())


model = json.loads(MODEL.read_text())
quartic_data = json.loads(QUARTIC.read_text())
pointing = json.loads(POINTING.read_text())
modular = json.loads(MODULAR.read_text())
audit = json.loads(AUDIT.read_text())
assert model["status"] == "PASS_EXACT_QQ_Q4O164_COMPACT_WEIERSTRASS_NORMALIZATION"
assert quartic_data["status"] == "PASS_EXACT_QQ_Q4O164_COMPACT_BINARY_QUARTIC"
assert pointing["status"] == "PASS_EXACT_QQ_Q4O164_C8_EQUATION_MARKING"
assert modular["status"] == "PASS_MOD131_Q4O164_ZERO_NODE_PARENT_DEGREE_TWO_RANK8_BRANCH"
assert audit["status"] == "PASS_EXACT_Q4O164_ZERO_ONE_NODE_PARENT_DEGREE_AUDIT"

RTQ = PolynomialRing(QQ, "t")
tq = RTQ.gen()
KQ = RTQ.fraction_field()
A_QQ = RTQ([QQ(value) for value in model["compact_model"]["A_coefficients_low_to_high"]])
B_QQ = RTQ([QQ(value) for value in model["compact_model"]["B_coefficients_low_to_high"]])


def reduce_mod_prime(value, field):
    value = QQ(value)
    return field(value.numerator()) / field(value.denominator())


F = GF(prime)
RTF = PolynomialRing(F, "t")
tf = RTF.gen()
A_F = RTF([reduce_mod_prime(value, F) for value in A_QQ])
B_F = RTF([reduce_mod_prime(value, F) for value in B_QQ])
seed = vector(F,
    modular["section"]["compact_x_coefficients_low_to_high"]
    + modular["section"]["compact_y_coefficients_low_to_high"]
)
assert len(seed) == 12

K = Qp(prime, prec=precision, type="capped-rel")
RT = PolynomialRing(K, "t")
t = RT.gen()
A = RT([K(value) for value in A_QQ])
B = RT([K(value) for value in B_QQ])


def polynomials(values, ring):
    return ring(list(values[:5])), ring(list(values[5:]))


def residual(values):
    X, Y = polynomials(values, RT)
    equation = Y**2 - X**3 - A * X - B
    return vector(K, [equation[index] for index in range(13)])


def jacobian(values, ring, surface_A):
    X, Y = polynomials(values, ring)
    x_derivative = -3 * X**2 - surface_A
    y_derivative = 2 * Y
    zero = ring.base_ring().zero()
    return matrix(ring.base_ring(), [[
        x_derivative[degree-shift]
        if 0 <= degree-shift <= x_derivative.degree() else zero
        for shift in range(5)
    ] + [
        y_derivative[degree-shift]
        if 0 <= degree-shift <= y_derivative.degree() else zero
        for shift in range(7)
    ] for degree in range(13)])


def minimum_valuation(values):
    nonzero = [value.valuation() for value in values if value]
    return min(nonzero) if nonzero else precision


def rational_reconstruct(value, digits):
    modulus = prime**digits
    residue = ZZ(value.lift()) % modulus
    return QQ(residue.rational_reconstruction(modulus))


X_F, Y_F = polynomials(seed, RTF)
assert Y_F**2 == X_F**3 + A_F * X_F + B_F
J_F = jacobian(seed, RTF, A_F)
rank = int(J_F.rank())
assert rank == 12
pivot_rows = list(map(int, J_F.transpose().pivots()))
assert pivot_rows == list(range(12))
determinant = int(matrix(F, [J_F.row(row) for row in pivot_rows]).det())
assert determinant == 130

values = vector(K, [K(value).add_bigoh(1) for value in seed])
known_precision = 1
iterations = []
while known_precision < precision:
    working_precision = min(2 * known_precision, precision)
    values = vector(K, [K(value.lift()).add_bigoh(working_precision) for value in values])
    full_residual = residual(values)
    chosen = vector(K, [full_residual[row] for row in pivot_rows])
    J = jacobian(values, RT, A)
    square = matrix(K, [J.row(row) for row in pivot_rows])
    correction = square.solve_right(-chosen)
    values += correction
    iterations.append({
        "working_precision_p_adic_digits": working_precision,
        "minimum_full_residual_valuation_after": int(minimum_valuation(residual(values))),
        "minimum_correction_valuation": int(minimum_valuation(correction)),
    })
    known_precision = working_precision

reconstruction_digits = precision - 10
reconstructed = [rational_reconstruct(value, reconstruction_digits) for value in values]
assert [reduce_mod_prime(value, F) for value in reconstructed] == list(seed)
X_QQ, Y_QQ = polynomials(reconstructed, RTQ)
assert Y_QQ**2 == X_QQ**3 + A_QQ * X_QQ + B_QQ

# Exact zero-node profile on the compact child.
RX = PolynomialRing(QQ, "x")
xvar = RX.gen()
supports = []
nodes = []
for fibre in model["compact_model"]["reducible_fibres"]:
    if fibre["support"] == "infinity":
        supports.append(None)
        cubic = xvar**3 + A_QQ[8] * xvar + B_QQ[12]
    else:
        support = QQ(fibre["support"])
        supports.append(support)
        cubic = xvar**3 + A_QQ(support) * xvar + B_QQ(support)
    repeated = cubic.gcd(cubic.derivative())
    assert repeated.degree() == 1
    nodes.append(-repeated[0] / repeated[1])
node_hits = []
for support, node in zip(supports, nodes):
    if support is None:
        node_hits.append(X_QQ[4] == node and Y_QQ[6] == 0)
    else:
        node_hits.append(X_QQ(support) == node and Y_QQ(support) == 0)
assert node_hits == [False, False, False, False]

# Invert the double-normalized C8-pointed quartic map exactly.
quartic = [
    RTQ([QQ(value) for value in row])
    for row in quartic_data["coefficients_in_T_low_to_high"]
]
raw_W0_record = pointing["selected_zero"]["quartic_ordinate"]
RU = PolynomialRing(QQ, "U")
raw_W0_numerator = RU([QQ(value) for value in raw_W0_record["numerator_coefficients_low_to_high"]])
raw_W0_denominator = RU([QQ(value) for value in raw_W0_record["denominator_coefficients_low_to_high"]])
child_c = QQ(model["exact_coordinate_change"]["c"])
child_s = QQ(model["exact_coordinate_change"]["s"])
parent_c = QQ(quartic_data["exact_coordinate_change"]["parent_c"])
W0 = KQ(RTQ(raw_W0_numerator(child_c * tq))) / (
    KQ(RTQ(raw_W0_denominator(child_c * tq))) * child_s * parent_c
)
e, d, c, b, a = map(KQ, quartic)
assert e == W0**2
a1 = d / W0
a2 = c - d**2 / (4 * W0**2)
a3 = 2 * W0 * b
b2 = a1**2 + 4 * a2
x_short, y_short = KQ(X_QQ), KQ(Y_QQ)
x_general = x_short / 9 - b2 / 12
y_general = y_short / 27 - (a1 * x_general + a3) / 2
V = KQ(2 * W0 * (x_general + a2) / y_general)
W = KQ((x_general * V**2 - d * V) / (2 * W0) - W0)
assert W**2 == sum(KQ(value) * V**index for index, value in enumerate(quartic))
parent_degree = max(V.numerator().degree(), V.denominator().degree())
assert parent_degree == 2

payload = {
    "schema": "elkies-k3.q4o164-zero-node-section-qq.v1",
    "status": "PASS_EXACT_QQ_Q4O164_ZERO_NODE_PARENT_DEGREE_TWO_RANK8_BRANCH",
    "prime": int(prime),
    "resolved_hensel": {
        "working_precision_p_adic_digits": precision,
        "variables": 12,
        "weierstrass_coefficient_equations": 13,
        "ordinary_mod131_jacobian_rank": rank,
        "selected_independent_equation_rows": pivot_rows,
        "selected_jacobian_determinant_mod131": determinant,
        "iterations": iterations,
    },
    "section": {
        "x_coefficients_low_to_high": [str(value) for value in X_QQ.list()],
        "y_coefficients_low_to_high": [str(value) for value in Y_QQ.list()],
        "maximum_x_rational_bits": max(map(coefficient_bits, X_QQ)),
        "maximum_y_rational_bits": max(map(coefficient_bits, Y_QQ)),
        "exact_compact_weierstrass_identity": True,
        "reduction_equals_mod131_seed": True,
        "exact_node_hits": node_hits,
    },
    "inverse_pointed_quartic": {
        "parent_base_numerator_coefficients_low_to_high": [str(value) for value in V.numerator().list()],
        "parent_base_denominator_coefficients_low_to_high": [str(value) for value in V.denominator().list()],
        "parent_base_degree": int(parent_degree),
        "exact_quartic_identity": True,
    },
    "marking_boundary": {
        "marked_NS_class_identified": False,
        "reason": (
            "Equation node hits and the degree of the inverse quartic base map do not by "
            "themselves attach the compact fibre supports to oriented physical root chains."
        ),
        "mod131_pool_relation": "P[5]+P[10] inside the stored rank-eight pool",
    },
    "method": {
        "large_Groebner_required": False,
        "node_incidence_rows_required": False,
        "runtime_seconds": time.monotonic() - started,
    },
    "proof_boundary": (
        "The displayed section is exact over QQ and has zero-node/inverse-parent-degree-two "
        "equation fingerprints. Its modular construction lies in the known rank-eight pool; "
        "the former marked-class and missing-direction claims are withdrawn."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {
            str(path.relative_to(ROOT)): sha256(path) for path in INPUTS
        },
    },
}
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "Q4O164ZERONODELIFT|rank={}|det={}|bits={}/{}|parent_degree={}|nodes={}|status={}|output={}".format(
        rank, determinant, payload["section"]["maximum_x_rational_bits"],
        payload["section"]["maximum_y_rational_bits"], parent_degree, node_hits,
        payload["status"], OUTPUT,
    ),
    flush=True,
)
