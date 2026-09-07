#!/usr/bin/env sage -python
"""Recover a q4/orbit164 zero-node rank-eight branch modulo 131.

Search sums and differences of the stored 128 pair-node sections.  Invert the
certified C8-pointed quartic map on every integral result and retain exactly
those with no reducible-node hit and old q4/o1584 base degree two.  The unique
answer is P[5]+P[10], hence it is visibly inside the searched rank-eight
subgroup.  This is a small group-law search, not an elimination or Groebner
calculation.
"""

import hashlib
import json
import time
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, factorial


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
RAW = LOCAL / "q4o1584-physical-q4o164-rr-qq.json"
COMPACT = LOCAL / "q4o164-compact-weierstrass-qq.json"
POINTING = LOCAL / "q4o164-c8-equation-marking-qq.json"
POOL = LOCAL / "q4o164-integral-sections-mod131.json"
OUTPUT = LOCAL / "q4o164-zero-node-section-mod131.json"
INPUTS = (RAW, COMPACT, POINTING, POOL)
PRIME = 131

started = time.monotonic()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


raw = json.loads(RAW.read_text())
compact = json.loads(COMPACT.read_text())
pointing = json.loads(POINTING.read_text())
pool = json.loads(POOL.read_text())
assert raw["status"] == "PASS_EXACT_QQ_PHYSICAL_Q4O164_2A3_2A1_RR_AND_JACOBIAN"
assert compact["status"] == "PASS_EXACT_QQ_Q4O164_COMPACT_WEIERSTRASS_NORMALIZATION"
assert pointing["status"] == "PASS_EXACT_QQ_Q4O164_C8_EQUATION_MARKING"
assert pool["status"] == "PASS_MOD131_Q4O164_PAIR_NODE_SECTION_SUBGROUP_RANK8"

F = GF(PRIME)
R = PolynomialRing(F, "U")
U = R.gen()
K = R.fraction_field()
RX = PolynomialRing(F, "x")
xvar = RX.gen()


def reduce_qq(value):
    value = QQ(value)
    return F(value.numerator()) / F(value.denominator())


def reduce_record(record):
    numerator = R([reduce_qq(value) for value in record["numerator_coefficients_low_to_high"]])
    denominator = R([reduce_qq(value) for value in record["denominator_coefficients_low_to_high"]])
    return K(numerator) / K(denominator)


A = R([reduce_qq(value) for value in raw["child"]["minimal_A_coefficients_low_to_high"]])
B = R([reduce_qq(value) for value in raw["child"]["minimal_B_coefficients_low_to_high"]])
E = EllipticCurve(K, [0, 0, 0, K(A), K(B)])
points = [
    E(
        K(R(record["x_coefficients_low_to_high"])),
        K(R(record["y_coefficients_low_to_high"])),
    )
    for record in pool["integral_subgroup"]["all_sections"]
]
assert len(points) == 128

# Reconstruct the C8-pointed generalized model used by the exact marking.
quartic_coefficients = [
    R([reduce_qq(value) for value in row])
    for row in raw["quartic"]["coefficients_in_T_low_to_high"]
]
t0 = reduce_qq(pointing["selected_zero"]["parent_base_support"])
W0 = reduce_record(pointing["selected_zero"]["quartic_ordinate"])
translated = [
    sum(
        K(coefficient) * F(factorial(index))
        / F(factorial(order) * factorial(index-order)) * K(t0)**(index-order)
        for index, coefficient in enumerate(quartic_coefficients)
        if index >= order
    )
    for order in range(5)
]
e, d, c, b, a = translated
assert e == W0**2
a1 = d / W0
a2 = c - d**2 / (4 * W0**2)
a3 = 2 * W0 * b
b2 = a1**2 + 4 * a2


def inverse_parent_base(point):
    x_short, y_short = K(point[0]), K(point[1])
    x_general = x_short / 9 - b2 / 12
    y_general = y_short / 27 - (a1 * x_general + a3) / 2
    if not y_general:
        return None
    z = 2 * W0 * (x_general + a2) / y_general
    T = K(t0) + z
    W = (x_general * z**2 - d * z) / (2 * W0) - W0
    quartic_value = sum(K(value) * T**index for index, value in enumerate(quartic_coefficients))
    if W**2 != quartic_value:
        raise ArithmeticError("inverse pointed map missed the quartic")
    return T, W


supports = [reduce_qq(record["support"]) for record in raw["child"]["finite_reducible_fibres"]]
finite_nodes = []
for support in supports:
    cubic = xvar**3 + A(support) * xvar + B(support)
    repeated = cubic.gcd(cubic.derivative())
    assert repeated.degree() == 1
    finite_nodes.append(-repeated[0] / repeated[1])
cubic_infinity = xvar**3 + A[8] * xvar + B[12]
repeated_infinity = cubic_infinity.gcd(cubic_infinity.derivative())
assert repeated_infinity.degree() == 1
node_infinity = -repeated_infinity[0] / repeated_infinity[1]


def is_integral_section(point):
    if point.is_zero():
        return False
    x, y = K(point[0]), K(point[1])
    return (
        x.denominator().degree() == 0 and y.denominator().degree() == 0
        and x.numerator().degree() <= 4 and y.numerator().degree() <= 6
    )


def node_hits(point):
    x, y = R(point[0]), R(point[1])
    hits = [
        x(support) == node and y(support) == 0
        for support, node in zip(supports, finite_nodes)
    ]
    hits.append(x[4] == node_infinity and y[6] == 0)
    return hits


candidates = {}
tests = 0
for left in range(len(points)):
    for right in range(left, len(points)):
        for operation, point in (("sum", points[left] + points[right]),
                                 ("difference", points[left] - points[right])):
            tests += 1
            if not is_integral_section(point):
                continue
            if any(node_hits(point)):
                continue
            inverse = inverse_parent_base(point)
            if inverse is None:
                continue
            parent_base, ordinate = inverse
            degree = max(parent_base.numerator().degree(), parent_base.denominator().degree())
            if degree != 2:
                continue
            key = (str(point[0]), str(point[1]))
            provenance = {
                "left_pool_index": left,
                "right_pool_index": right,
                "operation": operation,
            }
            if key in candidates:
                candidates[key]["pool_group_law_provenance"].append(provenance)
                continue
            candidates[key] = {
                "pool_group_law_provenance": [provenance],
                "raw_x_coefficients_low_to_high": [int(value) for value in R(point[0]).list()],
                "raw_y_coefficients_low_to_high": [int(value) for value in R(point[1]).list()],
                "raw_node_hits": node_hits(point),
                "parent_base": {
                    "numerator_coefficients_low_to_high": [int(value) for value in parent_base.numerator().list()],
                    "denominator_coefficients_low_to_high": [int(value) for value in parent_base.denominator().list()],
                    "degree": int(degree),
                },
                "inverse_quartic_identity": True,
            }

assert len(candidates) == 1
answer = next(iter(candidates.values()))
assert {
    "left_pool_index": 5, "right_pool_index": 10, "operation": "sum",
} in answer["pool_group_law_provenance"]

# Apply the already-certified compact child scaling.
child_c = reduce_qq(compact["exact_coordinate_change"]["c"])
child_s = reduce_qq(compact["exact_coordinate_change"]["s"])
x_raw = R(answer["raw_x_coefficients_low_to_high"])
y_raw = R(answer["raw_y_coefficients_low_to_high"])
x_compact = R(x_raw(child_c * U) / child_s**2)
y_compact = R(y_raw(child_c * U) / child_s**3)
A_compact = R([reduce_qq(value) for value in compact["compact_model"]["A_coefficients_low_to_high"]])
B_compact = R([reduce_qq(value) for value in compact["compact_model"]["B_coefficients_low_to_high"]])
assert y_compact**2 == x_compact**3 + A_compact * x_compact + B_compact
answer["compact_x_coefficients_low_to_high"] = [int(value) for value in x_compact.list()]
answer["compact_y_coefficients_low_to_high"] = [int(value) for value in y_compact.list()]
answer["compact_section_identity"] = True

payload = {
    "schema": "elkies-k3.q4o164-zero-node-section-mod131.v1",
    "status": "PASS_MOD131_Q4O164_ZERO_NODE_PARENT_DEGREE_TWO_RANK8_BRANCH",
    "prime": PRIME,
    "search": {
        "pool_sections": len(points),
        "pair_sum_difference_tests": tests,
        "unique_filtered_sections": len(candidates),
        "large_Groebner_required": False,
        "runtime_seconds": time.monotonic() - started,
    },
    "section": answer,
    "proof_boundary": (
        "This exact finite-field group-law search identifies a unique branch in the stored "
        "rank-eight pool with the zero-node and inverse-parent-degree-two fingerprints. "
        "Those coarse fingerprints do not identify a marked NS class."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {
            str(path.relative_to(ROOT)): sha256(path) for path in INPUTS
        },
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "Q4O164ZERONODE131|tests={}|combination={}+{}|compact={}/{}|status={}|output={}".format(
        tests, 5, 10,
        answer["compact_x_coefficients_low_to_high"], answer["compact_y_coefficients_low_to_high"],
        payload["status"], OUTPUT,
    ),
    flush=True,
)
