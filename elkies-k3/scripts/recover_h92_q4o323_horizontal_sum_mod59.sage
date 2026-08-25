#!/usr/bin/env sage -python
"""Recover q4/orbit323 simple-pole seeds as sums over GF(59).

In the compact q4/orbit208 chart the corrected horizontal decomposes into a
one-node polynomial section at t=0 and a two-node polynomial section at
t=1,infinity.  Exhaust the latter p^3 shell, add it to the already-certified
one-node modular shell, and retain sums with one simple pole and all three
node incidences.  No elimination or Groebner basis is used.
"""

import hashlib
import itertools
import json
import time
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
COMPACT = LOCAL / "q4o208-compact-weierstrass-qq.json"
ONE_NODE = LOCAL / "q4o208-q4o323-one-node-f1-mod59.json"
OUTPUT = LOCAL / "q4o208-q4o323-horizontal-sums-mod59.json"
PRIME = 59
started = time.monotonic()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


compact = json.loads(COMPACT.read_text())
one_node = json.loads(ONE_NODE.read_text())
assert compact["status"] == "PASS_EXACT_QQ_Q4O208_COMPACT_WEIERSTRASS_NORMALIZATION"
assert one_node["status"] == "PASS_MODP_Q4O323_ONE_NODE_ADJACENT_COMPONENT_SECTION_SCAN"

F = GF(PRIME)
R = PolynomialRing(F, "t")
t = R.gen()
K = R.fraction_field()
RX = PolynomialRing(F, "x")
xvar = RX.gen()


def reduce_qq(value):
    value = QQ(value)
    return F(value.numerator()) / F(value.denominator())


A = R([reduce_qq(value) for value in compact["compact_model"]["A_coefficients_low_to_high"]])
B = R([reduce_qq(value) for value in compact["compact_model"]["B_coefficients_low_to_high"]])
E = EllipticCurve(K, [0, 0, 0, K(A), K(B)])

nodes = []
nodal_reduction = []
for support in (F(0), F(1)):
    cubic = xvar**3 + A(support)*xvar + B(support)
    repeated = cubic.gcd(cubic.derivative())
    nodal_reduction.append(repeated.degree() == 1)
    nodes.append(-repeated[0]/repeated[1] if repeated.degree() == 1 else F(0))
cubic_infinity = xvar**3 + A[8]*xvar + B[12]
repeated_infinity = cubic_infinity.gcd(cubic_infinity.derivative())
assert repeated_infinity.degree() == 1
nodes.append(-repeated_infinity[0]/repeated_infinity[1])
nodal_reduction.append(True)
assert nodal_reduction == [True, False, True]
assert nodes[0] == 3


def transform_one_node(record):
    change = compact["exact_coordinate_change"]
    a, b, c, d, m = [reduce_qq(change[key]) for key in ("a", "b", "c", "d", "m")]
    N = a*t+b
    D = c*t+d
    X = R(sum(F(value)*N**degree*D**(4-degree)
              for degree, value in enumerate(record["x_coefficients_low_to_high"])) / m**2)
    Y = R(sum(F(value)*N**degree*D**(6-degree)
              for degree, value in enumerate(record["y_coefficients_low_to_high"])) / m**3)
    point = E(K(X), K(Y))
    assert Y**2 == X**3+A*X+B
    return point


one_node_points = []
for record in one_node["sections"]:
    point = transform_one_node(record)
    if point not in one_node_points:
        one_node_points.append(point)
assert len(one_node_points) == len(one_node["sections"])


def hits_node(point, index):
    if point.is_zero():
        return False
    x_coordinate, y_coordinate = map(K, point[:2])
    if index < 2:
        support = F(index)
        return (
            x_coordinate.denominator()(support) != 0
            and y_coordinate.denominator()(support) != 0
            and x_coordinate(support) == nodes[index]
            and y_coordinate(support) == 0
        )
    x_numerator, x_denominator = x_coordinate.numerator(), x_coordinate.denominator()
    y_numerator, y_denominator = y_coordinate.numerator(), y_coordinate.denominator()
    x_excess = x_numerator.degree() - x_denominator.degree()
    y_excess = y_numerator.degree() - y_denominator.degree()
    x_value = (
        x_numerator.leading_coefficient() / x_denominator.leading_coefficient()
        if x_excess == 4 else F.zero()
    )
    y_value = (
        y_numerator.leading_coefficient() / y_denominator.leading_coefficient()
        if y_excess == 6 else F.zero()
    )
    return x_excess <= 4 and y_excess <= 6 and x_value == nodes[2] and y_value == 0


# X(1)=the t=1 node and leading(X)=the infinity node leave three parameters.
two_node_points = []
tests = 0
square_x = 0
degree_six_roots = 0
incidence_1_infinity = 0
misses_t0 = 0
end_at_infinity = 0
for x0, x1, x2 in itertools.product(F, repeat=3):
    tests += 1
    x4 = nodes[2]
    x3 = nodes[1] - x0 - x1 - x2 - x4
    X = R([x0, x1, x2, x3, x4])
    rhs = X**3 + A*X + B
    if not rhs.is_square():
        continue
    square_x += 1
    root = rhs.sqrt()
    for Y in ({root, -root} if root else {root}):
        if Y.degree() > 6:
            continue
        degree_six_roots += 1
        point = E(K(X), K(Y))
        if not hits_node(point, 1) or not hits_node(point, 2):
            continue
        incidence_1_infinity += 1
        if hits_node(point, 0):
            continue
        misses_t0 += 1
        if not hits_node(2*point, 2):
            continue
        end_at_infinity += 1
        two_node_points.append(point)


def rational_record(value):
    value = K(value)
    numerator = value.numerator()
    denominator = value.denominator()
    return {
        "numerator_coefficients_low_to_high": [int(entry) for entry in numerator.list()],
        "denominator_coefficients_low_to_high": [int(entry) for entry in denominator.list()],
        "degrees_numerator_denominator": [int(numerator.degree()), int(denominator.degree())],
    }


sums = {}
for one_index, one in enumerate(one_node_points):
    for two_index, two in enumerate(two_node_points):
        point = one + two
        if point.is_zero():
            continue
        x_coordinate, y_coordinate = map(K, point[:2])
        if x_coordinate.denominator().degree() != 2 or y_coordinate.denominator().degree() != 3:
            continue
        if not all(hits_node(point, index) for index in range(3)):
            continue
        if not hits_node(2*point, 0):
            continue
        key = (tuple(x_coordinate.numerator().list()), tuple(x_coordinate.denominator().list()),
               tuple(y_coordinate.numerator().list()), tuple(y_coordinate.denominator().list()))
        sums[key] = {
            "one_node_index": one_index,
            "two_node_index": two_index,
            "x": rational_record(x_coordinate),
            "y": rational_record(y_coordinate),
            "node_hits_at_compact_0_1_infinity": [hits_node(point, index) for index in range(3)],
            "double_node_hits_at_compact_0_1_infinity": [hits_node(2*point, index) for index in range(3)],
        }

payload = {
    "schema": "elkies-k3.h3-q4o208-q4o323-horizontal-sums-mod59.v1",
    "status": "PASS_MOD59_Q4O323_SIMPLE_POLE_SUM_SEEDS",
    "prime": PRIME,
    "search": {
        "two_node_polynomial_x_tests": tests,
        "square_x": square_x,
        "degree_at_most_six_roots": degree_six_roots,
        "t1_infinity_incidence_sections": incidence_1_infinity,
        "also_misses_t0": misses_t0,
        "also_end_component_at_infinity": end_at_infinity,
        "one_node_sections": len(one_node_points),
        "two_node_sections": len(two_node_points),
        "simple_pole_sum_seeds": len(sums),
        "large_Groebner_required": False,
        "runtime_seconds": time.monotonic()-started,
    },
    "two_node_sections": [
        {"x": rational_record(point[0]), "y": rational_record(point[1])}
        for point in two_node_points
    ],
    "horizontal_sum_seeds": list(sums.values()),
    "proof_boundary": (
        "This exhaustive mod-59 construction supplies simple-pole section seeds as sums "
        "of the one-node and complementary two-node polynomial shells. The t=1 I4 has "
        "cuspidal bad reduction and is used only as an incidence condition. Exact lifting "
        "and marked q4/orbit323 identification remain separate gates."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in (COMPACT, ONE_NODE)],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in (COMPACT, ONE_NODE)},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "Q4O323SUM59|tests={}|two_node={}|sums={}|seconds={:.3f}|status={}|output={}".format(
        tests, len(two_node_points), len(sums), payload["search"]["runtime_seconds"],
        payload["status"], OUTPUT,
    ), flush=True,
)
