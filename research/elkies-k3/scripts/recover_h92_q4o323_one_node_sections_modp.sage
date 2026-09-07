#!/usr/bin/env sage -python
"""Recover the missing one-node polynomial shell on q4/orbit208 modulo p.

The corrected q4/orbit323 horizontal class is a Mordell--Weil sum of two
sections disjoint from the current zero.  One summand meets two I4 nodes and
is already present in ``recover_h92_q4o323_horizontal_mod131.sage``.  The
other meets only one I4 node.  This script exhausts the resulting p^4
polynomial-x chart, tests the Weierstrass right hand side for being a square,
and retains the adjacent-component branches at the selected I4 fibre.  No
Groebner basis or characteristic-zero inference is used.
"""

import argparse
import hashlib
import itertools
import json
import time
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, matrix


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
MODEL = LOCAL / "q24-2a5-physical-q4o208-rr-qq.json"
MARKING = LOCAL / "q24-2a5-physical-q4o208-equation-marking-qq.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, default=59)
parser.add_argument("--fibre-index", type=int, default=1, choices=(0, 1, 2))
parser.add_argument("--output", type=Path)
args = parser.parse_args()
PRIME = int(args.prime)
SELECTED = int(args.fibre_index)
OUTPUT = args.output or LOCAL / f"q4o208-q4o323-one-node-f{SELECTED}-mod{PRIME}.json"
if not OUTPUT.is_absolute():
    OUTPUT = ROOT / OUTPUT

started = time.monotonic()
model = json.loads(MODEL.read_text())
marking = json.loads(MARKING.read_text())
assert model["status"] == "PASS_EXACT_QQ_PHYSICAL_Q4O208_3A3_RR_AND_JACOBIAN"
assert marking["status"] == "PASS_EXACT_QQ_PHYSICAL_Q4O208_C5_EQUATION_MARKING"

F = GF(PRIME)
R = PolynomialRing(F, "u")
u = R.gen()
K = R.fraction_field()
RX = PolynomialRing(F, "x")
xvar = RX.gen()


def reduce_qq(value):
    value = QQ(value)
    return F(value.numerator()) / F(value.denominator())


A = R([reduce_qq(value) for value in model["child"]["minimal_A_coefficients_low_to_high"]])
B = R([reduce_qq(value) for value in model["child"]["minimal_B_coefficients_low_to_high"]])
E = EllipticCurve(K, [0, 0, 0, K(A), K(B)])
fibre_records = list(marking["physical_fibres"].items())
supports = [reduce_qq(record["support"]) for unused, record in fibre_records]
identity_indices = [int(record["identity_component_index"]) for unused, record in fibre_records]
assert len(set(supports)) == 3
discriminant = -F(16) * (4 * A**3 + 27 * B**2)
assert discriminant.degree() == 24
assert all(((u - support)**4).divides(discriminant) for support in supports)

nodes = []
nodal_reduction = []
for support in supports:
    cubic = xvar**3 + A(support) * xvar + B(support)
    repeated = cubic.gcd(cubic.derivative())
    nodal_reduction.append(repeated.degree() == 1)
    if repeated.degree() == 1:
        nodes.append(-repeated[0] / repeated[1])
    else:
        # A bad auxiliary reduction may merge the node to the cusp x=0.
        # It is deliberately not used as an exclusion filter below.
        assert cubic == xvar**3
        nodes.append(F.zero())
assert nodal_reduction[SELECTED]


def node_hit(point, index):
    if point.is_zero():
        return False
    x_coordinate, y_coordinate = K(point[0]), K(point[1])
    support = supports[index]
    return (
        x_coordinate.denominator()(support) != 0
        and y_coordinate.denominator()(support) != 0
        and x_coordinate(support) == nodes[index]
        and y_coordinate(support) == 0
    )


def constrained_x(coefficients):
    """Degree-at-most-four X with X(selected support)=selected node."""
    support = supports[SELECTED]
    first = list(coefficients)
    if support == 0:
        return R([nodes[SELECTED]] + first)
    last = (
        nodes[SELECTED]
        - sum(first[degree] * support**degree for degree in range(4))
    ) / support**4
    return R(first + [last])


def coefficient_jacobian(X, Y):
    dx = -3 * X**2 - A
    dy = 2 * Y
    return matrix(F, 13, 12, lambda degree, variable: (
        (dx * u**variable)[degree]
        if variable < 5
        else (dy * u**(variable - 5))[degree]
    ))


def incidence_rows():
    answer = matrix(F, 2, 12)
    support = supports[SELECTED]
    for degree in range(5):
        answer[0, degree] = support**degree
    for degree in range(7):
        answer[1, 5 + degree] = support**degree
    return answer


candidates = {}
tests = 0
for coefficients in itertools.product(F, repeat=4):
    tests += 1
    X = constrained_x(coefficients)
    rhs = X**3 + A * X + B
    if not rhs.is_square():
        continue
    root = rhs.sqrt()
    for Y in ({root, -root} if root else {root}):
        if Y.degree() > 6:
            continue
        point = E(K(X), K(Y))
        hits = [node_hit(point, index) for index in range(3)]
        if not hits[SELECTED]:
            continue
        if any(hits[index] for index in range(3) if index != SELECTED and nodal_reduction[index]):
            continue
        double_hits = [node_hit(2 * point, index) for index in range(3)]
        # Physical index 2 is adjacent to the physical identity index 3.
        # Its relative Z/4 label is +/-1, so doubling still meets the node.
        if not double_hits[SELECTED]:
            continue
        jacobian = coefficient_jacobian(X, Y)
        augmented = jacobian.stack(incidence_rows())
        if augmented.rank() != 12:
            continue
        key = (tuple(X.list()), tuple(Y.list()))
        candidates[key] = {
            "x_coefficients_low_to_high": [int(value) for value in X.list()],
            "y_coefficients_low_to_high": [int(value) for value in Y.list()],
            "node_hits": hits,
            "double_node_hits": double_hits,
            "ordinary_weierstrass_jacobian_rank": int(jacobian.rank()),
            "resolved_node_augmented_jacobian_rank": int(augmented.rank()),
        }

answers = list(candidates.values())
payload = {
    "schema": "elkies-k3.h3-q4o208-q4o323-one-node-sections-modp.v1",
    "status": "PASS_MODP_Q4O323_ONE_NODE_ADJACENT_COMPONENT_SECTION_SCAN",
    "prime": PRIME,
    "selected_fibre": {
        "index": SELECTED,
        "label": fibre_records[SELECTED][0],
        "support": int(supports[SELECTED]),
        "node": int(nodes[SELECTED]),
        "physical_identity_component_index": identity_indices[SELECTED],
        "required_physical_component_index": 2,
        "required_relative_component_distance_mod4": 1,
        "nodal_reduction_at_all_physical_fibres": nodal_reduction,
    },
    "search": {
        "constrained_polynomial_x_tests": tests,
        "unique_filtered_sections": len(answers),
        "large_Groebner_required": False,
        "runtime_seconds": time.monotonic() - started,
    },
    "sections": answers,
    "proof_boundary": (
        "This is an exhaustive finite-field scan of the selected-node, adjacent-I4-component "
        "polynomial chart. Auxiliary fibres with bad nodal reduction are not used as exclusion "
        "filters. Exact QQ lifting, physical component orientation, and the "
        "q4/orbit323 horizontal group-sum identification remain separate gates."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in (MODEL, MARKING)],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (MODEL, MARKING)
        },
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "Q4O323ONENODE|prime={}|fibre={}|tests={}|candidates={}|seconds={:.3f}|status={}|output={}".format(
        PRIME, fibre_records[SELECTED][0], tests, len(answers),
        payload["search"]["runtime_seconds"], payload["status"], OUTPUT,
    ),
    flush=True,
)
