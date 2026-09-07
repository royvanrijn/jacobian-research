#!/usr/bin/env sage -python
"""Recover primitive odd-I4 polynomial sections on q4/orbit164 modulo p.

The physical norm-four shell shows that every P.O=0 section outside the
known rank-eight hyperplane either misses all reducible nodes or meets an end
component of one of the two I4 fibres.  This script treats the latter case.
It fixes the nodal value of a degree-four x polynomial, exhausts the remaining
p^4 coefficients, tests whether x^3+A*x+B is a polynomial square, and retains
exactly-one-node branches whose double still meets that I4 node.  The last
condition distinguishes component labels 1/3 from the opposite label 2.

This is a bounded univariate polynomial-square scan.  It uses no Groebner
basis and makes no characteristic-zero or marked-lattice claim.
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
MODEL = LOCAL / "q4o164-compact-weierstrass-qq.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, default=41)
parser.add_argument("--fibre", choices=("finite", "infinity"), required=True)
parser.add_argument("--output", type=Path)
args = parser.parse_args()
PRIME = int(args.prime)
OUTPUT = args.output or LOCAL / f"q4o164-odd-{args.fibre}-i4-sections-mod{PRIME}.json"
if not OUTPUT.is_absolute():
    OUTPUT = ROOT / OUTPUT

started = time.monotonic()
model = json.loads(MODEL.read_text())
assert model["status"] == "PASS_EXACT_QQ_Q4O164_COMPACT_WEIERSTRASS_NORMALIZATION"

F = GF(PRIME)
R = PolynomialRing(F, "t")
t = R.gen()
K = R.fraction_field()
RX = PolynomialRing(F, "x")
xvar = RX.gen()


def reduce_qq(value):
    value = QQ(value)
    return F(value.numerator()) / F(value.denominator())


A = R([reduce_qq(value) for value in model["compact_model"]["A_coefficients_low_to_high"]])
B = R([reduce_qq(value) for value in model["compact_model"]["B_coefficients_low_to_high"]])
assert (A.degree(), B.degree()) == (8, 12)
E = EllipticCurve(K, [0, 0, 0, K(A), K(B)])

supports = []
nodes = []
orders = []
labels = []
for index, fibre in enumerate(model["compact_model"]["reducible_fibres"]):
    if fibre["support"] == "infinity":
        support = None
        cubic = xvar**3 + A[8] * xvar + B[12]
        label = "infinity_I4"
    else:
        support = reduce_qq(fibre["support"])
        cubic = xvar**3 + A(support) * xvar + B(support)
        label = f"finite_{index}_{fibre['kodaira']}"
    repeated = cubic.gcd(cubic.derivative())
    assert repeated.degree() == 1
    supports.append(support)
    nodes.append(-repeated[0] / repeated[1])
    orders.append(int(fibre["kodaira"][1:]))
    labels.append(label)

selected = next(
    index for index, (support, order) in enumerate(zip(supports, orders))
    if order == 4 and ((support is None) == (args.fibre == "infinity"))
)


def node_hit(point, index):
    if point.is_zero():
        return False
    x_coordinate, y_coordinate = K(point[0]), K(point[1])
    support = supports[index]
    if support is None:
        return (
            x_coordinate.denominator().degree() == 0
            and y_coordinate.denominator().degree() == 0
            and x_coordinate.numerator()[4] == nodes[index]
            and y_coordinate.numerator()[6] == 0
        )
    return (
        x_coordinate.denominator()(support) != 0
        and y_coordinate.denominator()(support) != 0
        and x_coordinate(support) == nodes[index]
        and y_coordinate(support) == 0
    )


def constrained_x(coefficients):
    support = supports[selected]
    node = nodes[selected]
    if support is None:
        return R(list(coefficients) + [node])
    if support == 0:
        return R([node] + list(coefficients))
    first = list(coefficients)
    last = (node - sum(first[degree] * support**degree for degree in range(4))) / support**4
    return R(first + [last])


def coefficient_jacobian(X, Y):
    dx = -3 * X**2 - A
    dy = 2 * Y
    return matrix(F, 13, 12, lambda degree, variable: (
        (dx * t**variable)[degree]
        if variable < 5
        else (dy * t**(variable - 5))[degree]
    ))


def incidence_rows():
    answer = matrix(F, 2, 12)
    support = supports[selected]
    if support is None:
        answer[0, 4] = 1
        answer[1, 11] = 1
    else:
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
        hits = [node_hit(point, index) for index in range(4)]
        if sum(hits) != 1 or not hits[selected]:
            continue
        double_hits = [node_hit(2 * point, index) for index in range(4)]
        if not double_hits[selected]:
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
    "schema": "elkies-k3.q4o164-odd-i4-sections-modp.v1",
    "status": "PASS_MODP_Q4O164_ODD_I4_SECTION_SCAN",
    "prime": PRIME,
    "selected_fibre": {
        "index": selected,
        "label": labels[selected],
        "node": int(nodes[selected]),
    },
    "search": {
        "constrained_polynomial_x_tests": tests,
        "unique_filtered_sections": len(answers),
        "large_Groebner_required": False,
        "runtime_seconds": time.monotonic() - started,
    },
    "sections": answers,
    "proof_boundary": (
        "Exhaustive finite-field polynomial-square scan for exactly-one-node sections "
        "on an odd component of the selected I4 fibre. Exact QQ lifting, independence, "
        "and marked identification remain separate gates."
    ),
    "inputs": {
        "paths": [str(MODEL.relative_to(ROOT))],
        "sha256": {str(MODEL.relative_to(ROOT)): hashlib.sha256(MODEL.read_bytes()).hexdigest()},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "Q4O164ODDI4|prime={}|fibre={}|tests={}|candidates={}|status={}|output={}".format(
        PRIME, labels[selected], tests, len(answers), payload["status"], OUTPUT,
    ),
    flush=True,
)
