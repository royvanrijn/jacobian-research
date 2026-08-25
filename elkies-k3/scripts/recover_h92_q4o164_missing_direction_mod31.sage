#!/usr/bin/env sage -python
"""Recover a low-pole q4/orbit164 section in the ninth MW direction mod 31.

The exact lattice shell contains a unique P.O=0 class with ninth MW
coordinate +1, one hit at the infinity I4 node, and q4/orbit1584 parent
degree three.  Impose the infinity-node value on a polynomial x-coordinate,
enumerate the remaining 31^4 possibilities, and retain polynomial squares on
the compact Weierstrass model with that exact parent-degree fingerprint.

This is a bounded univariate polynomial-square scan; no Groebner basis or
multivariate elimination is used.
"""

import argparse
import hashlib
import itertools
import json
import time
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, matrix


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
MODEL = LOCAL / "q4o164-compact-weierstrass-qq.json"
QUARTIC = LOCAL / "q4o164-compact-binary-quartic-qq.json"
POINTING = LOCAL / "q4o164-c8-equation-marking-qq.json"
INPUTS = (MODEL, QUARTIC, POINTING)

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, default=31)
parser.add_argument("--output", type=Path)
args = parser.parse_args()
PRIME = int(args.prime)
OUTPUT = args.output or LOCAL / f"q4o164-missing-direction-parent-degree3-mod{PRIME}.json"
if not OUTPUT.is_absolute():
    OUTPUT = ROOT / OUTPUT

started = time.monotonic()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


model = json.loads(MODEL.read_text())
quartic_data = json.loads(QUARTIC.read_text())
pointing = json.loads(POINTING.read_text())
assert model["status"] == "PASS_EXACT_QQ_Q4O164_COMPACT_WEIERSTRASS_NORMALIZATION"
assert quartic_data["status"] == "PASS_EXACT_QQ_Q4O164_COMPACT_BINARY_QUARTIC"
assert pointing["status"] == "PASS_EXACT_QQ_Q4O164_C8_EQUATION_MARKING"

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
assert A.degree() == 8 and B.degree() == 12

# The infinity I4 node is read in the compact line-bundle trivialization.
cubic_infinity = xvar**3 + A[8] * xvar + B[12]
repeated_infinity = cubic_infinity.gcd(cubic_infinity.derivative())
assert repeated_infinity.degree() == 1
node_infinity = -repeated_infinity[0] / repeated_infinity[1]

# Reconstruct the inverse pointed-quartic map directly on the doubly compact
# model.  Its translated quartic variable is the old q4/orbit1584 base V.
quartic = [
    R([reduce_qq(value) for value in row])
    for row in quartic_data["coefficients_in_T_low_to_high"]
]
e, d, c, b, a = map(K, quartic)
# The stored raw zero ordinate has bad reduction at 31, while its doubly
# normalized square e has good reduction.  Either square-root sign changes
# only the pointed orientation; the old-base degree fingerprint is invariant.
assert R(e).is_square()
W0 = K(R(e).sqrt())
a1 = d / W0
a2 = c - d**2 / (4 * W0**2)
a3 = 2 * W0 * b
b2 = a1**2 + 4 * a2


def inverse_parent_base(X, Y):
    x_general = K(X) / 9 - b2 / 12
    y_general = K(Y) / 27 - (a1 * x_general + a3) / 2
    if not y_general:
        return None
    V = K(2 * W0 * (x_general + a2) / y_general)
    W = K((x_general * V**2 - d * V) / (2 * W0) - W0)
    if W**2 != sum(K(value) * V**index for index, value in enumerate(quartic)):
        raise ArithmeticError("inverse pointed map missed the compact quartic")
    return V


def coefficient_jacobian(X, Y):
    dx = -3 * X**2 - A
    dy = 2 * Y
    return matrix(F, 13, 12, lambda degree, variable: (
        (dx * t**variable)[degree]
        if variable < 5
        else (dy * t**(variable - 5))[degree]
    ))


# Exact finite supports/nodes are retained only as a fingerprint filter.
supports = []
nodes = []
for fibre in model["compact_model"]["reducible_fibres"][:-1]:
    support = reduce_qq(fibre["support"])
    cubic = xvar**3 + A(support) * xvar + B(support)
    repeated = cubic.gcd(cubic.derivative())
    assert repeated.degree() == 1
    supports.append(support)
    nodes.append(-repeated[0] / repeated[1])


def node_hits(X, Y):
    return [
        X(support) == node and Y(support) == 0
        for support, node in zip(supports, nodes)
    ] + [X[4] == node_infinity and Y[6] == 0]


candidates = {}
tests = 0
for coefficients in itertools.product(F, repeat=4):
    tests += 1
    X = R(list(coefficients) + [node_infinity])
    rhs = X**3 + A * X + B
    if not rhs.is_square():
        continue
    root = rhs.sqrt()
    for Y in ({root, -root} if root else {root}):
        if Y.degree() > 6:
            continue
        hits = node_hits(X, Y)
        if hits != [False, False, False, True]:
            continue
        V = inverse_parent_base(X, Y)
        if V is None:
            continue
        parent_degree = max(V.numerator().degree(), V.denominator().degree())
        if parent_degree != 3:
            continue
        J = coefficient_jacobian(X, Y)
        # Fixing the resolved infinity-node branch contributes the two literal
        # incidence rows x_4=node and y_6=0.
        incidence = matrix(F, 2, 12)
        incidence[0, 4] = 1
        incidence[1, 11] = 1
        augmented = J.stack(incidence)
        if augmented.rank() != 12:
            continue
        key = (tuple(X.list()), tuple(Y.list()))
        candidates[key] = {
            "x_coefficients_low_to_high": [int(value) for value in X.list()],
            "y_coefficients_low_to_high": [int(value) for value in Y.list()],
            "node_hits": hits,
            "parent_base": {
                "numerator_coefficients_low_to_high": [int(value) for value in V.numerator().list()],
                "denominator_coefficients_low_to_high": [int(value) for value in V.denominator().list()],
                "degree": int(parent_degree),
            },
            "ordinary_weierstrass_jacobian_rank": int(J.rank()),
            "resolved_node_augmented_jacobian_rank": int(augmented.rank()),
        }

assert 1 <= len(candidates) <= 8
answers = list(candidates.values())
assert all(answer["resolved_node_augmented_jacobian_rank"] == 12 for answer in answers)

payload = {
    "schema": "elkies-k3.q4o164-missing-direction-parent-degree3-modp.v1",
    "status": "PASS_MODP_Q4O164_MISSING_DIRECTION_PARENT_DEGREE_THREE_CANDIDATES",
    "prime": PRIME,
    "search": {
        "constrained_polynomial_x_tests": tests,
        "unique_filtered_sections": len(candidates),
        "infinity_node_x_leading_coefficient": int(node_infinity),
        "large_Groebner_required": False,
        "runtime_seconds": time.monotonic() - started,
    },
    "sections": answers,
    "marked_fingerprint": {
        "P_dot_O": 0,
        "ninth_MW_coordinate": 1,
        "node_hits": [False, False, False, True],
        "q4o1584_parent_degree": 3,
        "uniqueness_source": "complete exact q4/o164 norm-four lattice shell",
    },
    "proof_boundary": (
        "The exhaustive mod-31 polynomial-square scan recovers every branch with the "
        "declared lattice fingerprint and a full-rank resolved node chart. Exact QQ lifting, "
        "marked identification, and the q8 horizontal group-law word remain separate gates."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in INPUTS},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "Q4O164MISSING31|tests={}|candidates={}|seeds={}|parent_degree=3|status={}|output={}".format(
        tests, len(answers), [
            [answer["x_coefficients_low_to_high"], answer["y_coefficients_low_to_high"]]
            for answer in answers
        ], payload["status"], OUTPUT,
    ),
    flush=True,
)
