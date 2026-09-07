#!/usr/bin/env sage -python
"""Exhaust q4/orbit164 zero-node polynomial sections modulo a small prime.

Enumerate all degree-at-most-four x polynomials, test whether x^3+A*x+B is
a polynomial square of degree at most six, and retain sections avoiding every
reducible-fibre node.  This is a bounded univariate square test (p^5 cases),
not a Groebner calculation.  The output is a finite-field construction aid;
QQ lifting and MW independence remain separate gates.
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

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, default=19)
parser.add_argument("--output", type=Path)
args = parser.parse_args()
PRIME = int(args.prime)
OUTPUT = args.output or LOCAL / f"q4o164-zero-node-sections-mod{PRIME}.json"
if not OUTPUT.is_absolute():
    OUTPUT = ROOT / OUTPUT

started = time.monotonic()
model = json.loads(MODEL.read_text())
assert model["status"] == "PASS_EXACT_QQ_Q4O164_COMPACT_WEIERSTRASS_NORMALIZATION"

F = GF(PRIME)
R = PolynomialRing(F, "t")
t = R.gen()
RX = PolynomialRing(F, "x")
xvar = RX.gen()


def reduce_qq(value):
    value = QQ(value)
    return F(value.numerator()) / F(value.denominator())


A = R([reduce_qq(value) for value in model["compact_model"]["A_coefficients_low_to_high"]])
B = R([reduce_qq(value) for value in model["compact_model"]["B_coefficients_low_to_high"]])
assert (A.degree(), B.degree()) == (8, 12)

supports = []
nodes = []
for fibre in model["compact_model"]["reducible_fibres"]:
    if fibre["support"] == "infinity":
        support = None
        cubic = xvar**3 + A[8] * xvar + B[12]
    else:
        support = reduce_qq(fibre["support"])
        cubic = xvar**3 + A(support) * xvar + B(support)
    repeated = cubic.gcd(cubic.derivative())
    assert repeated.degree() == 1
    supports.append(support)
    nodes.append(-repeated[0] / repeated[1])


def node_hits(X, Y):
    answers = []
    for support, node in zip(supports, nodes):
        answers.append(
            X[4] == node and Y[6] == 0
            if support is None else X(support) == node and Y(support) == 0
        )
    return answers


def coefficient_jacobian(X, Y):
    dx = -3 * X**2 - A
    dy = 2 * Y
    return matrix(F, 13, 12, lambda degree, variable: (
        (dx * t**variable)[degree]
        if variable < 5
        else (dy * t**(variable - 5))[degree]
    ))


candidates = {}
tests = 0
leading_rejected = 0
for coefficients in itertools.product(F, repeat=5):
    tests += 1
    leading_rhs = coefficients[4]**3 + A[8] * coefficients[4] + B[12]
    if not leading_rhs.is_square():
        leading_rejected += 1
        continue
    X = R(coefficients)
    rhs = X**3 + A * X + B
    if not rhs.is_square():
        continue
    root = rhs.sqrt()
    for Y in ({root, -root} if root else {root}):
        if Y.degree() > 6:
            continue
        hits = node_hits(X, Y)
        if any(hits):
            continue
        jacobian = coefficient_jacobian(X, Y)
        if jacobian.rank() != 12:
            continue
        key = (tuple(X.list()), tuple(Y.list()))
        candidates[key] = {
            "x_coefficients_low_to_high": [int(value) for value in X.list()],
            "y_coefficients_low_to_high": [int(value) for value in Y.list()],
            "node_hits": hits,
            "ordinary_weierstrass_jacobian_rank": int(jacobian.rank()),
        }

answers = list(candidates.values())
payload = {
    "schema": "elkies-k3.q4o164-zero-node-sections-modp.v1",
    "status": "PASS_MODP_Q4O164_ZERO_NODE_SECTION_SCAN",
    "prime": PRIME,
    "search": {
        "polynomial_x_tests": tests,
        "leading_coefficient_nonsquare_rejections": leading_rejected,
        "unique_regular_zero_node_sections": len(answers),
        "large_Groebner_required": False,
        "runtime_seconds": time.monotonic() - started,
    },
    "sections": answers,
    "proof_boundary": (
        "Exhaustive finite-field polynomial-square scan for regular P.O=0 sections "
        "avoiding every reducible node. QQ lifting and MW independence are separate gates."
    ),
    "inputs": {
        "paths": [str(MODEL.relative_to(ROOT))],
        "sha256": {str(MODEL.relative_to(ROOT)): hashlib.sha256(MODEL.read_bytes()).hexdigest()},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "Q4O164ZERONODESCAN|prime={}|tests={}|candidates={}|status={}|output={}".format(
        PRIME, tests, len(answers), payload["status"], OUTPUT,
    ),
    flush=True,
)
