#!/usr/bin/env sage -python
"""Verify and package the exact pole-reduced H92 P2 lift."""

from sage.all import GF, QQ, PolynomialRing

import argparse
import hashlib
import json
from importlib.machinery import SourceFileLoader
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-p2-hensel.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-p2-lift.json"
ANCHOR = ROOT / "elkies-k3/scripts/verify_h3_noncm_q6_source_anchor.sage"
H92 = ROOT / "artifacts/local/humbert-inputs/92/igusa92.txt"
CANDIDATE = ROOT / "artifacts/generated-results/h92-p2-candidate-mod-100003-500.json"


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

record = json.loads(args.input.read_text())
assert record["complete"]
ring = PolynomialRing(QQ, "t")
t = ring.gen()
Z = ring([QQ(value) for value in record["Z"]])
X = ring([QQ(value) for value in record["X"]])
Y = ring([QQ(value) for value in record["Y"]])
assert Z.degree() == 21 and Z.leading_coefficient() == 1
assert X.degree() == 46 and Y.degree() == 69
assert X.gcd(Z) == 1 and Y.gcd(Z) == 1

anchor = SourceFileLoader("h92_p2_lift_anchor", str(ANCHOR)).load_module()
r, s = anchor.EXPECTED_H92
_, formulas = anchor.parse_h92(H92)
A1, A, B1, B, B2 = tuple(QQ(value(r, s)) for value in formulas)
a = A1*t**3 + A*t**4
b = B1*t**5 + B*t**6 + B2*t**7
assert Y**2 == X**3 + a*X*Z**4 + b*Z**6


def reduce_polynomial(value, field):
    return value.parent().change_ring(field)([
        field(coefficient.numerator()) / field(coefficient.denominator())
        for coefficient in value.list()
    ])


candidate = json.loads(CANDIDATE.read_text())
prime = candidate["prime"]
field = GF(prime)
mod_ring = PolynomialRing(field, "t")
tm = mod_ring.gen()
zm = reduce_polynomial(Z, field)
xm = reduce_polynomial(X, field)
ym = reduce_polynomial(Y, field)
am = field(A1.numerator()) / field(A1.denominator()) * tm**3 + field(A.numerator()) / field(A.denominator()) * tm**4
bm = (
    field(B1.numerator()) / field(B1.denominator()) * tm**5
    + field(B.numerator()) / field(B.denominator()) * tm**6
    + field(B2.numerator()) / field(B2.denominator()) * tm**7
)
matches = 0
opposite_matches = 0
for target, candidate_x, candidate_y in candidate["values"][:100]:
    value = field(target)
    curve = __import__("sage.all", fromlist=["EllipticCurve"]).EllipticCurve(
        [0, 0, 0, am(value), bm(value)]
    )
    point = curve(xm(value) / zm(value)**2, ym(value) / zm(value)**3)
    doubled = 2*point
    assert doubled.xy()[0] == field(candidate_x)
    if doubled.xy()[1] == field(candidate_y):
        matches += 1
    elif -doubled.xy()[1] == field(candidate_y):
        opposite_matches += 1
    else:
        raise AssertionError("candidate y-coordinate is not a doubled section value")
assert matches + opposite_matches == 100
orientation = "P2" if matches else "-P2"


payload = {
    "schema": "elkies-k3.h92-p2-exact-lift.v1",
    "status": "PASS_EXACT_H92_P2",
    "input": {
        "path": str(args.input.relative_to(ROOT)),
        "sha256": hashlib.sha256(args.input.read_bytes()).hexdigest(),
        "prime": record["prime"],
        "precision": record["precision"],
    },
    "coordinate": {
        "base": "t",
        "x": {"numerator": [str(value) for value in X.list()], "denominator": [str(value) for value in (Z**2).list()]},
        "y": {"numerator": [str(value) for value in Y.list()], "denominator": [str(value) for value in (Z**3).list()]},
        "Z": [str(value) for value in Z.list()],
        "degrees": {
            "Z": int(Z.degree()),
            "x": [int(X.degree()), int((Z**2).degree())],
            "y": [int(Y.degree()), int((Z**3).degree())],
        },
    },
    "exact_weierstrass_identity": True,
    "canonical_divisor_class": "3D-29H",
    "modular_double_check": {
        "candidate": str(CANDIDATE.relative_to(ROOT)),
        "prime": int(prime),
        "fibers": int(100),
        "x_matches": int(100),
        "y_matches": int(matches),
        "y_opposite_matches": int(opposite_matches),
        "orientation": orientation,
    },
}
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92P2VERIFY|x_degrees=46,42|y_degrees=69,63|"
    f"double_check_fibers=100|orientation={orientation}|status=PASS_EXACT",
    flush=True,
)
