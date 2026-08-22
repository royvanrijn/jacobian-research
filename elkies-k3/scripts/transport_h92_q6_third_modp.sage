#!/usr/bin/env sage -python
"""Diagnose naive old-group-law transport of the third H3 q=6 divisor.

This is not a point search.  The Weyl-aware divisor certificate fixes the
third source word as ``22*(-P1)-P2`` *after root and vertical corrections are
removed*.  This script evaluates only that generic-fibre word over
``GF(p)(u)``.  Its q=6 pencil degree is deliberately not one, proving that a
correct compiler must transport the full divisor class and its resolved
vertical correction rather than merely evaluate an MW word.
"""

import argparse
import json
from importlib.machinery import SourceFileLoader
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
ANCHOR = ROOT / "elkies-k3/scripts/verify_h3_noncm_q6_source_anchor.sage"
H92 = ROOT / "artifacts/local/humbert-inputs/92/igusa92.txt"
P1 = ROOT / "artifacts/generated-results/elkies-k3-h92-p1-lift.json"
P2 = ROOT / "artifacts/generated-results/elkies-k3-h92-p2-hensel-100003-p1024.json"
RR = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-global-rr.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-third-transport-mod-100003.json"


def coefficient(field, value):
    value = QQ(value)
    return field(ZZ(value.numerator())) / field(ZZ(value.denominator()))


def polynomial(ring, field, values):
    return ring([coefficient(field, value) for value in values])


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, default=100003)
parser.add_argument("--p1", type=Path, default=P1)
parser.add_argument("--p2", type=Path, default=P2)
parser.add_argument("--rr", type=Path, default=RR)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--minus-p1-multiple", type=int, default=22)
parser.add_argument(
    "--p2-coordinate-multiple", type=int, default=1,
    help="multiple of the reconstructed coordinate, whose H3-frame sign is -P2",
)
args = parser.parse_args()

prime = ZZ(args.prime)
if not prime.is_prime() or prime in (2, 3):
    raise ValueError("prime must be good and different from 2,3")
finite = GF(prime)
u_ring = PolynomialRing(finite, "u")
u = u_ring.gen()
field = u_ring.fraction_field()

anchor = SourceFileLoader("h92_q6_mod_transport_anchor", str(ANCHOR)).load_module()
r, s = anchor.EXPECTED_H92
_, formulas = anchor.parse_h92(H92)
A1, A, B1, B, B2 = (coefficient(finite, value(r, s)) for value in formulas)
old_a = field(A1) / u**3 + field(A) / u**4
old_b = field(B1) / u**5 + field(B) / u**6 + field(B2) / u**7
curve = EllipticCurve(field, [0, 0, 0, old_a, old_b])

p1_payload = json.loads(args.p1.read_text())
xp1 = field(polynomial(u_ring, finite, p1_payload["x_entrance_base"]["numerator_coefficients"])) / field(polynomial(u_ring, finite, p1_payload["x_entrance_base"]["denominator_coefficients"]))
yp1 = field(polynomial(u_ring, finite, p1_payload["y_entrance_base"]["numerator_coefficients"])) / field(polynomial(u_ring, finite, p1_payload["y_entrance_base"]["denominator_coefficients"]))
p1 = curve(xp1, yp1)

p2_payload = json.loads(args.p2.read_text())
if p2_payload["schema"] != "elkies-k3.h92-p2-hensel-lift.v1" or not p2_payload["complete"]:
    raise ValueError("P2 input must be a complete Hensel reconstruction")


def reciprocal_polynomial(values):
    return sum(coefficient(finite, value) / u**index for index, value in enumerate(values))


z2 = reciprocal_polynomial(p2_payload["Z"])
xp2 = reciprocal_polynomial(p2_payload["X"]) / z2**2
yp2 = reciprocal_polynomial(p2_payload["Y"]) / z2**3
p2 = curve(xp2, yp2)
assert p1 in curve and p2 in curve
print("H92Q6THIRDMOD|stage=inputs|prime={}".format(prime), flush=True)

# The default word comes from certify_h3_q6_weyl_section_transport.sage.
word_a = ZZ(args.minus_p1_multiple)
word_b = ZZ(args.p2_coordinate_multiple)
# The exact Hensel coordinate is -P2 in the pinned H3 frame: its q=6 pencil
# degree is 50, whereas the frame P2 class has D.P2=44.  Thus the source
# projection 22*(-P1)-P2 is 22*(-p1)+p2 in these coordinates.
third = word_a * (-p1) + word_b * p2
if third.is_zero():
    raise ArithmeticError("transported third word vanished modulo p")
x3, y3 = third.xy()
print(
    "H92Q6THIRDMOD|stage=word|x_degree={}|y_degree={}".format(
        max(x3.numerator().degree(), x3.denominator().degree()),
        max(y3.numerator().degree(), y3.denominator().degree()),
    ),
    flush=True,
)

rr = json.loads(args.rr.read_text())
h = polynomial(u_ring, finite, p1_payload["structured_denominator"]["Z4_coefficients"])


def coefficient_pair(entry):
    numerator_a = polynomial(u_ring, finite, entry["A_coefficients_low_to_high"])
    numerator_b = polynomial(u_ring, finite, entry["B_coefficients_low_to_high"])
    return field(numerator_a) / h**2, field(numerator_b) / h


(a0, b0), (a1, b1) = tuple(coefficient_pair(entry) for entry in rr["kernel"]["sections"])
if x3 == xp1 and y3 == -yp1:
    # The marked chord has its allowed pole at -P1.
    parameter = b1 / b0
elif x3 == xp1 and y3 == yp1:
    m = (3 * xp1**2 + old_a) / (2 * yp1)
    parameter = (a1 + b1 * m) / (a0 + b0 * m)
else:
    m = (y3 - yp1) / (x3 - xp1)
    parameter = (a1 + b1 * m) / (a0 + b0 * m)
parameter_degree = max(parameter.numerator().degree(), parameter.denominator().degree())
# This is a fixed regression value for the recorded H3 word at the chosen
# good prime.  Degree one here would signal that a root/vertical correction
# had accidentally been folded into the group-law calculation.
negative_guard = (word_a, word_b) == (22, 1)
if negative_guard and parameter_degree != 4769:
    raise ArithmeticError("unexpected naive-word pencil degree {}".format(parameter_degree))

payload = {
    "schema": "elkies-k3.h92-q6-third-transport-modp.v1",
    "status": (
        "PASS_DETECTS_MISSING_VERTICAL_CORRECTION"
        if negative_guard else "PASS_DETERMINISTIC_MW_PROJECTION_EVALUATION"
    ),
    "prime": int(prime),
    "coordinate_orientation_in_h3_frame": "reconstructed_p2=-P2",
    "word": "{}*(-P1)-P2".format(word_a) if word_b == 1 else "{}*(-P1){}{}*(reconstructed_-P2)".format(word_a, "+" if word_b >= 0 else "", word_b),
    "old_word_coordinate_degrees": {
        "x": [int(x3.numerator().degree()), int(x3.denominator().degree())],
        "y": [int(y3.numerator().degree()), int(y3.denominator().degree())],
    },
    "new_parameter": {
        "numerator_coefficients_low_to_high": [int(value) for value in parameter.numerator().list()],
        "denominator_coefficients_low_to_high": [int(value) for value in parameter.denominator().list()],
        "degree": int(parameter_degree),
    },
    "boundary": (
        "The old MW word omits its resolved vertical correction: its q=6 pencil degree is 4769 rather than one. This artifact is a negative regression guard, not a child-section reconstruction."
        if negative_guard else
        "This is a deterministic generic-fibre MW evaluation only; it does not identify a transported child divisor."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92Q6THIRDMOD|parameter_degree={}|status={}".format(
        parameter_degree,
        "PASS_DETECTS_MISSING_VERTICAL_CORRECTION" if negative_guard else "PASS_DETERMINISTIC_MW_PROJECTION_EVALUATION",
    ),
    flush=True,
)
