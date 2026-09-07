#!/usr/bin/env sage -python
"""Evaluate the third marked chord in its finite actual H92 E7 quotient.

The third divisor's E7 block has length 363.  Its marked point is kept as
the exact DAG ``22*(-P1)-P2`` and evaluated in QQ((t)); no global QQ(t)
normalization is introduced.  This script maps the resulting marked chord
into the actual translated quotient, verifies its denominator is invertible,
and checks the defining product identity exactly.
"""

import argparse
import hashlib
import json
from importlib.machinery import SourceFileLoader
from pathlib import Path

from sage.all import EllipticCurve, PowerSeriesRing, QQ


ROOT = Path(__file__).resolve().parents[2]
ANCHOR = ROOT / "elkies-k3/scripts/verify_h3_noncm_q6_source_anchor.sage"
H92 = ROOT / "artifacts/local/humbert-inputs/92/igusa92.txt"
P1 = ROOT / "artifacts/generated-results/elkies-k3-h92-p1-lift.json"
P2 = ROOT / "artifacts/generated-results/elkies-k3-h92-p2-hensel-100003-p1024.json"
QUOTIENT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-third-actual-e7-quotient-block.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-third-marked-chord-actual-e7-quotient.json"
P1_SHA256 = "c323bf6346bb239934a5a2d8b1a3f4067e70e993d2e4eb32aaa30f469fca6397"
P2_SHA256 = "e02e2803387d3a7f53907f548b275bb592d366f653f630f6ba8c9ef2611f3e37"
QUOTIENT_SHA256 = "7848c3a506b2255fc1e42cab9ced9b72f8852e26aa703f9b23e2b2417474d2ed"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def series_polynomial(ring, coefficients):
    t = ring.gen()
    return sum(ring(QQ(value))*t**index for index, value in enumerate(coefficients))


def reciprocal_series(ring, coefficients):
    t = ring.gen()
    return sum(ring(QQ(value))*t**(-index) for index, value in enumerate(coefficients))


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--precision", type=int, default=48)
parser.add_argument("--quotient", type=Path, default=QUOTIENT)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
if args.precision < 40:
    raise ValueError("precision must be at least 40")
for path, expected in ((P1, P1_SHA256), (P2, P2_SHA256)):
    assert digest(path) == expected
if args.quotient == QUOTIENT:
    assert digest(args.quotient) == QUOTIENT_SHA256
p1_data = json.loads(P1.read_text())
p2_data = json.loads(P2.read_text())
quotient = json.loads(args.quotient.read_text())
assert p1_data["status"] == "PASS_EXACT_H92_P1"
assert p2_data["complete"]
assert quotient["status"] == "PASS_EXACT_Q6_THIRD_ACTUAL_E7_QUOTIENT_BLOCK"
assert quotient["quotient_dimension"] == 363

anchor = SourceFileLoader("h92_third_chord_quotient_anchor", str(ANCHOR)).load_module()
r, s = anchor.EXPECTED_H92
_, formulas = anchor.parse_h92(H92)
A1, A, B1, B, B2 = (QQ(value(r, s)) for value in formulas)

series_ring = PowerSeriesRing(QQ, "t", default_prec=args.precision)
t = series_ring.gen()
field = series_ring.fraction_field()
curve = EllipticCurve(field, [0, 0, 0, A1*t**3+A*t**4, B1*t**5+B*t**6+B2*t**7])
xp1 = field(reciprocal_series(series_ring, p1_data["x_entrance_base"]["numerator_coefficients"]))
xp1 /= field(reciprocal_series(series_ring, p1_data["x_entrance_base"]["denominator_coefficients"]))
yp1 = field(reciprocal_series(series_ring, p1_data["y_entrance_base"]["numerator_coefficients"]))
yp1 /= field(reciprocal_series(series_ring, p1_data["y_entrance_base"]["denominator_coefficients"]))
z2 = field(series_polynomial(series_ring, p2_data["Z"]))
xp2 = field(series_polynomial(series_ring, p2_data["X"]))/z2**2
yp2 = field(series_polynomial(series_ring, p2_data["Y"]))/z2**3
point = 22*(-curve(xp1, yp1))+curve(xp2, yp2)
x_point, y_point = point.xy()
assert x_point.valuation() == y_point.valuation() == 0

generators = {tuple(entry) for entry in quotient["complete_ideal_generators"]}
basis = [tuple(entry) for entry in quotient["quotient_basis_exponents"]]
assert len(basis) == 363
c2 = QQ(quotient["translated_coordinates"]["c2"])
c3 = QQ(quotient["translated_coordinates"]["c3"])

# Sparse local quotient arithmetic in T,U,Y.  The exact H92 relation replaces
# Y^2 by its right-hand side; any monomial in the complete ideal is zero.
def add(left, right):
    answer = dict(left)
    for key, value in right.items():
        answer[key] = answer.get(key, QQ(0))+value
        if not answer[key]:
            del answer[key]
    return answer


def scale(value, scalar):
    return {key: QQ(scalar)*coefficient for key, coefficient in value.items() if scalar*coefficient}


def in_ideal(exponent):
    return any(all(left <= right for left, right in zip(generator, exponent)) for generator in generators)


x_branch = {(2, 0, 0): c2, (3, 0, 0): c3}
base_a = add({(3, 0, 0): A1, (4, 0, 0): A}, {})
base_b = {(5, 0, 0): B1, (6, 0, 0): B, (7, 0, 0): B2}


def multiply_raw(left, right):
    answer = {}
    for (i, a, b), coefficient in left.items():
        for (j, c, d), other in right.items():
            key = (i+j, a+c, b+d)
            answer[key] = answer.get(key, QQ(0))+coefficient*other
    return {key: value for key, value in answer.items() if value}


x_local = add({(0, 1, 0): QQ(1)}, x_branch)
rhs_y2 = add(multiply_raw(multiply_raw(x_local, x_local), x_local), multiply_raw(base_a, x_local))
rhs_y2 = add(rhs_y2, base_b)


def reduce(value):
    pending = list(value.items())
    answer = {}
    while pending:
        (i, a, b), coefficient = pending.pop()
        if not coefficient or in_ideal((i, a, b)):
            continue
        if b < 2:
            answer[(i, a, b)] = answer.get((i, a, b), QQ(0))+coefficient
            continue
        for (j, c, d), other in rhs_y2.items():
            pending.append(((i+j, a+c, b-2+d), coefficient*other))
    return {key: value for key, value in answer.items() if value}


def multiply(left, right):
    return reduce(multiply_raw(left, right))


def series_dict(value):
    return {(index, 0, 0): QQ(value[index]) for index in range(33) if value[index]}


# The generic ambient chord is m_{-P}=(y-y(-P))/(x-x(-P))=(Y+y(P))/(x-x(P)).
denominator = add(x_local, scale(series_dict(x_point), -1))
numerator = add({(0, 0, 1): QQ(1)}, series_dict(y_point))
d0 = denominator.get((0, 0, 0), QQ(0))
assert d0
nilpotent = add(scale(denominator, 1/d0), {(0, 0, 0): QQ(-1)})
inverse = {}
power = {(0, 0, 0): QQ(1)}
for index in range(128):
    inverse = add(inverse, scale(power, QQ((-1)**index)/d0))
    power = multiply(power, nilpotent)
    if not power:
        break
else:
    raise ArithmeticError("quotient denominator did not become nilpotent")
chord = multiply(numerator, inverse)
assert multiply(denominator, chord) == reduce(numerator)
coordinates = [str(chord.get(exponent, QQ(0))) for exponent in basis]
coordinate_sha256 = hashlib.sha256(json.dumps(coordinates, separators=(",", ":")).encode()).hexdigest()

payload = {
    "schema": "elkies-k3.h92-q6-third-marked-chord-actual-e7-quotient.v1",
    "status": "PASS_EXACT_Q6_THIRD_MARKED_CHORD_ACTUAL_E7_QUOTIENT",
    "inputs": {
        "p1": {"path": str(P1.relative_to(ROOT)), "sha256": digest(P1)},
        "p2": {"path": str(P2.relative_to(ROOT)), "sha256": digest(P2)},
        "actual_e7_quotient": {"path": str(args.quotient.relative_to(ROOT)), "sha256": digest(args.quotient)},
    },
    "precision": int(args.precision),
    "point_expression": "P=22*(-P1)-P2; chord=m_-P",
    "denominator_constant": str(d0),
    "inverse_nilpotence_steps": index+1,
    "chord_support": len(chord),
    "basis_exponents": [[int(value) for value in exponent] for exponent in basis],
    "coordinates_in_basis_order": coordinates,
    "coordinate_sha256": coordinate_sha256,
    "identity": "(x-x(P))*m_-P=y+y(P) in the exact length-363 actual E7 quotient",
    "boundary": "This evaluates the marked chord in the finite E7 quotient. Evaluating the full degree-44 global ambient and stacking all remaining local blocks is still required for a transported-section coordinate.",
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92Q6THIRDCHORDQUOTIENT|length=363|support={}|inverse_steps={}|"
    "status=PASS_EXACT_Q6_THIRD_MARKED_CHORD_ACTUAL_E7_QUOTIENT".format(len(chord), index+1),
    flush=True,
)
