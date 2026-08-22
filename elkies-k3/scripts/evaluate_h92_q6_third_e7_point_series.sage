#!/usr/bin/env sage -python
"""Evaluate P=22*(-P1)-P2 as an exact truncated point at the H92 E7 base.

The third marked point is evaluated from its expression DAG over the Laurent
series field QQ((t)), with t the exact finite H92 base parameter. This avoids
global rational-function normalization and provides the local data needed for
the actual E7 chart pullbacks.
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
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-third-e7-point-series.json"
P1_SHA256 = "c323bf6346bb239934a5a2d8b1a3f4067e70e993d2e4eb32aaa30f469fca6397"
P2_SHA256 = "e02e2803387d3a7f53907f548b275bb592d366f653f630f6ba8c9ef2611f3e37"

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def series_polynomial(ring, coefficients):
    t = ring.gen()
    return sum(ring(QQ(value)) * t**index for index, value in enumerate(coefficients))

def reciprocal_series(ring, coefficients):
    t = ring.gen()
    return sum(ring(QQ(value)) * t**(-index) for index, value in enumerate(coefficients))

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--precision", type=int, default=96)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
if args.precision < 32:
    raise ValueError("precision must be at least 32")
assert digest(P1) == P1_SHA256 and digest(P2) == P2_SHA256
p1_data = json.loads(P1.read_text())
p2_data = json.loads(P2.read_text())
assert p1_data["status"] == "PASS_EXACT_H92_P1"
assert p2_data["complete"] and p2_data["schema"] == "elkies-k3.h92-p2-hensel-lift.v1"

anchor = SourceFileLoader("h92_third_series_anchor", str(ANCHOR)).load_module()
r, s = anchor.EXPECTED_H92
_, formulas = anchor.parse_h92(H92)
A1, A, B1, B, B2 = (QQ(value(r, s)) for value in formulas)

series_ring = PowerSeriesRing(QQ, "t", default_prec=args.precision)
t = series_ring.gen()
field = series_ring.fraction_field()
old_a = field(A1 * t**3 + A * t**4)
old_b = field(B1 * t**5 + B * t**6 + B2 * t**7)
curve = EllipticCurve(field, [0, 0, 0, old_a, old_b])

# P1 is stored in u=1/t, while P2 is stored directly in t.
xp1 = field(reciprocal_series(series_ring, p1_data["x_entrance_base"]["numerator_coefficients"]))
xp1 /= field(reciprocal_series(series_ring, p1_data["x_entrance_base"]["denominator_coefficients"]))
yp1 = field(reciprocal_series(series_ring, p1_data["y_entrance_base"]["numerator_coefficients"]))
yp1 /= field(reciprocal_series(series_ring, p1_data["y_entrance_base"]["denominator_coefficients"]))
z2 = field(series_polynomial(series_ring, p2_data["Z"]))
xp2 = field(series_polynomial(series_ring, p2_data["X"])) / z2**2
yp2 = field(series_polynomial(series_ring, p2_data["Y"])) / z2**3
p1 = curve(xp1, yp1)
p2 = curve(xp2, yp2)
assert p1 in curve and p2 in curve

point = 22 * (-p1) + p2
assert not point.is_zero()
x_point, y_point = point.xy()
assert x_point.valuation() == 0 and y_point.valuation() == 0
assert x_point[0] and y_point[0] and y_point[0]**2 == x_point[0]**3

payload = {
    "schema": "elkies-k3.h92-q6-third-e7-point-series.v1",
    "status": "PASS_EXACT_Q6_THIRD_E7_SERIES_POINT",
    "precision": int(args.precision),
    "point_expression": "22*(-P1)-P2; reconstructed coordinate is -P2",
    "valuations": {"x": int(x_point.valuation()), "y": int(y_point.valuation())},
    "specialization": {
        "x_nonzero": bool(x_point[0]), "y_nonzero": bool(y_point[0]),
        "smooth_affine_cubic_check": bool(y_point[0]**2 == x_point[0]**3),
    },
    "compiler_instruction": "The chord at -P has numerator and denominator units at the E7 exceptional locus after the actual chart pullbacks; evaluate its higher jets locally from this DAG/series representation.",
    "boundary": "This is a local series evaluation of the marked point, not a complete chart quotient or a global child-section coordinate.",
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print("H92Q6THIRDSERIES|precision={}|valuations=0,0|status=PASS_EXACT_Q6_THIRD_E7_SERIES_POINT".format(args.precision), flush=True)
