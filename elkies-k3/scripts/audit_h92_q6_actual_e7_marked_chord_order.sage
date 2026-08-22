#!/usr/bin/env sage -python
"""Audit the raw q=6 marked chord on the actual E7_2--E7_5 chart.

The old marked chord is ``m=(y-y(P1))/(x-x(P1))``.  This script substitutes
the exact P1 Laurent functions and the actual H92 blow-up pullback

    t=Z*U*(Z-A1/B1),
    x=Z^2*U^2*(Z-A1/B1),
    y=Z^3*U^2*(Z-A1/B1)*Y.

It computes Z-orders in the rational chart function field, i.e. at the
generic point of E7_5 (Z=0).  This is intentionally a regression audit: it
does not infer a chart from the III* type and does not construct a quotient.
"""

import argparse
import hashlib
import json
from importlib.machinery import SourceFileLoader
from pathlib import Path

from sage.all import PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[2]
ANCHOR = ROOT / "elkies-k3/scripts/verify_h3_noncm_q6_source_anchor.sage"
H92 = ROOT / "artifacts/local/humbert-inputs/92/igusa92.txt"
P1 = ROOT / "artifacts/generated-results/elkies-k3-h92-p1-lift.json"
PULLBACKS = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-actual-e7-chart-pullbacks.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-actual-e7-marked-chord-order-audit.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial(ring, coefficients):
    return ring([QQ(value) for value in coefficients])


def invert_base(rational_u):
    """Write an exact QQ(u) function as a QQ(t) function for t=1/u."""
    numerator = rational_u.numerator()
    denominator = rational_u.denominator()
    t_ring = PolynomialRing(QQ, "t")
    t = t_ring.gen()
    t_field = t_ring.fraction_field()
    numerator_reverse = t_ring(list(reversed(numerator.list())))
    denominator_reverse = t_ring(list(reversed(denominator.list())))
    return t_field(
        t**(denominator.degree()-numerator.degree())
        * numerator_reverse/denominator_reverse
    )


def variable_order(value, ring, variable_index):
    value = ring.fraction_field()(value)

    def order(polynomial_value):
        polynomial_value = ring(polynomial_value)
        return min(exponent[variable_index] for exponent, coefficient in polynomial_value.dict().items() if coefficient)

    return order(value.numerator())-order(value.denominator())


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--p1", type=Path, default=P1)
parser.add_argument("--pullbacks", type=Path, default=PULLBACKS)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

p1 = json.loads(args.p1.read_text())
pullbacks = json.loads(args.pullbacks.read_text())
assert p1["status"] == "PASS_EXACT_H92_P1"
assert pullbacks["status"] == "PASS_EXACT_H92_E7_CHART_PULLBACKS"

anchor = SourceFileLoader("h92_q6_marked_chord_audit_anchor", str(ANCHOR)).load_module()
r, s = anchor.EXPECTED_H92
_, formulas = anchor.parse_h92(H92)
A1, A, B1, B, B2 = (QQ(value(r, s)) for value in formulas)

u_ring = PolynomialRing(QQ, "u")
u_field = u_ring.fraction_field()
x_p_u = u_field(polynomial(u_ring, p1["x_entrance_base"]["numerator_coefficients"]))
x_p_u /= u_field(polynomial(u_ring, p1["x_entrance_base"]["denominator_coefficients"]))
y_p_u = u_field(polynomial(u_ring, p1["y_entrance_base"]["numerator_coefficients"]))
y_p_u /= u_field(polynomial(u_ring, p1["y_entrance_base"]["denominator_coefficients"]))
x_p_t = invert_base(x_p_u)
y_p_t = invert_base(y_p_u)

ring = PolynomialRing(QQ, names=("Z", "U", "Y"))
Z, U, Y = ring.gens()
field = ring.fraction_field()
second = -A1/B1
t = field(Z*U*(Z+second))
x = field(Z**2*U**2*(Z+second))
y = field(Z**3*U**2*(Z+second)*Y)
x_p = field(x_p_t(t))
y_p = field(y_p_t(t))
chord = (y-y_p)/(x-x_p)

orders = {
    "t": variable_order(t, ring, 0),
    "x_minus_xP": variable_order(x-x_p, ring, 0),
    "y_minus_yP": variable_order(y-y_p, ring, 0),
    "m": variable_order(chord, ring, 0),
    "m_over_t": variable_order(chord/t, ring, 0),
}
assert orders == {"t": 1, "x_minus_xP": 3, "y_minus_yP": 3, "m": 0, "m_over_t": -1}

payload = {
    "schema": "elkies-k3.h92-q6-actual-e7-marked-chord-order-audit.v1",
    "status": "REJECTS_Q6_MARKED_E7_FRAME_AS_STATED",
    "inputs": {
        "h92_source": {"path": str(H92.relative_to(ROOT)), "sha256": digest(H92)},
        "p1": {"path": str(args.p1.relative_to(ROOT)), "sha256": digest(args.p1)},
        "actual_chart_pullbacks": {"path": str(args.pullbacks.relative_to(ROOT)), "sha256": digest(args.pullbacks)},
    },
    "chart": "E7_2--E7_5",
    "generic_component": "E7_5=(Z=0)",
    "raw_chord": "m=(y-y(P1))/(x-x(P1))",
    "exact_Z_orders": orders,
    "conclusion": (
        "The raw marked chord has ord_Z(m/t)=-1 at the generic point of E7_5. "
        "Therefore a frame asserting m/t is generically regular on this chart "
        "cannot be used without an additional Z-line-bundle trivialization."
    ),
    "boundary": (
        "This audits only the raw chord in one actual chart. It does not supply "
        "the corrected q6/q8 E7 module, a quotient condition, or a new pencil."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "H92Q6MARKEDE7AUDIT|ord_t=1|ord_dx=3|ord_dy=3|ord_m_over_t=-1|"
    "status=REJECTS_Q6_MARKED_E7_FRAME_AS_STATED",
    flush=True,
)
