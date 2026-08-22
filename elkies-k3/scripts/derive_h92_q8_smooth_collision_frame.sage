#!/usr/bin/env sage -python
"""Derive the regular q=8 chord frame at the four smooth P1.O collisions.

Set ``p=y(P1)/x(P1)``, ``q=(m-p)/h`` and ``X=h^2*x``.  The chord relation
``y=y(P1)+m*(x-x(P1))`` and the old Weierstrass equation have the marked
factor ``x-x(P1)``.  After substituting ``m=p+h*q`` and ``x=X/h^2``, clearing
``h^6``, and removing that marked factor, the residual relation is a monic
quadratic in X.  Its coefficients are regular at h=0.

Thus ``q`` and ``X`` give the actual smooth-collision algebra in which the
degree-18 generic basis is represented by

    1,q,...,q^9, X,Xq,...,Xq^7.

This script establishes the exact regular algebra/frame needed for a q=8
smooth quotient.  Identifying the precise q=8 line-bundle lattice inside it
and forming its principal-part condition matrix remain separate steps.
"""

import argparse
import hashlib
import json
from importlib.machinery import SourceFileLoader
from pathlib import Path

from sage.all import PolynomialRing, QQ, gcd


ROOT = Path(__file__).resolve().parents[2]
ANCHOR = ROOT / "elkies-k3/scripts/verify_h3_noncm_q6_source_anchor.sage"
H92 = ROOT / "artifacts/local/humbert-inputs/92/igusa92.txt"
P1 = ROOT / "artifacts/generated-results/elkies-k3-h92-p1-lift.json"
SMOOTH = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-smooth-po-module.json"
GENERIC = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-generic-rr-ambient.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-smooth-collision-frame.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial(ring, coefficients):
    return ring([QQ(value) for value in coefficients])


def serialize_q_polynomial(value):
    return [str(coefficient) for coefficient in value.list()]


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--p1", type=Path, default=P1)
parser.add_argument("--smooth", type=Path, default=SMOOTH)
parser.add_argument("--generic", type=Path, default=GENERIC)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

p1 = json.loads(args.p1.read_text())
smooth = json.loads(args.smooth.read_text())
generic = json.loads(args.generic.read_text())
assert p1["status"] == "PASS_EXACT_H92_P1"
assert smooth["status"] == "PASS_EXACT_SMOOTH_PO_CHORD_MODULE"
assert generic["status"] == "PASS_EXACT_Q8_GENERIC_RR_AMBIENT"
assert generic["dimension"] == 18

anchor = SourceFileLoader("h92_q8_smooth_frame_anchor", str(ANCHOR)).load_module()
r, s = anchor.EXPECTED_H92
_, formulas = anchor.parse_h92(H92)
A1, A, B1, B, B2 = (QQ(value(r, s)) for value in formulas)

u_ring = PolynomialRing(QQ, "u")
u = u_ring.gen()
u_field = u_ring.fraction_field()
h = polynomial(u_ring, p1["structured_denominator"]["Z4_coefficients"])
assert h.degree() == 4 and gcd(h, h.derivative()) == 1
x_p = u_field(polynomial(u_ring, p1["x_entrance_base"]["numerator_coefficients"]))
x_p /= u_field(polynomial(u_ring, p1["x_entrance_base"]["denominator_coefficients"]))
y_p = u_field(polynomial(u_ring, p1["y_entrance_base"]["numerator_coefficients"]))
y_p /= u_field(polynomial(u_ring, p1["y_entrance_base"]["denominator_coefficients"]))
old_a = A1/u**3+A/u**4
old_b = B1/u**5+B/u**6+B2/u**7
assert y_p**2 == x_p**3+old_a*x_p+old_b

q_ring = PolynomialRing(u_field, "q")
q = q_ring.gen()
x_ring = PolynomialRing(q_ring, "X")
X = x_ring.gen()
p = y_p/x_p
m = q_ring(p)+q_ring(h)*q
x = x_ring(X)/x_ring(h**2)
y = x_ring(y_p)+x_ring(m)*(x-x_ring(x_p))
relation = y**2-x**3-x_ring(old_a)*x-x_ring(old_b)

# h^6 clears the x substitution.  The marked x=x(P1) factor is then
# removed exactly, leaving the quadratic relation in the collision frame.
scaled = x_ring(h**6*relation)
marked_factor = X-x_ring(h**2*x_p)
quadratic, remainder = scaled.quo_rem(marked_factor)
assert not remainder
assert quadratic.degree() == 2 and quadratic.leading_coefficient() == -1
quadratic = x_ring(-quadratic)
assert quadratic.leading_coefficient() == 1

regularity = []
for x_power, coefficient in enumerate(quadratic.list()):
    coefficient = q_ring(coefficient)
    q_data = []
    for q_power, value in enumerate(coefficient.list()):
        value = u_field(value)
        numerator = u_ring(value.numerator())
        denominator = u_ring(value.denominator())
        assert gcd(h, denominator) == 1
        valuation = "infinity" if numerator == 0 else int(numerator.valuation(h))
        q_data.append({
            "q_power": q_power,
            "numerator_h_valuation": valuation,
            "denominator_h_unit": True,
        })
    regularity.append({"X_power": x_power, "q_coefficients": q_data})

frame = (
    [{"X_power": 0, "q_power": power} for power in range(10)]
    + [{"X_power": 1, "q_power": power} for power in range(8)]
)
assert len(frame) == 18

payload = {
    "schema": "elkies-k3.h92-q8-smooth-collision-frame.v1",
    "status": "PASS_EXACT_Q8_SMOOTH_COLLISION_FRAME",
    "inputs": {
        "h92_source": {"path": str(H92.relative_to(ROOT)), "sha256": digest(H92)},
        "p1": {"path": str(args.p1.relative_to(ROOT)), "sha256": digest(args.p1)},
        "q6_smooth_module": {"path": str(args.smooth.relative_to(ROOT)), "sha256": digest(args.smooth)},
        "q8_generic_ambient": {"path": str(args.generic.relative_to(ROOT)), "sha256": digest(args.generic)},
    },
    "coordinates": {
        "q": "(m-y(P1)/x(P1))/h",
        "X": "h^2*x",
        "collision_polynomial": str(h),
    },
    "quadratic_relation": {
        "equation": "X^2 + alpha(q)*X + beta(q)=0",
        "alpha_coefficients_low_to_high": serialize_q_polynomial(quadratic[1]),
        "beta_coefficients_low_to_high": serialize_q_polynomial(quadratic[0]),
        "regularity_at_h": regularity,
        "derivation": "h^6*Weierstrass(y(P1)+m*(x-x(P1)),X/h^2)/(X-h^2*x(P1))",
    },
    "regular_degree_18_frame": frame,
    "conclusion": (
        "The marked chord quadratic becomes monic in X=h^2*x with all "
        "coefficients regular at the four smooth h=0 fibres. The listed 18 "
        "q/X monomials are the resulting exact local algebra frame."
    ),
    "boundary": (
        "This derives the regular collision algebra but does not identify the "
        "complete q8 line-bundle lattice, form a principal-part condition "
        "matrix, or claim a global kernel or pencil."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "H92Q8SMOOTHFRAME|quadratic_in_X=1|frame=18|h_regular=1|"
    "status=PASS_EXACT_Q8_SMOOTH_COLLISION_FRAME",
    flush=True,
)
