#!/usr/bin/env sage -python
"""Certify the actual q8 local module frame at the E7_4--E7_3 node.

This is a chart-level prerequisite for a q8 node quotient, not a Kodaira
normal-form substitution.  On the actual H92 edge chart, the old coordinates
are

    t=Z^3 U^2,  x=Z^4 U^3,  y=Z^6 U^4 Y,

and the surface equation has the form ``Y^2-U*H(Z,U)=0`` with ``H`` a unit.
The P1 entrance functions have orders two and three in t.  Consequently
``x-x(P1)`` and ``y-y(P1)`` are respectively a unit times x and y in the
completed local ring.  The q6 marked module is therefore exactly ``t*R`` at
this node, so the q8 gluing condition is the actual regularity test

    (Z^4*Y^6)*f/t^9 in R.

It only supplies this one local trivialization.  Expanding its finite
principal-part quotient, and all other nodes/overlaps, remains separate.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import PolynomialRing, QQ, sage_eval


ROOT = Path(__file__).resolve().parents[2]
P1 = ROOT / "artifacts/generated-results/elkies-k3-h92-p1-lift.json"
PULLBACKS = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-actual-e7-chart-pullbacks.json"
GLUING = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-actual-e7-gluing.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-e7-4-3-node-frame.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial(ring, coefficients):
    return ring([QQ(value) for value in coefficients])


def invert_base(rational_u):
    """Rewrite QQ(u) in QQ(t), centered at the old E7 fibre t=0."""
    numerator = rational_u.numerator()
    denominator = rational_u.denominator()
    t_ring = PolynomialRing(QQ, "t")
    t = t_ring.gen()
    field = t_ring.fraction_field()
    return field(
        t**(denominator.degree()-numerator.degree())
        * t_ring(list(reversed(numerator.list())))
        / t_ring(list(reversed(denominator.list())))
    )


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--p1", type=Path, default=P1)
parser.add_argument("--pullbacks", type=Path, default=PULLBACKS)
parser.add_argument("--gluing", type=Path, default=GLUING)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

p1 = json.loads(args.p1.read_text())
pullbacks = json.loads(args.pullbacks.read_text())
gluing = json.loads(args.gluing.read_text())
assert p1["status"] == "PASS_EXACT_H92_P1"
assert pullbacks["status"] == "PASS_EXACT_H92_E7_CHART_PULLBACKS"
assert gluing["status"] == "PASS_EXACT_Q8_ACTUAL_E7_GLUING"

chart = next(item for item in pullbacks["charts"] if item["name"] == "E7_4--E7_3")
edge = next(item for item in gluing["actual_edge_chart_gluing"] if item["name"] == "E7_4--E7_3")
assert chart["surface_equation"] == edge["actual_h92_surface_equation"]
assert edge["w_cartier_equation"] == "Z^4*Y^6"

ring = PolynomialRing(QQ, names=("Z", "U", "Y"))
Z, U, Y = ring.gens()
surface = ring(sage_eval(chart["surface_equation"], locals={"Z": Z, "U": U, "Y": Y}))
t_value = ring(sage_eval(chart["old_coordinate_pullback"]["t"], locals={"Z": Z, "U": U, "Y": Y}))
x_value = ring(sage_eval(chart["old_coordinate_pullback"]["x"], locals={"Z": Z, "U": U, "Y": Y}))
y_value = ring(sage_eval(chart["old_coordinate_pullback"]["y"], locals={"Z": Z, "U": U, "Y": Y}))
assert (t_value, x_value, y_value) == (Z**3*U**2, Z**4*U**3, Z**6*U**4*Y)

# The actual surface, rather than a type-III* diagram, supplies the completed
# regular parameters.  Its equation is exactly Y^2-U*H with H(0)=1.
quotient, remainder = (Y**2-surface).quo_rem(U)
assert not remainder
unit_h = quotient
assert unit_h(0, 0, 0) == 1
assert surface == Y**2-U*unit_h
assert surface.derivative(U)(0, 0, 0) == -1

u_ring = PolynomialRing(QQ, "u")
u_field = u_ring.fraction_field()
x_p = u_field(polynomial(u_ring, p1["x_entrance_base"]["numerator_coefficients"]))
x_p /= u_field(polynomial(u_ring, p1["x_entrance_base"]["denominator_coefficients"]))
y_p = u_field(polynomial(u_ring, p1["y_entrance_base"]["numerator_coefficients"]))
y_p /= u_field(polynomial(u_ring, p1["y_entrance_base"]["denominator_coefficients"]))
x_p_t = invert_base(x_p)
y_p_t = invert_base(y_p)
t_ring = x_p_t.parent()
t = t_ring.gen()
assert x_p_t.valuation() == 2
assert y_p_t.valuation() == 3
x_entrance_unit = x_p_t/t**2
y_entrance_unit = y_p_t/t**3
assert x_entrance_unit(0) != 0
assert y_entrance_unit(0) != 0

# In the completed local ring R=QQ[[Z,Y]]/(surface), U=Y^2/H.  The following
# two ratios show that P1 corrections are in the maximal ideal after dividing
# by the leading old coordinate.  They are displayed as exact formal terms;
# H is a unit by the actual derivative check above.
x_p_over_x = "Z^2*U*(xP(t)/t^2)"
y_p_over_y = "Z^3*Y^3*(yP(t)/t^3)/H(Z,U)^2"
assert x_value == Z**4*U**3
assert t_value**2/x_value == Z**2*U
assert t_value**3/y_value == Z**3*U**2/Y

payload = {
    "schema": "elkies-k3.h92-q8-e7-4-3-node-frame.v1",
    "status": "PASS_EXACT_Q8_E7_4_3_NODE_FRAME",
    "inputs": {
        "p1": {"path": str(args.p1.relative_to(ROOT)), "sha256": digest(args.p1)},
        "actual_pullbacks": {"path": str(args.pullbacks.relative_to(ROOT)), "sha256": digest(args.pullbacks)},
        "q8_gluing": {"path": str(args.gluing.relative_to(ROOT)), "sha256": digest(args.gluing)},
    },
    "chart": {
        "name": chart["name"],
        "surface_equation": str(surface),
        "actual_unit_factor_H": str(unit_h),
        "completed_relation": "U=Y^2/H(Z,U), with H(0)=1",
        "old_pullbacks": {"t": str(t_value), "x": str(x_value), "y": str(y_value)},
    },
    "P1_entrance": {
        "ord_t_xP": 2,
        "ord_t_yP": 3,
        "xP_over_t2_at_zero": str(x_entrance_unit(0)),
        "yP_over_t3_at_zero": str(y_entrance_unit(0)),
    },
    "unit_comparisons": {
        "xP_over_x": x_p_over_x,
        "yP_over_y": y_p_over_y,
        "conclusion": "Both displayed ratios lie in the completed maximal ideal, so x-xP1=x*unit and y-yP1=y*unit.",
    },
    "module_frame": {
        "m": "(y-yP1)/(x-xP1)=Z^2*U*Y*unit",
        "q6_fractional_module": "t*(x-xP1,y-yP1)/(x-xP1)=t*(1,m)=t*R",
        "q8_twist_cartier_equation": "g=Z^4*Y^6",
        "q8_condition": "g*f belongs to (t*R)^9, equivalently g*f/t^9 belongs to R",
        "term_template": "u^i*x^a*m^b/h(u)^k maps to g*t^(4*k-i-9)*x^a*m^b/h_reverse(t)^k",
    },
    "boundary": (
        "This certifies the actual q6/q8 local module frame at one E7 node. "
        "It does not expand the finite two-variable principal-part quotient, "
        "cover the other nodes or overlaps, assemble the q8 matrix, or claim "
        "a pencil or child surface."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "H92Q8E743NODEFRAME|q6_generator=t|q8_factor=Z4Y6|"
    "status=PASS_EXACT_Q8_E7_4_3_NODE_FRAME",
    flush=True,
)
