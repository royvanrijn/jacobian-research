#!/usr/bin/env sage -python
"""Certify the cancellation-sensitive q8 frame at the H92 E7_2--E7_5 node.

Unlike the five unmarked nodes, ``x(P1)`` has the same leading order as x
here.  The exact H92 second-U chart nevertheless shows that x-x(P1) has the
factor Z^3 U^2 (Z-A1/B1) times a unit, while y-y(P1) has that factor times Y
times a unit.  Thus m is regular and the q6 module is still t*R at the *node*.
This must not be confused with the distinct marked smooth point -P1 on E7_5,
where the q6 frame is nontrivial.
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
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-e7-2-5-node-frame.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial(ring, coefficients):
    return ring([QQ(value) for value in coefficients])


def invert_base(rational_u):
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

chart = next(item for item in pullbacks["charts"] if item["name"] == "E7_2--E7_5")
edge = next(item for item in gluing["actual_edge_chart_gluing"] if item["name"] == "E7_2--E7_5")
assert edge["w_cartier_equation"] == "Z^6*Y^5"

ring = PolynomialRing(QQ, names=("Z", "U", "Y"))
Z, U, Y = ring.gens()
surface = ring(sage_eval(chart["surface_equation"], locals={"Z": Z, "U": U, "Y": Y}))
t_value = ring(sage_eval(chart["old_coordinate_pullback"]["t"], locals={"Z": Z, "U": U, "Y": Y}))
x_value = ring(sage_eval(chart["old_coordinate_pullback"]["x"], locals={"Z": Z, "U": U, "Y": Y}))
y_value = ring(sage_eval(chart["old_coordinate_pullback"]["y"], locals={"Z": Z, "U": U, "Y": Y}))
assert surface.derivative(U)(0, 0, 0) != 0
unit_h, remainder = (Y**2-surface).quo_rem(U)
assert not remainder and unit_h(0, 0, 0) != 0

u_ring = PolynomialRing(QQ, "u")
u_field = u_ring.fraction_field()
x_p = u_field(polynomial(u_ring, p1["x_entrance_base"]["numerator_coefficients"]))
x_p /= u_field(polynomial(u_ring, p1["x_entrance_base"]["denominator_coefficients"]))
y_p = u_field(polynomial(u_ring, p1["y_entrance_base"]["numerator_coefficients"]))
y_p /= u_field(polynomial(u_ring, p1["y_entrance_base"]["denominator_coefficients"]))
x_p_t = invert_base(x_p)
y_p_t = invert_base(y_p)
t = x_p_t.parent().gen()
assert x_p_t.valuation() == 2 and y_p_t.valuation() == 3
c2 = QQ((x_p_t/t**2)(0))
assert c2 != 0

# The actual chart writes C=Z-A1/B1.  The leading P1 cancellation is exact:
# x-c2*t^2=-c2*Z^3*U^2*C, because c2*(-A1/B1)=1.
linear_factor = t_value/(Z*U)
assert linear_factor(0, 0, 0) != 0
assert x_value == Z**2*U**2*linear_factor
assert c2*linear_factor(0, 0, 0) == 1
leading_difference = ring(x_value-c2*t_value**2)
assert leading_difference == -c2*Z**3*U**2*linear_factor
x_tail = x_p_t/t**2-c2
assert x_tail.valuation() >= 1
# x(P1)=c2*t^2+t^3*r(t), hence the displayed leading difference remains a
# unit after adding the tail: its correction has one extra U factor.
assert t_value**3/(Z**3*U**2*linear_factor) == U*linear_factor**2
assert y_value == Z**3*U**2*linear_factor*Y
assert t_value**3/(Z**3*U**2*linear_factor) == U*linear_factor**2

payload = {
    "schema": "elkies-k3.h92-q8-e7-2-5-node-frame.v1",
    "status": "PASS_EXACT_Q8_E7_2_5_NODE_FRAME",
    "inputs": {
        "p1": {"path": str(args.p1.relative_to(ROOT)), "sha256": digest(args.p1)},
        "actual_pullbacks": {"path": str(args.pullbacks.relative_to(ROOT)), "sha256": digest(args.pullbacks)},
        "q8_gluing": {"path": str(args.gluing.relative_to(ROOT)), "sha256": digest(args.gluing)},
    },
    "chart": {
        "name": chart["name"], "surface_equation": str(surface),
        "completed_relation": "U=Y^2/unit", "linear_factor": str(linear_factor),
        "pullbacks": {"t": str(t_value), "x": str(x_value), "y": str(y_value)},
    },
    "P1_cancellation": {
        "c2": str(c2), "leading_x_difference": str(leading_difference),
        "xP_tail": "x(P1)-c2*t^2 is divisible by t^3",
        "conclusion": "x-x(P1)=Z^3*U^2*linear_factor*unit",
    },
    "module_frame": {
        "y_difference": "y-y(P1)=Z^3*U^2*linear_factor*Y*unit",
        "m": "(y-y(P1))/(x-x(P1))=Y*unit",
        "q6_fractional_module": "t*(x-x(P1),y-y(P1))/(x-x(P1))=t*(1,m)=t*R",
        "q8_twist_cartier_equation": "g=Z^6*Y^5",
        "q8_condition": "Z^6*Y^5*f/t^9 belongs to R",
    },
    "boundary": (
        "This is the E7_2--E7_5 node frame only. The marked smooth point "
        "-P1 has its separate nontrivial frame, and no node quotient, overlap, "
        "global q8 matrix, or pencil is asserted."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "H92Q8E725NODEFRAME|q6_generator=t|q8_factor=Z6Y5|"
    "status=PASS_EXACT_Q8_E7_2_5_NODE_FRAME",
    flush=True,
)
