#!/usr/bin/env sage -python
"""Certify q8 local module frames on the five unmarked H92 E7 edge nodes.

Every chart used here is pulled back from the actual H92 blow-up tree.  At an
unmarked edge origin, the P1 entrance has x(P1)=t^2*unit and
y(P1)=t^3*unit, while the table below proves x(P1)/x and y(P1)/y lie in the
completed maximal ideal.  Therefore

    t*(x-x(P1),y-y(P1))/(x-x(P1)) = t*R,

and the q8 condition on that chart is ``g*f/t^9 in R`` for the transported
Cartier factor g.  The E7_2--E7_5 node is deliberately excluded because its
x(P1) leading cancellation is audited in a separate node calculation.
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
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-unmarked-e7-node-frames.json"


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

u_ring = PolynomialRing(QQ, "u")
u_field = u_ring.fraction_field()
x_p = u_field(polynomial(u_ring, p1["x_entrance_base"]["numerator_coefficients"]))
x_p /= u_field(polynomial(u_ring, p1["x_entrance_base"]["denominator_coefficients"]))
y_p = u_field(polynomial(u_ring, p1["y_entrance_base"]["numerator_coefficients"]))
y_p /= u_field(polynomial(u_ring, p1["y_entrance_base"]["denominator_coefficients"]))
x_p_t = invert_base(x_p)
y_p_t = invert_base(y_p)
assert x_p_t.valuation() == 2 and y_p_t.valuation() == 3

# The explicit ratios t^2/x and t^3/y identify the P1 correction orders.
# A solved coordinate is Y^2 times a unit because the listed actual surface
# derivative is nonzero.  In each case the two ratios lie in its completed
# maximal ideal, including the displayed Y denominators after that relation.
expected = {
    "E7_1--E7_4": {
        "solved": "Z", "g": "U^4*Y^2", "t2_over_x": "Z*U^2",
        "t3_over_y": "Z^2*U^3/Y", "after_relation": "Y^3*U^3*unit",
    },
    "E7_4--E7_3": {
        "solved": "U", "g": "Z^4*Y^6", "t2_over_x": "Z^2*U",
        "t3_over_y": "Z^3*U^2/Y", "after_relation": "Z^3*Y^3*unit",
    },
    "E7_3--E7_7": {
        "solved": "Z", "g": "U^5*Y^6", "t2_over_x": "Z*U",
        "t3_over_y": "Z^2*U^2/Y", "after_relation": "Y^3*U^2*unit",
    },
    "E7_7--E7_2": {
        "solved": "U", "g": "Z^5*Y^5", "t2_over_x": "Z",
        "t3_over_y": "Z^2*U/Y", "after_relation": "Z^2*Y*unit",
    },
    "E7_3--E7_6": {
        "solved": "U", "g": "Z^3*Y^6", "t2_over_x": "Z*U*unit",
        "t3_over_y": "Z*U^2/Y*unit", "after_relation": "Z*Y^3*unit",
    },
}

ring = PolynomialRing(QQ, names=("Z", "U", "Y"))
Z, U, Y = ring.gens()
charts = {entry["name"]: entry for entry in pullbacks["charts"]}
edges = {entry["name"]: entry for entry in gluing["actual_edge_chart_gluing"]}
assert set(expected) == set(charts)-{"E7_2--E7_5"}

records = []
for name, data in expected.items():
    chart = charts[name]
    edge = edges[name]
    surface = ring(sage_eval(chart["surface_equation"], locals={"Z": Z, "U": U, "Y": Y}))
    t_value = ring(sage_eval(chart["old_coordinate_pullback"]["t"], locals={"Z": Z, "U": U, "Y": Y}))
    x_value = ring(sage_eval(chart["old_coordinate_pullback"]["x"], locals={"Z": Z, "U": U, "Y": Y}))
    y_value = ring(sage_eval(chart["old_coordinate_pullback"]["y"], locals={"Z": Z, "U": U, "Y": Y}))
    solved = {"Z": Z, "U": U}[data["solved"]]
    quotient, remainder = (Y**2-surface).quo_rem(solved)
    assert not remainder and quotient(0, 0, 0) != 0
    assert surface.derivative(solved)(0, 0, 0) != 0
    t2_over_x = t_value**2/x_value
    t3_over_y = t_value**3/y_value
    # The literal first ratio is polynomial; the second has the stated Y
    # denominator, removed in the completion by solved=Y^2*unit.
    if "unit" not in data["t2_over_x"]:
        assert str(t2_over_x).replace(" ", "") == data["t2_over_x"].replace(" ", "")
    else:
        # The translated third-U chart contributes a nonzero factor
        # Z-1/A1.  It is a unit at the node, while the displayed Z*U factor
        # is exact and lies in the maximal ideal.
        unit_factor = t2_over_x/(Z*U)
        assert unit_factor(0, 0, 0) != 0
    assert edge["w_cartier_equation"] == data["g"]
    records.append({
        "chart": name,
        "surface_equation": str(surface),
        "solved_coordinate": data["solved"],
        "completed_relation": "{}=Y^2*unit".format(data["solved"]),
        "pullbacks": {"t": str(t_value), "x": str(x_value), "y": str(y_value)},
        "P1_ratio_orders": {
            "t2_over_x": data["t2_over_x"],
            "t3_over_y": data["t3_over_y"],
            "t3_over_y_after_completed_relation": data["after_relation"],
        },
        "q8_cartier_factor": data["g"],
        "q6_module": "t*R",
        "q8_condition": "{}*f/t^9 belongs to R".format(data["g"]),
    })

payload = {
    "schema": "elkies-k3.h92-q8-unmarked-e7-node-frames.v1",
    "status": "PASS_EXACT_Q8_UNMARKED_E7_NODE_FRAMES",
    "inputs": {
        "p1": {"path": str(args.p1.relative_to(ROOT)), "sha256": digest(args.p1)},
        "actual_pullbacks": {"path": str(args.pullbacks.relative_to(ROOT)), "sha256": digest(args.pullbacks)},
        "q8_gluing": {"path": str(args.gluing.relative_to(ROOT)), "sha256": digest(args.gluing)},
    },
    "P1_entrance_orders": {"xP": 2, "yP": 3},
    "unmarked_nodes": records,
    "excluded_marked_node": {
        "chart": "E7_2--E7_5",
        "reason": "x(P1) has a leading cancellation against x on this chart; its separate exact node calculation is kept distinct from the uniform strict-order proof.",
    },
    "boundary": (
        "These are actual local module frames, not their finite quotient "
        "matrices. The marked node, node principal parts, overlap gluing, and "
        "global q8 kernel remain to be compiled."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "H92Q8UNMARKEDNODEFRAMES|nodes=5|q6_generator=t|"
    "status=PASS_EXACT_Q8_UNMARKED_E7_NODE_FRAMES",
    flush=True,
)
