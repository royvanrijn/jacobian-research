#!/usr/bin/env sage -python
"""Derive the principal-frame clearing identity at the H92 E7_4--E7_3 node.

The actual node frame has q6 module t*R, hence q8 regularity is g*f/t^9 in R.
This script rewrites the marked chord in the completed chart as

    m = Z^2*U*Y * B/A,

where A and B are explicit units.  It then clears only these units and the
reversed h denominator.  A finite endpoint linear combination is regular at
the node exactly when its displayed polynomial numerator is divisible by one
power t^T in the actual local ring.  This is the structured input for the
finite two-variable quotient; no ten-generator singular-model ideal is used.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import PolynomialRing, QQ, sage_eval


ROOT = Path(__file__).resolve().parents[2]
P1 = ROOT / "artifacts/generated-results/elkies-k3-h92-p1-lift.json"
PULLBACKS = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-actual-e7-chart-pullbacks.json"
AMBIENT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-endpoint-rr-ambient.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-e7-4-3-principal-node-clearing.json"


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
parser.add_argument("--ambient", type=Path, default=AMBIENT)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

p1 = json.loads(args.p1.read_text())
pullbacks = json.loads(args.pullbacks.read_text())
ambient = json.loads(args.ambient.read_text())
assert p1["status"] == "PASS_EXACT_H92_P1"
assert pullbacks["status"] == "PASS_EXACT_H92_E7_CHART_PULLBACKS"
assert ambient["status"] in (
    "PASS_EXACT_Q8_ENDPOINT_RR_AMBIENT", "PASS_EXACT_Q8_ENLARGED_ENDPOINT_RR_AMBIENT",
)

chart = next(item for item in pullbacks["charts"] if item["name"] == "E7_4--E7_3")
ring = PolynomialRing(QQ, names=("Z", "U", "Y"))
Z, U, Y = ring.gens()
surface = ring(sage_eval(chart["surface_equation"], locals={"Z": Z, "U": U, "Y": Y}))
t_value = ring(sage_eval(chart["old_coordinate_pullback"]["t"], locals={"Z": Z, "U": U, "Y": Y}))
x_value = ring(sage_eval(chart["old_coordinate_pullback"]["x"], locals={"Z": Z, "U": U, "Y": Y}))
y_value = ring(sage_eval(chart["old_coordinate_pullback"]["y"], locals={"Z": Z, "U": U, "Y": Y}))
assert (t_value, x_value, y_value) == (Z**3*U**2, Z**4*U**3, Z**6*U**4*Y)
unit_h, remainder = (Y**2-surface).quo_rem(U)
assert not remainder and unit_h(0, 0, 0) == 1

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
r = x_p_t/t**2
s = y_p_t/t**3
assert r.valuation() == 0 and s.valuation() == 0
r_num, r_den = t_ring(r.numerator()), t_ring(r.denominator())
s_num, s_den = t_ring(s.numerator()), t_ring(s.denominator())
assert r_den(0) and s_den(0)
h_u = polynomial(u_ring, p1["structured_denominator"]["Z4_coefficients"])
h_reverse = t_ring(list(reversed(h_u.list())))
assert h_reverse(0) != 0

r_num_chart, r_den_chart = ring(r_num(t_value)), ring(r_den(t_value))
s_num_chart, s_den_chart = ring(s_num(t_value)), ring(s_den(t_value))
h_chart = ring(h_reverse(t_value))
# x(P1)/x=Z^2*U*r and y(P1)/y=Z^3*U*Y*s/H, using U*H=Y^2.
a_numerator = r_den_chart-Z**2*U*r_num_chart
b_numerator = unit_h*s_den_chart-Z**3*U*Y*s_num_chart
assert a_numerator(0, 0, 0) != 0
assert b_numerator(0, 0, 0) != 0
assert h_chart(0, 0, 0) != 0

max_h_power = max(int(entry["h_power"]) for entry in ambient["ambient_basis"])
common_t_power = 9+max(
    int(entry["u_power"])-4*int(entry["h_power"])
    for entry in ambient["ambient_basis"]
)
assert common_t_power >= 9
term_records = []
for index, entry in enumerate(ambient["ambient_basis"]):
    a, b = int(entry["x_power"]), int(entry["m_power"])
    i, k = int(entry["u_power"]), int(entry["h_power"])
    t_exponent = common_t_power+4*k-i-9
    assert t_exponent >= 0
    term_records.append({
        "basis_index": index,
        "t_exponent_after_common_clearing": t_exponent,
        "numerator_template": (
            "(Z^4*Y^6)*t^{}*x^{}*(Z^2*U*Y)^{}*B^{}*r_den^{}*"
            "H^{}*s_den^{}*A^{}*h_reverse^{}".format(
                t_exponent, a, b, b, b, 9-b, 9-b, 9-b, max_h_power-k
            )
        ),
    })

payload = {
    "schema": "elkies-k3.h92-q8-e7-4-3-principal-node-clearing.v1",
    "status": "PASS_EXACT_Q8_E7_4_3_PRINCIPAL_NODE_CLEARING",
    "inputs": {
        "p1": {"path": str(args.p1.relative_to(ROOT)), "sha256": digest(args.p1)},
        "actual_pullbacks": {"path": str(args.pullbacks.relative_to(ROOT)), "sha256": digest(args.pullbacks)},
        "endpoint_ambient": {"path": str(args.ambient.relative_to(ROOT)), "sha256": digest(args.ambient)},
    },
    "chart": {
        "name": "E7_4--E7_3", "surface_equation": str(surface),
        "completed_relation": "U=Y^2/H with H(0)=1", "H": str(unit_h),
        "pullbacks": {"t": str(t_value), "x": str(x_value), "y": str(y_value)},
    },
    "unit_chord_factors": {
        "A": "1-Z^2*U*(xP/t^2)=A_num/r_den",
        "A_numerator": str(a_numerator),
        "B": "1-Z^3*U*Y*(yP/t^3)/H=B_num/(H*s_den)",
        "B_numerator": str(b_numerator),
        "m": "Z^2*U*Y*B_num*r_den/(H*s_den*A_num)",
        "all_denominators_are_units_at_node": True,
    },
    "common_clearing": {
        "T": common_t_power,
        "K": max_h_power,
        "unit_multiplier": "A_num^9*H^9*s_den^9*h_reverse^K",
        "equivalence": (
            "For a coefficient combination of the displayed term numerators, "
            "g*f/t^9 is regular iff its common-cleared numerator belongs to t^T*R."
        ),
    },
    "term_numerator_templates": term_records,
    "boundary": (
        "This is an exact principal-frame clearing identity for one actual "
        "node. It does not yet calculate the finite quotient R/(t^T), process "
        "the other nodes or overlaps, or produce a q8 kernel or pencil."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "H92Q8E743PRINCIPALCLEARING|ambient={}|T={}|K={}|"
    "status=PASS_EXACT_Q8_E7_4_3_PRINCIPAL_NODE_CLEARING".format(
        ambient["ambient_dimension"], common_t_power, max_h_power,
    ),
    flush=True,
)
