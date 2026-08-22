#!/usr/bin/env sage -python
"""Derive exact principal clearing identities on all six H92 E7 nodes.

At every resolved E7 edge node the q6 marked module has the local frame
``t*R`` (with the cancellation-sensitive E7_2--E7_5 case certified
separately).  Thus the q8 condition is ``g*f/t^9 in R``.  This script derives
the chord in each *actual* blow-up chart, without guessing it from component
labels, and clears its unit denominators uniformly for a finite q8 ambient.

It only packages the six exact local principal conditions.  Their finite
images, their overlap compatibility, and their common kernel are deliberately
left to the compiler's local-condition evaluator.
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
FRAMES = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-unmarked-e7-node-frames.json"
MARKED_FRAME = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-e7-2-5-node-frame.json"
AMBIENT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-endpoint-rr-ambient.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-e7-node-principal-clearings.json"


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


def common_monomial_exponents(value):
    """Return the maximal monomial dividing a nonzero Z,U,Y polynomial."""
    terms = list(value.dict())
    assert terms
    return tuple(min(exponent[index] for exponent in terms) for index in range(3))


def monomial_string(exponents):
    names = ("Z", "U", "Y")
    pieces = []
    for name, exponent in zip(names, exponents):
        if exponent == 1:
            pieces.append(name)
        elif exponent:
            pieces.append("{}^{}".format(name, exponent))
    return "*".join(pieces) if pieces else "1"


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--p1", type=Path, default=P1)
parser.add_argument("--pullbacks", type=Path, default=PULLBACKS)
parser.add_argument("--gluing", type=Path, default=GLUING)
parser.add_argument("--frames", type=Path, default=FRAMES)
parser.add_argument("--marked-frame", type=Path, default=MARKED_FRAME)
parser.add_argument("--ambient", type=Path, default=AMBIENT)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

p1 = json.loads(args.p1.read_text())
pullbacks = json.loads(args.pullbacks.read_text())
gluing = json.loads(args.gluing.read_text())
frames = json.loads(args.frames.read_text())
marked_frame = json.loads(args.marked_frame.read_text())
ambient = json.loads(args.ambient.read_text())
assert p1["status"] == "PASS_EXACT_H92_P1"
assert pullbacks["status"] == "PASS_EXACT_H92_E7_CHART_PULLBACKS"
assert gluing["status"] == "PASS_EXACT_Q8_ACTUAL_E7_GLUING"
assert frames["status"] == "PASS_EXACT_Q8_UNMARKED_E7_NODE_FRAMES"
assert marked_frame["status"] == "PASS_EXACT_Q8_E7_2_5_NODE_FRAME"
assert ambient["status"] in {
    "PASS_EXACT_Q8_ENDPOINT_RR_AMBIENT",
    "PASS_EXACT_Q8_ENLARGED_ENDPOINT_RR_AMBIENT",
}

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
r = x_p_t / t**2
s = y_p_t / t**3
assert r.valuation() == 0 and s.valuation() == 0
r_num, r_den = t_ring(r.numerator()), t_ring(r.denominator())
s_num, s_den = t_ring(s.numerator()), t_ring(s.denominator())
assert r_den(0) and s_den(0)
h_reverse = t_ring(list(reversed(polynomial(
    u_ring, p1["structured_denominator"]["Z4_coefficients"]
).list())))
assert h_reverse(0)

K = max(int(entry["h_power"]) for entry in ambient["ambient_basis"])
T = 9 + max(
    int(entry["u_power"])-4*int(entry["h_power"])
    for entry in ambient["ambient_basis"]
)
assert T >= 9

ring = PolynomialRing(QQ, names=("Z", "U", "Y"))
Z, U, Y = ring.gens()
chart_by_name = {entry["name"]: entry for entry in pullbacks["charts"]}
edge_by_name = {entry["name"]: entry for entry in gluing["actual_edge_chart_gluing"]}
assert set(chart_by_name) == set(edge_by_name)

records = []
for name in sorted(chart_by_name):
    chart = chart_by_name[name]
    edge = edge_by_name[name]
    surface = ring(sage_eval(chart["surface_equation"], locals={"Z": Z, "U": U, "Y": Y}))
    t_value = ring(sage_eval(chart["old_coordinate_pullback"]["t"], locals={"Z": Z, "U": U, "Y": Y}))
    x_value = ring(sage_eval(chart["old_coordinate_pullback"]["x"], locals={"Z": Z, "U": U, "Y": Y}))
    y_value = ring(sage_eval(chart["old_coordinate_pullback"]["y"], locals={"Z": Z, "U": U, "Y": Y}))
    # Lx=(Mx*A)/r_den and Ly=(My*B)/s_den exactly.  A, r_den,
    # s_den are units; B is allowed to vanish and records the local chord.
    numerator_x = x_value*ring(r_den(t_value))-t_value**2*ring(r_num(t_value))
    numerator_y = y_value*ring(s_den(t_value))-t_value**3*ring(s_num(t_value))
    mx = common_monomial_exponents(numerator_x)
    my = common_monomial_exponents(numerator_y)
    m_exponents = tuple(my[index]-mx[index] for index in range(3))
    assert all(exponent >= 0 for exponent in m_exponents)
    mx_value = ring.monomial(*mx)
    my_value = ring.monomial(*my)
    A = numerator_x // mx_value
    B = numerator_y // my_value
    assert numerator_x == mx_value*A and numerator_y == my_value*B
    assert A(0, 0, 0) != 0
    assert ring(r_den(t_value))(0, 0, 0) != 0
    assert ring(s_den(t_value))(0, 0, 0) != 0
    assert ring(h_reverse(t_value))(0, 0, 0) != 0
    term_records = []
    for index, entry in enumerate(ambient["ambient_basis"]):
        a, b = int(entry["x_power"]), int(entry["m_power"])
        i, k = int(entry["u_power"]), int(entry["h_power"])
        exponent = T+4*k-i-9
        assert exponent >= 0 and 0 <= b <= 9
        term_records.append({
            "basis_index": index,
            "t_exponent_after_common_clearing": exponent,
            "numerator_template": (
                "({g})*t^{e}*x^{a}*({m})^{b}*B^{b}*r_den^{b}*"
                "A^{ab}*s_den^{ab}*h_reverse^{hk}"
            ).format(
                g=edge["w_cartier_equation"], e=exponent, a=a,
                m=monomial_string(m_exponents), b=b, ab=9-b, hk=K-k,
            ),
        })
    records.append({
        "chart": name,
        "surface_equation": str(surface),
        "pullbacks": {"t": str(t_value), "x": str(x_value), "y": str(y_value)},
        "q8_cartier_factor": edge["w_cartier_equation"],
        "chord": {
            "Lx": "({})*A/r_den".format(monomial_string(mx)),
            "Ly": "({})*B/s_den".format(monomial_string(my)),
            "m": "({})*B*r_den/(A*s_den)".format(monomial_string(m_exponents)),
            "A_is_unit": True,
            "r_den_is_unit": True,
            "s_den_is_unit": True,
        },
        "common_clearing": {
            "unit_multiplier": "A^9*s_den^9*h_reverse^{}".format(K),
            "equivalence": (
                "A combination of the displayed term numerators is in "
                "(t^{}) in this actual chart iff g*f/t^9 is regular."
            ).format(T),
        },
        "term_numerator_templates": term_records,
    })

payload = {
    "schema": "elkies-k3.h92-q8-e7-node-principal-clearings.v1",
    "status": "PASS_EXACT_Q8_E7_NODE_PRINCIPAL_CLEARINGS",
    "inputs": {
        "p1": {"path": str(args.p1.relative_to(ROOT)), "sha256": digest(args.p1)},
        "actual_pullbacks": {"path": str(args.pullbacks.relative_to(ROOT)), "sha256": digest(args.pullbacks)},
        "q8_gluing": {"path": str(args.gluing.relative_to(ROOT)), "sha256": digest(args.gluing)},
        "unmarked_node_frames": {"path": str(args.frames.relative_to(ROOT)), "sha256": digest(args.frames)},
        "marked_node_frame": {"path": str(args.marked_frame.relative_to(ROOT)), "sha256": digest(args.marked_frame)},
        "endpoint_ambient": {"path": str(args.ambient.relative_to(ROOT)), "sha256": digest(args.ambient)},
    },
    "common_parameters": {"T": T, "K": K, "ambient_dimension": len(ambient["ambient_basis"])},
    "nodes": records,
    "boundary": (
        "These six identities are exact local chart conditions. They do not "
        "supply their local quotient matrices, Cech overlap relations, a "
        "characteristic-zero q8 kernel, h0(D), a pencil, or a child model."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "H92Q8E7NODECLEARINGS|nodes={}|ambient={}|T={}|K={}|"
    "status=PASS_EXACT_Q8_E7_NODE_PRINCIPAL_CLEARINGS".format(
        len(records), len(ambient["ambient_basis"]), T, K,
    ),
    flush=True,
)
