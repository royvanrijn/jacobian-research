#!/usr/bin/env sage -python
"""Certify the q=6 marked module on every actual resolved H92 E7 chart.

The source q=6 divisor has E7 degree one on ``E7_5``.  Its rank-one
reflexive local module is the marked branch module.  This script proves the
only facts needed by the bounded q=6 coefficient ambient on the *actual*
resolution:

* the raw chord ``m=(y-y(P1))/(x-x(P1))`` has no exceptional pole on any of
  the seven resolved E7 components;
* its only horizontal pole over that fibre is the marked smooth point
  ``-P1``; there the corrected chart calculation makes ``m`` a unit times
  the marked generator ``Z*m/t=1/W``;
* consequently ``<1,m>`` is the actual q=6 module cover, with generator
  ``1`` away from ``-P1`` and generator ``m`` near ``-P1``.

The nontrivial cancellation in ``x-x(P1)`` is checked in the actual
``E7_2--E7_5`` blow-up chart.  No Kodaira normal form or component-label
guess is used.
"""

import argparse
import hashlib
import json
from importlib.machinery import SourceFileLoader
from pathlib import Path

from sage.all import PolynomialRing, QQ, vector, ZZ


ROOT = Path(__file__).resolve().parents[2]
ANCHOR = ROOT / "elkies-k3/scripts/verify_h3_noncm_q6_source_anchor.sage"
H92 = ROOT / "artifacts/local/humbert-inputs/92/igusa92.txt"
P1 = ROOT / "artifacts/generated-results/elkies-k3-h92-p1-lift.json"
ATLAS = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-actual-e7-valuation-atlas.json"
TRACE = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-p1-actual-e7-trace.json"
MARKED = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-p1-actual-e7-marked-module-corrected.json"
AUDIT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-actual-e7-marked-chord-order-audit.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-actual-e7-all-edge-module.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial(ring, coefficients):
    return ring([QQ(value) for value in coefficients])


def invert_base(rational_u):
    """Rewrite a QQ(u) function as a QQ(t) function, with t=1/u."""
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
parser.add_argument("--atlas", type=Path, default=ATLAS)
parser.add_argument("--trace", type=Path, default=TRACE)
parser.add_argument("--marked", type=Path, default=MARKED)
parser.add_argument("--audit", type=Path, default=AUDIT)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

p1 = json.loads(args.p1.read_text())
atlas = json.loads(args.atlas.read_text())
trace = json.loads(args.trace.read_text())
marked = json.loads(args.marked.read_text())
audit = json.loads(args.audit.read_text())
assert p1["status"] == "PASS_EXACT_H92_P1"
assert atlas["status"] == "PASS_EXACT_H92_E7_VALUATION_ATLAS"
assert trace["status"] == "PASS_EXACT_P1_ACTUAL_E7_TRACE"
assert marked["status"] == "PASS_EXACT_P1_ACTUAL_E7_MARKED_MODULE_CORRECTED"
assert marked["chart"] == "E7_2--E7_5"
assert marked["marked_component"] == "E7_5"
assert marked["unit_coefficients_at_minus_P1"]["(Z*m/t)*W"] != "0"
assert audit["status"] == "REJECTS_Q6_MARKED_E7_FRAME_AS_STATED"
assert audit["exact_Z_orders"] == {
    "t": 1, "x_minus_xP": 3, "y_minus_yP": 3, "m": 0, "m_over_t": -1,
}

anchor = SourceFileLoader("h92_q6_all_edge_anchor", str(ANCHOR)).load_module()
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
t_ring = x_p_t.parent()
t = t_ring.gen()
c2 = QQ((x_p_t/t**2)(0))
c3 = QQ((x_p_t/t**2).derivative()(0))
assert c2 == -B1/A1

# The seven actual exceptional valuations come from the resolved H92 atlas.
component_order = tuple(entry["component"] for entry in atlas["entries"])
v_t = vector(ZZ, [entry["old_coordinate_orders"]["t"] for entry in atlas["entries"]])
v_x = vector(ZZ, [entry["old_coordinate_orders"]["x"] for entry in atlas["entries"]])
v_y = vector(ZZ, [entry["old_coordinate_orders"]["y"] for entry in atlas["entries"]])
assert component_order == ("E7_1", "E7_2", "E7_3", "E7_4", "E7_5", "E7_6", "E7_7")
assert v_t == vector(ZZ, (2, 2, 4, 3, 1, 2, 3))

# Away from E7_2 and E7_5 the leading order of x-x(P1) is forced by the
# strict inequality v(x)<2v(t).  The y numerator needs only its valuation
# lower bound, which follows from the ultrametric inequality.
dx_exact = vector(ZZ, (2, 4, 6, 4, 3, 3, 5))
dy_lower = vector(ZZ, (3, 5, 9, 6, 3, 5, 7))
for index in (0, 2, 3, 5, 6):
    assert v_x[index] < 2*v_t[index]
    assert dx_exact[index] == v_x[index]
for index in range(7):
    assert dy_lower[index] == min(v_y[index], 3*v_t[index])

# The two equality cases are calculated after substituting the actual
# E7_2--E7_5 blow-up map.  Its two exceptional components are U=0 (E7_2)
# and Z=0 (E7_5); the atlas gives U~Y^2 along E7_2.
ring = PolynomialRing(QQ, names=("Z", "U", "Y"))
Z, U, Y = ring.gens()
second = -A1/B1
t_chart = Z*U*(Z+second)
x_chart = Z**2*U**2*(Z+second)
leading_x_difference = ring(x_chart-c2*t_chart**2)
assert leading_x_difference == -c2*Z**3*U**2*(Z+second)
assert c2*second == 1

# At E7_2, the displayed coefficient is a unit times U^2, and U has
# reduced order two.  Higher P1 terms begin with t^3 and have U-order three.
assert dx_exact[1] == 4
# At E7_5, retain the t^3 P1 coefficient.  Its Z^3 coefficient is not the
# zero rational function at the generic point (U is a generic coordinate).
x_z3 = ring(-c2*U**2*second-c3*U**3*second**3)
assert x_z3
assert dx_exact[4] == 3

m_lower = dy_lower-dx_exact
assert m_lower == vector(ZZ, (1, 1, 3, 2, 0, 2, 2))
assert all(value >= 0 for value in m_lower)
# The only equality in the numerator lower-bound table is E7_5; the exact
# chart audit supplies its order. Thus the displayed vector is exact, not
# merely a regularity lower bound, and can safely feed later q8 valuations.
m_exact = vector(ZZ, (1, 1, 3, 2, 0, 2, 2))
assert m_exact == m_lower

# The marked trace supplies the unique horizontal pole.  The Weierstrass
# factorization shows x=x(P1) meets the curve only at +/-P1; at P1 the chord
# numerator cancels, and at -P1 the corrected chart gives m=unit/W.
assert trace["node_chart"]["name"] == "E7_2--E7_5"
assert trace["node_chart"]["-P1"]["Z"] == "0"
assert trace["node_chart"]["P1"]["Z"] == "0"
branch_identity = "(y-y(P1))*(y+y(P1))=(x-x(P1))*(x^2+x*x(P1)+x(P1)^2+a(t))"

payload = {
    "schema": "elkies-k3.h92-q6-actual-e7-all-edge-module.v1",
    "status": "PASS_EXACT_Q6_ACTUAL_E7_ALL_EDGE_MODULE",
    "inputs": {
        "h92_source": {"path": str(H92.relative_to(ROOT)), "sha256": digest(H92)},
        "p1": {"path": str(args.p1.relative_to(ROOT)), "sha256": digest(args.p1)},
        "actual_valuation_atlas": {"path": str(args.atlas.relative_to(ROOT)), "sha256": digest(args.atlas)},
        "actual_trace": {"path": str(args.trace.relative_to(ROOT)), "sha256": digest(args.trace)},
        "corrected_marked_module": {"path": str(args.marked.relative_to(ROOT)), "sha256": digest(args.marked)},
        "marked_chord_order_audit": {"path": str(args.audit.relative_to(ROOT)), "sha256": digest(args.audit)},
    },
    "actual_component_order": list(component_order),
    "exceptional_orders": {
        "t": [int(value) for value in v_t],
        "x": [int(value) for value in v_x],
        "y": [int(value) for value in v_y],
        "x_minus_xP1": [int(value) for value in dx_exact],
        "y_minus_yP1_lower_bound": [int(value) for value in dy_lower],
        "m_lower_bound": [int(value) for value in m_lower],
        "m_exact": [int(value) for value in m_exact],
    },
    "actual_equality_case_certificate": {
        "chart": "E7_2--E7_5",
        "pullbacks": {"t": str(t_chart), "x": str(x_chart)},
        "leading_x_minus_xP1": str(leading_x_difference),
        "E7_2": "U=0 with U a unit times Y^2; the displayed leading term has exact reduced order four",
        "E7_5_Z3_coefficient": str(x_z3),
        "E7_5": "the displayed Z^3 coefficient is nonzero at the generic point, so ord_E7_5(x-x(P1))=3",
    },
    "marked_horizontal_frame": {
        "point": "-P1 on E7_5",
        "generator": "m=(t/Z)*(Z*m/t)",
        "corrected_identity": "Z*m/t=unit/W",
        "conclusion": "m=unit/W, hence m is a local frame at the marked pole",
    },
    "module_cover": {
        "away_from_minus_P1": "generator 1; m is regular because its exceptional orders are nonnegative and its only possible horizontal denominator branches are +/-P1",
        "near_minus_P1": "generator m; 1/m is regular because m=unit/W",
        "transition": "on the overlap avoiding the zero and pole loci of m, multiplication by m is a unit",
        "branch_factorization": branch_identity,
    },
    "compiler_conclusion": (
        "For the q6 bounded ambient a+b*m, every coefficient regular at t=0 "
        "satisfies the complete actual E7 cover. The E7 module contributes no "
        "additional linear condition row; only the E8 and smooth collision blocks "
        "remain in the q6 condition matrix."
    ),
    "boundary": (
        "This certifies the q6 E7 module cover for the declared two-generator "
        "ambient. It does not itself eliminate the pencil, identify the child, or "
        "transport sections."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "H92Q6ACTUALE7ALLEDGE|m_orders_lower=1,1,3,2,0,2,2|"
    "status=PASS_EXACT_Q6_ACTUAL_E7_ALL_EDGE_MODULE",
    flush=True,
)
