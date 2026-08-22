#!/usr/bin/env sage -python
"""Follow the marked q8 section through the two IV* entrance blow-ups.

This is deliberately smaller than a full IV* resolution.  It starts from the
exact q6-child Jacobian and the certified q8 marking, translates the finite
IV* place to ``u=0``, and proves that the marked section enters the smooth
locus of the second ordinary u-chart.  This supplies the exact local jet from
which a resolved component module can be constructed; it does not choose the
remaining E6 component labels or construct the global q8 pencil.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[2]
CHILD = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-jacobian.json"
MARKING = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-marking.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-ivstar-entrance.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial(ring, coefficients):
    return ring([QQ(value) for value in coefficients])


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--child", type=Path, default=CHILD)
parser.add_argument("--marking", type=Path, default=MARKING)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
for name in ("child", "marking", "output"):
    setattr(args, name, getattr(args, name).resolve())

child = json.loads(args.child.read_text())
marking = json.loads(args.marking.read_text())
assert child["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert marking["status"] == "PASS_EXACT_Q6_CHILD_Q8_MARKING"

T_ring = PolynomialRing(QQ, "T")
T = T_ring.gen()
F_T = T_ring.fraction_field()
model = child["minimal_short_weierstrass"]
A = polynomial(T_ring, model["A_coefficients_low_to_high"])
B = polynomial(T_ring, model["B_coefficients_low_to_high"])
section = marking["selected_q8"]["relative_child_section_standard_jacobian_coordinates"]
x_section = F_T(polynomial(T_ring, section["x_numerator_coefficients_low_to_high"])) / F_T(
    polynomial(T_ring, section["x_denominator_coefficients_low_to_high"])
)
y_section = F_T(polynomial(T_ring, section["y_numerator_coefficients_low_to_high"])) / F_T(
    polynomial(T_ring, section["y_denominator_coefficients_low_to_high"])
)
assert y_section**2 == x_section**3 + A * x_section + B

iv_fibre = next(item for item in child["finite_fibres"] if item["kodaira"] == "IV*")
factor = T_ring(iv_fibre["factor"])
assert factor.degree() == 1 and tuple(iv_fibre["minimal_orders"]) == (3, 4, 8)
base_point = -factor[0] / factor[1]

u_ring = PolynomialRing(QQ, "u")
u = u_ring.gen()
F_u = u_ring.fraction_field()


def translate(value):
    """Translate a QQ(T) rational function to QQ(u), T=base_point+u."""
    numerator = u_ring(value.numerator()(base_point + u))
    denominator = u_ring(value.denominator()(base_point + u))
    return F_u(numerator) / F_u(denominator)


A_u = u_ring(A(base_point + u))
B_u = u_ring(B(base_point + u))
a_u, remainder = A_u.quo_rem(u**3)
assert not remainder and a_u(0)
b_u, remainder = B_u.quo_rem(u**4)
assert not remainder and b_u(0)
x_u = translate(x_section)
y_u = translate(y_section)
assert x_u.valuation() == 2 and y_u.valuation() == 2
X_section = x_u / u**2
Y_section = y_u / u**2
assert X_section.denominator()(0) and Y_section.denominator()(0)
assert Y_section(0) and Y_section(0)**2 == b_u(0)

R = PolynomialRing(QQ, names=("u", "x", "y"))
u_chart, x_chart, y_chart = R.gens()
a = R(a_u(u_chart))
b = R(b_u(u_chart))
f0 = y_chart**2 - x_chart**3 - u_chart**3 * a * x_chart - u_chart**4 * b


def strict(value, substitutions, exceptional, multiplicity):
    transformed = R(value(*substitutions))
    quotient, remainder = transformed.quo_rem(exceptional**multiplicity)
    assert not remainder
    return quotient


# The marked section has x,y orders two, so it stays at the origin in the
# first u-chart and becomes a smooth point on the second u-chart's exceptional
# divisor.  These are actual substitutions in the q6-child model, not a
# Kodaira-label proxy.
f1 = strict(f0, (u_chart, u_chart * x_chart, u_chart * y_chart), u_chart, 2)
f2 = strict(f1, (u_chart, u_chart * x_chart, u_chart * y_chart), u_chart, 2)
assert f1 == y_chart**2 - u_chart*x_chart**3 - u_chart**2*a*x_chart - u_chart**2*b
assert f2 == y_chart**2 - u_chart**2*x_chart**3 - u_chart*a*x_chart - b
assert f2(u, X_section, Y_section) == 0
entrance_point = (QQ(0), QQ(X_section(0)), QQ(Y_section(0)))
assert f2(*entrance_point) == 0
assert f2.derivative(y_chart)(*entrance_point) == 2 * entrance_point[2]
assert f2.derivative(y_chart)(*entrance_point) != 0

# On this actual second chart the exceptional fibre splits over QQ into the
# two branches y=+c and y=-c, c=Y_section(0).  Both are smooth everywhere in
# this chart because d(f2)/dy is nonzero there.  The marked section chooses
# one branch exactly; subsequent blow-ups elsewhere in the IV* resolution do
# not change its germ at this point.
c = entrance_point[2]
R_fibre = PolynomialRing(QQ, names=("x", "y"))
x_fibre, y_fibre = R_fibre.gens()
fibre_restriction = R_fibre(f2(0, x_fibre, y_fibre))
assert fibre_restriction == y_fibre**2 - c**2
assert fibre_restriction == (y_fibre - c) * (y_fibre + c)
branch_data = {
    "exceptional_fibre_restriction": "y^2-c^2=(y-c)(y+c), c=Y_section(0)",
    "marked_branch": "y-c=0",
    "other_branch": "y+c=0",
    "both_branches_smooth_in_second_chart": True,
    "marked_branch_derivative_in_y": str(2 * c),
    "resolution_consequence": (
        "The marked germ is already smooth after the second u-chart blow-up; "
        "further IV* resolution centres do not alter that marked point."
    ),
}

payload = {
    "schema": "elkies-k3.h92-q6-child-q8-ivstar-entrance.v1",
    "status": "PASS_EXACT_Q6_CHILD_Q8_IVSTAR_ENTRANCE",
    "inputs": {
        "child_jacobian": {"path": str(args.child.relative_to(ROOT)), "sha256": digest(args.child)},
        "q8_marking": {"path": str(args.marking.relative_to(ROOT)), "sha256": digest(args.marking)},
    },
    "ivstar": {
        "factor": str(factor),
        "local_base": "u=T-(-factor[0]/factor[1])",
        "integral_germ": "y^2=x^3+u^3*a(u)*x+u^4*b(u), with a(0)b(0)!=0",
        "section_orders": {"x": 2, "y": 2},
        "two_u_chart_transforms": {
            "first": str(f1),
            "second": str(f2),
            "maps": [
                "(u,x,y)=(u_1,u_1*x_1,u_1*y_1)",
                "(u_1,x_1,y_1)=(u_2,u_2*x_2,u_2*y_2)",
            ],
        },
        "marked_entrance": {
            "point": "u=0, x=X_section(0), y=Y_section(0)",
            "proof": "Y_section(0)^2=b(0) and d(f2)/dy=2*Y_section(0) is nonzero",
            "smooth": True,
        },
        "second_chart_exceptional_branches": branch_data,
    },
    "boundary": (
        "This follows the q8 marked section into a genuine smooth IV* blow-up chart. "
        "It does not complete the IV* resolution, identify the E6 component label, "
        "derive the II*/IV* finite modules, or construct a global q8 pencil."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92Q6CHILDQ8IV|orders=2,2|u_charts=2|entrance=smooth|"
    "status=PASS_EXACT_Q6_CHILD_Q8_IVSTAR_ENTRANCE",
    flush=True,
)
