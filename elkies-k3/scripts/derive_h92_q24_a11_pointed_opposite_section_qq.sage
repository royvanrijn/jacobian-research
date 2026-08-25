#!/usr/bin/env sage -python
"""Derive the opposite pointed-quartic section on the exact A11 Jacobian.

The resolved D12-to-A11 binary quartic is pointed at ``(V,W)=(alpha,q)``.
That point maps to the Jacobian zero, but the conjugate point over the same
old-base value, ``(alpha,-q)``, is another rational point over ``QQ(T)``.
This script evaluates the degree-one pointed-quartic map by its exact limit at
that point and audits its pole order and I12 component.  It uses only
univariate rational-function arithmetic; there is no section ansatz or
Groebner calculation.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import PolynomialRing, PowerSeriesRing, QQ


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--output",
    type=Path,
    default=LOCAL / "q24-a11-pointed-opposite-section-qq.json",
)
args = parser.parse_args()

A11 = LOCAL / "q24-d12-to-a11-orbit42-resolved-rr-qq.json"
PARENT = LOCAL / "q24-d13-to-d12-component-valuation-qq.json"
for path in (A11, PARENT):
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")

a11 = json.loads(A11.read_text())
parent = json.loads(PARENT.read_text())
assert a11["status"] == "PASS_EXACT_Q24_D12_Q6_A11_COMPONENT_VALUATION_RR"
assert parent["status"] == "PASS_EXACT_Q24_D13_TO_D12_COMPONENT_VALUATION_RR"

TQ = PolynomialRing(QQ, "T")
T = TQ.gen()
KT = TQ.fraction_field()
VQ = PolynomialRing(KT, "V")
V = VQ.gen()

coefficients = [
    KT(TQ(text))
    for text in a11["quartic"]["coefficients_in_T_low_to_high"]
]
assert len(coefficients) == 5
quartic = sum(coefficients[index] * V**index for index in range(5))

i8 = next(
    row for row in parent["child"]["finite_fibres"]
    if row["kodaira"] == "I8*"
)
V0Q = PolynomialRing(QQ, "V")
i8_factor = V0Q(str(i8["factor"]))
assert i8_factor.degree() == 1
alpha = -KT(i8_factor[0]) / KT(i8_factor[1])

q_squared = KT(quartic(alpha))
if not q_squared.is_square():
    raise ArithmeticError("the pointed quartic value is not a square")
q = KT(q_squared.sqrt())
assert q**2 == q_squared and q != 0

e, d0, c0, b0, a0 = coefficients
a = a0
b = b0 + 4 * alpha * a0
c = c0 + 3 * alpha * b0 + 6 * alpha**2 * a0
d = d0 + 2 * alpha * c0 + 3 * alpha**2 * b0 + 4 * alpha**3 * a0

# Expansion of the opposite branch W=-q+w1*u+w2*u^2+w3*u^3+... .
w1 = -d / (2 * q)
w2 = (w1**2 - c) / (2 * q)
w3 = (2 * w1 * w2 - b) / (2 * q)
xg = 2 * q * w2
yg = 4 * q**2 * w3

a1 = d / q
a2 = c - d**2 / (4 * q**2)
a3 = 2 * q * b
b2 = a1**2 + 4 * a2
x = KT(9 * (xg + b2 / 12))
y = KT(27 * (yg + (a1 * xg + a3) / 2))

A = TQ([QQ(value) for value in a11["child"]["minimal_A_coefficients_low_to_high"]])
B = TQ([QQ(value) for value in a11["child"]["minimal_B_coefficients_low_to_high"]])
assert y**2 == x**3 + KT(A) * x + KT(B)


def power_root(poly, exponent):
    poly = TQ(poly)
    leading = QQ(poly.leading_coefficient())
    if not leading.is_square() and exponent % 2 == 0:
        raise ArithmeticError("denominator leading coefficient is not a square")
    answer = TQ.one()
    for factor, multiplicity in (poly / leading).factor():
        if int(multiplicity) % exponent:
            raise ArithmeticError("section denominator is not a perfect power")
        answer *= factor.monic() ** (int(multiplicity) // exponent)
    return answer.monic()


Zx = power_root(x.denominator(), 2)
Zy = power_root(y.denominator(), 3)
if Zx != Zy:
    raise ArithmeticError("x and y denominators give different section Z")
Z = Zx
X = TQ(x * Z**2)
Y = TQ(y * Z**3)
assert Y**2 == X**3 + A * X * Z**4 + B * Z**6

pole_order = int(Z.degree())
degree_bounds = {
    "X": int(X.degree()),
    "Y": int(Y.degree()),
    "Z": pole_order,
}

# Exact formal audit at the unique I12 fibre.  For a split I_n fibre the
# component depth up to negation is min(v(x-center),v(y),n/2).
delta_rows = a11["child"]["discriminant_factorization"]
i12_factor = TQ(next(row["factor"] for row in delta_rows if int(row["multiplicity"]) == 12))
assert i12_factor.degree() == 1
beta = -QQ(i12_factor[0]) / QQ(i12_factor[1])
SQ = PowerSeriesRing(QQ, "s", default_prec=15)
s = SQ.gen()
A_series = SQ(A(T + beta))
B_series = SQ(B(T + beta))
center = SQ(-3 * B_series[0] / (2 * A_series[0]))
for unused in range(7):
    center = (center + (-A_series / 3) / center) / 2
assert (center**2 + A_series / 3).valuation() >= 14
g_center = center**3 + A_series * center + B_series
assert g_center.valuation() == 12

x_series = SQ(x(T + beta))
y_series = SQ(y(T + beta))
vx = int((x_series - center).valuation())
vy = int(y_series.valuation())
component_depth = min(vx, vy, 6)
correction = QQ(component_depth * (12 - component_depth)) / 12
height = QQ(4 + 2 * pole_order) - correction

payload = {
    "schema": "elkies-k3.h3-q24-a11-pointed-opposite-section-qq.v1",
    "status": "PASS_EXACT_A11_POINTED_OPPOSITE_SECTION_QQ",
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in (A11, PARENT)],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (A11, PARENT)
        },
    },
    "pointed_quartic": {
        "old_base_value": str(alpha),
        "point_used_as_zero": [str(alpha), str(q)],
        "opposite_point": [str(alpha), str(-q)],
        "limit_formula_verified": True,
    },
    "section": {
        "X_coefficients_low_to_high": [str(value) for value in X.list()],
        "Y_coefficients_low_to_high": [str(value) for value in Y.list()],
        "Z_coefficients_low_to_high": [str(value) for value in Z.list()],
        "degrees_X_Y_Z": [degree_bounds["X"], degree_bounds["Y"], pole_order],
        "exact_weierstrass_identity": True,
    },
    "height_profile": {
        "P_dot_O": pole_order,
        "I12_root": str(beta),
        "valuation_x_minus_center": vx,
        "valuation_y": vy,
        "I12_component_up_to_negation": [component_depth, 12 - component_depth],
        "local_correction": str(correction),
        "height": str(height),
    },
    "method": (
        "Exact limit of the pointed binary-quartic isomorphism at the opposite "
        "rational point over the same old-base value; univariate QQ(T) arithmetic only."
    ),
    "proof_boundary": (
        "This constructs and profiles one exact A11 section. Its six-coordinate "
        "MW marking and membership in the target bridge coset still require exact "
        "lattice/equation identification."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print(
    "A11POINTEDOPPOSITE|degrees={},{},{}|component={}|correction={}|height={}|status={}".format(
        degree_bounds["X"], degree_bounds["Y"], pole_order,
        component_depth, correction, height, payload["status"],
    ),
    flush=True,
)
print(f"OUTPUT|{args.output.resolve()}", flush=True)
