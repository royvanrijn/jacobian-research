#!/usr/bin/env sage -python
"""Derive the opposite pointed-quartic section on the exact q24 D12 model.

The D13-to-D12 quartic is pointed above the old I9* base value.  Evaluate the
degree-one pointed-quartic map at the opposite ordinate by its exact limit,
then apply the certified raw-to-minimal scaling.  This is univariate QQ(V)
arithmetic only and is intended to test a missing-parent-MW shortcut before
any section solve.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--output",
    type=Path,
    default=LOCAL / "q24-d12-pointed-opposite-section-qq.json",
)
args = parser.parse_args()

Q24 = LOCAL / "q24-d13-to-d12-component-valuation-qq.json"
PARENT = LOCAL / "q8-corrected2cover-qq-child.json"
for path in (Q24, PARENT):
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")

q24 = json.loads(Q24.read_text())
parent = json.loads(PARENT.read_text())
assert q24["status"] == "PASS_EXACT_Q24_D13_TO_D12_COMPONENT_VALUATION_RR"
assert parent["status"] == "PASS_EXACT_CORRECTED_Q8_D13_CHILD"

VQ = PolynomialRing(QQ, "V")
V = VQ.gen()
KV = VQ.fraction_field()
UQ = PolynomialRing(KV, "U")
U = UQ.gen()

coefficients = [
    KV(text) for text in q24["quartic"]["coefficients_in_U_low_to_high"]
]
quartic = sum(coefficients[index] * U**index for index in range(5))
assert quartic.degree() == 4

old_child = parent["child"]
i9 = next(row for row in old_child["finite_fibres"] if row["kodaira"] == "I9*")
TUQ = PolynomialRing(QQ, "U")
i9_factor = TUQ(str(i9["factor"]))
assert i9_factor.degree() == 1
alpha = -KV(i9_factor[0]) / KV(i9_factor[1])

q_squared = KV(quartic(alpha))
if not q_squared.is_square():
    raise ArithmeticError("the q24 pointed value is not a square")
q = KV(q_squared.sqrt())
assert q**2 == q_squared and q != 0

e, d0, c0, b0, a0 = coefficients
a = a0
b = b0 + 4 * alpha * a0
c = c0 + 3 * alpha * b0 + 6 * alpha**2 * a0
d = d0 + 2 * alpha * c0 + 3 * alpha**2 * b0 + 4 * alpha**3 * a0
w1 = -d / (2 * q)
w2 = (w1**2 - c) / (2 * q)
w3 = (2 * w1 * w2 - b) / (2 * q)
xg = 2 * q * w2
yg = 4 * q**2 * w3
a1 = d / q
a2 = c - d**2 / (4 * q**2)
a3 = 2 * q * b
b2 = a1**2 + 4 * a2
x_raw = KV(9 * (xg + b2 / 12))
y_raw = KV(27 * (yg + (a1 * xg + a3) / 2))

raw_A = KV(q24["jacobian_raw"]["A"])
raw_B = KV(q24["jacobian_raw"]["B"])
assert y_raw**2 == x_raw**3 + raw_A * x_raw + raw_B
A = VQ([QQ(value) for value in q24["child"]["minimal_A_coefficients_low_to_high"]])
B = VQ([QQ(value) for value in q24["child"]["minimal_B_coefficients_low_to_high"]])
cA = raw_A / KV(A)
cB = raw_B / KV(B)
scale2 = cB / cA
if not cB.is_square():
    raise ArithmeticError("raw/minimal y scale is not a square")
scale3 = KV(cB.sqrt())
assert scale2**2 == cA and scale3**2 == cB
x = KV(x_raw / scale2)
y = KV(y_raw / scale3)
assert y**2 == x**3 + KV(A) * x + KV(B)


def power_root(poly, exponent):
    poly = VQ(poly)
    leading = QQ(poly.leading_coefficient())
    if exponent % 2 == 0 and not leading.is_square():
        raise ArithmeticError("denominator leading coefficient is not a square")
    answer = VQ.one()
    for factor, multiplicity in (poly / leading).factor():
        if int(multiplicity) % exponent:
            raise ArithmeticError("section denominator is not a perfect power")
        answer *= factor.monic() ** (int(multiplicity) // exponent)
    return answer.monic()


Z = power_root(x.denominator(), 2)
if Z != power_root(y.denominator(), 3):
    raise ArithmeticError("x/y denominators give different Z")
X = VQ(x * Z**2)
Y = VQ(y * Z**3)
assert Y**2 == X**3 + A * X * Z**4 + B * Z**6

payload = {
    "schema": "elkies-k3.h3-q24-d12-pointed-opposite-section-qq.v1",
    "status": "PASS_EXACT_Q24_D12_POINTED_OPPOSITE_SECTION_QQ",
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in (Q24, PARENT)],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (Q24, PARENT)
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
        "degrees_X_Y_Z": [int(X.degree()), int(Y.degree()), int(Z.degree())],
        "P_dot_O": int(Z.degree()),
        "exact_weierstrass_identity": True,
    },
    "method": (
        "Exact opposite-point limit on the q24 pointed quartic followed by the "
        "certified raw-to-minimal Jacobian scaling."
    ),
    "proof_boundary": (
        "This constructs one exact D12 section. Its D12 MW vector and whether it "
        "supplies the missing parent coordinate require a marking fingerprint."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "Q24D12POINTEDOPPOSITE|degrees={},{},{}|PO={}|status={}".format(
        X.degree(), Y.degree(), Z.degree(), Z.degree(), payload["status"]
    ),
    flush=True,
)
print(f"OUTPUT|{args.output.resolve()}", flush=True)
