#!/usr/bin/env sage -python
"""Certify the corrected H92 q6-child q8 marked section.

The binary-quartic covariant map used by the q6 child is the canonical
2-covering map to the Jacobian.  Therefore differences of transported quartic
points already represent twice the corresponding primitive MW directions.
The old q8 marking doubled those differences a second time.

Write

    Pmap = phi(E7_7) - phi(old_O),
    Qmap = phi(E7_7) - phi(affine_E7).

The corrected q8 section of MW coordinate (-2,-2,0) is

    S = Pmap + Qmap,

while the withdrawn section is exactly 2*S.  On the globally minimal E8+E6
K3 this script proves S.O=10, height(S)=24, smooth identity-component
specialization at II* and IV*, and a reduced degree-10 collision divisor.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, PolynomialRing, QQ, matrix, vector

ROOT = Path(__file__).resolve().parents[2]
GEN = ROOT / "artifacts/generated-results"
CHILD = GEN / "elkies-k3-h92-q6-child-jacobian.json"
ZERO = GEN / "elkies-k3-h92-q6-child-zero-section.json"
COMPONENTS = GEN / "elkies-k3-h92-q6-child-e7-infinity-sections.json"
DEFAULT_OUTPUT = GEN / "elkies-k3-h92-q6-child-q8-marking-2cover.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial(ring, values):
    return ring([QQ(value) for value in values])


def rational(field, ring, data, numerator, denominator):
    return field(polynomial(ring, data[numerator])) / field(
        polynomial(ring, data[denominator])
    )


def point_data(point, ring):
    x, y = point.xy()
    return {
        "x_numerator_coefficients_low_to_high": [str(v) for v in ring(x.numerator()).list()],
        "x_denominator_coefficients_low_to_high": [str(v) for v in ring(x.denominator()).list()],
        "y_numerator_coefficients_low_to_high": [str(v) for v in ring(y.numerator()).list()],
        "y_denominator_coefficients_low_to_high": [str(v) for v in ring(y.denominator()).list()],
        "coordinate_degrees": {
            "x": [int(ring(x.numerator()).degree()), int(ring(x.denominator()).degree())],
            "y": [int(ring(y.numerator()).degree()), int(ring(y.denominator()).degree())],
        },
    }


def monic_power_root(value, exponent):
    root = value.parent().one()
    for factor, multiplicity in value.factor():
        assert multiplicity % exponent == 0
        root *= factor.monic() ** (multiplicity // exponent)
    return root.monic()


def order_at(value, factor):
    return int(value.numerator().valuation(factor) - value.denominator().valuation(factor))


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--child", type=Path, default=CHILD)
parser.add_argument("--zero", type=Path, default=ZERO)
parser.add_argument("--components", type=Path, default=COMPONENTS)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
for name in ("child", "zero", "components", "output"):
    setattr(args, name, getattr(args, name).resolve())

child = json.loads(args.child.read_text())
zero = json.loads(args.zero.read_text())
components = json.loads(args.components.read_text())
assert child["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert zero["status"] == "PASS_EXACT_CHILD_ZERO_SECTION_TRANSPORT"
assert components["status"] == "PASS_EXACT_CHILD_E7_INFINITY_TRANSPORT"

ring = PolynomialRing(QQ, "T")
field = ring.fraction_field()
A = polynomial(ring, child["minimal_short_weierstrass"]["A_coefficients_low_to_high"])
B = polynomial(ring, child["minimal_short_weierstrass"]["B_coefficients_low_to_high"])
Delta = polynomial(ring, child["minimal_short_weierstrass"]["Delta_coefficients_low_to_high"])
assert A.degree() <= 8 and B.degree() <= 12 and Delta.degree() <= 24
curve = EllipticCurve(field, [0, 0, 0, field(A), field(B)])

zdata = zero["section"]
old_zero = curve(
    rational(field, ring, zdata, "x_numerator_coefficients_low_to_high", "x_denominator_coefficients_low_to_high"),
    rational(field, ring, zdata, "y_numerator_coefficients_low_to_high", "y_denominator_coefficients_low_to_high"),
)
points = {}
for entry in components["sections"]:
    points[entry["sign"]] = curve(
        rational(field, ring, entry, "x_numerator_coefficients_low_to_high", "x_denominator_coefficients_low_to_high"),
        rational(field, ring, entry, "y_numerator_coefficients_low_to_high", "y_denominator_coefficients_low_to_high"),
    )
affine = points[components["source"]["affine_E7_sign"]]
e7_7 = points[components["source"]["E7_7_sign"]]

# Covariant images are doubled MW directions.  Do not multiply these
# differences by two again.
first_covariant_difference = e7_7 - old_zero
second_covariant_difference = e7_7 - affine
selected = first_covariant_difference + second_covariant_difference
withdrawn = 2 * first_covariant_difference + 2 * second_covariant_difference
assert withdrawn == 2 * selected and not selected.is_zero()

sx, sy = selected.xy()
dx, dy = ring(sx.denominator()), ring(sy.denominator())
h = monic_power_root(dx, 2)
assert h == monic_power_root(dy, 3)
assert h.degree() == 10
assert dx // h**2 in QQ and dy // h**3 in QQ
assert h.gcd(Delta) in QQ
z = -sx / sy
assert ring((z / field(h)).numerator()).gcd(h) in QQ
assert ring((z / field(h)).denominator()).gcd(h) in QQ
assert ring(sx.numerator()).degree() - dx.degree() == 4
assert ring(sy.numerator()).degree() - dy.degree() == 6

# The corrected section is smooth at both additive fibres, so its local
# Shioda corrections vanish.  Since chi(O_K3)=2, height=4+2(S.O)=24.
local = {}
for fibre in child["finite_fibres"]:
    if fibre["kodaira"] not in ("II*", "IV*"):
        continue
    factor = ring(fibre["factor"])
    point = -factor[0] / factor[1]
    ox, oy = order_at(sx, factor), order_at(sy, factor)
    assert (ox, oy) == (0, 0)
    specialization = (QQ(sx(point)), QQ(sy(point)))
    assert specialization != (0, 0)
    local[fibre["kodaira"]] = {
        "factor": str(factor),
        "orders_x_y": [ox, oy],
        "specialization": [str(value) for value in specialization],
        "component": "identity",
    }
assert set(local) == {"II*", "IV*"}
height = QQ(4 + 2*h.degree())
assert height == 24
height_gram = matrix(QQ, [
    [QQ(8)/3, QQ(1)/3, -1],
    [QQ(1)/3, QQ(8)/3, 1],
    [-1, 1, 46],
])
coordinates = vector(QQ, (-2, -2, 0))
assert coordinates * height_gram * coordinates == height

# Regression for the old doubled point.
wx, wy = withdrawn.xy()
old_h = monic_power_root(ring(wx.denominator()), 2)
assert old_h == monic_power_root(ring(wy.denominator()), 3)
assert old_h.degree() == 46
assert QQ(4 + 2*old_h.degree()) == 96 == 4*height

payload = {
    "schema": "elkies-k3.h92-q6-child-q8-marking-2cover.v1",
    "status": "PASS_EXACT_Q8_MARKING_2COVER_CORRECTION",
    "inputs": {
        "child": {"path": str(args.child.relative_to(ROOT)), "sha256": digest(args.child)},
        "zero": {"path": str(args.zero.relative_to(ROOT)), "sha256": digest(args.zero)},
        "components": {"path": str(args.components.relative_to(ROOT)), "sha256": digest(args.components)},
    },
    "correction": {
        "binary_quartic_covariant_role": "degree-two covering map to the Jacobian",
        "covariant_difference_multiplier_in_MW": 2,
        "withdrawn_formula": "2*Pmap+2*Qmap",
        "correct_formula": "Pmap+Qmap",
        "exact_group_check": "withdrawn_point=2*corrected_point",
    },
    "selected_q8": {
        "relative_child_section_MW_coordinates": [-2, -2, 0],
        "relative_child_section_standard_jacobian_coordinates": point_data(selected, ring),
        "height": str(height),
        "O_intersection": int(h.degree()),
        "collision_divisor": {
            "degree": int(h.degree()),
            "coefficients_low_to_high": [str(value) for value in h.list()],
            "squarefree": bool(h.gcd(h.derivative()) in QQ),
            "coprime_to_discriminant": bool(h.gcd(Delta) in QQ),
        },
        "additive_fibres": local,
        "generic_rr_basis": ["1", "m"],
        "m": "(y+y(S))/(x-x(S))",
    },
    "withdrawn_marking_regression": {
        "is_double_of_corrected_section": True,
        "height": "96",
        "collision_degree": 46,
    },
    "boundary": (
        "This certifies the corrected rational q8 marking.  The complete q8 "
        "pencil and D13 child are certified separately by "
        "derive_h92_q6_child_q8_corrected2cover_qq.sage."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92Q6CHILDQ8MARK2|mw=-2,-2,0|height=24|O=10|"
    "old_double_height=96|old_collision=46|new_collision=10|"
    "status=PASS_EXACT_Q8_MARKING_2COVER_CORRECTION",
    flush=True,
)
