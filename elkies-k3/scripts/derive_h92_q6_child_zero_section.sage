#!/usr/bin/env sage -python
"""Transport the old zero section to the first H3 q=6 child Jacobian.

The first exact q=6 pencil is represented by a binary quartic in the old
base coordinate ``u``.  The old zero section has new-base value ``T=b1/b0``;
after cancelling its common ``u^3`` factor, this gives a linear formula
``u=u_O(T)``.  The resulting quartic value is a square in ``QQ(T)``.

For a binary quartic ``f`` write ``H=(f_xx*f_zz-f_xz^2)/3`` and
``G=f_x*H_z-f_z*H_x``.  The classical identity

    G^2 = -16/3 H^3 + 256 I H f^2 - 1024/3 J f^3

gives the exact Jacobian map

    X=-3H/(4f),  Y=9G*w/(32f^2),  w^2=f,

to ``Y^2=X^3-27I*X-27J``.  We replay this map at ``u_O(T)`` and then the
minimalizing fourth/sixth-power change used by the child certificate.  The
output is one genuine child section; it does not claim a full MW basis or a
Shioda Gram replay.
"""

import argparse
import hashlib
import json
from importlib.machinery import SourceFileLoader
from pathlib import Path

from sage.all import PolynomialRing, QQ, prod


ROOT = Path(__file__).resolve().parents[2]
RR = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-global-rr.json"
CHILD = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-jacobian.json"
SECTION = ROOT / "artifacts/generated-results/elkies-k3-h92-p1-lift.json"
SECTION_SHA256 = "c323bf6346bb239934a5a2d8b1a3f4067e70e993d2e4eb32aaa30f469fca6397"
ANCHOR = ROOT / "elkies-k3/scripts/verify_h3_noncm_q6_source_anchor.sage"
H92 = ROOT / "artifacts/local/humbert-inputs/92/igusa92.txt"
CORE = ROOT / "elkies-k3/scripts/elliptic_neighbor_compiler.sage"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-zero-section.json"

exec(compile(CORE.read_text(), str(CORE), "exec"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial(ring, coefficients):
    return ring([QQ(value) for value in coefficients])


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--rr", type=Path, default=RR)
parser.add_argument("--child", type=Path, default=CHILD)
parser.add_argument("--section", type=Path, default=SECTION)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument(
    "--component-output", type=Path,
    help="write the two transported points at binary-quartic infinity",
)
args = parser.parse_args()

rr = json.loads(args.rr.read_text())
child = json.loads(args.child.read_text())
assert rr["status"] == "PASS_EXACT_GLOBAL_RR_KERNEL"
assert child["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert digest(args.section) == SECTION_SHA256
section = json.loads(args.section.read_text())

u_ring = PolynomialRing(QQ, "u")
u = u_ring.gen()
u_field = u_ring.fraction_field()
h = polynomial(u_ring, section["structured_denominator"]["Z4_coefficients"])


def coefficient_pair(entry):
    A = polynomial(u_ring, entry["A_coefficients_low_to_high"])
    B = polynomial(u_ring, entry["B_coefficients_low_to_high"])
    return u_field(A) / u_field(h**2), u_field(B) / u_field(h)


(a0, b0), (a1, b1) = tuple(
    coefficient_pair(entry) for entry in rr["kernel"]["sections"]
)
B0 = u_ring(b0 * h)
B1 = u_ring(b1 * h)
assert B0.valuation() == 3 and B1.valuation() >= 3

T_ring = PolynomialRing(QQ, "T")
T = T_ring.gen()
T_field = T_ring.fraction_field()
old_base_ring = PolynomialRing(T_field, "u")
uu = old_base_ring.gen()
old_base_field = old_base_ring.fraction_field()


def lift_polynomial(value):
    return old_base_ring([T_field(coefficient) for coefficient in value.list()])


def transport(value):
    return old_base_field(
        old_base_ring([T_field(coefficient) for coefficient in value.numerator().list()])
    ) / old_base_field(
        old_base_ring([T_field(coefficient) for coefficient in value.denominator().list()])
    )


# At old O, m has a simple pole, so T=(a1+b1*m)/(a0+b0*m) specializes to b1/b0.
zero_line = (lift_polynomial(B1) - T * lift_polynomial(B0)) // uu**3
assert zero_line.degree() == 1
u_zero = -zero_line[0] / zero_line[1]

x_p = u_field(polynomial(
    u_ring, section["x_entrance_base"]["numerator_coefficients"]
)) / u_field(polynomial(u_ring, section["x_entrance_base"]["denominator_coefficients"]))
y_p = u_field(polynomial(
    u_ring, section["y_entrance_base"]["numerator_coefficients"]
)) / u_field(polynomial(u_ring, section["y_entrance_base"]["denominator_coefficients"]))
anchor = SourceFileLoader("h92_q6_child_zero_anchor", str(ANCHOR)).load_module()
h92_ring, h92_formulas = anchor.parse_h92(H92)
r92, s92 = anchor.EXPECTED_H92
A1, A, B1_old, B_old, B2_old = tuple(QQ(value(r92, s92)) for value in h92_formulas)
old_a = A1 / u**3 + A / u**4

m = (transport(a1) - T * transport(a0)) / (T * transport(b0) - transport(b1))
radicand = (
    m**4 - 6 * transport(x_p) * m**2 + 8 * transport(y_p) * m
    - 3 * transport(x_p)**2 - 4 * transport(old_a)
)
numerator = old_base_ring(radicand.numerator())
denominator = old_base_ring(radicand.denominator())
quartic, square_factor = squarefree_binary_quartic(radicand, old_base_ring)
assert quartic.degree() == 4
quartic_at_zero = T_field(quartic(u_zero))
assert quartic_at_zero.is_square()
w_zero = quartic_at_zero.sqrt()

# Classical binary-quartic covariants, with coefficients in QQ(T).
binary_ring = PolynomialRing(T_field, names=("x", "z"))
x, z = binary_ring.gens()
binary_quartic = sum(
    binary_ring(quartic[index]) * x**index * z**(4 - index)
    for index in range(5)
)
H = (
    binary_quartic.derivative(x, 2) * binary_quartic.derivative(z, 2)
    - binary_quartic.derivative(x).derivative(z)**2
) / 3
G = (
    binary_quartic.derivative(x) * H.derivative(z)
    - binary_quartic.derivative(z) * H.derivative(x)
)
I = 12 * quartic[4] * quartic[0] - 3 * quartic[3] * quartic[1] + quartic[2]**2
J = (
    72 * quartic[4] * quartic[2] * quartic[0]
    + 9 * quartic[3] * quartic[2] * quartic[1]
    - 27 * quartic[4] * quartic[1]**2
    - 27 * quartic[3]**2 * quartic[0]
    - 2 * quartic[2]**3
)
assert G**2 == -QQ(16) / 3 * H**3 + 256 * I * H * binary_quartic**2 - QQ(1024) / 3 * J * binary_quartic**3
H_zero = T_field(H(x=u_zero, z=1))
G_zero = T_field(G(x=u_zero, z=1))
raw_x = -QQ(3) / 4 * H_zero / quartic_at_zero
raw_y = QQ(9) / 32 * G_zero * w_zero / quartic_at_zero**2
assert raw_y**2 == raw_x**3 - 27 * I * raw_x - 27 * J

minimal_A = T_field([
    QQ(value) for value in child["minimal_short_weierstrass"]["A_coefficients_low_to_high"]
])
minimal_B = T_field([
    QQ(value) for value in child["minimal_short_weierstrass"]["B_coefficients_low_to_high"]
])
fourth_power_ratio = minimal_A / (-27 * I)
minimalizing_unit = T_field(1)
remaining_constant = fourth_power_ratio
for factor, exponent in fourth_power_ratio.factor():
    assert exponent % 4 == 0
    minimalizing_unit *= factor**(exponent // 4)
    remaining_constant /= factor**exponent
assert remaining_constant in QQ and QQ(remaining_constant).nth_root(4)**4 == remaining_constant
minimalizing_unit *= QQ(remaining_constant).nth_root(4)
child_x = minimalizing_unit**2 * raw_x
child_y = minimalizing_unit**3 * raw_y
assert child_y**2 == child_x**3 + minimal_A * child_x + minimal_B

child_x_numerator = T_ring(child_x.numerator())
child_x_denominator = T_ring(child_x.denominator())
child_y_numerator = T_ring(child_y.numerator())
child_y_denominator = T_ring(child_y.denominator())
child_z = child_x_denominator.sqrt()
assert child_z.degree() == 4
assert child_x_denominator == child_z**2
assert child_y_denominator == (
    child_y_denominator.leading_coefficient() / child_z.leading_coefficient()**3
) * child_z**3

# The old III* fibre is the old-base infinity of the binary quartic.  Its
# affine E7 component and its E7_7 exceptional component both have D-degree
# one, so their two generic points are the two signs above binary-quartic
# infinity.  Their resolved-chart sign assignment is a separate certificate;
# here we transport both exact points without guessing a Kodaira chart.
infinity_value = T_field(quartic[4])
assert infinity_value.is_square()
infinity_root = infinity_value.sqrt()


def infinity_coefficient(value, power):
    """Return the coefficient of u^-power in a rational function at infinity."""
    numerator = old_base_ring(value.numerator())
    denominator = old_base_ring(value.denominator())
    assert numerator.degree() + power == denominator.degree()
    return T_field(numerator.leading_coefficient() / denominator.leading_coefficient())


# In the first x-chart of the actual E7 resolution, the old affine component
# is Z=0, U=r^2, Y=r^3.  Hence m=r, while b0 vanishes at infinity.  The
# pencil limit gives r as a rational function of T.  The discriminant square
# root is 2*U-m^2=r^2 at leading order; normalizing it by square_factor gives
# the sign of the binary-quartic point at infinity.
a0_infinity = infinity_coefficient(transport(a0), 0)
a1_infinity = infinity_coefficient(transport(a1), 0)
b1_infinity = infinity_coefficient(transport(b1), 0)
assert infinity_coefficient(transport(b0), 1)
affine_r = (T * a0_infinity - a1_infinity) / b1_infinity
square_factor_infinity = infinity_coefficient(square_factor, 2)
affine_w_infinity = affine_r**2 / square_factor_infinity
assert affine_w_infinity**2 == infinity_value
if affine_w_infinity == infinity_root:
    affine_sign = "plus"
elif affine_w_infinity == -infinity_root:
    affine_sign = "minus"
else:
    raise ArithmeticError("affine E7 chart did not select a quartic-infinity sign")


def transport_binary_quartic_point(x_value, z_value, w_value):
    value = T_field(binary_quartic(x=x_value, z=z_value))
    assert w_value**2 == value
    h_value = T_field(H(x=x_value, z=z_value))
    g_value = T_field(G(x=x_value, z=z_value))
    raw_x_value = -QQ(3) / 4 * h_value / value
    raw_y_value = QQ(9) / 32 * g_value * w_value / value**2
    result_x = minimalizing_unit**2 * raw_x_value
    result_y = minimalizing_unit**3 * raw_y_value
    assert result_y**2 == result_x**3 + minimal_A * result_x + minimal_B
    return result_x, result_y


infinity_sections = []
for sign, w_value in (("plus", infinity_root), ("minus", -infinity_root)):
    infinity_x, infinity_y = transport_binary_quartic_point(T_field(1), T_field(0), w_value)
    infinity_sections.append({
        "binary_quartic_point": "[1:0:{}*sqrt(leading_coefficient)]".format(
            "+" if sign == "plus" else "-"
        ),
        "sign": sign,
        "x_numerator_coefficients_low_to_high": [str(value) for value in T_ring(infinity_x.numerator()).list()],
        "x_denominator_coefficients_low_to_high": [str(value) for value in T_ring(infinity_x.denominator()).list()],
        "y_numerator_coefficients_low_to_high": [str(value) for value in T_ring(infinity_y.numerator()).list()],
        "y_denominator_coefficients_low_to_high": [str(value) for value in T_ring(infinity_y.denominator()).list()],
    })
assert infinity_sections[0] != infinity_sections[1]

if args.component_output:
    component_payload = {
        "schema": "elkies-k3.h92-q6-child-e7-infinity-sections.v1",
        "status": "PASS_EXACT_CHILD_E7_INFINITY_TRANSPORT",
        "inputs": {
            "global_rr": {"path": str(args.rr.relative_to(ROOT)), "sha256": digest(args.rr)},
            "child_jacobian": {"path": str(args.child.relative_to(ROOT)), "sha256": digest(args.child)},
        },
        "source": {
            "old_base": "u=infinity (the III* fibre)",
            "candidate_curves": ["old E7_7", "old affine E7"],
            "affine_E7_chart": "Z=0, U=r^2, Y=r^3 in the first x-chart of the E7 resolution",
            "affine_E7_sign": affine_sign,
            "E7_7_sign": "minus" if affine_sign == "plus" else "plus",
            "boundary": "The affine sign is fixed by its explicit E7 chart and the normalized discriminant square root. The complementary infinity point is E7_7 after the E7 resolution graph certificate, not by a Kodaira-symbol inference.",
        },
        "sections": infinity_sections,
    }
    args.component_output.parent.mkdir(parents=True, exist_ok=True)
    args.component_output.write_text(json.dumps(component_payload, indent=2, sort_keys=True) + "\n")
    print("H92Q6CHILDE7|binary_infinity_points=2|status=PASS_EXACT_CHILD_E7_INFINITY_TRANSPORT", flush=True)

payload = {
    "schema": "elkies-k3.h92-q6-child-zero-section.v1",
    "status": "PASS_EXACT_CHILD_ZERO_SECTION_TRANSPORT",
    "inputs": {
        "global_rr": {"path": str(args.rr.relative_to(ROOT)), "sha256": digest(args.rr)},
        "child_jacobian": {"path": str(args.child.relative_to(ROOT)), "sha256": digest(args.child)},
        "marked_section": {"path": str(args.section.relative_to(ROOT)), "sha256": SECTION_SHA256},
    },
    "old_zero_map": {
        "new_base_value": "T=b1/b0 at old O",
        "u_of_T": str(u_zero),
        "quartic_value_is_square": True,
    },
    "covariant_map": {
        "X": "-3*H/(4*f)",
        "Y": "9*G*w/(32*f^2)",
        "identity": "G^2=-16/3*H^3+256*I*H*f^2-1024/3*J*f^3",
    },
    "section": {
        "x_numerator_coefficients_low_to_high": [str(value) for value in child_x_numerator.list()],
        "x_denominator_coefficients_low_to_high": [str(value) for value in child_x_denominator.list()],
        "y_numerator_coefficients_low_to_high": [str(value) for value in child_y_numerator.list()],
        "y_denominator_coefficients_low_to_high": [str(value) for value in child_y_denominator.list()],
        "denominator_root_degree": int(child_z.degree()),
        "denominator_profile": "x=NX/Z4^2, y=NY/Z4^3",
    },
    "boundary": "This is one transported child section. It does not identify a saturated MW basis, its height, or the next q=8 divisor function.",
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92Q6CHILDZERO|u_of_T_degree=1|quartic_square=1|Z4_degree=4|"
    "status=PASS_EXACT_CHILD_ZERO_SECTION_TRANSPORT",
    flush=True,
)
