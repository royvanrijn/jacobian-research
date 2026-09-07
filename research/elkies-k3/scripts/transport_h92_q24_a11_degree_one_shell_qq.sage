#!/usr/bin/env sage -python
"""Transport the degree-one old curves to the exact A11 Jacobian.

The resolved orbit42 pencil has two exact identity-shell curves whose degree
over the new base is one (equation-shell indices 7 and 17), and one exact
spinor-shell curve of degree one (spinor index 0).  Restricting the exact
pencil to any of these curves gives a Mobius function T(V), hence an exact
rational point (V(T), W(T)) on the binary quartic.  The classical quartic
covariants then give an exact point on

    y^2 = x^3 - 27 I(T) x - 27 J(T).

This is a characteristic-zero construction; it uses no section ansatz,
Groebner basis, modular search, or Hensel lift.  The covariant map has degree
two, so this script does not by itself identify either output with a chosen
MW vector in the pinned A11 zero marking.  It provides the exact geometric
foothold needed for that pointed marking step.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--output",
    type=Path,
    default=LOCAL / "q24-a11-degree-one-shell-covariants-qq.json",
)
args = parser.parse_args()

RR_PATH = LOCAL / "q24-d12-to-a11-orbit42-resolved-rr-qq.json"
PARENT_PATH = LOCAL / "q24-d13-to-d12-component-valuation-qq.json"
CANDIDATES_PATH = LOCAL / "q24-orbit42-exact-section-candidates-qq.json"
ZERO_PATH = LOCAL / "q24-orbit42-rational-zero-pole-sections-qq.json"
SPINOR_PATH = LOCAL / "q24-orbit42-spinor-zero-pole-sections-qq.json"
MODEL_PATH = LOCAL / "q24-orbit42-zero-pole-seeds-mod-43.json"
INPUTS = (RR_PATH, PARENT_PATH, CANDIDATES_PATH, ZERO_PATH, SPINOR_PATH, MODEL_PATH)
for path in INPUTS:
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")

rr = json.loads(RR_PATH.read_text())
parent = json.loads(PARENT_PATH.read_text())
candidates = json.loads(CANDIDATES_PATH.read_text())
zero = json.loads(ZERO_PATH.read_text())
spinor = json.loads(SPINOR_PATH.read_text())
model = json.loads(MODEL_PATH.read_text())["exact_model"]

assert rr["status"] == "PASS_EXACT_Q24_D12_Q6_A11_COMPONENT_VALUATION_RR"
assert parent["status"] == "PASS_EXACT_Q24_D13_TO_D12_COMPONENT_VALUATION_RR"
assert candidates["status"] == "PASS_EXACT_Q42_ORBIT42_SECTION_CANDIDATES_QQ"
assert zero["status"] == "PASS_EXACT_Q42_RATIONAL_ZERO_POLE_SECTIONS_QQ"
assert spinor["status"] == "PASS_EXACT_Q42_SPINOR_ZERO_POLE_SECTIONS_QQ"

UQ = PolynomialRing(QQ, "u")
u = UQ.gen()
KU = UQ.fraction_field()
VQ = PolynomialRing(QQ, "V")
V = VQ.gen()
KV = VQ.fraction_field()
TQ = PolynomialRing(QQ, "T")
T = TQ.gen()
KT = TQ.fraction_field()
WV = PolynomialRing(KT, "W")
W = WV.gen()


def evaluate_u(value, argument):
    value = KU(value)
    return KV(VQ(value.numerator())(argument)) / KV(
        VQ(value.denominator())(argument)
    )


def rational_square_root(value, polynomial_ring):
    """Return the exact rational-function square root, or fail."""

    value = polynomial_ring.fraction_field()(value)
    numerator = polynomial_ring(value.numerator())
    denominator = polynomial_ring(value.denominator())
    answer = value.parent().one()
    for polynomial, direction in ((numerator, 1), (denominator, -1)):
        if not polynomial:
            return value.parent().zero()
        leading = QQ(polynomial.leading_coefficient())
        if not leading.is_square():
            raise ArithmeticError("rational-function square has nonsquare leading coefficient")
        root = value.parent()(leading.sqrt())
        for factor, multiplicity in (polynomial / leading).factor():
            if int(multiplicity) % 2:
                raise ArithmeticError("rational-function value is not a square")
            root *= value.parent()(factor.monic()) ** (int(multiplicity) // 2)
        answer = answer * root if direction == 1 else answer / root
    if answer**2 != value:
        raise ArithmeticError("rational square-root reconstruction failed")
    return answer


def normalized_rational_record(value, polynomial_ring):
    value = polynomial_ring.fraction_field()(value)
    numerator = polynomial_ring(value.numerator())
    denominator = polynomial_ring(value.denominator())
    scale = denominator.leading_coefficient()
    numerator /= scale
    denominator /= scale
    return {
        "numerator_coefficients_low_to_high": [str(item) for item in numerator.list()],
        "denominator_coefficients_low_to_high": [str(item) for item in denominator.list()],
        "numerator_degree": int(numerator.degree()),
        "denominator_degree": int(denominator.degree()),
    }


u_of_V = KV(rr["coordinate_change"]["u_of_V"])
x_scale = KV(rr["coordinate_change"]["x_scale"])
y_scale = KV(rr["coordinate_change"]["y_scale"])
A_parent = VQ(
    [QQ(value) for value in parent["child"]["minimal_A_coefficients_low_to_high"]]
)
B_parent = VQ(
    [QQ(value) for value in parent["child"]["minimal_B_coefficients_low_to_high"]]
)

selected = candidates["candidates"][rr["selected_section"]["candidate_index"]]
X_shell = UQ([QQ(value) for value in selected["X_coefficients_low_to_high"]])
Y_shell = UQ([QQ(value) for value in selected["Y_coefficients_low_to_high"]])
Z_shell = UQ([QQ(value) for value in selected["Z_coefficients_low_to_high"]])
x_parent = x_scale * evaluate_u(KU(X_shell) / KU(Z_shell**2), u_of_V)
y_parent = y_scale * evaluate_u(KU(Y_shell) / KU(Z_shell**3), u_of_V)
E_parent = EllipticCurve(KV, [0, 0, 0, KV(A_parent), KV(B_parent)])
P_parent = E_parent(x_parent, y_parent)

Z_denominator = VQ(x_parent.denominator()).monic()
Z_parent = VQ.one()
for factor, multiplicity in Z_denominator.factor():
    if int(multiplicity) % 2:
        raise ArithmeticError("selected parent-section denominator is not a square")
    Z_parent *= factor.monic() ** (int(multiplicity) // 2)
X_parent = VQ(x_parent * Z_parent**2)
Y_parent = VQ(y_parent * Z_parent**3)
assert E_parent(KV(X_parent) / KV(Z_parent**2), KV(Y_parent) / KV(Z_parent**3)) == P_parent

alpha = QQ(model["I8star_root"])
collision_modulus = Z_parent**2
X_inverse = VQ(X_parent.inverse_mod(collision_modulus))
rr_pairs = []
for BB in (VQ.one(), V):
    AA = VQ((BB * Y_parent * X_inverse) % collision_modulus)
    AA -= AA(alpha) / Z_parent(alpha) ** 2 * Z_parent**2
    assert AA(alpha) == 0
    assert (AA * X_parent - BB * Y_parent) % collision_modulus == 0
    rr_pairs.append((AA, BB))

(AA0, BB0), (AA1, BB1) = rr_pairs
a0, b0 = KV(AA0) / KV(Z_parent**2), KV(BB0) / KV(Z_parent)
a1, b1 = KV(AA1) / KV(Z_parent**2), KV(BB1) / KV(Z_parent)

quartic = WV(
    [KT(TQ(value)) for value in rr["quartic"]["coefficients_in_T_low_to_high"]]
)
if quartic.degree() != 4:
    raise ArithmeticError("stored exact quartic is not degree four")
e, d, c, b, a = quartic.list()
I = 12 * a * e - 3 * b * d + c**2
J = 72 * a * c * e + 9 * b * c * d - 27 * a * d**2 - 27 * b**2 * e - 2 * c**3
A_child = KT(-27 * I)
B_child = KT(-27 * J)
assert A_child == KT(TQ([QQ(value) for value in rr["child"]["minimal_A_coefficients_low_to_high"]]))
assert B_child == KT(TQ([QQ(value) for value in rr["child"]["minimal_B_coefficients_low_to_high"]]))
E_child = EllipticCurve(KT, [0, 0, 0, A_child, B_child])


def covariants_at(x_value):
    g0 = b**2 / 16 - a * c / 6
    g1 = b * c / 12 - a * d / 2
    g2 = c**2 / 12 - b * d / 8 - a * e
    g3 = c * d / 12 - b * e / 2
    g4 = d**2 / 16 - c * e / 6
    ux = 4 * a * x_value**3 + 3 * b * x_value**2 + 2 * c * x_value + d
    uy = b * x_value**3 + 2 * c * x_value**2 + 3 * d * x_value + 4 * e
    g = g0 * x_value**4 + g1 * x_value**3 + g2 * x_value**2 + g3 * x_value + g4
    gx = 4 * g0 * x_value**3 + 3 * g1 * x_value**2 + 2 * g2 * x_value + g3
    gy = g1 * x_value**3 + 2 * g2 * x_value**2 + 3 * g3 * x_value + 4 * g4
    h = (ux * gy - uy * gx) / 8
    return g, h


targets = (
    ("identity", 7, zero["sections"][7]),
    ("identity", 17, zero["sections"][17]),
    ("spinor", 0, spinor["sections"][0]),
)
rows = []
for shell_kind, section_index, section in targets:
    x_shell_curve = UQ([QQ(value) for value in section["x_coefficients_low_to_high"]])
    y_shell_curve = UQ([QQ(value) for value in section["y_coefficients_low_to_high"]])
    x_curve = x_scale * evaluate_u(KU(x_shell_curve), u_of_V)
    y_curve = y_scale * evaluate_u(KU(y_shell_curve), u_of_V)
    assert y_curve**2 == x_curve**3 + KV(A_parent) * x_curve + KV(B_parent)

    chord = (y_curve + y_parent) / (x_curve - x_parent)
    new_base = (a1 + b1 * chord) / (a0 + b0 * chord)
    numerator = VQ(new_base.numerator())
    denominator = VQ(new_base.denominator())
    common = numerator.gcd(denominator)
    numerator //= common
    denominator //= common
    if max(numerator.degree(), denominator.degree()) != 1:
        raise ArithmeticError(f"shell index {section_index} is not degree one")

    inverse_equation = WV(KT(T) * WV(denominator) - WV(numerator))
    if inverse_equation.degree() != 1:
        raise ArithmeticError("Mobius inverse did not remain linear")
    V_of_T = KT(-inverse_equation[0] / inverse_equation[1])
    if KT(new_base.numerator()(V_of_T)) / KT(new_base.denominator()(V_of_T)) != KT(T):
        raise ArithmeticError("Mobius inversion failed")

    quartic_value = KT(quartic(V_of_T))
    ordinate = rational_square_root(quartic_value, TQ)
    g_value, h_value = covariants_at(V_of_T)
    x_jacobian = KT(36 * g_value / ordinate**2)
    y_jacobian = KT(108 * h_value / ordinate**3)
    point = E_child(x_jacobian, y_jacobian)
    assert point[1] ** 2 == point[0] ** 3 + A_child * point[0] + B_child

    # Point the quartic at (V_of_T,+ordinate).  The opposite ordinate at the
    # same V is the residual chord curve.  Unlike the covariant image, its
    # image below is obtained through a degree-one pointed isomorphism.
    shifted = WV(quartic(W + V_of_T))
    shifted_e, shifted_d, shifted_c, shifted_b, shifted_a = shifted.list()
    if shifted_e != ordinate**2:
        raise ArithmeticError("pointed quartic constant term is not the chosen square")
    a1_pointed = shifted_d / ordinate
    a2_pointed = shifted_c - shifted_d**2 / (4 * ordinate**2)
    a3_pointed = 2 * ordinate * shifted_b
    a4_pointed = -4 * ordinate**2 * shifted_a
    a6_pointed = a2_pointed * a4_pointed
    b2_pointed = a1_pointed**2 + 4 * a2_pointed
    b4_pointed = 2 * a4_pointed + a1_pointed * a3_pointed
    b6_pointed = a3_pointed**2 + 4 * a6_pointed
    c4_pointed = b2_pointed**2 - 24 * b4_pointed
    c6_pointed = -b2_pointed**3 + 36 * b2_pointed * b4_pointed - 216 * b6_pointed
    short_A = -c4_pointed / 48
    short_B = -c6_pointed / 864
    if A_child != 81 * short_A or B_child != 729 * short_B:
        raise ArithmeticError("pointed quartic and invariant Jacobian normalizations disagree")

    x_general = -a2_pointed
    y_general = a1_pointed * a2_pointed - a3_pointed
    if (
        y_general**2 + a1_pointed * x_general * y_general + a3_pointed * y_general
        != x_general**3
        + a2_pointed * x_general**2
        + a4_pointed * x_general
        + a6_pointed
    ):
        raise ArithmeticError("opposite point missed the pointed generalized model")
    x_short = x_general + b2_pointed / 12
    y_short = y_general + (a1_pointed * x_general + a3_pointed) / 2
    x_opposite = KT(9 * x_short)
    y_opposite = KT(27 * y_short)
    opposite_point = E_child(x_opposite, y_opposite)
    assert opposite_point[1] ** 2 == opposite_point[0] ** 3 + A_child * opposite_point[0] + B_child

    rows.append(
        {
            "shell_kind": shell_kind,
            "equation_shell_index": section_index,
            "old_section_pair_index": (
                int(section["pair_index"]) if "pair_index" in section else None
            ),
            "old_section_sign": int(section["sign"]),
            "new_base_degree": 1,
            "T_of_V": normalized_rational_record(new_base, VQ),
            "V_of_T": normalized_rational_record(V_of_T, TQ),
            "quartic_ordinate": normalized_rational_record(ordinate, TQ),
            "covariant_point": {
                "x": normalized_rational_record(x_jacobian, TQ),
                "y": normalized_rational_record(y_jacobian, TQ),
            },
            "pointed_opposite_section": {
                "chosen_quartic_zero_ordinate_sign": 1,
                "opposite_ordinate_sign": -1,
                "map_degree": 1,
                "equals_negative_covariant_point": bool(
                    x_opposite == x_jacobian and y_opposite == -y_jacobian
                ),
            },
            "exact_quartic_identity": True,
            "exact_child_identity": True,
            "exact_pointed_opposite_identity": True,
        }
    )

payload = {
    "schema": "elkies-k3.h3-q24-a11-degree-one-curve-covariants-qq.v2",
    "status": "PASS_EXACT_A11_DEGREE_ONE_CURVE_COVARIANTS_QQ",
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in INPUTS
        },
    },
    "degree_one_shell_indices": [7, 17],
    "degree_one_spinor_indices": [0],
    "points": rows,
    "method": (
        "exact restriction of the resolved RR pencil, Mobius inversion, exact quartic "
        "square root, classical degree-two covariant, and degree-one pointed-quartic "
        "transport of the opposite residual chord curve"
    ),
    "proof_boundary": (
        "All three covariant images and their pointed opposite sections are exact "
        "characteristic-zero points on the certified A11 Jacobian. The covariant has "
        "degree two, while each opposite section uses a degree-one pointed quartic "
        "isomorphism with its corresponding shell curve as zero. Identifying the "
        "resulting zero translations with the pinned A11 MW marking is the next step."
    ),
}

args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "A11SHELLQQ|identity_indices=7,17|spinor_indices=0|covariant_points=3|"
    "pointed_opposites=3|"
    "method=mobius+quartic_covariant+pointed_opposite|"
    f"status={payload['status']}",
    flush=True,
)
print(f"OUTPUT|{args.output}", flush=True)
