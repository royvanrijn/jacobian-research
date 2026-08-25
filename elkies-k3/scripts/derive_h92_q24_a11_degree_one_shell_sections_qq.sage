#!/usr/bin/env sage -python
"""Point the exact A11 quartic at its fixed equation zero.

The degree-one shell curves S7 and S17 already give exact rational points on
the D12-to-A11 binary quartic.  The covariant map doubles those points.  This
script instead applies the degree-one pointed-quartic isomorphism based at the
fixed old-I8* point used for the A11 equation.  It therefore recovers the two
shell sections themselves, without a section ansatz, Groebner basis, or
Hensel lift.

Both ordinate signs are checked exactly.  Reduction at the pinned good prime
selects the sign agreeing with the independently reconstructed marked trace.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--covariants",
    type=Path,
    default=LOCAL / "q24-a11-degree-one-shell-covariants-qq.json",
)
parser.add_argument(
    "--marked-reduction",
    type=Path,
    default=LOCAL / "q24-a11-pinned-zero-section-mod100003.json",
)
parser.add_argument(
    "--output",
    type=Path,
    default=LOCAL / "q24-a11-degree-one-shell-sections-qq.json",
)
args = parser.parse_args()

A11 = LOCAL / "q24-d12-to-a11-orbit42-resolved-rr-qq.json"
PARENT = LOCAL / "q24-d13-to-d12-component-valuation-qq.json"
COVARIANTS = args.covariants.resolve()
MARKED = args.marked_reduction.resolve()
INPUTS = (A11, PARENT, COVARIANTS, MARKED)
for path in INPUTS:
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")

a11 = json.loads(A11.read_text())
parent = json.loads(PARENT.read_text())
covariants = json.loads(COVARIANTS.read_text())
marked = json.loads(MARKED.read_text())
assert a11["status"] == "PASS_EXACT_Q24_D12_Q6_A11_COMPONENT_VALUATION_RR"
assert parent["status"] == "PASS_EXACT_Q24_D13_TO_D12_COMPONENT_VALUATION_RR"
assert covariants["status"] == "PASS_EXACT_A11_DEGREE_ONE_CURVE_COVARIANTS_QQ"
assert marked["status"] == "PASS_Q24_A11_PINNED_ZERO_SECTION_RECONSTRUCTION_MODP"

TQ = PolynomialRing(QQ, "T")
T = TQ.gen()
KT = TQ.fraction_field()
VQ = PolynomialRing(KT, "V")
V = VQ.gen()


def rational_from_record(record):
    numerator = TQ([QQ(value) for value in record["numerator_coefficients_low_to_high"]])
    denominator = TQ(
        [QQ(value) for value in record["denominator_coefficients_low_to_high"]]
    )
    return KT(numerator) / KT(denominator)


def power_root(poly, exponent):
    poly = TQ(poly)
    leading = QQ(poly.leading_coefficient())
    if exponent % 2 == 0 and not leading.is_square():
        raise ArithmeticError("denominator leading coefficient is not a square")
    answer = TQ.one()
    for factor, multiplicity in (poly / leading).factor():
        if int(multiplicity) % exponent:
            raise ArithmeticError("section denominator is not a perfect power")
        answer *= factor.monic() ** (int(multiplicity) // exponent)
    return answer.monic()


def normalized_section(point, A, B):
    x, y = map(KT, point.xy())
    Zx = power_root(x.denominator(), 2)
    Zy = power_root(y.denominator(), 3)
    if Zx != Zy:
        raise ArithmeticError("x and y denominators give different section Z")
    Z = Zx
    X = TQ(x * Z**2)
    Y = TQ(y * Z**3)
    if Y**2 != X**3 + A * X * Z**4 + B * Z**6:
        raise ArithmeticError("normalized section missed the exact A11 equation")
    return X, Y, Z


def section_record(section):
    X, Y, Z = section
    return {
        "X_coefficients_low_to_high": [str(value) for value in X.list()],
        "Y_coefficients_low_to_high": [str(value) for value in Y.list()],
        "Z_coefficients_low_to_high": [str(value) for value in Z.list()],
        "degrees_X_Y_Z": [int(X.degree()), int(Y.degree()), int(Z.degree())],
        "exact_weierstrass_identity": True,
    }


coefficients = [KT(TQ(value)) for value in a11["quartic"]["coefficients_in_T_low_to_high"]]
quartic = sum(coefficients[index] * V**index for index in range(5))
e0, d0, c0, b0, a0 = coefficients

i8 = next(row for row in parent["child"]["finite_fibres"] if row["kodaira"] == "I8*")
V0Q = PolynomialRing(QQ, "V")
i8_factor = V0Q(str(i8["factor"]))
if i8_factor.degree() != 1:
    raise ArithmeticError("old I8* factor is not linear")
alpha = -KT(i8_factor[0]) / KT(i8_factor[1])
q_squared = KT(quartic(alpha))
if not q_squared.is_square():
    raise ArithmeticError("fixed pointed quartic value is not a square")
q = KT(q_squared.sqrt())

# Coefficients after V=alpha+u.
a = a0
b = b0 + 4 * alpha * a0
c = c0 + 3 * alpha * b0 + 6 * alpha**2 * a0
d = d0 + 2 * alpha * c0 + 3 * alpha**2 * b0 + 4 * alpha**3 * a0
a1 = d / q
a2 = c - d**2 / (4 * q**2)
a3 = 2 * q * b
a4 = -4 * q**2 * a
b2 = a1**2 + 4 * a2

A = TQ([QQ(value) for value in a11["child"]["minimal_A_coefficients_low_to_high"]])
B = TQ([QQ(value) for value in a11["child"]["minimal_B_coefficients_low_to_high"]])
E = EllipticCurve(KT, [0, 0, 0, KT(A), KT(B)])


def fixed_pointed_image(v_value, w_value):
    """The degree-one quartic map taking (alpha,q) to infinity."""

    u = v_value - alpha
    if not u:
        raise ArithmeticError("shell point coincides identically with the fixed origin")
    x_general = (2 * q * (w_value + q) + d * u) / u**2
    y_general = (
        4 * q**2 * (w_value + q)
        + 2 * q * d * u
        + (2 * q * c - d**2 / (2 * q)) * u**2
    ) / u**3
    if (
        y_general**2 + a1 * x_general * y_general + a3 * y_general
        != x_general**3 + a2 * x_general**2 + a4 * x_general + a2 * a4
    ):
        raise ArithmeticError("pointed image missed the generalized model")
    x = KT(9 * (x_general + b2 / 12))
    y = KT(27 * (y_general + (a1 * x_general + a3) / 2))
    return E(x, y)


p = ZZ(marked["prime"])
F = GF(p)
RF = PolynomialRing(F, "T")
KF = RF.fraction_field()
AF = RF([F(value.numerator()) / F(value.denominator()) for value in A.list()])
BF = RF([F(value.numerator()) / F(value.denominator()) for value in B.list()])
EF = EllipticCurve(KF, [0, 0, 0, KF(AF), KF(BF)])


def reduce_section(section):
    X, Y, Z = section
    Xp = RF([F(value.numerator()) / F(value.denominator()) for value in X.list()])
    Yp = RF([F(value.numerator()) / F(value.denominator()) for value in Y.list()])
    Zp = RF([F(value.numerator()) / F(value.denominator()) for value in Z.list()])
    return EF(KF(Xp) / KF(Zp**2), KF(Yp) / KF(Zp**3))


def marked_point(index):
    row = marked["selected_trace_sections"][f"S{index}"]
    X = RF(row["X_coefficients_low_to_high"])
    Y = RF(row["Y_coefficients_low_to_high"])
    Z = RF(row["Z_coefficients_low_to_high"])
    return EF(KF(X) / KF(Z**2), KF(Y) / KF(Z**3))


rows = []
for source in covariants["points"]:
    index = int(source["equation_shell_index"])
    if source["shell_kind"] != "identity" or index not in (7, 17):
        continue
    v_value = rational_from_record(source["V_of_T"])
    ordinate = rational_from_record(source["quartic_ordinate"])
    if ordinate**2 != KT(quartic(v_value)):
        raise ArithmeticError(f"stored S{index} ordinate missed the quartic")
    candidates = []
    for sign in (1, -1):
        point = fixed_pointed_image(v_value, sign * ordinate)
        for equation_sign in (1, -1):
            signed_point = equation_sign * point
            section = normalized_section(signed_point, A, B)
            candidates.append(
                (sign, equation_sign, section, reduce_section(section))
            )
    target = marked_point(index)
    selected = [
        (sign, equation_sign, section)
        for sign, equation_sign, section, reduction in candidates
        if reduction == target
    ]
    if len(selected) != 1:
        raise ArithmeticError(f"S{index} has {len(selected)} marked sign matches modulo {p}")
    selected_sign, equation_sign, selected_section = selected[0]
    rows.append(
        {
            "equation_shell_index": index,
            "selected_quartic_ordinate_sign": int(selected_sign),
            "pointed_image_equation_sign": int(equation_sign),
            "section": section_record(selected_section),
            "exact_source_quartic_identity": True,
            "exact_pointed_generalized_identity": True,
            "exact_child_identity": True,
            "marked_reduction_prime": int(p),
            "exact_reduction_matches_marked_trace": True,
        }
    )

if [row["equation_shell_index"] for row in rows] != [7, 17]:
    raise ArithmeticError("did not recover both degree-one identity shells")

payload = {
    "schema": "elkies-k3.h3-q24-a11-degree-one-shell-sections-qq.v1",
    "status": "PASS_EXACT_MARKED_A11_DEGREE_ONE_SHELL_SECTIONS_QQ",
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in INPUTS
        },
    },
    "fixed_quartic_origin": [str(alpha), str(q)],
    "sections": rows,
    "method": (
        "exact degree-one pointed-quartic isomorphism based at the fixed old-I8* "
        "point, with sign selected by exact reduction at one independently marked good prime"
    ),
    "large_Groebner_required": False,
    "proof_boundary": (
        "This exactly constructs and marks the A11 equation sections S7 and S17. "
        "It does not yet construct the degree-three S5 trace or the full q8 H0."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print(
    "A11DEG1SECTIONS|{}|status={}".format(
        "|".join(
            "S{}=sign{},degrees{}".format(
                row["equation_shell_index"],
                "{}/{}".format(
                    row["selected_quartic_ordinate_sign"],
                    row["pointed_image_equation_sign"],
                ),
                tuple(row["section"]["degrees_X_Y_Z"]),
            )
            for row in rows
        ),
        payload["status"],
    ),
    flush=True,
)
print(f"OUTPUT|{args.output.resolve()}", flush=True)
