#!/usr/bin/env sage -python
"""Transport the lifted q4/o1584 degree-one curves to q4/orbit164 over QQ.

Restrict the certified q4/o164 chord pencil to both signs of the four lifted
q4/o1584 curves.  Degree-one restrictions are inverted by a Mobius transform,
giving exact points on the stored binary quartic.  The quartic is then pointed
at the already-certified C8 zero and converted by the standard degree-one
pointed-quartic map to the stored short q4/o164 Jacobian.  This preserves the
primitive zero marking; no section ansatz, elimination, or Groebner basis is
used on the child.
"""

import hashlib
import json
import time
from pathlib import Path

from sage.all import EllipticCurve, PolynomialRing, QQ, ZZ, factorial


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
PARENT = LOCAL / "q4o208-physical-q4o1584-rr-qq.json"
PARENT_MARKING = LOCAL / "q4o1584-second-affine-equation-marking-qq.json"
COMPACT = LOCAL / "q4o1584-compact-weierstrass-qq.json"
SECTIONS = LOCAL / "q4o1584-degree1-all-node-sections-qq.json"
Q164 = LOCAL / "q4o1584-physical-q4o164-rr-qq.json"
C8 = LOCAL / "q4o164-c8-equation-marking-qq.json"
Q164_COMPACT = LOCAL / "q4o164-compact-weierstrass-qq.json"
OUTPUT = LOCAL / "q4o1584-degree1-sections-to-q4o164-qq.json"
INPUTS = (PARENT, PARENT_MARKING, COMPACT, SECTIONS, Q164, C8, Q164_COMPACT)

started = time.monotonic()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rational_bits(values):
    entries = [QQ(value) for value in values]
    return max(
        max(abs(ZZ(value.numerator())).nbits(), ZZ(value.denominator()).nbits())
        for value in entries
    )


def record(value, polynomial_ring):
    value = polynomial_ring.fraction_field()(value)
    numerator = polynomial_ring(value.numerator())
    denominator = polynomial_ring(value.denominator())
    common = numerator.gcd(denominator)
    numerator //= common
    denominator //= common
    scale = denominator.leading_coefficient()
    numerator /= scale
    denominator /= scale
    return {
        "numerator_coefficients_low_to_high": [str(entry) for entry in numerator.list()],
        "denominator_coefficients_low_to_high": [str(entry) for entry in denominator.list()],
        "degrees_numerator_denominator": [int(numerator.degree()), int(denominator.degree())],
    }


parent = json.loads(PARENT.read_text())
parent_marking = json.loads(PARENT_MARKING.read_text())
compact = json.loads(COMPACT.read_text())
sections = json.loads(SECTIONS.read_text())
q164 = json.loads(Q164.read_text())
c8 = json.loads(C8.read_text())
q164_compact = json.loads(Q164_COMPACT.read_text())
assert parent["status"] == "PASS_EXACT_QQ_PHYSICAL_Q4O1584_D4_A3_3A1_RR_AND_JACOBIAN"
assert parent_marking["status"] == "PASS_EXACT_QQ_Q4O1584_SECOND_AFFINE_POINTING_AND_C0_MARKING"
assert compact["status"] == "PASS_EXACT_QQ_Q4O1584_COMPACT_WEIERSTRASS_NORMALIZATION"
assert sections["status"] == "PASS_EXACT_QQ_Q4O1584_FOUR_DEGREE1_ALL_NODE_SECTION_PAIRS"
assert q164["status"] == "PASS_EXACT_QQ_PHYSICAL_Q4O164_2A3_2A1_RR_AND_JACOBIAN"
assert c8["status"] == "PASS_EXACT_QQ_Q4O164_C8_EQUATION_MARKING"
assert q164_compact["status"] == "PASS_EXACT_QQ_Q4O164_COMPACT_WEIERSTRASS_NORMALIZATION"

RV = PolynomialRing(QQ, "V")
V = RV.gen()
KV = RV.fraction_field()
RU = PolynomialRing(QQ, "U")
U = RU.gen()
KU = RU.fraction_field()
RT = PolynomialRing(QQ, "T")
T = RT.gen()

A_parent = RT([QQ(value) for value in parent["child"]["minimal_A_coefficients_low_to_high"]])
B_parent = RT([QQ(value) for value in parent["child"]["minimal_B_coefficients_low_to_high"]])


def read_parent_polynomial(point_record):
    numerator = RT([QQ(value) for value in point_record["numerator_coefficients_low_to_high"]])
    denominator = RT([QQ(value) for value in point_record["denominator_coefficients_low_to_high"]])
    quotient = numerator / denominator
    assert quotient.denominator() == 1
    return RT(quotient)


c0 = parent_marking["old_A11_component_0_on_second_affine_pointed_child"]
PX = read_parent_polynomial(c0["x"])
PY = read_parent_polynomial(c0["y"])
assert PY**2 == PX**3 + A_parent * PX + B_parent

change = compact["exact_coordinate_change"]
base_a = QQ(change["a"])
base_c = QQ(change["c"])
scale_m = QQ(change["m"])
T_of_V = base_a + base_c * V

linear_T = RT([QQ(value) for value in q164["resolved_RR"]["interpolating_linear_polynomial"]])
denominator_T = RT([QQ(value) for value in q164["resolved_RR"]["common_denominator"]])
quartic_coefficients = [
    RU([QQ(value) for value in coefficients])
    for coefficients in q164["quartic"]["coefficients_in_T_low_to_high"]
]


def evaluate_polynomial(polynomial, argument):
    target = argument.parent()
    return sum(
        target(QQ(coefficient)) * argument**index
        for index, coefficient in enumerate(polynomial)
    )


def curve_coordinates(section, sign):
    Z = RV([QQ(value) for value in section["Z_coefficients_low_to_high"]])
    X = RV([QQ(value) for value in section["X_coefficients_low_to_high"]])
    Y = sign * RV([QQ(value) for value in section["Y_coefficients_low_to_high"]])
    x_compact = KV(X) / KV(Z**2)
    y_compact = KV(Y) / KV(Z**3)
    x_parent = scale_m**2 * x_compact
    y_parent = scale_m**3 * y_compact
    assert y_parent**2 == x_parent**3 + evaluate_polynomial(A_parent, T_of_V) * x_parent + evaluate_polynomial(B_parent, T_of_V)
    return x_parent, y_parent


restrictions = []
for pair_index, section in enumerate(sections["resolved_hensel"]["sections"]):
    for sign in (1, -1):
        x_parent, y_parent = curve_coordinates(section, sign)
        px = evaluate_polynomial(PX, T_of_V)
        py = evaluate_polynomial(PY, T_of_V)
        chord = (y_parent + py) / (x_parent - px)
        new_base = (
            chord - evaluate_polynomial(linear_T, T_of_V)
        ) / evaluate_polynomial(denominator_T, T_of_V)
        numerator = RV(new_base.numerator())
        denominator = RV(new_base.denominator())
        common = numerator.gcd(denominator)
        numerator //= common
        denominator //= common
        degree = max(int(numerator.degree()), int(denominator.degree()))
        W = (
            2 * x_parent + px - chord**2
        ) / evaluate_polynomial(denominator_T, T_of_V)

        # Verify the surface-to-quartic map before selecting degree-one curves.
        quartic_on_curve = sum(
            evaluate_polynomial(coefficient, new_base) * KV(T_of_V)**index
            for index, coefficient in enumerate(quartic_coefficients)
        )
        assert W**2 == quartic_on_curve
        restrictions.append({
            "pair_index": pair_index,
            "sign": sign,
            "new_base_degree": degree,
            "new_base": record(new_base, RV),
            "quartic_identity_before_inversion": True,
            "_numerator": numerator,
            "_denominator": denominator,
            "_W": W,
        })
        print(
            f"Q4O164D1TRANSPORT|pair={pair_index}|sign={sign}|degree={degree}|"
            f"elapsed={time.monotonic()-started:.3f}",
            flush=True,
        )

degree_one = [row for row in restrictions if row["new_base_degree"] == 1]
assert len(degree_one) == 2

# The certified q4/o164 zero C8 is the point (T,W)=(t0,W0(U)).
t0 = QQ(c8["selected_zero"]["parent_base_support"])
W0_record = c8["selected_zero"]["quartic_ordinate"]
W0 = KU(RU([QQ(value) for value in W0_record["numerator_coefficients_low_to_high"]])) / KU(
    RU([QQ(value) for value in W0_record["denominator_coefficients_low_to_high"]])
)

quartic_KU = [KU(value) for value in quartic_coefficients]
translated = []
for order in range(5):
    translated.append(sum(
        coefficient * QQ(factorial(index)) / QQ(factorial(order) * factorial(index-order)) * KU(t0)**(index-order)
        for index, coefficient in enumerate(quartic_KU)
        if index >= order
    ))
e, d, c, b, a = translated
assert e == W0**2
a1 = d / W0
a2 = c - d**2 / (4 * W0**2)
a3 = 2 * W0 * b
a4 = -4 * W0**2 * a
a6 = a2 * a4
b2 = a1**2 + 4 * a2
b4 = 2 * a4 + a1 * a3
b6 = a3**2 + 4 * a6
c4 = b2**2 - 24 * b4
c6 = -b2**3 + 36 * b2 * b4 - 216 * b6
A_child = KU(RU([QQ(value) for value in q164["child"]["minimal_A_coefficients_low_to_high"]]))
B_child = KU(RU([QQ(value) for value in q164["child"]["minimal_B_coefficients_low_to_high"]]))
assert 81 * (-c4 / 48) == A_child
assert 729 * (-c6 / 864) == B_child
E_child = EllipticCurve(KU, [0, 0, 0, A_child, B_child])
child_c = QQ(q164_compact["exact_coordinate_change"]["c"])
child_s = QQ(q164_compact["exact_coordinate_change"]["s"])
A_child_compact = RU([QQ(value) for value in q164_compact["compact_model"]["A_coefficients_low_to_high"]])
B_child_compact = RU([QQ(value) for value in q164_compact["compact_model"]["B_coefficients_low_to_high"]])

transported = []
for row in degree_one:
    numerator = row.pop("_numerator")
    denominator = row.pop("_denominator")
    W_of_V = row.pop("_W")
    inverse = RU([KU(-numerator[0]), KU(denominator[0])])
    # Solve U*denominator(V)-numerator(V)=0 for V over QQ(U).
    constant = KU(U * denominator[0] - numerator[0])
    linear_coefficient = KU(U * denominator[1] - numerator[1])
    V_of_U = -constant / linear_coefficient
    assert evaluate_polynomial(numerator, V_of_U) / evaluate_polynomial(denominator, V_of_U) == KU(U)
    T_curve = KU(base_a + base_c * V_of_U)
    W_curve = KU(RV(W_of_V.numerator())(V_of_U)) / KU(RV(W_of_V.denominator())(V_of_U))
    quartic_value = sum(coefficient * T_curve**index for index, coefficient in enumerate(quartic_KU))
    assert W_curve**2 == quartic_value

    local_x = T_curve - t0
    if not local_x:
        raise ArithmeticError("degree-one curve coincides with the C8 zero")
    x_general = (2 * W0 * (W_curve + W0) + d * local_x) / local_x**2
    y_general = (
        4 * W0**2 * (W_curve + W0)
        + 2 * W0 * (d * local_x + c * local_x**2)
        - d**2 * local_x**2 / (2 * W0)
    ) / local_x**3
    assert (
        y_general**2 + a1 * x_general * y_general + a3 * y_general
        == x_general**3 + a2 * x_general**2 + a4 * x_general + a6
    )
    x_short = KU(9 * (x_general + b2 / 12))
    y_short = KU(27 * (y_general + (a1 * x_general + a3) / 2))
    point = E_child(x_short, y_short)
    assert point[1]**2 == point[0]**3 + A_child * point[0] + B_child
    x_compact = KU(RU(x_short.numerator())(child_c * U)) / (
        KU(RU(x_short.denominator())(child_c * U)) * child_s**2
    )
    y_compact = KU(RU(y_short.numerator())(child_c * U)) / (
        KU(RU(y_short.denominator())(child_c * U)) * child_s**3
    )
    assert y_compact**2 == x_compact**3 + KU(A_child_compact) * x_compact + KU(B_child_compact)
    row.update({
        "old_base_as_function_of_q4o164_base": record(T_curve, RU),
        "quartic_ordinate_as_function_of_q4o164_base": record(W_curve, RU),
        "compact_child_x": record(x_compact, RU),
        "compact_child_y": record(y_compact, RU),
        "pointed_map_degree": 1,
        "selected_zero": "old_A11_component_8",
        "exact_pointed_generalized_identity": True,
        "exact_short_child_identity": True,
        "maximum_compact_child_coordinate_rational_bits": rational_bits(
            list(x_compact.numerator()) + list(x_compact.denominator())
            + list(y_compact.numerator()) + list(y_compact.denominator())
        ),
    })
    transported.append(row)

for row in restrictions:
    row.pop("_numerator", None)
    row.pop("_denominator", None)
    row.pop("_W", None)

payload = {
    "schema": "elkies-k3.q4o1584-degree1-sections-to-q4o164-qq.v1",
    "status": "PASS_EXACT_QQ_Q4O164_TWO_PRIMITIVE_ONE_NODE_SECTIONS",
    "restriction_degrees": [
        {key: row[key] for key in ("pair_index", "sign", "new_base_degree")}
        for row in restrictions
    ],
    "degree_one_count": len(degree_one),
    "degree_one_sections": transported,
    "coefficient_growth": {
        "maximum_compact_child_coordinate_rational_bits": max(
            row["maximum_compact_child_coordinate_rational_bits"] for row in transported
        ),
    },
    "method": {
        "large_Groebner_required": False,
        "child_section_ansatz_required": False,
        "construction": (
            "exact chord-pencil restriction, Mobius inversion, and degree-one quartic "
            "pointing at the certified C8 zero"
        ),
        "runtime_seconds": time.monotonic() - started,
    },
    "proof_boundary": (
        "Both displayed points are exact primitive zero-marked sections on the certified "
        "q4/orbit164 short model. Their precise NS lattice labels and their independence "
        "from the existing rank-eight pair-node subgroup remain separate checks."
    ),
    "next_required": (
        "Match the two sections to the enumerated one-node classes and select the sign that "
        "extends the exact q4/orbit164 section subgroup to rank nine."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {
            str(path.relative_to(ROOT)): sha256(path) for path in INPUTS
        },
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "Q4O164D1TRANSPORT|degree_one={}|bits={}|status={}|output={}".format(
        len(transported), payload["coefficient_growth"]["maximum_compact_child_coordinate_rational_bits"],
        payload["status"], OUTPUT,
    ),
    flush=True,
)
