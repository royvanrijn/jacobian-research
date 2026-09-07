#!/usr/bin/env sage -python
"""Point the physical q4/orbit1584 quartic at second_I6_affine over QQ.

The selected zero is already an exact section of the q4/orbit208 parent.
Evaluating the q4/orbit1584 RR quotient on that section gives a Mobius
function of the old base, so one exact inversion supplies a point on the
stored quartic.  Translating the binary quartic at that point gives the
pointed generalized Weierstrass model and its exact standard short model.

The two signs above the constant old_A11_component_0 support are attached
and distinguished by their exact specialization at the three finite I2
fibres: C0 hits one nonidentity node, while C4 hits two.  No Groebner basis
or global discriminant factorization is used.
"""

import hashlib
import json
import time
from pathlib import Path

from sage.all import EllipticCurve, PolynomialRing, QQ, factorial


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
Q208 = LOCAL / "q24-2a5-physical-q4o208-equation-marking-qq.json"
Q1584 = LOCAL / "q4o208-physical-q4o1584-rr-qq.json"
OUTPUT = LOCAL / "q4o1584-second-affine-equation-marking-qq.json"
INPUTS = (Q208, Q1584)

started = time.monotonic()


def log(stage, **fields):
    tail = "|".join(f"{key}={value}" for key, value in fields.items())
    print(
        f"Q4O1584MARKQQ|stage={stage}|elapsed={time.monotonic()-started:.3f}"
        + (f"|{tail}" if tail else ""),
        flush=True,
    )


q208 = json.loads(Q208.read_text())
q1584 = json.loads(Q1584.read_text())
assert q208["status"] == "PASS_EXACT_QQ_PHYSICAL_Q4O208_C5_EQUATION_MARKING"
assert q1584["status"] == "PASS_EXACT_QQ_PHYSICAL_Q4O1584_D4_A3_3A1_RR_AND_JACOBIAN"

RT = PolynomialRing(QQ, "T")
T = RT.gen()
KT = RT.fraction_field()
RV = PolynomialRing(QQ, "V")
V = RV.gen()
KV = RV.fraction_field()
RL = PolynomialRing(KV, "L")
L = RL.gen()


def read_rational(record, ring, field):
    numerator = ring([QQ(value) for value in record["numerator_coefficients_low_to_high"]])
    denominator = ring([QQ(value) for value in record["denominator_coefficients_low_to_high"]])
    return field(numerator) / field(denominator)


def record(value, ring, field):
    value = field(value)
    numerator = ring(value.numerator())
    denominator = ring(value.denominator())
    common = numerator.gcd(denominator)
    numerator //= common
    denominator //= common
    if denominator.leading_coefficient() < 0:
        numerator, denominator = -numerator, -denominator
    return {
        "numerator_coefficients_low_to_high": [str(entry) for entry in numerator.list()],
        "denominator_coefficients_low_to_high": [str(entry) for entry in denominator.list()],
        "degrees_numerator_denominator": [int(numerator.degree()), int(denominator.degree())],
    }


def rational_bits(values):
    entries = []
    for value in values:
        value = KV(value)
        entries.extend(QQ(entry) for entry in value.numerator().list())
        entries.extend(QQ(entry) for entry in value.denominator().list())
    return max(
        max(abs(entry.numerator()).nbits(), entry.denominator().nbits())
        for entry in entries
    )


first = q208["first_I6_affine_component_on_C5_pointed_child"]
second = q208["second_I6_affine_component_on_C5_pointed_child"]
PX = read_rational(first["x"], RT, KT)
PY = read_rational(first["y"], RT, KT)
QX = read_rational(second["x"], RT, KT)
QY = read_rational(second["y"], RT, KT)
s = QQ(q1584["resolved_RR"]["support"])
c0 = QQ(q1584["resolved_RR"]["double_branch_c0"])
c1 = QQ(q1584["resolved_RR"]["first_jet_c1"])

chord = KT((QY + PY) / (QX - PX))
new_base = KT((chord - c0 - c1 * (T - s)) / (T - s) ** 2)
base_num = RT(new_base.numerator())
base_den = RT(new_base.denominator())
if (base_num.degree(), base_den.degree()) != (1, 1):
    raise ArithmeticError("second affine does not give a Mobius q4/o1584 base map")
n0, n1 = QQ(base_num[0]), QQ(base_num[1])
d0, d1 = QQ(base_den[0]), QQ(base_den[1])
T_of_V = KV((KV(n0) - KV(V) * KV(d0)) / (KV(V) * KV(d1) - KV(n1)))
if KV(base_num(T_of_V) / base_den(T_of_V)) != KV(V):
    raise ArithmeticError("second-affine Mobius inversion failed")


def substitute_t(value):
    value = KT(value)
    return KV(RT(value.numerator())(T_of_V) / RT(value.denominator())(T_of_V))


L0 = KV(T_of_V - s)
W0 = substitute_t((2 * QX + PX - chord**2) / (T - s) ** 2)
quartic_coefficients = [
    RV([QQ(entry) for entry in coefficients])
    for coefficients in q1584["quartic"]["coefficients_in_L_low_to_high"]
]
quartic = RL([KV(value) for value in quartic_coefficients])
if W0**2 != quartic(L0):
    raise ArithmeticError("second affine misses the exact q4/o1584 quartic")
log(
    "ZERO", selected="second_I6_affine", base_degrees="1/1",
    L0_degrees=f"{L0.numerator().degree()}/{L0.denominator().degree()}",
)

# Translate L=L0+z and point the quartic at (0,W0).
translated = [
    KV(quartic.derivative(order)(L0) / factorial(order))
    for order in range(5)
]
e, d, c, b, a = translated
if e != W0**2:
    raise ArithmeticError("translated quartic constant is not the selected square")
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
pointed_A = -c4 / 48
pointed_B = -c6 / 864
A_short = KV(RV([QQ(value) for value in q1584["child"]["minimal_A_coefficients_low_to_high"]]))
B_short = KV(RV([QQ(value) for value in q1584["child"]["minimal_B_coefficients_low_to_high"]]))
if 81 * pointed_A != A_short or 729 * pointed_B != B_short:
    raise ArithmeticError("second-affine pointed model misses the stored short Jacobian")
E = EllipticCurve(KV, [0, 0, 0, A_short, B_short])
log("POINTED_JACOBIAN", A=True, B=True, short_scaling="9/27")


def square_root_rational(value):
    value = KV(value)
    numerator = RV(value.numerator())
    denominator = RV(value.denominator())
    nlc = QQ(numerator.leading_coefficient())
    dlc = QQ(denominator.leading_coefficient())
    if not nlc.is_square() or not dlc.is_square():
        raise ArithmeticError("rational square has nonsquare leading coefficient")
    nroot = RV(nlc.sqrt())
    droot = RV(dlc.sqrt())
    for factor, exponent in numerator.monic().factor():
        if exponent % 2:
            raise ArithmeticError("quartic numerator is not a square")
        nroot *= factor ** (exponent // 2)
    for factor, exponent in denominator.monic().factor():
        if exponent % 2:
            raise ArithmeticError("quartic denominator is not a square")
        droot *= factor ** (exponent // 2)
    answer = KV(nroot / droot)
    if answer**2 != value:
        raise ArithmeticError("rational square-root reconstruction failed")
    return answer


def pointed_coordinates(old_l, ordinate):
    z = KV(old_l - L0)
    if not z:
        if ordinate == W0:
            raise ArithmeticError("selected quartic origin maps to infinity")
        if ordinate != -W0:
            raise ArithmeticError("unexpected ordinate above selected old-base coordinate")
        # Exact limit at the opposite point above z=0.
        xg = -a2
        yg = a1 * a2 - a3
    else:
        xg = (2 * W0 * (ordinate + W0) + d * z) / z**2
        yg = (
            4 * W0**2 * (ordinate + W0) + 2 * W0 * d * z
            + (2 * W0 * c - d**2 / (2 * W0)) * z**2
        ) / z**3
    x_short = KV(9 * (xg + b2 / 12))
    y_short = KV(27 * (yg + (a1 * xg + a3) / 2))
    return E(x_short, y_short)


c0_support = QQ(q208["physical_fibres"]["first_old_I6_I4"]["support"])
c0_old_l = KV(c0_support - s)
c0_ordinate = square_root_rational(quartic(c0_old_l))
c0_pair = [
    pointed_coordinates(c0_old_l, c0_ordinate),
    pointed_coordinates(c0_old_l, -c0_ordinate),
]
if c0_pair[0] == c0_pair[1]:
    raise ArithmeticError("constant-support q4/o1584 pair collapsed")


def map_parent_section(section_record):
    sx = read_rational(section_record["x"], RT, KT)
    sy = read_rational(section_record["y"], RT, KT)
    section_chord = KT((sy + PY) / (sx - PX))
    section_base = KT(
        (section_chord - c0 - c1 * (T - s)) / (T - s) ** 2
    )
    numerator = RT(section_base.numerator())
    denominator = RT(section_base.denominator())
    if (numerator.degree(), denominator.degree()) != (1, 1):
        raise ArithmeticError("parent section is not degree one for q4/o1584")
    sn0, sn1 = QQ(numerator[0]), QQ(numerator[1])
    sd0, sd1 = QQ(denominator[0]), QQ(denominator[1])
    old_t = KV(
        (KV(sn0) - KV(V) * KV(sd0))
        / (KV(V) * KV(sd1) - KV(sn1))
    )

    def at_old_t(value):
        value = KT(value)
        return KV(
            RT(value.numerator())(old_t) / RT(value.denominator())(old_t)
        )

    old_l = KV(old_t - s)
    ordinate = at_old_t(
        (2 * sx + PX - section_chord**2) / (T - s) ** 2
    )
    if ordinate**2 != quartic(old_l):
        raise ArithmeticError("parent section misses q4/o1584 quartic")
    return pointed_coordinates(old_l, ordinate)


c7_point = map_parent_section(q208["old_A11_component_7_on_C5_pointed_child"])
if c0_pair[0] + c0_pair[1] != c7_point:
    raise ArithmeticError("constant-support pair does not sum to C7")

i2_supports = [
    QQ(item["support"])
    for item in q1584["child"]["finite_reducible_fibres"]
    if item["kodaira"] == "I2"
]
if len(i2_supports) != 3:
    raise ArithmeticError("stored q4/o1584 model does not have three finite I2 fibres")


def hits_i2_node(point, support):
    x_coordinate, y_coordinate = KV(point[0]), KV(point[1])
    if x_coordinate.denominator()(support) == 0 or y_coordinate.denominator()(support) == 0:
        return False
    node = QQ(-3 * B_short(support) / (2 * A_short(support)))
    if node**3 + A_short(support) * node + B_short(support):
        raise ArithmeticError("stored I2 support has the wrong cubic node")
    return x_coordinate(support) == node and y_coordinate(support) == 0


i2_profiles = [
    [bool(hits_i2_node(point, support)) for support in i2_supports]
    for point in c0_pair
]
i2_hit_counts = [sum(profile) for profile in i2_profiles]
if sorted(i2_hit_counts) != [1, 2]:
    raise ArithmeticError(
        f"constant-support pair has unexpected I2 profiles: {i2_profiles}"
    )
c0_index = i2_hit_counts.index(1)
c0_point = c0_pair[c0_index]
c4_point = c0_pair[1 - c0_index]
log(
    "C0_PAIR", support_mod103=26,
    selected_index=c0_index,
    I2_hit_counts=",".join(map(str, i2_hit_counts)),
    x_degrees=";".join(
        f"{point[0].numerator().degree()}/{point[0].denominator().degree()}"
        for point in c0_pair
    ),
)

payload = {
    "schema": "elkies-k3.q4o1584-second-affine-equation-marking-qq.v1",
    "status": "PASS_EXACT_QQ_Q4O1584_SECOND_AFFINE_POINTING_AND_C0_MARKING",
    "selected_zero": {
        "label": "second_I6_affine_component",
        "q4o1584_base_map": record(
            KV(base_num(T_of_V) / base_den(T_of_V)), RV, KV
        ),
        "old_base_coordinate": record(L0, RV, KV),
        "quartic_ordinate": record(W0, RV, KV),
        "exact_quartic_identity": True,
        "maps_to": "point at infinity on the pointed generalized Weierstrass model",
    },
    "pointed_jacobian": {
        "standard_short_scaling": "x_short=9*x_pointed_short, y_short=27*y_pointed_short",
        "exact_A_identity": True,
        "exact_B_identity": True,
    },
    "old_A11_component_7_on_second_affine_pointed_child": {
        "x": record(c7_point[0], RV, KV),
        "y": record(c7_point[1], RV, KV),
        "exact_child_identity": True,
    },
    "old_A11_component_0_on_second_affine_pointed_child": {
        "selected_candidate_index": c0_index,
        "x": record(c0_point[0], RV, KV),
        "y": record(c0_point[1], RV, KV),
        "finite_I2_nonidentity_profile": i2_profiles[c0_index],
        "finite_I2_nonidentity_count": 1,
        "exact_child_identity": True,
    },
    "old_A11_component_4_on_second_affine_pointed_child": {
        "x": record(c4_point[0], RV, KV),
        "y": record(c4_point[1], RV, KV),
        "finite_I2_nonidentity_profile": i2_profiles[1 - c0_index],
        "finite_I2_nonidentity_count": 2,
        "exact_child_identity": True,
    },
    "constant_support_candidate_pair": [
        {
            "x": record(point[0], RV, KV),
            "y": record(point[1], RV, KV),
            "exact_child_identity": True,
        }
        for point in c0_pair
    ],
    "coefficient_growth": {
        "selected_zero_maximum_rational_bits": rational_bits([L0, W0]),
        "C0_pair_maximum_rational_bits": rational_bits(
            [coordinate for point in c0_pair for coordinate in (point[0], point[1])]
        ),
    },
    "method": {
        "large_Groebner_required": False,
        "global_discriminant_factorization_required": False,
        "runtime_seconds": time.monotonic() - started,
    },
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in INPUTS
        },
    },
    "proof_boundary": (
        "The selected second-I6 affine zero is exact over QQ(V), its pointed quartic "
        "has exactly the stored q4/o1584 short Jacobian, and the two constant-support "
        "sections above the old_A11_component_0 support are exact. The labeled C7 "
        "section is their exact group-law sum. Exact specialization at the three I2 "
        "fibres distinguishes old_A11_component_0 (one nonidentity hit) from "
        "old_A11_component_4 (two nonidentity hits). This completes the two "
        "equation inputs needed to compile q4/o164."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
log(
    "DONE", zero="second_I6_affine", C0_pair=2,
    status=payload["status"], output=OUTPUT,
)
