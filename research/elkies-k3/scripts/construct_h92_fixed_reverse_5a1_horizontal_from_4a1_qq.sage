#!/usr/bin/env sage
"""HISTORICAL_DIAGNOSTIC: rejected direct q52 section transport over QQ.

The stored 3A1 sections are multisections, not degree-one curves, for the
physical reflected q114 pencil.  Consequently the asserted direct pointed
transport below is invalid and deliberately fails before producing an
artifact.  The active replacement is
``construct_h92_fixed_reverse_5a1_abel_word_mod131.sage``, which fibrewise
Abel-reduces the multisections before taking the cancellation-heavy words.
This file is retained only to document the rejected shortcut.
"""

import hashlib
import json
import time
from pathlib import Path

from sage.all import PolynomialRing, QQ, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
SOURCE_SURFACE = LOCAL / "fixed-reverse-3a1-rr-qq.json"
SOURCE_CURVES = LOCAL / "fixed-reverse-4a1-horizontal-from-3a1-qq.json"
CURRENT_RR = LOCAL / "fixed-reverse-4a1-rr-qq.json"
CURRENT_POINTING = LOCAL / "fixed-reverse-4a1-pointing-qq.json"
Q52_AUDIT = LOCAL / "fixed-reverse-5a1-physical-nef-audit.json"
MANIFEST = LOCAL / "h3-r17-backward-exact-lift-manifest.json"
OUTPUT = LOCAL / "fixed-reverse-5a1-horizontal-from-4a1-qq.json"


def read_json(path):
    return json.loads(path.read_text())


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coeffs(poly):
    return [str(value) for value in poly.list()]


def maximum_bits(polys):
    answer = 0
    for poly in polys:
        for value in poly:
            value = QQ(value)
            answer = max(
                answer,
                abs(ZZ(value.numerator())).nbits(),
                ZZ(value.denominator()).nbits(),
            )
    return int(answer)


started = time.monotonic()
source_surface = read_json(SOURCE_SURFACE)
source_curves = read_json(SOURCE_CURVES)
current_rr = read_json(CURRENT_RR)
current_pointing = read_json(CURRENT_POINTING)
q52_audit = read_json(Q52_AUDIT)
manifest = read_json(MANIFEST)
assert source_surface["status"] == "PASS_EXACT_QQ_FIXED_REVERSE_3A1_RR_JACOBIAN"
assert source_curves["status"] == "PASS_EXACT_QQ_FIXED_REVERSE_4A1_HORIZONTAL_ON_3A1"
assert current_rr["status"] == "PASS_EXACT_QQ_FIXED_REVERSE_4A1_RR_JACOBIAN"
assert current_pointing["status"] == "PASS_EXACT_QQ_FIXED_REVERSE_4A1_POINTING"
assert q52_audit["status"] == "PASS_EXACT_QQ_FIXED_REVERSE_5A1_PHYSICAL_NEF"

# -------------------------------------------------------------------------
# Exact lattice words in the physical 4A1 Mordell--Weil tail.
# -------------------------------------------------------------------------
T114 = matrix(ZZ, manifest["forward_steps"][-4]["transition"])
source_basis_records = source_curves["exact_3A1_MW_basis"]
source_class_records = source_curves["candidate_construction"]
assert len(source_basis_records) == len(source_class_records) == 14
source_classes = [
    vector(ZZ, record["class_in_3A1_coordinates"])
    for record in source_class_records
]
physical_4a1_classes = [vector(ZZ, value * T114) for value in source_classes]
assert all(value[1] == 1 for value in physical_4a1_classes)
tail_matrix = matrix(ZZ, [list(value[6:]) for value in physical_4a1_classes])
assert tail_matrix.rank() == 13
selected_indices = next(
    tuple(index for index in range(14) if index != omitted)
    for omitted in range(14)
    if abs(matrix(ZZ, [tail_matrix.row(index) for index in range(14) if index != omitted]).det()) == 1
)
selected_tails = matrix(ZZ, [tail_matrix.row(index) for index in selected_indices])
assert abs(selected_tails.det()) == 1

physical_fibre = vector(
    ZZ, q52_audit["physical_component_reduction"]["reduced_fibre_in_4A1_coordinates"]
)
target_section = vector(ZZ, [
    16, 1, -3, -7, 4, 0, 2, -3, -2, 0, 0, 0, 1, -1, 3, 0, 1, 0, 0,
])
old_zero = vector(ZZ, [-1, 1] + [0] * 17)
vertical_residual = physical_fibre - old_zero - target_section
assert list(vertical_residual) == [-6, 0, 0, 1] + [0] * 15
assert physical_fibre == old_zero + target_section + vertical_residual

root_classes = [
    vector(ZZ, values)
    for values in q52_audit["full_marked_transport"][
        "physical_5A1_root_classes_in_4A1_coordinates"
    ]
]
horizontal_root_classes = [
    root if root[1] == 1 else -root
    for root in root_classes if abs(root[1]) == 1
]
vertical_root_classes = [root for root in root_classes if root[1] == 0]
assert len(horizontal_root_classes) == 3 and len(vertical_root_classes) == 2


def word_for_tail(target):
    word = vector(ZZ, target) * selected_tails.inverse()
    assert word * selected_tails == vector(ZZ, target)
    return vector(ZZ, word)


target_word = word_for_tail(target_section[6:])
root_words = [word_for_tail(root[6:]) for root in horizontal_root_classes]

# -------------------------------------------------------------------------
# Exact physical q114 pointed-quartic transport 3A1 -> 4A1.
# -------------------------------------------------------------------------
RS = PolynomialRing(QQ, "s")
s = RS.gen()
KS = RS.fraction_field()
A3 = RS(source_surface["child"]["minimal_A_coefficients_low_to_high"])
B3 = RS(source_surface["child"]["minimal_B_coefficients_low_to_high"])
horizontal = source_curves["section"]
X = RS(horizontal["x_numerator_coefficients_low_to_high"])
Y = RS(horizontal["y_numerator_coefficients_low_to_high"])
Z = RS(horizontal["Z_coefficients_low_to_high"])
Hx = KS(X) / KS(Z ** 2)
Hy = KS(Y) / KS(Z ** 3)
basis_pairs = current_rr["smooth_RR"]["basis_pairs"]
AA0 = RS(basis_pairs[0]["AA_coefficients_low_to_high"])
BB0 = RS(basis_pairs[0]["BB_coefficients_low_to_high"])
AA1 = RS(basis_pairs[1]["AA_coefficients_low_to_high"])
BB1 = RS(basis_pairs[1]["BB_coefficients_low_to_high"])

RT = PolynomialRing(QQ, "t")
t = RT.gen()
KT = RT.fraction_field()
U = PolynomialRing(RT, "s")
ss = U.gen()
quartic = sum((
    RT(values) * ss ** degree
    for degree, values in enumerate(
        current_rr["binary_quartic"]["coefficients_in_old_u_low_to_high"]
    )
), U.zero())
square_factor = sum((
    RT(values) * ss ** degree
    for degree, values in enumerate(
        current_rr["binary_quartic"]["square_factor_coefficients_in_old_u_low_to_high"]
    )
), U.zero())


def rational_from_record(record):
    return KT(RT(record["numerator_coefficients_low_to_high"])) / KT(
        RT(record["denominator_coefficients_low_to_high"])
    )


alpha = KT(current_pointing["fixed_zero"]["old_I2_support"])
q_origin = rational_from_record(current_pointing["fixed_zero"]["quartic_ordinate"])


def point_from_record(record):
    x = KS(RS(record["x_numerator_coefficients_low_to_high"])) / KS(
        RS(record["x_denominator_coefficients_low_to_high"])
    )
    y = KS(RS(record["y_numerator_coefficients_low_to_high"])) / KS(
        RS(record["y_denominator_coefficients_low_to_high"])
    )
    assert y ** 2 == x ** 3 + KS(A3) * x + KS(B3)
    return x, y


def restrict_to_quartic(P):
    px, py = P
    slope = (py + Hy) / (px - Hx)
    base = -(KS(AA0) + KS(BB0 * Z) * slope) / (
        KS(AA1) + KS(BB1 * Z) * slope
    )
    numerator, denominator = RS(base.numerator()), RS(base.denominator())
    assert max(numerator.degree(), denominator.degree()) == 1
    old_base = KT(numerator[0] - t * denominator[0]) / KT(
        t * denominator[1] - numerator[1]
    )

    def eval_rt_at_ks(poly):
        answer = KS.zero()
        for coefficient in reversed(RT(poly).list()):
            answer = answer * base + KS(coefficient)
        return answer

    def eval_bivariate_at_ks(poly):
        answer = KS.zero()
        for coefficient in reversed(U(poly).list()):
            answer = answer * KS(s) + eval_rt_at_ks(coefficient)
        return answer

    bb_value = KS(BB0) + base * KS(BB1)
    ordinate_s = (
        bb_value ** 2 * (2 * px + Hx - slope ** 2)
        / eval_bivariate_at_ks(square_factor)
    )
    assert ordinate_s ** 2 == eval_bivariate_at_ks(quartic)

    def eval_ks_at_kt(value):
        value = KS(value)

        def evaluate(poly):
            answer = KT.zero()
            for coefficient in reversed(RS(poly).list()):
                answer = answer * old_base + KT(coefficient)
            return answer

        return evaluate(value.numerator()) / evaluate(value.denominator())

    return old_base, eval_ks_at_kt(ordinate_s)


def eval_u_at_alpha(poly):
    return KT(U(poly)(U(alpha)))


e = eval_u_at_alpha(quartic)
d = eval_u_at_alpha(quartic.derivative())
c = eval_u_at_alpha(quartic.derivative(2)) / 2
b = eval_u_at_alpha(quartic.derivative(3)) / 6
a = eval_u_at_alpha(quartic.derivative(4)) / 24
assert e == q_origin ** 2
a1 = d / q_origin
a2 = c - d ** 2 / (4 * q_origin ** 2)
a3 = 2 * q_origin * b
b2 = a1 ** 2 + 4 * a2
A4 = RT(current_rr["child"]["minimal_A_coefficients_low_to_high"])
B4 = RT(current_rr["child"]["minimal_B_coefficients_low_to_high"])


def pointed_image(P):
    old_base, ordinate = restrict_to_quartic(P)
    relative = old_base - alpha
    assert relative
    x_general = (
        2 * q_origin * (ordinate + q_origin) + d * relative
    ) / relative ** 2
    y_general = (
        4 * q_origin ** 2 * (ordinate + q_origin)
        + 2 * q_origin * d * relative
        + (2 * q_origin * c - d ** 2 / (2 * q_origin)) * relative ** 2
    ) / relative ** 3
    x = KT(9 * (x_general + b2 / 12))
    y = KT(27 * (y_general + (a1 * x_general + a3) / 2))
    assert y ** 2 == x ** 3 + KT(A4) * x + KT(B4)
    return x, y


source_points = [point_from_record(record["section"]) for record in source_basis_records]
transported_points = [pointed_image(source_points[index]) for index in selected_indices]


def neg(P):
    return None if P is None else (P[0], -P[1])


def add(P, Q):
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2:
        if y1 == -y2:
            return None
        slope = (3 * x1 ** 2 + KT(A4)) / (2 * y1)
    else:
        slope = (y2 - y1) / (x2 - x1)
    x3 = slope ** 2 - x1 - x2
    return x3, slope * (x1 - x3) - y1


def multiply(P, coefficient):
    coefficient = ZZ(coefficient)
    if coefficient < 0:
        return neg(multiply(P, -coefficient))
    answer = None
    addend = P
    while coefficient:
        if coefficient & 1:
            answer = add(answer, addend)
        coefficient >>= 1
        if coefficient:
            addend = add(addend, addend)
    return answer


def compose(word):
    summands = [
        multiply(P, coefficient)
        for P, coefficient in zip(transported_points, word) if coefficient
    ]
    while len(summands) > 1:
        summands = [
            add(summands[index], summands[index + 1])
            if index + 1 < len(summands) else summands[index]
            for index in range(0, len(summands), 2)
        ]
    point = summands[0]
    assert point[1] ** 2 == point[0] ** 3 + KT(A4) * point[0] + KT(B4)
    return point


target_point = compose(target_word)
root_points = [compose(word) for word in root_words]


def normalized_record(point, expected_p_dot_o=None):
    x, y = point
    xn, xd = RT(x.numerator()), RT(x.denominator())
    yn, yd = RT(y.numerator()), RT(y.denominator())
    assert xd.is_square() and yd == xd.sqrt() ** 3
    z = xd.sqrt()
    if expected_p_dot_o is not None:
        assert z.degree() <= expected_p_dot_o
    return {
        "x_numerator_coefficients_low_to_high": coeffs(xn),
        "x_denominator_coefficients_low_to_high": coeffs(xd),
        "y_numerator_coefficients_low_to_high": coeffs(yn),
        "y_denominator_coefficients_low_to_high": coeffs(yd),
        "Z_coefficients_low_to_high": coeffs(z),
        "degrees_X_Y_Z": [int(xn.degree()), int(yn.degree()), int(z.degree())],
        "P_dot_O": int(expected_p_dot_o if expected_p_dot_o is not None else z.degree()),
        "literal_weierstrass_identity": True,
        "maximum_rational_bits": maximum_bits((xn, xd, yn, yd)),
    }


target_record = normalized_record(target_point, 15)
root_records = [normalized_record(point) for point in root_points]
transported_records = [normalized_record(point) for point in transported_points]
payload = {
    "schema": "elkies-k3.fixed-reverse-5a1-horizontal-from-4a1-qq.v1",
    "status": "PASS_EXACT_QQ_FIXED_REVERSE_5A1_HORIZONTAL_ON_4A1",
    "fixed_edge": {
        "forward": "5A1/MW12 --q4 orbit52--> 4A1/MW13",
        "physical_reverse_fibre_in_4A1_coordinates": list(map(int, physical_fibre)),
        "identity": "F_5A1=O+P-C2-6F_4A1",
        "target_section_in_4A1_coordinates": list(map(int, target_section)),
        "vertical_residual_in_4A1_coordinates": list(map(int, vertical_residual)),
        "fibre_twist": -6,
        "P_dot_O": 15,
        "selected_basis_indices": list(map(int, selected_indices)),
        "target_word": list(map(int, target_word)),
    },
    "section": target_record,
    "effective_5A1_horizontal_roots_on_4A1_source": [
        {
            "class_in_4A1_coordinates": list(map(int, root_class)),
            "basis_word": list(map(int, word)),
            "section": record,
        }
        for root_class, word, record in zip(
            horizontal_root_classes, root_words, root_records
        )
    ],
    "remaining_vertical_5A1_roots": {
        "count": len(vertical_root_classes),
        "classes_in_4A1_coordinates": [
            list(map(int, root)) for root in vertical_root_classes
        ],
    },
    "exact_4A1_MW_basis": [
        {
            "source_basis_index": int(index),
            "MW_tail": list(map(int, physical_4a1_classes[index][6:])),
            "section": record,
        }
        for index, record in zip(selected_indices, transported_records)
    ],
    "prescribed_5A1_zero": {
        "class_in_4A1_coordinates": q52_audit["full_marked_transport"][
            "physical_zero_in_4A1_coordinates"
        ],
        "identification": "fourth old 4A1 nonidentity I2 component",
    },
    "method": {
        "degree_one_physical_q114_pointed_transport": True,
        "unimodular_4A1_MW_tail_basis": True,
        "exact_integral_words": True,
        "exact_group_law": True,
        "groebner_or_surface_elimination": False,
        "runtime_seconds": time.monotonic() - started,
    },
    "inputs": {
        str(path.relative_to(ROOT)): sha256(path)
        for path in (
            SOURCE_SURFACE, SOURCE_CURVES, CURRENT_RR, CURRENT_POINTING,
            Q52_AUDIT, MANIFEST,
        )
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "FIXEDREVERSE5A1HORIZONTAL|basis=13|PO=15|degrees={}|roots=3|vertical=2|"
    "bits={}|seconds={:.3f}|status={}|output={}".format(
        target_record["degrees_X_Y_Z"], target_record["maximum_rational_bits"],
        payload["method"]["runtime_seconds"], payload["status"], OUTPUT,
    ),
    flush=True,
)
