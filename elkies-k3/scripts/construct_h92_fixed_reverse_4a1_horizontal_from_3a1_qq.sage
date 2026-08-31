#!/usr/bin/env sage
"""Construct the fixed 4A1 reverse horizontal on the exact 3A1 equation.

Enumerate short vectors in twice the exact 2A1 MW height lattice and retain
the 77 sections which have degree one for the 3A1 pencil.  A q114-adapted
unimodular subset gives short exact words for the marked section and the two
horizontal roots.  Construct those points from the exact 2A1 MW basis, apply
the third pointed-quartic map, and finish with exact 3A1 group law.  The
vertical correction 3*C1+2*C2 lowers P.O from 34 to 21.
"""

import hashlib
import itertools
import json
import time
from pathlib import Path

from sage.all import IntegralLattice, PolynomialRing, QQ, ZZ, block_diagonal_matrix, ceil, floor, identity_matrix, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
MANIFEST = LOCAL / "h3-r17-backward-exact-lift-manifest.json"
BRIDGE = LOCAL / "q24-equation-d13-to-pinned-r17.json"
TWO_A1_RR = LOCAL / "fixed-reverse-2a1-rr-qq.json"
TWO_A1_GENERATORS = LOCAL / "fixed-reverse-3a1-horizontal-from-2a1-qq.json"
THREE_A1_RR = LOCAL / "fixed-reverse-3a1-rr-qq.json"
THREE_A1_POINTING = LOCAL / "fixed-reverse-3a1-pointing-qq.json"
OUTPUT = LOCAL / "fixed-reverse-4a1-horizontal-from-3a1-qq.json"


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
            answer = max(answer, abs(ZZ(value.numerator())).nbits(), ZZ(value.denominator()).nbits())
    return int(answer)


started = time.monotonic()
manifest = read_json(MANIFEST)
bridge = read_json(BRIDGE)
two_a1_rr = read_json(TWO_A1_RR)
two_a1_generators = read_json(TWO_A1_GENERATORS)
three_a1_rr = read_json(THREE_A1_RR)
three_a1_pointing = read_json(THREE_A1_POINTING)
assert two_a1_rr["status"] == "PASS_EXACT_QQ_FIXED_REVERSE_2A1_RR_JACOBIAN"
assert two_a1_generators["status"] == "PASS_EXACT_QQ_FIXED_REVERSE_3A1_HORIZONTAL_ON_2A1"
assert three_a1_rr["status"] == "PASS_EXACT_QQ_FIXED_REVERSE_3A1_RR_JACOBIAN"
assert three_a1_pointing["status"] == "PASS_EXACT_QQ_FIXED_REVERSE_3A1_POINTING"

# -------------------------------------------------------------------------
# Lattice enumeration and sparse integral construction words.
# -------------------------------------------------------------------------
G1 = matrix(ZZ, bridge["final_a1_frame"])
step_2a1_to_a1 = manifest["forward_steps"][-2]
step_3a1_to_2a1 = manifest["forward_steps"][-3]
step_4a1_to_3a1 = manifest["forward_steps"][-4]
assert [step["orbit"] for step in (step_2a1_to_a1, step_3a1_to_2a1, step_4a1_to_3a1)] == [981, 498, 114]
T981 = matrix(ZZ, step_2a1_to_a1["transition"])
T498 = matrix(ZZ, step_3a1_to_2a1["transition"])
T114 = matrix(ZZ, step_4a1_to_3a1["transition"])
T981_inverse = T981.inverse().change_ring(ZZ)
T498_inverse = T498.inverse().change_ring(ZZ)
T114_inverse = T114.inverse().change_ring(ZZ)
E19 = identity_matrix(ZZ, 19)

Ugram = matrix(ZZ, [[0, 1], [1, 0]])
Q1 = block_diagonal_matrix(Ugram, -G1)
Q2 = T981_inverse * Q1 * T981_inverse.transpose()
assert Q2[:2, :2] == Ugram and not Q2[:2, 2:]
G2 = (-Q2[2:, 2:]).change_ring(ZZ)
Rroot = G2[:2, :2]
assert Rroot == 2 * identity_matrix(ZZ, 2)
coupling = G2[:2, 2:]
H2 = G2[2:, 2:] - coupling.transpose() * Rroot.inverse() * coupling
H2_twice = (2 * H2).change_ring(ZZ)
Q3 = T498_inverse * Q2 * T498_inverse.transpose()
assert Q3[:2, :2] == Ugram and not Q3[:2, 2:]
G3 = (-Q3[2:, 2:]).change_ring(ZZ)
R3 = G3[:3, :3]
C3 = G3[:3, 3:]
H3 = G3[3:, 3:] - C3.transpose() * R3.inverse() * C3

candidates = []
for shell in IntegralLattice(H2_twice).short_vectors(13):
    for tail_value in shell:
        z = vector(ZZ, tail_value)
        raw = vector(ZZ, [0, 0] + list(z))
        dual = vector(QQ, (raw * G2)[:2]) * Rroot.inverse()
        root_choices = [
            sorted({ZZ(floor(-value)), ZZ(ceil(-value))})
            for value in dual
        ]
        for root_coordinates in itertools.product(*root_choices):
            pframe = vector(ZZ, list(root_coordinates) + list(z))
            norm = pframe * G2 * pframe
            if norm < 4 or (norm - 4) % 2:
                continue
            p_dot_o = ZZ((norm - 4) // 2)
            section_2a1 = vector(ZZ, [p_dot_o + 1, 1] + list(pframe))
            section_3a1 = vector(ZZ, section_2a1 * T498)
            if section_3a1[1] == 1:
                candidates.append({
                    "2A1_tail": z,
                    "2A1_class": section_2a1,
                    "3A1_class": section_3a1,
                })
assert len(candidates) == 77
candidate_tails = matrix(ZZ, [list(row["3A1_class"][5:]) for row in candidates])
assert candidate_tails.rank() == 14
assert candidate_tails.elementary_divisors()[:14] == [1] * 14

basis_candidate_indices = [0, 2, 4, 5, 6, 8, 10, 11, 12, 14, 15, 16, 26, 32]
basis_candidate_tails = matrix(ZZ, [
    list(candidates[index]["3A1_class"][5:])
    for index in basis_candidate_indices
])
assert abs(basis_candidate_tails.det()) == 1


def basis_solution(target):
    answer = vector(ZZ, target) * basis_candidate_tails.inverse()
    assert answer * basis_candidate_tails == vector(ZZ, target)
    return vector(ZZ, answer)


fibre_4a1_in_3a1 = vector(ZZ, E19.row(0) * T114_inverse)
assert list(fibre_4a1_in_3a1) == [18, 2, 14, -3, -1, -1, 2, 2, -2, -3, 3, -3, -2, 3, -3, -3, 3, -3, 0]
old_zero = vector(ZZ, [-1, 1] + [0] * 17)
old_fibre = vector(ZZ, E19.row(0))
vertical_correction = 3 * vector(ZZ, E19.row(2)) + 2 * vector(ZZ, E19.row(3))
target_section_3a1 = fibre_4a1_in_3a1 - old_zero - vertical_correction + 3 * old_fibre
assert list(target_section_3a1) == [22, 1, 11, -5, -1, -1, 2, 2, -2, -3, 3, -3, -2, 3, -3, -3, 3, -3, 0]
assert fibre_4a1_in_3a1 == old_zero + target_section_3a1 + vertical_correction - 3 * old_fibre

horizontal_word = basis_solution(target_section_3a1[5:])
effective_new_root_classes = [
    -vector(ZZ, E19.row(root_index) * T114_inverse)
    for root_index in (3, 4)
]
root_words = [basis_solution(root_class[5:]) for root_class in effective_new_root_classes]
assert [[(basis_candidate_indices[i], int(c)) for i, c in enumerate(word) if c] for word in [horizontal_word] + root_words] == [
    [(0, 2), (2, 1), (5, 1), (6, -1), (8, 2), (12, 1), (14, -2), (15, 3), (16, -4), (26, -3), (32, -1)],
    [(5, 1), (6, -1), (8, 1), (11, 1), (14, -1), (15, 2), (16, -2), (26, -2)],
    [(0, 1), (8, 1), (14, -1), (15, 1), (16, -1), (26, -1)],
]
needed_candidate_indices = basis_candidate_indices

# -------------------------------------------------------------------------
# Exact 2A1 generator pool and group law.
# -------------------------------------------------------------------------
RS = PolynomialRing(QQ, "s")
s = RS.gen()
KS = RS.fraction_field()
A2 = RS(two_a1_rr["child"]["minimal_A_coefficients_low_to_high"])
B2 = RS(two_a1_rr["child"]["minimal_B_coefficients_low_to_high"])


def point_from_record(record, ring, field, Acurve, Bcurve):
    x = field(ring(record["x_numerator_coefficients_low_to_high"])) / field(ring(record["x_denominator_coefficients_low_to_high"]))
    y = field(ring(record["y_numerator_coefficients_low_to_high"])) / field(ring(record["y_denominator_coefficients_low_to_high"]))
    assert y ** 2 == x ** 3 + field(Acurve) * x + field(Bcurve)
    return x, y


def neg(P):
    return None if P is None else (P[0], -P[1])


def add(P, Q, Acurve, Bcurve):
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2:
        if y1 == -y2:
            return None
        slope = (3 * x1 ** 2 + Acurve) / (2 * y1)
    else:
        slope = (y2 - y1) / (x2 - x1)
    x3 = slope ** 2 - x1 - x2
    y3 = slope * (x1 - x3) - y1
    return x3, y3


def multiply(P, coefficient, Acurve, Bcurve):
    coefficient = ZZ(coefficient)
    if coefficient < 0:
        return neg(multiply(P, -coefficient, Acurve, Bcurve))
    answer = None
    addend = P
    while coefficient:
        if coefficient & 1:
            answer = add(answer, addend, Acurve, Bcurve)
        coefficient >>= 1
        if coefficient:
            addend = add(addend, addend, Acurve, Bcurve)
    return answer


def compose(points, word, Acurve, Bcurve):
    summands = [multiply(P, c, Acurve, Bcurve) for P, c in zip(points, word) if c]
    while len(summands) > 1:
        summands = [
            add(summands[i], summands[i + 1], Acurve, Bcurve)
            if i + 1 < len(summands) else summands[i]
            for i in range(0, len(summands), 2)
        ]
    return summands[0]


generator_rows = two_a1_generators["exact_2A1_MW_basis"]
assert len(generator_rows) == 15
generator_tails = matrix(ZZ, [row["MW_tail"] for row in generator_rows])
assert abs(generator_tails.det()) == 1
generator_points = [
    point_from_record(row["section"], RS, KS, A2, B2)
    for row in generator_rows
]


def generator_word_for_2a1_tail(target):
    answer = vector(ZZ, target) * generator_tails.inverse()
    answer = vector(ZZ, answer)
    assert answer * generator_tails == vector(ZZ, target)
    return answer


candidate_2a1_words = {}
candidate_2a1_points = {}
for index in needed_candidate_indices:
    word = generator_word_for_2a1_tail(candidates[index]["2A1_tail"])
    assert sum(abs(value) for value in word) <= 15
    candidate_2a1_words[index] = word
    candidate_2a1_points[index] = compose(generator_points, word, KS(A2), KS(B2))

# -------------------------------------------------------------------------
# Third pointed-quartic transport 2A1 -> 3A1.
# -------------------------------------------------------------------------
horizontal_3a1 = two_a1_generators["section"]
X = RS(horizontal_3a1["x_numerator_coefficients_low_to_high"])
Y = RS(horizontal_3a1["y_numerator_coefficients_low_to_high"])
Z = RS(horizontal_3a1["Z_coefficients_low_to_high"])
Hx = KS(X) / KS(Z ** 2)
Hy = KS(Y) / KS(Z ** 3)
basis_pairs = three_a1_rr["smooth_RR"]["basis_pairs"]
AA0 = RS(basis_pairs[0]["AA_coefficients_low_to_high"])
BB0 = RS(basis_pairs[0]["BB_coefficients_low_to_high"])
AA1 = RS(basis_pairs[1]["AA_coefficients_low_to_high"])
BB1 = RS(basis_pairs[1]["BB_coefficients_low_to_high"])

RT = PolynomialRing(QQ, "t")
t = RT.gen()
KT = RT.fraction_field()
U = PolynomialRing(RT, "s")
ss = U.gen()
quartic = sum(
    (RT(values) * ss ** degree for degree, values in enumerate(three_a1_rr["binary_quartic"]["coefficients_in_old_u_low_to_high"])),
    U.zero(),
)
square_factor = sum(
    (RT(values) * ss ** degree for degree, values in enumerate(three_a1_rr["binary_quartic"]["square_factor_coefficients_in_old_u_low_to_high"])),
    U.zero(),
)


def rational_from_record(record, ring, field):
    return field(ring(record["numerator_coefficients_low_to_high"])) / field(ring(record["denominator_coefficients_low_to_high"]))


alpha = KT(three_a1_pointing["fixed_zero"]["old_I2_support"])
q = rational_from_record(three_a1_pointing["fixed_zero"]["quartic_ordinate"], RT, KT)


def restrict_to_quartic(P):
    px, py = P
    slope = (py + Hy) / (px - Hx)
    base = -(KS(AA0) + KS(BB0 * Z) * slope) / (KS(AA1) + KS(BB1 * Z) * slope)
    numerator, denominator = RS(base.numerator()), RS(base.denominator())
    assert max(numerator.degree(), denominator.degree()) == 1
    old_base = KT(numerator[0] - t * denominator[0]) / KT(t * denominator[1] - numerator[1])

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
    ordinate_s = bb_value ** 2 * (2 * px + Hx - slope ** 2) / eval_bivariate_at_ks(square_factor)
    assert ordinate_s ** 2 == eval_bivariate_at_ks(quartic)

    def eval_ks_at_kt(value):
        value = KS(value)

        def evaluate(poly):
            answer = KT.zero()
            for coefficient in reversed(RS(poly).list()):
                answer = answer * old_base + KT(coefficient)
            return answer

        return evaluate(value.numerator()) / evaluate(value.denominator())

    ordinate = eval_ks_at_kt(ordinate_s)
    return old_base, ordinate


def eval_u_at_alpha(poly):
    return KT(U(poly)(U(alpha)))


e = eval_u_at_alpha(quartic)
d = eval_u_at_alpha(quartic.derivative())
c = eval_u_at_alpha(quartic.derivative(2)) / 2
b = eval_u_at_alpha(quartic.derivative(3)) / 6
a = eval_u_at_alpha(quartic.derivative(4)) / 24
assert e == q ** 2
a1 = d / q
a2 = c - d ** 2 / (4 * q ** 2)
a3 = 2 * q * b
a4 = -4 * q ** 2 * a
b2 = a1 ** 2 + 4 * a2
A3 = RT(three_a1_rr["child"]["minimal_A_coefficients_low_to_high"])
B3 = RT(three_a1_rr["child"]["minimal_B_coefficients_low_to_high"])


def pointed_image(P):
    old_base, ordinate = restrict_to_quartic(P)
    relative = old_base - alpha
    assert relative
    x_general = (2 * q * (ordinate + q) + d * relative) / relative ** 2
    y_general = (
        4 * q ** 2 * (ordinate + q) + 2 * q * d * relative
        + (2 * q * c - d ** 2 / (2 * q)) * relative ** 2
    ) / relative ** 3
    x = KT(9 * (x_general + b2 / 12))
    y = KT(27 * (y_general + (a1 * x_general + a3) / 2))
    assert y ** 2 == x ** 3 + KT(A3) * x + KT(B3)
    return x, y


candidate_3a1_points = {
    index: pointed_image(candidate_2a1_points[index])
    for index in needed_candidate_indices
}


def compose_candidate_word(word):
    summands = []
    for basis_position, coefficient in enumerate(word):
        if not coefficient:
            continue
        candidate_index = basis_candidate_indices[basis_position]
        point = multiply(candidate_3a1_points[candidate_index], coefficient, KT(A3), KT(B3))
        tail = coefficient * vector(ZZ, candidates[candidate_index]["3A1_class"][5:])
        summands.append((point, tail, f"{int(coefficient)}*{candidate_index}"))
    schedule = []
    maximum_intermediate_height = max(QQ(tail * H3 * tail) for unused, tail, unused_label in summands)
    while len(summands) > 1:
        unused_height, left, right = min(
            (QQ((summands[i][1] + summands[j][1]) * H3 * (summands[i][1] + summands[j][1])), i, j)
            for i in range(len(summands)) for j in range(i + 1, len(summands))
        )
        point_left, tail_left, label_left = summands[left]
        point_right, tail_right, label_right = summands[right]
        new_point = add(point_left, point_right, KT(A3), KT(B3))
        new_tail = tail_left + tail_right
        new_height = QQ(new_tail * H3 * new_tail)
        maximum_intermediate_height = max(maximum_intermediate_height, new_height)
        schedule.append({
            "left": label_left, "right": label_right,
            "result_height": str(new_height),
        })
        new_label = f"({label_left}+{label_right})"
        summands = [
            item for index, item in enumerate(summands)
            if index not in (left, right)
        ] + [(new_point, new_tail, new_label)]
    point = summands[0][0]
    assert point[1] ** 2 == point[0] ** 3 + KT(A3) * point[0] + KT(B3)
    return point, schedule, maximum_intermediate_height


target_point, target_schedule, target_maximum_height = compose_candidate_word(horizontal_word)
root_constructions = [compose_candidate_word(word) for word in root_words]
effective_root_points = [record[0] for record in root_constructions]


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
        "finite_Z_degree": int(z.degree()),
        "P_dot_O": int(z.degree() if expected_p_dot_o is None else expected_p_dot_o),
        "literal_weierstrass_identity": True,
        "maximum_rational_bits": maximum_bits((xn, xd, yn, yd)),
    }


target_record = normalized_record(target_point, 21)
assert target_record["degrees_X_Y_Z"] == [46, 69, 21]
root_records = [
    normalized_record(point, expected)
    for point, expected in zip(effective_root_points, (9, 4))
]
candidate_2a1_records = [
    normalized_record(candidate_2a1_points[index], int(candidates[index]["2A1_class"][0] - 1))
    for index in basis_candidate_indices
]
three_a1_basis_points = [candidate_3a1_points[index] for index in basis_candidate_indices]
three_a1_basis_records = [
    normalized_record(point, int(candidates[index]["3A1_class"][0] - 1))
    for index, point in zip(basis_candidate_indices, three_a1_basis_points)
]

payload = {
    "schema": "elkies-k3.fixed-reverse-4a1-horizontal-from-3a1-qq.v1",
    "status": "PASS_EXACT_QQ_FIXED_REVERSE_4A1_HORIZONTAL_ON_3A1",
    "fixed_edge": {
        "forward": "4A1/MW13 --q4 orbit114--> 3A1/MW14",
        "reverse_fibre_in_3A1_coordinates": list(map(int, fibre_4a1_in_3a1)),
        "identity": "F_4A1=O+P+3*C1+2*C2-3*F_3A1",
        "target_section_in_3A1_coordinates": list(map(int, target_section_3a1)),
        "vertical_correction_in_3A1_coordinates": list(map(int, vertical_correction)),
        "fibre_twist": -3,
        "P_dot_O": 21,
        "short_2A1_degree_one_count": len(candidates),
        "3A1_tail_rank": int(candidate_tails.rank()),
        "selected_sparse_combination": [
            [int(basis_candidate_indices[i]), int(horizontal_word[i])]
            for i in range(len(horizontal_word)) if horizontal_word[i]
        ],
    },
    "candidate_construction": [
        {
            "candidate_index": int(index),
            "class_in_2A1_coordinates": list(map(int, candidates[index]["2A1_class"])),
            "class_in_3A1_coordinates": list(map(int, candidates[index]["3A1_class"])),
            "2A1_generator_word": list(map(int, candidate_2a1_words[index])),
            "2A1_section": source_record,
            "3A1_section": child_record,
            "exact_2A1_group_law": True,
            "exact_degree_one_third_pointed_transport": True,
        }
        for index, source_record, child_record in zip(
            basis_candidate_indices, candidate_2a1_records, three_a1_basis_records
        )
    ],
    "section": target_record,
    "effective_4A1_horizontal_roots_on_3A1_source": [
        {
            "class_in_3A1_coordinates": list(map(int, root_class)),
            "candidate_word": [
                [int(basis_candidate_indices[i]), int(word[i])]
                for i in range(len(word)) if word[i]
            ],
            "section": record,
        }
        | {
            "addition_schedule": construction[1],
            "maximum_intermediate_MW_height": str(construction[2]),
        }
        for root_class, word, record, construction in zip(
            effective_new_root_classes, root_words, root_records, root_constructions
        )
    ],
    "exact_3A1_MW_basis": [
        {
            "basis_index": int(basis_index),
            "candidate_index": int(candidate_index),
            "MW_tail": list(map(int, candidates[candidate_index]["3A1_class"][5:])),
            "2A1_generator_word": list(map(int, candidate_2a1_words[candidate_index])),
            "section": record,
        }
        for basis_index, candidate_index, record in zip(
            range(14), basis_candidate_indices, three_a1_basis_records
        )
    ],
    "prescribed_4A1_zero": {
        "class_in_3A1_coordinates": list(map(int, E19.row(4))),
        "old_I2_support_source": "fixed-reverse-3a1-pointing-qq.json effective_horizontal_components[1]",
    },
    "remaining_vertical_4A1_roots": {
        "classes_in_3A1_coordinates": [
            list(map(int, vector(ZZ, E19.row(index) * T114_inverse)))
            for index in (2, 5)
        ],
        "count": 2,
    },
    "method": {
        "2A1_height_lattice_short_vector_bound": 12,
        "q114_adapted_unimodular_3A1_basis_determinant": int(basis_candidate_tails.det()),
        "target_addition_schedule": target_schedule,
        "target_maximum_intermediate_MW_height": str(target_maximum_height),
        "target_final_MW_height": "91/2",
        "two_exact_group_laws": True,
        "third_degree_one_pointed_transport": True,
        "groebner_or_elimination": False,
        "runtime_seconds": time.monotonic() - started,
    },
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in (MANIFEST, BRIDGE, TWO_A1_RR, TWO_A1_GENERATORS, THREE_A1_RR, THREE_A1_POINTING)],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in (MANIFEST, BRIDGE, TWO_A1_RR, TWO_A1_GENERATORS, THREE_A1_RR, THREE_A1_POINTING)},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "FIXEDREVERSE4A1HORIZONTAL|short_2A1_degree1=77|tail_rank=14|candidates=14|"
    "degrees=(46,69,21)|PdotO=21|bits={}|seconds={:.3f}|status={}|output={}".format(
        target_record["maximum_rational_bits"], payload["method"]["runtime_seconds"],
        payload["status"], OUTPUT,
    ), flush=True,
)
