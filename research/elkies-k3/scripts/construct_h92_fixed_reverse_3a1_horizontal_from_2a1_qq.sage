#!/usr/bin/env sage
"""Construct the fixed 3A1 reverse horizontal on the exact 2A1 equation.

Enumerate short vectors in twice the exact A1 MW height lattice, restore the
A1 root correction, and retain sections having degree one for the 2A1 pencil.
Seventy-eight such sections through doubled height norm twelve span the full
rank-15 2A1 MW tail lattice.  Build the prescribed orbit498 horizontal and
the two new horizontal root curves by integral Smith words, exact A1 group
law, the second degree-one pointed-quartic map, and exact 2A1 group law.
"""

import hashlib
import json
import time
from pathlib import Path

from sage.all import IntegralLattice, PolynomialRing, QQ, ZZ, ceil, floor, identity_matrix, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
MANIFEST = LOCAL / "h3-r17-backward-exact-lift-manifest.json"
BRIDGE = LOCAL / "q24-equation-d13-to-pinned-r17.json"
A1_MODEL = LOCAL / "fixed-final-a1-reverse-rr-qq.json"
A1_GENERATORS = LOCAL / "fixed-reverse-2a1-horizontal-from-a1-qq.json"
TWO_A1_RR = LOCAL / "fixed-reverse-2a1-rr-qq.json"
TWO_A1_POINTING = LOCAL / "fixed-reverse-2a1-pointing-qq.json"
OUTPUT = LOCAL / "fixed-reverse-3a1-horizontal-from-2a1-qq.json"


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
a1_model = read_json(A1_MODEL)
a1_generators_artifact = read_json(A1_GENERATORS)
two_a1_rr = read_json(TWO_A1_RR)
two_a1_pointing = read_json(TWO_A1_POINTING)
assert a1_generators_artifact["status"] == "PASS_EXACT_QQ_FIXED_REVERSE_2A1_HORIZONTAL_ON_A1"
assert two_a1_rr["status"] == "PASS_EXACT_QQ_FIXED_REVERSE_2A1_RR_JACOBIAN"
assert two_a1_pointing["status"] == "PASS_EXACT_QQ_FIXED_REVERSE_2A1_POINTING"

# -------------------------------------------------------------------------
# Lattice enumeration and sparse integral construction words.
# -------------------------------------------------------------------------
G1 = matrix(ZZ, bridge["final_a1_frame"])
Rroot = G1[:1, :1]
coupling = G1[:1, 1:]
H1 = G1[1:, 1:] - coupling.transpose() * Rroot.inverse() * coupling
H1_twice = (2 * H1).change_ring(ZZ)
assert H1_twice.det() == 31064064

step_2a1_to_a1 = manifest["forward_steps"][-2]
step_3a1_to_2a1 = manifest["forward_steps"][-3]
assert step_2a1_to_a1["orbit"] == 981 and step_3a1_to_2a1["orbit"] == 498
T981 = matrix(ZZ, step_2a1_to_a1["transition"])
T498 = matrix(ZZ, step_3a1_to_2a1["transition"])
T498_inverse = T498.inverse().change_ring(ZZ)
E19 = identity_matrix(ZZ, 19)

candidates = []
for shell in IntegralLattice(H1_twice).short_vectors(13):
    for tail_value in shell:
        z = vector(ZZ, tail_value)
        raw = vector(ZZ, [0] + list(z))
        dual = QQ((raw * G1)[0]) / QQ(G1[0, 0])
        for root_coordinate in sorted({ZZ(floor(-dual)), ZZ(ceil(-dual))}):
            # The set above is {-floor(dual),-ceil(dual)}; retain both closest
            # integers because the A1 discriminant coset can have two lifts.
            pframe = vector(ZZ, [root_coordinate] + list(z))
            norm = pframe * G1 * pframe
            if norm < 4 or (norm - 4) % 2:
                continue
            p_dot_o = ZZ((norm - 4) // 2)
            section_a1 = vector(ZZ, [p_dot_o + 1, 1] + list(pframe))
            section_2a1 = vector(ZZ, section_a1 * T981)
            if section_2a1[1] == 1:
                candidates.append({
                    "A1_tail": z,
                    "A1_class": section_a1,
                    "2A1_class": section_2a1,
                })
assert len(candidates) == 78
candidate_tails = matrix(ZZ, [list(row["2A1_class"][4:]) for row in candidates])
assert candidate_tails.rank() == 15
assert candidate_tails.elementary_divisors()[:15] == [1] * 15


def smith_solution(target):
    linear = candidate_tails.transpose()
    smith, left, right = linear.smith_form()
    transformed = left * vector(ZZ, target)
    coordinates = vector(ZZ, [0] * len(candidates))
    for index in range(15):
        assert smith[index, index] == 1
        coordinates[index] = transformed[index]
    answer = vector(ZZ, right * coordinates)
    assert answer * candidate_tails == vector(ZZ, target)
    return answer


fibre_3a1_in_2a1 = vector(ZZ, E19.row(0) * T498_inverse)
assert list(fibre_3a1_in_2a1[:2]) == [4, 2]
target_section_2a1 = vector(ZZ, [7, 1] + list(fibre_3a1_in_2a1[2:]))
horizontal_word = smith_solution(target_section_2a1[4:])
assert [(i, int(c)) for i, c in enumerate(horizontal_word) if c] == [
    (50, -1), (56, 1), (68, -1), (70, -2), (71, -1),
    (72, 1), (73, -1), (74, 2), (75, 3), (77, 1),
]

effective_new_root_classes = [
    -vector(ZZ, E19.row(root_index) * T498_inverse)
    for root_index in (2, 4)
]
root_words = [smith_solution(root_class[4:]) for root_class in effective_new_root_classes]
assert [[(i, int(c)) for i, c in enumerate(word) if c] for word in root_words] == [
    [(61, -1), (67, 1), (68, -1), (71, -1), (72, 1), (74, 1), (75, 1), (76, -1)],
    [(46, -1), (56, -1), (61, -1), (67, 3), (68, -3), (69, 2),
     (70, 1), (71, -1), (72, 2), (73, -1), (74, 2), (75, 2), (76, -4)],
]
# These already-needed short sections contain a unimodular MW basis.  Keeping
# the basis in this form avoids reconstructing long unit-vector Smith words on
# the 2A1 equation (the latter produce much larger intermediate functions).
basis_candidate_indices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 15, 22]
basis_candidate_tails = matrix(ZZ, [
    list(candidates[index]["2A1_class"][4:])
    for index in basis_candidate_indices
])
assert basis_candidate_tails.det() == 1
# A generic 2A1 basis need not map onto the saturated 3A1 MW lattice.  This
# subset of the same short candidate pool has determinant one after T498 and
# therefore gives the next source equation a saturated rank-14 basis without
# any division or saturation reconstruction.
three_a1_basis_candidate_indices = [71, 1, 70, 72, 75, 46, 74, 4, 15, 61, 73, 56, 9, 77]
three_a1_basis_tails = matrix(ZZ, [
    list((candidates[index]["2A1_class"] * T498)[5:])
    for index in three_a1_basis_candidate_indices
])
assert three_a1_basis_tails.det() == 1
needed_candidate_indices = sorted(set(
    index
    for word in [horizontal_word] + root_words
    for index, coefficient in enumerate(word)
    if coefficient
) | set(basis_candidate_indices) | set(three_a1_basis_candidate_indices))
assert needed_candidate_indices == [
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 15, 22,
    46, 50, 56, 61, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77,
]

# -------------------------------------------------------------------------
# Exact A1 generator pool and group law.
# -------------------------------------------------------------------------
RS = PolynomialRing(QQ, "s")
s = RS.gen()
KS = RS.fraction_field()
A1 = RS(a1_model["child"]["minimal_A_coefficients_low_to_high"])
B1 = RS(a1_model["child"]["minimal_B_coefficients_low_to_high"])


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
    assert y3 ** 2 == x3 ** 3 + Acurve * x3 + Bcurve
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


generator_rows = a1_generators_artifact["construction_sections"]
assert len(generator_rows) == 22
generator_tails = matrix(ZZ, [row["class_in_A1_coordinates"][3:] for row in generator_rows])
assert generator_tails.rank() == 16 and generator_tails.elementary_divisors()[:16] == [1] * 16
generator_points = [
    point_from_record(row["A1_section"], RS, KS, A1, B1)
    for row in generator_rows
]


def generator_word_for_a1_tail(target):
    linear = generator_tails.transpose()
    smith, left, right = linear.smith_form()
    transformed = left * vector(ZZ, target)
    coordinates = vector(ZZ, [0] * len(generator_rows))
    for index in range(16):
        assert smith[index, index] == 1
        coordinates[index] = transformed[index]
    answer = vector(ZZ, right * coordinates)
    assert answer * generator_tails == vector(ZZ, target)
    return answer


candidate_a1_words = {}
candidate_a1_points = {}
for index in needed_candidate_indices:
    word = generator_word_for_a1_tail(candidates[index]["A1_tail"])
    candidate_a1_words[index] = word
    point = compose(generator_points, word, KS(A1), KS(B1))
    candidate_a1_points[index] = point

# -------------------------------------------------------------------------
# Second pointed-quartic transport A1 -> 2A1.
# -------------------------------------------------------------------------
horizontal_2a1 = a1_generators_artifact["section"]
X = RS(horizontal_2a1["x_numerator_coefficients_low_to_high"])
Y = RS(horizontal_2a1["y_numerator_coefficients_low_to_high"])
Z = RS(horizontal_2a1["Z_coefficients_low_to_high"])
Hx = KS(X) / KS(Z ** 2)
Hy = KS(Y) / KS(Z ** 3)
basis_pairs = two_a1_rr["smooth_RR"]["basis_pairs"]
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
    (RT(values) * ss ** degree for degree, values in enumerate(two_a1_rr["binary_quartic"]["coefficients_in_old_u_low_to_high"])),
    U.zero(),
)
square_factor = sum(
    (RT(values) * ss ** degree for degree, values in enumerate(two_a1_rr["binary_quartic"]["square_factor_coefficients_in_old_u_low_to_high"])),
    U.zero(),
)


def rational_from_record(record, ring, field):
    return field(ring(record["numerator_coefficients_low_to_high"])) / field(ring(record["denominator_coefficients_low_to_high"]))


alpha = KT(two_a1_pointing["fixed_zero"]["old_I2_support"])
q = rational_from_record(two_a1_pointing["fixed_zero"]["quartic_ordinate"], RT, KT)


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
A2 = RT(two_a1_rr["child"]["minimal_A_coefficients_low_to_high"])
B2 = RT(two_a1_rr["child"]["minimal_B_coefficients_low_to_high"])


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
    assert y ** 2 == x ** 3 + KT(A2) * x + KT(B2)
    return x, y


candidate_2a1_points = {
    index: pointed_image(candidate_a1_points[index])
    for index in needed_candidate_indices
}


def compose_candidate_word(word):
    indices = [index for index, coefficient in enumerate(word) if coefficient]
    return compose(
        [candidate_2a1_points[index] for index in indices],
        vector(ZZ, [word[index] for index in indices]),
        KT(A2), KT(B2),
    )


target_point = compose_candidate_word(horizontal_word)
effective_root_points = [compose_candidate_word(word) for word in root_words]


def normalized_record(point, expected_p_dot_o=None):
    x, y = point
    xn, xd = RT(x.numerator()), RT(x.denominator())
    yn, yd = RT(y.numerator()), RT(y.denominator())
    assert xd.is_square() and yd == xd.sqrt() ** 3
    z = xd.sqrt()
    if expected_p_dot_o is not None:
        assert z.degree() == expected_p_dot_o
    return {
        "x_numerator_coefficients_low_to_high": coeffs(xn),
        "x_denominator_coefficients_low_to_high": coeffs(xd),
        "y_numerator_coefficients_low_to_high": coeffs(yn),
        "y_denominator_coefficients_low_to_high": coeffs(yd),
        "Z_coefficients_low_to_high": coeffs(z),
        "degrees_X_Y_Z": [int(xn.degree()), int(yn.degree()), int(z.degree())],
        "P_dot_O": int(z.degree()),
        "literal_weierstrass_identity": True,
        "maximum_rational_bits": maximum_bits((xn, xd, yn, yd)),
    }


target_record = normalized_record(target_point, 6)
assert target_record["degrees_X_Y_Z"] == [16, 24, 6]
root_records = [normalized_record(point, 0) for point in effective_root_points]
two_a1_basis_points = [candidate_2a1_points[index] for index in basis_candidate_indices]
two_a1_basis_records = [normalized_record(point) for point in two_a1_basis_points]
three_a1_basis_source_records = [
    normalized_record(candidate_2a1_points[index])
    for index in three_a1_basis_candidate_indices
]

payload = {
    "schema": "elkies-k3.fixed-reverse-3a1-horizontal-from-2a1-qq.v1",
    "status": "PASS_EXACT_QQ_FIXED_REVERSE_3A1_HORIZONTAL_ON_2A1",
    "fixed_edge": {
        "forward": "3A1/MW14 --q4 orbit498--> 2A1/MW15",
        "reverse_fibre_in_2A1_coordinates": list(map(int, fibre_3a1_in_2a1)),
        "identity": "F_3A1=O+P-2F",
        "target_section_in_2A1_coordinates": list(map(int, target_section_2a1)),
        "P_dot_O": 6,
        "short_A1_degree_one_count": len(candidates),
        "2A1_tail_rank": int(candidate_tails.rank()),
        "selected_sparse_combination": [[int(i), int(horizontal_word[i])] for i in range(len(horizontal_word)) if horizontal_word[i]],
    },
    "candidate_construction": [
        {
            "candidate_index": int(index),
            "class_in_A1_coordinates": list(map(int, candidates[index]["A1_class"])),
            "class_in_2A1_coordinates": list(map(int, candidates[index]["2A1_class"])),
            "A1_generator_word": list(map(int, candidate_a1_words[index])),
            "exact_A1_group_law": True,
            "exact_degree_one_second_pointed_transport": True,
        }
        for index in needed_candidate_indices
    ],
    "section": target_record,
    "effective_3A1_horizontal_roots_on_2A1_source": [
        {
            "class_in_2A1_coordinates": list(map(int, root_class)),
            "candidate_word": [[int(i), int(word[i])] for i in range(len(word)) if word[i]],
            "section": record,
        }
        for root_class, word, record in zip(effective_new_root_classes, root_words, root_records)
    ],
    "exact_2A1_MW_basis": [
        {
            "basis_index": int(basis_index),
            "candidate_index": int(candidate_index),
            "MW_tail": list(map(int, candidates[candidate_index]["2A1_class"][4:])),
            "A1_generator_word": list(map(int, candidate_a1_words[candidate_index])),
            "section": record,
        }
        for basis_index, candidate_index, record in zip(
            range(15), basis_candidate_indices, two_a1_basis_records
        )
    ],
    "exact_2A1_sections_for_3A1_MW_basis": [
        {
            "basis_index": int(basis_index),
            "candidate_index": int(candidate_index),
            "class_in_2A1_coordinates": list(map(int, candidates[candidate_index]["2A1_class"])),
            "class_in_3A1_coordinates": list(map(int, candidates[candidate_index]["2A1_class"] * T498)),
            "3A1_MW_tail": list(map(int, (candidates[candidate_index]["2A1_class"] * T498)[5:])),
            "A1_generator_word": list(map(int, candidate_a1_words[candidate_index])),
            "2A1_section": record,
        }
        for basis_index, candidate_index, record in zip(
            range(14), three_a1_basis_candidate_indices, three_a1_basis_source_records
        )
    ],
    "remaining_vertical_3A1_root": {
        "class_in_2A1_coordinates": list(map(int, vector(ZZ, E19.row(3) * T498_inverse))),
        "role": "old second 2A1 component; identified by the remaining child I2 support after the two horizontal roots are contracted",
    },
    "method": {
        "A1_height_lattice_short_vector_bound": 12,
        "integral_smith_words": True,
        "direct_unimodular_2A1_basis_determinant": int(basis_candidate_tails.det()),
        "direct_unimodular_3A1_basis_determinant": int(three_a1_basis_tails.det()),
        "two_exact_group_laws": True,
        "second_degree_one_pointed_transport": True,
        "groebner_or_elimination": False,
        "runtime_seconds": time.monotonic() - started,
    },
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in (MANIFEST, BRIDGE, A1_MODEL, A1_GENERATORS, TWO_A1_RR, TWO_A1_POINTING)],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in (MANIFEST, BRIDGE, A1_MODEL, A1_GENERATORS, TWO_A1_RR, TWO_A1_POINTING)},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "FIXEDREVERSE3A1HORIZONTAL|short_A1_degree1=78|tail_rank=15|candidates=30|"
    "degrees=(16,24,6)|PdotO=6|bits={}|seconds={:.3f}|status={}|output={}".format(
        target_record["maximum_rational_bits"], payload["method"]["runtime_seconds"],
        payload["status"], OUTPUT,
    ), flush=True,
)
