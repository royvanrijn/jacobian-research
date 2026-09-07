#!/usr/bin/env sage
"""Construct the fixed 2A1 reverse horizontal on the exact A1 equation.

The target A1 section is a degree-15 curve in the q12 rootless presentation,
so transporting it directly is unattractive.  Instead enumerate the pinned
R17 norm-four sections meeting the reverse A1 fibre once.  Twenty-one occur;
their A1 Mordell--Weil tails span rank 15 and the desired tail has a sparse
integral expression in seven of them.  Construct those seven curves by exact
rootless group law, apply the degree-one pointed-quartic map, and add them on
the exact A1 Jacobian.

No Groebner basis, elimination, or nonlinear section solve is used.
"""

import hashlib
import json
import time
from pathlib import Path

from sage.all import IntegralLattice, PolynomialRing, QQ, ZZ, identity_matrix, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
PINNED = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
BRIDGE = LOCAL / "q24-equation-d13-to-pinned-r17.json"
MANIFEST = LOCAL / "h3-r17-backward-exact-lift-manifest.json"
ROOTLESS_MODEL = LOCAL / "q12o5867-smooth-rr-qq.json"
ROOTLESS_BASIS = LOCAL / "q12o5867-rootless-selected-basis-qq.json"
ROOTLESS_HEIGHTS = LOCAL / "q12o5867-rootless-height-basis-qq.json"
A1_CURVES = LOCAL / "fixed-final-a1-horizontal-from-q12-endpoint-qq.json"
A1_RR = LOCAL / "fixed-final-a1-reverse-rr-qq.json"
OUTPUT = LOCAL / "fixed-reverse-2a1-horizontal-from-a1-qq.json"


def read_json(path):
    return json.loads(path.read_text())


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_matrix(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


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
bridge = read_json(BRIDGE)
manifest = read_json(MANIFEST)
rootless_model = read_json(ROOTLESS_MODEL)
rootless_basis = read_json(ROOTLESS_BASIS)
rootless_heights = read_json(ROOTLESS_HEIGHTS)
a1_curves = read_json(A1_CURVES)
a1_rr = read_json(A1_RR)
assert a1_rr["status"] == "PASS_EXACT_QQ_FIXED_FINAL_A1_REVERSE_RR_JACOBIAN"

Gp = load_matrix(PINNED)
A1_to_pinned = matrix(ZZ, bridge["final_a1_to_pinned_r17_transition"])
C = matrix(ZZ, rootless_heights["basis_to_pinned_rank17"])
assert abs(A1_to_pinned.det()) == 1 and C.det() == -1

# -------------------------------------------------------------------------
# Exact lattice plan for the reverse 2A1 fibre and its horizontal section.
# -------------------------------------------------------------------------
step = manifest["forward_steps"][-2]
assert step["parent"] == "2A1/MW15" and step["child"] == "A1/MW16" and step["orbit"] == 981
T = matrix(ZZ, step["transition"])
E19 = identity_matrix(ZZ, 19)
fibre_2a1_in_a1 = vector(ZZ, E19.row(0) * T.inverse().change_ring(ZZ))
assert list(fibre_2a1_in_a1[:2]) == [4, 2]
target_section_a1 = vector(ZZ, [7, 1] + list(fibre_2a1_in_a1[2:]))

fibre_a1_pinned = vector(ZZ, a1_curves["fixed_edge"]["parent_fibre_in_pinned_rootless_coordinates"])
wd = vector(ZZ, fibre_a1_pinned[2:])
norm4 = [vector(ZZ, value) for value in IntegralLattice(Gp).short_vectors(5)[4]]
degree_one = [value for value in norm4 if 6 - value * Gp * wd == 1]
assert len(norm4) == 2622 and len(degree_one) == 21

a1_tails = []
for value in degree_one:
    section_pinned = vector(ZZ, [1, 1] + list(value))
    section_a1 = vector(ZZ, section_pinned * A1_to_pinned)
    assert section_a1[1] == 1
    a1_tails.append(list(section_a1[3:]))
tail_matrix = matrix(ZZ, a1_tails)
assert tail_matrix.rank() == 15
target_tail = vector(ZZ, target_section_a1[3:])
assert target_tail in tail_matrix.row_module()

# Smith solve for c*tail_matrix=target_tail, with free Smith coordinates zero.
linear = tail_matrix.transpose()
smith, left, right = linear.smith_form()
transformed_target = left * target_tail
z = vector(ZZ, [0] * len(degree_one))
for index in range(linear.rank()):
    assert smith[index, index] == 1
    z[index] = transformed_target[index]
combination = vector(ZZ, right * z)
assert combination * tail_matrix == target_tail
assert [(i, int(c)) for i, c in enumerate(combination) if c] == [
    (1, 1), (7, -1), (9, 1), (12, -1), (13, 1), (15, 1), (18, -2),
]

# One norm-six rootless section completes the rank-15 norm-four pool to the
# full A1 MW rank.  Its q12-basis word has L1 norm three.  This lets us also
# construct the two effective roots of the new 2A1 frame for the marking gate.
norm6_degree_one = [
    value for value in [vector(ZZ, v) for v in IntegralLattice(Gp).short_vectors(7)[6]]
    if 8 - value * Gp * wd == 1
]
assert len(norm6_degree_one) == 68
extender = norm6_degree_one[3]
assert list(extender * C.inverse().change_ring(ZZ)) == [0, 0, 0, 0, 0, -1, 0, 0, -1, 0, 0, 0, 0, -1, 0, 0, 0]
full_vectors = degree_one + [extender]
full_tails = a1_tails + [
    list((vector(ZZ, [2, 1] + list(extender)) * A1_to_pinned)[3:])
]
full_tail_matrix = matrix(ZZ, full_tails)
assert full_tail_matrix.rank() == 16


def integral_tail_solution(target):
    linear_full = full_tail_matrix.transpose()
    smith_full, left_full, right_full = linear_full.smith_form()
    transformed = left_full * vector(ZZ, target)
    coordinates = vector(ZZ, [0] * len(full_vectors))
    for smith_index in range(16):
        assert smith_full[smith_index, smith_index] == 1
        coordinates[smith_index] = transformed[smith_index]
    answer = vector(ZZ, right_full * coordinates)
    assert answer * full_tail_matrix == vector(ZZ, target)
    return answer


effective_root_classes_a1 = [
    -vector(ZZ, E19.row(root_index) * T.inverse().change_ring(ZZ))
    for root_index in (2, 3)
]
root_combinations = [integral_tail_solution(root[3:]) for root in effective_root_classes_a1]
assert [[(i, int(c)) for i, c in enumerate(word) if c] for word in root_combinations] == [
    [(1, 1), (14, 1), (15, 1), (17, -1), (19, -1), (21, -1)],
    [(7, -1), (9, 1), (13, 1), (14, -1), (17, 1), (18, -1)],
]

# -------------------------------------------------------------------------
# Exact group law on the rootless equation.
# -------------------------------------------------------------------------
R = PolynomialRing(QQ, "u")
u = R.gen()
K = R.fraction_field()
A0 = R(rootless_model["child"]["minimal_A_coefficients_low_to_high"])
B0 = R(rootless_model["child"]["minimal_B_coefficients_low_to_high"])


def source_basis_point(record):
    section = record["section"]
    x = K(R(section["x_coefficients_low_to_high"]))
    y = K(R(section["y_coefficients_low_to_high"]))
    assert y ** 2 == x ** 3 + K(A0) * x + K(B0)
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


rootless_points = [source_basis_point(record) for record in rootless_basis["sections"]]
selected_indices = [index for index, coefficient in enumerate(combination) if coefficient]
# Materialize the complete 22-section generating pool.  Its A1 MW tails span
# Z^16, so later reverse stages can construct arbitrary A1 sections without
# returning to modular section solving.
construction_indices = list(range(len(full_vectors)))
source_words = {}
source_points = {}
for index in construction_indices:
    source_word = vector(ZZ, full_vectors[index] * C.inverse().change_ring(ZZ))
    assert source_word * C == full_vectors[index]
    source_words[index] = source_word
    source_points[index] = compose(rootless_points, source_word, K(A0), K(B0))
    source_p_dot_o = (full_vectors[index] * Gp * full_vectors[index] - 4) // 2
    assert source_points[index][0].denominator().degree() == 2 * source_p_dot_o
    assert source_points[index][1].denominator().degree() == 3 * source_p_dot_o

# -------------------------------------------------------------------------
# Degree-one pointed-quartic transport to the exact A1 Jacobian.
# -------------------------------------------------------------------------
horizontal = a1_curves["section"]
X = R(horizontal["x_numerator_coefficients_low_to_high"])
Y = R(horizontal["y_numerator_coefficients_low_to_high"])
Z = R(horizontal["Z_coefficients_low_to_high"])
Hx = K(X) / K(Z ** 2)
Hy = K(Y) / K(Z ** 3)
basis_pairs = a1_rr["smooth_RR"]["basis_pairs"]
AA0 = R(basis_pairs[0]["AA_coefficients_low_to_high"])
BB0 = R(basis_pairs[0]["BB_coefficients_low_to_high"])
AA1 = R(basis_pairs[1]["AA_coefficients_low_to_high"])
BB1 = R(basis_pairs[1]["BB_coefficients_low_to_high"])

S = PolynomialRing(QQ, "s")
s = S.gen()
KS = S.fraction_field()
U = PolynomialRing(S, "u")
uu = U.gen()
quartic = sum(
    (S(values) * uu ** degree for degree, values in enumerate(a1_rr["binary_quartic"]["coefficients_in_old_u_low_to_high"])),
    U.zero(),
)
square_factor = sum(
    (S(values) * uu ** degree for degree, values in enumerate(a1_rr["binary_quartic"]["square_factor_coefficients_in_old_u_low_to_high"])),
    U.zero(),
)


def restriction(P):
    qx, qy = P
    slope = (qy + Hy) / (qx - Hx)
    f0 = K(AA0) + K(BB0 * Z) * slope
    f1 = K(AA1) + K(BB1 * Z) * slope
    base = -f0 / f1

    def eval_s_at_k(poly):
        answer = K.zero()
        for coefficient in reversed(S(poly).list()):
            answer = answer * base + K(coefficient)
        return answer

    def eval_bivariate(poly):
        answer = K.zero()
        for coefficient in reversed(U(poly).list()):
            answer = answer * K(u) + eval_s_at_k(coefficient)
        return answer

    bb_value = K(BB0) + base * K(BB1)
    ordinate = bb_value ** 2 * (2 * qx + Hx - slope ** 2) / eval_bivariate(square_factor)
    assert ordinate ** 2 == eval_bivariate(quartic)
    return base, ordinate


def invert_degree_one(base):
    numerator = R(base.numerator())
    denominator = R(base.denominator())
    assert max(numerator.degree(), denominator.degree()) == 1
    return KS(numerator[0] - s * denominator[0]) / KS(s * denominator[1] - numerator[1])


def eval_k_at_ks(value, old_base):
    value = K(value)

    def evaluate(poly):
        answer = KS.zero()
        for coefficient in reversed(R(poly).list()):
            answer = answer * old_base + KS(coefficient)
        return answer

    return evaluate(value.numerator()) / evaluate(value.denominator())


def eval_u_poly(poly, old_base):
    answer = KS.zero()
    for coefficient in reversed(U(poly).list()):
        answer = answer * old_base + KS(S(coefficient))
    return answer


zero_record = a1_curves["fixed_A1_zero_on_rootless_source"]
zero_point = (
    K(R(zero_record["x_numerator_coefficients_low_to_high"])) / K(R(zero_record["x_denominator_coefficients_low_to_high"])),
    K(R(zero_record["y_numerator_coefficients_low_to_high"])) / K(R(zero_record["y_denominator_coefficients_low_to_high"])),
)
zero_base_u, zero_ordinate_u = restriction(zero_point)
alpha = invert_degree_one(zero_base_u)
q = eval_k_at_ks(zero_ordinate_u, alpha)
quartic_coefficients = [KS(value) for value in quartic.list()]
e, d0, c0, b0, a0 = quartic_coefficients
d = eval_u_poly(quartic.derivative(), alpha)
c = eval_u_poly(quartic.derivative(2), alpha) / 2
b = eval_u_poly(quartic.derivative(3), alpha) / 6
a = eval_u_poly(quartic.derivative(4), alpha) / 24
assert eval_u_poly(quartic, alpha) == q ** 2
a1 = d / q
a2 = c - d ** 2 / (4 * q ** 2)
a3 = 2 * q * b
a4 = -4 * q ** 2 * a
b2 = a1 ** 2 + 4 * a2
A1 = S(a1_rr["child"]["minimal_A_coefficients_low_to_high"])
B1 = S(a1_rr["child"]["minimal_B_coefficients_low_to_high"])


def pointed_image(P):
    base_u, ordinate_u = restriction(P)
    old_base = invert_degree_one(base_u)
    ordinate = eval_k_at_ks(ordinate_u, old_base)
    assert ordinate ** 2 == eval_u_poly(quartic, old_base)
    relative = old_base - alpha
    assert relative
    x_general = (2 * q * (ordinate + q) + d * relative) / relative ** 2
    y_general = (
        4 * q ** 2 * (ordinate + q)
        + 2 * q * d * relative
        + (2 * q * c - d ** 2 / (2 * q)) * relative ** 2
    ) / relative ** 3
    x = KS(9 * (x_general + b2 / 12))
    y = KS(27 * (y_general + (a1 * x_general + a3) / 2))
    assert y ** 2 == x ** 3 + KS(A1) * x + KS(B1)
    return x, y


a1_points = {index: pointed_image(source_points[index]) for index in construction_indices}
ordered_points = [a1_points[index] for index in selected_indices]
ordered_word = vector(ZZ, [combination[index] for index in selected_indices])
target = compose(ordered_points, ordered_word, KS(A1), KS(B1))
x, y = target
assert y ** 2 == x ** 3 + KS(A1) * x + KS(B1)
xn, xd = S(x.numerator()), S(x.denominator())
yn, yd = S(y.numerator()), S(y.denominator())
assert xd.is_square() and yd == xd.sqrt() ** 3
target_Z = xd.sqrt()
assert target_Z.degree() == 6
assert (xn.degree(), xd.degree(), yn.degree(), yd.degree()) == (16, 12, 24, 18)


def normalized_record(point, expected_p_dot_o):
    px, py = point
    pxn, pxd = S(px.numerator()), S(px.denominator())
    pyn, pyd = S(py.numerator()), S(py.denominator())
    assert pxd.is_square() and pyd == pxd.sqrt() ** 3
    pz = pxd.sqrt()
    assert pz.degree() == expected_p_dot_o
    return {
        "x_numerator_coefficients_low_to_high": coeffs(pxn),
        "x_denominator_coefficients_low_to_high": coeffs(pxd),
        "y_numerator_coefficients_low_to_high": coeffs(pyn),
        "y_denominator_coefficients_low_to_high": coeffs(pyd),
        "Z_coefficients_low_to_high": coeffs(pz),
        "degrees_X_Y_Z": [int(pxn.degree()), int(pyn.degree()), int(pz.degree())],
        "P_dot_O": int(expected_p_dot_o),
        "literal_weierstrass_identity": True,
        "maximum_rational_bits": maximum_bits((pxn, pxd, pyn, pyd)),
    }


effective_root_points = []
for root_word, expected_p_dot_o in zip(root_combinations, (1, 0)):
    indices = [index for index, coefficient in enumerate(root_word) if coefficient]
    point = compose(
        [a1_points[index] for index in indices],
        vector(ZZ, [root_word[index] for index in indices]),
        KS(A1), KS(B1),
    )
    effective_root_points.append(normalized_record(point, expected_p_dot_o))

a1_generator_classes = []
a1_generator_records = []
for index in construction_indices:
    pinned_a = 1 if full_vectors[index] * Gp * full_vectors[index] == 4 else 2
    a1_class = vector(ZZ, [pinned_a, 1] + list(full_vectors[index])) * A1_to_pinned
    assert a1_class[1] == 1
    a1_generator_classes.append(a1_class)
    a1_generator_records.append(normalized_record(a1_points[index], int(a1_class[0] - 1)))

payload = {
    "schema": "elkies-k3.fixed-reverse-2a1-horizontal-from-a1-qq.v1",
    "status": "PASS_EXACT_QQ_FIXED_REVERSE_2A1_HORIZONTAL_ON_A1",
    "fixed_edge": {
        "forward": "2A1/MW15 --q4 orbit981--> A1/MW16",
        "reverse_fibre_in_A1_coordinates": list(map(int, fibre_2a1_in_a1)),
        "identity": "F_2A1=O+P-2F",
        "target_section_in_A1_coordinates": list(map(int, target_section_a1)),
        "P_dot_O": 6,
        "degree_one_norm4_count": len(degree_one),
        "degree_one_tail_rank": int(tail_matrix.rank()),
        "selected_sparse_combination": [[int(index), int(combination[index])] for index in selected_indices],
    },
    "construction_sections": [
        {
            "degree_one_index": int(index),
            "pinned_vector": list(map(int, full_vectors[index])),
            "pinned_norm": int(full_vectors[index] * Gp * full_vectors[index]),
            "rootless_basis_word": list(map(int, source_words[index])),
            "A1_combination_coefficient": int(combination[index]) if index < len(combination) else 0,
            "exact_rootless_section": True,
            "rootless_P_dot_O": int((full_vectors[index] * Gp * full_vectors[index] - 4) // 2),
            "exact_degree_one_pointed_transport": True,
            "class_in_A1_coordinates": list(map(int, a1_generator_classes[index])),
            "A1_section": a1_generator_records[index],
        }
        for index in construction_indices
    ],
    "effective_2A1_roots_on_A1_source": [
        {
            "class_in_A1_coordinates": list(map(int, root_class)),
            "combination_in_degree_one_pool": [[int(index), int(root_word[index])] for index in range(len(root_word)) if root_word[index]],
            "section": section_record,
        }
        for root_class, root_word, section_record in zip(
            effective_root_classes_a1, root_combinations, effective_root_points,
        )
    ],
    "section": {
        "x_numerator_coefficients_low_to_high": coeffs(xn),
        "x_denominator_coefficients_low_to_high": coeffs(xd),
        "y_numerator_coefficients_low_to_high": coeffs(yn),
        "y_denominator_coefficients_low_to_high": coeffs(yd),
        "Z_coefficients_low_to_high": coeffs(target_Z),
        "degrees_X_Y_Z": [int(xn.degree()), int(yn.degree()), int(target_Z.degree())],
        "P_dot_O": 6,
        "literal_weierstrass_identity": True,
        "maximum_rational_bits": maximum_bits((xn, xd, yn, yd)),
    },
    "method": {
        "norm4_enumeration": True,
        "integral_smith_solve": True,
        "exact_rootless_and_A1_group_law": True,
        "degree_one_pointed_quartic_transport": True,
        "groebner_or_elimination": False,
        "runtime_seconds": time.monotonic() - started,
    },
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in (PINNED, BRIDGE, MANIFEST, ROOTLESS_MODEL, ROOTLESS_BASIS, ROOTLESS_HEIGHTS, A1_CURVES, A1_RR)],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in (PINNED, BRIDGE, MANIFEST, ROOTLESS_MODEL, ROOTLESS_BASIS, ROOTLESS_HEIGHTS, A1_CURVES, A1_RR)},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "FIXEDREVERSE2A1HORIZONTAL|norm4_degree1=21|tail_rank=15|summands=7|"
    "degrees=(16,24,6)|PdotO=6|bits={}|seconds={:.3f}|status={}|output={}".format(
        payload["section"]["maximum_rational_bits"], payload["method"]["runtime_seconds"],
        payload["status"], OUTPUT,
    ), flush=True,
)
