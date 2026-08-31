#!/usr/bin/env sage -python
"""Recover the physical q4/o323 horizontal section modulo 131.

The q4/o208 child has three I4 fibres and MW rank eight.  Polynomial
sections meeting two prescribed nonidentity components are found by a
three-parameter interpolation search: no Groebner basis is used.  The
16-section intersection of two such shells is matched to the exact NS shell
by component labels and all pairwise section intersections.  The two exact
constant-old-base points over the second old I6 support distinguish C7 from
the second affine component by the exact group relation

    second_affine = first_affine - C7.

Three overlapping 131^3 shells test whether polynomial sections span the
target q4/o323 horizontal class modulo the trivial lattice.  They do not: the
polynomial subgroup has rank seven in the rank-eight Mordell--Weil group, and
the simple-pole target supplies the missing direction.  This is an exact
finite-field/lattice obstruction, not a characteristic-zero equation lift.
"""

import hashlib
import json
import time
from collections import Counter, defaultdict
from itertools import permutations, product
from pathlib import Path

from sage.all import (
    EllipticCurve, GF, PolynomialRing, QQ, ZZ, block_diagonal_matrix,
    factorial, matrix, pari, vector,
)


ROOT = Path(__file__).resolve().parents[2]
if not (ROOT / "MATH_STATUS.json").exists():
    ROOT = Path.cwd()
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"
RR_PATH = LOCAL / "q24-2a5-physical-q4o208-rr-qq.json"
PARENT_PATH = LOCAL / "q24-a11-to-2a5-q8-resolved-rr-qq.json"
EQUATION_MARKING_PATH = LOCAL / "q24-2a5-physical-q4o208-equation-marking-qq.json"
SOURCE_MARKING_PATH = GENERATED / "elkies-k3-h3-a5a5-physical-component-chamber-marking.json"
ROUTE_PATH = GENERATED / "elkies-k3-h3-a5a5-physical-q4o208-to-pinned-r17-certificate.json"
SUFFIX_AUDIT_PATH = GENERATED / "elkies-k3-h3-q4o208-canonical-suffix-physical-nef-audit.json"
COMPACT_PATH = LOCAL / "q4o208-compact-weierstrass-qq.json"
EXACT_LIFTS_PATH = LOCAL / "q4o208-q4o323-horizontal-resolved-qq.json"
SUMMANDS_PATH = LOCAL / "q4o208-q4o323-horizontal-via-summands-qq.json"
OUTPUT = LOCAL / "q4o208-physical-q4o323-horizontal-mod131.json"
INPUTS = (
    RR_PATH, PARENT_PATH, EQUATION_MARKING_PATH, SOURCE_MARKING_PATH,
    ROUTE_PATH, SUFFIX_AUDIT_PATH, COMPACT_PATH, EXACT_LIFTS_PATH, SUMMANDS_PATH,
)

started = time.monotonic()


def log(stage, **fields):
    tail = "|".join(f"{key}={value}" for key, value in fields.items())
    print(
        f"Q4O323HMOD131|stage={stage}|elapsed={time.monotonic()-started:.3f}"
        + (f"|{tail}" if tail else ""),
        flush=True,
    )


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


rr = json.loads(RR_PATH.read_text())
parent = json.loads(PARENT_PATH.read_text())
equation_marking = json.loads(EQUATION_MARKING_PATH.read_text())
source_marking = json.loads(SOURCE_MARKING_PATH.read_text())
route = json.loads(ROUTE_PATH.read_text())
suffix_audit = json.loads(SUFFIX_AUDIT_PATH.read_text())
compact = json.loads(COMPACT_PATH.read_text())
exact_lifts = json.loads(EXACT_LIFTS_PATH.read_text())
summands = json.loads(SUMMANDS_PATH.read_text())
assert rr["status"] == "PASS_EXACT_QQ_PHYSICAL_Q4O208_3A3_RR_AND_JACOBIAN"
assert parent["status"] == "PASS_EXACT_Q24_A11_Q8_2A5_RESOLVED_RR"
assert equation_marking["status"] == "PASS_EXACT_QQ_PHYSICAL_Q4O208_C5_EQUATION_MARKING"
assert source_marking["status"] == "PASS_EXACT_A5A5_PHYSICAL_COMPONENT_CHAMBER_MARKING"
assert route["status"] == "PASS_EXACT_PHYSICAL_Q4O208_3A3_TO_PINNED_R17"
assert suffix_audit["status"] == "PASS_EXACT_Q4O208_CANONICAL_SUFFIX_PHYSICAL_WALL_CORRECTION"
assert compact["status"] == "PASS_EXACT_QQ_Q4O208_COMPACT_WEIERSTRASS_NORMALIZATION"
assert exact_lifts["status"] == "PASS_EXACT_QQ_Q4O323_RESOLVED_SIMPLE_POLE_HORIZONTAL"
assert summands["status"] == "PASS_EXACT_QQ_Q4O323_SUM_REPRODUCES_BRANCH33_INVERSE_TARGET_EXCLUDED"


# Exact constant-old-base points over the second old I6 support.
RUQ = PolynomialRing(QQ, "U")
UQ = RUQ.gen()
KUQ = RUQ.fraction_field()
RTQ = PolynomialRing(KUQ, "T")
TQ = RTQ.gen()
R0 = PolynomialRing(QQ, "T")


def rational_function(record):
    numerator = RUQ([QQ(value) for value in record["numerator_coefficients_low_to_high"]])
    denominator = RUQ([QQ(value) for value in record["denominator_coefficients_low_to_high"]])
    return KUQ(numerator) / KUQ(denominator)


i6_roots = []
for record in parent["child"]["discriminant_factorization"]:
    if int(record["multiplicity"]) != 6:
        continue
    factor = R0(record["factor"])
    i6_roots.append(QQ(-factor[0] / factor[1]))
assert len(i6_roots) == 2
Fp103 = GF(103)
roots_by_reduction = {int(Fp103(value)): value for value in i6_roots}
beta, gamma = roots_by_reduction[89], roots_by_reduction[68]

quartic = RTQ([
    KUQ(value) for value in rr["quartic"]["coefficients_in_old_T_low_to_high"]
])
jet = [
    KUQ(quartic.derivative(order)(KUQ(beta)) / factorial(order))
    for order in range(5)
]
e, d, c, b, a = jet
w0 = rational_function(equation_marking["resolved_C5_slice"]["quartic_ordinate_on_C5"])
assert e == w0**2

kernel = matrix(QQ, [
    [QQ(value) for value in row] for row in rr["resolved_RR"]["kernel_basis"]
])
b0, b1 = kernel[0, 3], kernel[1, 3]
denominator = KUQ(UQ) * KUQ(b0) - KUQ(b1)
a1 = d / w0
a2 = c - d**2 / (4 * w0**2)
a3 = 2 * w0 * b
a4 = -4 * w0**2 * a
a6 = a2 * a4
b2 = a1**2 + 4 * a2
A_exact = KUQ(RUQ([QQ(value) for value in rr["child"]["minimal_A_coefficients_low_to_high"]]))
B_exact = KUQ(RUQ([QQ(value) for value in rr["child"]["minimal_B_coefficients_low_to_high"]]))
E_exact = EllipticCurve(KUQ, [0, 0, 0, A_exact, B_exact])


def square_root_rational(value):
    numerator = RUQ(value.numerator())
    denominator_value = RUQ(value.denominator())
    numerator_lc = QQ(numerator.leading_coefficient())
    denominator_lc = QQ(denominator_value.leading_coefficient())
    assert numerator_lc.is_square() and denominator_lc.is_square()
    numerator_root = RUQ(numerator_lc.sqrt())
    denominator_root = RUQ(denominator_lc.sqrt())
    for factor, exponent in numerator.monic().factor():
        assert exponent % 2 == 0
        numerator_root *= factor ** (exponent // 2)
    for factor, exponent in denominator_value.monic().factor():
        assert exponent % 2 == 0
        denominator_root *= factor ** (exponent // 2)
    root = KUQ(numerator_root / denominator_root)
    assert root**2 == value
    return root


gamma_ordinate = square_root_rational(KUQ(quartic(KUQ(gamma))))
z_gamma = KUQ(gamma - beta)


def pointed_map(ordinate):
    x_general = (2 * w0 * (ordinate + w0) + d * z_gamma) / z_gamma**2
    y_general = (
        4 * w0**2 * (ordinate + w0) + 2 * w0 * d * z_gamma
        + (2 * w0 * c - d**2 / (2 * w0)) * z_gamma**2
    ) / z_gamma**3
    x_short = KUQ(denominator**4 * 9 * (x_general + b2 / 12))
    y_short = KUQ(
        denominator**6 * 27 * (y_general + (a1 * x_general + a3) / 2)
    )
    point = E_exact(x_short, y_short)
    assert y_short**2 == x_short**3 + A_exact * x_short + B_exact
    return point


gamma_points = [pointed_map(gamma_ordinate), pointed_map(-gamma_ordinate)]
first_affine_record = equation_marking["first_I6_affine_component_on_C5_pointed_child"]
first_affine = E_exact(
    rational_function(first_affine_record["x"]),
    rational_function(first_affine_record["y"]),
)
assert gamma_points[0] + gamma_points[1] == first_affine
log("EXACT_GAMMA", relation="gamma_point_0+gamma_point_1=first_affine")


# Good-reduction finite-field model and interpolation shells.
finite = GF(131)
RU = PolynomialRing(finite, "U")
U = RU.gen()
KU = RU.fraction_field()
RX = PolynomialRing(finite, "x")
x_variable = RX.gen()


def reduce_qq(value):
    value = QQ(value)
    return finite(value.numerator()) / finite(value.denominator())


def reduce_exact(value):
    numerator = RU([reduce_qq(coefficient) for coefficient in RUQ(value.numerator()).list()])
    denominator_value = RU([
        reduce_qq(coefficient) for coefficient in RUQ(value.denominator()).list()
    ])
    return KU(numerator) / KU(denominator_value)


A = RU([reduce_qq(value) for value in rr["child"]["minimal_A_coefficients_low_to_high"]])
B = RU([reduce_qq(value) for value in rr["child"]["minimal_B_coefficients_low_to_high"]])
E = EllipticCurve(KU, [0, 0, 0, KU(A), KU(B)])
supports = [
    reduce_qq(record["support"]) for record in equation_marking["physical_fibres"].values()
]
identity_indices = [
    int(record["identity_component_index"])
    for record in equation_marking["physical_fibres"].values()
]
assert len(set(supports)) == 3
discriminant = -finite(16) * (4 * A**3 + 27 * B**2)
assert all(((U - support)**4).divides(discriminant) for support in supports)
nodes = {}
for support in supports:
    cubic = x_variable**3 + A(support) * x_variable + B(support)
    repeated = cubic.gcd(cubic.derivative())
    assert repeated.degree() == 1
    nodes[support] = -repeated[0] / repeated[1]


def polynomial_square_roots(polynomial):
    if polynomial == 0:
        return (RU.zero(),)
    shift = next(value for value in finite if polynomial(value) != 0)
    shifted = polynomial(U + shift)
    constant = shifted[0]
    if not constant.is_square():
        return ()
    roots = []
    for first in constant.sqrt(all=True):
        coefficients = [first]
        for degree in range(1, 7):
            known = sum(
                (coefficients[left] * coefficients[degree-left]
                 for left in range(1, degree)),
                finite.zero(),
            )
            coefficients.append((shifted[degree] - known) / (2 * first))
        candidate_shifted = RU(coefficients)
        if candidate_shifted**2 == shifted:
            roots.append(candidate_shifted(U - shift))
    return tuple(roots)


def search_pair(pair):
    pair_supports = [supports[index] for index in pair]
    interpolation = RU.lagrange_polynomial([
        (support, nodes[support]) for support in pair_supports
    ])
    product_value = RU.one()
    for support in pair_supports:
        product_value *= U - support
    vanishing = product_value
    answers = []
    search_started = time.monotonic()
    for constant in finite:
        for linear in finite:
            for quadratic in finite:
                X = interpolation + (constant + linear * U + quadratic * U**2) * vanishing
                for Y in polynomial_square_roots(X**3 + A * X + B):
                    answers.append(E(KU(X), KU(Y)))
    log(
        "PAIR_SEARCH", pair=f"{pair[0]}{pair[1]}", tests=131**3,
        sections=len(answers), seconds=f"{time.monotonic()-search_started:.3f}",
    )
    return answers


def hits_node(point, support):
    if point.is_zero():
        return False
    x_coordinate, y_coordinate = point[0], point[1]
    return (
        x_coordinate.denominator()(support) != 0
        and x_coordinate(support) == nodes[support]
        and y_coordinate(support) == 0
    )


shell_01 = search_pair((0, 1))
shell_02 = search_pair((0, 2))
shell_12 = search_pair((1, 2))
anchors = [
    point for point in shell_01
    if all(hits_node(point, support) for support in supports)
]
assert len(anchors) == 16
assert set(anchors) == {
    point for point in shell_02
    if all(hits_node(point, support) for support in supports)
}

references = {
    support: next(point for point in anchors if hits_node(2 * point, support))
    for support in supports
}


def raw_component_label(point, support):
    answers = [
        multiplier for multiplier in range(4)
        if not hits_node(point - multiplier * references[support], support)
    ]
    assert len(answers) == 1
    return answers[0]


def component_profile(point, orientation=(1, 1, 1)):
    return tuple(
        (identity_indices[index]
         + orientation[index] * raw_component_label(point, support)) % 4
        for index, support in enumerate(supports)
    )


def section_intersection(left, right):
    difference = left - right
    if difference.is_zero():
        return -2
    x_coordinate = difference[0]
    numerator_degree = int(x_coordinate.numerator().degree())
    denominator_degree = int(x_coordinate.denominator().degree())
    assert denominator_degree % 2 == 0
    infinity_excess = max(0, numerator_degree - denominator_degree - 4)
    assert infinity_excess % 2 == 0
    return denominator_degree // 2 + infinity_excess // 2


# Exact norm-four section shell in the marked NS lattice.
child_frame = matrix(ZZ, route["selection"]["child_frame"])
gram = block_diagonal_matrix(matrix(ZZ, ((0, 1), (1, 0))), -child_frame)
cycles = [
    [vector(ZZ, value) for value in record["components_in_cycle_order"]]
    for record in equation_marking["physical_fibres"].values()
]
half_short = matrix(ZZ, pari(child_frame).qfminim(4)[2]).transpose().rows()
lattice_sections = []
for half_vector in half_short:
    for short_vector in (vector(ZZ, half_vector), -vector(ZZ, half_vector)):
        if short_vector * child_frame * short_vector != 4:
            continue
        section = vector(ZZ, [1, 1] + list(short_vector))
        hits = [
            [int(section * gram * component) for component in cycle]
            for cycle in cycles
        ]
        if not all(sorted(row) == [0, 0, 0, 1] for row in hits):
            continue
        profile = tuple(row.index(1) for row in hits)
        lattice_sections.append((section, profile))
assert len(lattice_sections) == 112
lattice_anchors = [
    (section, profile) for section, profile in lattice_sections
    if all(profile[index] != identity_indices[index] for index in range(3))
]
assert len(lattice_anchors) == 16

# There are two component-orientation choices, related by global inversion.
# Test the small profile permutations exactly, then use the exact C7 point.
anchor_intersections = matrix(
    ZZ, 16, 16,
    lambda left, right: section_intersection(anchors[left], anchors[right]),
)
lattice_intersections = matrix(
    ZZ, 16, 16,
    lambda left, right: int(lattice_anchors[left][0] * gram * lattice_anchors[right][0]),
)
anchor_solutions = []
for orientation in product((1, 3), repeat=3):
    point_profiles = [component_profile(point, orientation) for point in anchors]
    lattice_profiles = [profile for unused, profile in lattice_anchors]
    if Counter(point_profiles) != Counter(lattice_profiles):
        continue
    point_groups = defaultdict(list)
    lattice_groups = defaultdict(list)
    for index, profile in enumerate(point_profiles):
        point_groups[profile].append(index)
    for index, profile in enumerate(lattice_profiles):
        lattice_groups[profile].append(index)
    groups = sorted(point_groups)
    for choices in product(*[
        tuple(permutations(lattice_groups[group])) for group in groups
    ]):
        assignment = {
            point_index: lattice_index
            for group, permutation in zip(groups, choices)
            for point_index, lattice_index in zip(point_groups[group], permutation)
        }
        if all(
            anchor_intersections[left, right]
            == lattice_intersections[assignment[left], assignment[right]]
            for left in range(16) for right in range(16)
        ):
            anchor_solutions.append((orientation, assignment))
assert len(anchor_solutions) == 2

log("ANCHOR_MATCH", solutions=len(anchor_solutions), ambiguity="global_inversion")


# Reduce every exact simple-pole lift into the already marked mod-131 shell.
# The compact coordinate t and this script's q4/o208 coordinate U satisfy
# U=(a*t+b)/(c*t+d), with the stated weighted Weierstrass transform.
change = compact["exact_coordinate_change"]
change_a, change_b, change_c, change_d, change_m = [
    reduce_qq(change[key]) for key in ("a", "b", "c", "d", "m")
]
compact_t = KU((change_b-change_d*U)/(change_c*U-change_a))
compact_denominator = KU(change_c*compact_t+change_d)


def evaluate_compact_polynomial(values):
    return sum(
        (reduce_qq(value)*compact_t**degree for degree, value in enumerate(values)),
        KU.zero(),
    )


def reduced_exact_branch(record, sign):
    Zc = evaluate_compact_polynomial(record["Z_coefficients_low_to_high"])
    Xc = evaluate_compact_polynomial(record["X_coefficients_low_to_high"])
    Yc = evaluate_compact_polynomial(record["Y_coefficients_low_to_high"])
    x_old = KU(change_m**2*Xc/(Zc**2*compact_denominator**4))
    y_old = KU(sign*change_m**3*Yc/(Zc**3*compact_denominator**6))
    point = E(x_old, y_old)
    return point


exact_polynomial_points = []
for node_type, records in summands["exact_QQ_polynomial_sections"].items():
    for record in records:
        Xc = evaluate_compact_polynomial(record["X_coefficients_low_to_high"])
        Yc = evaluate_compact_polynomial(record["Y_coefficients_low_to_high"])
        point = E(
            KU(change_m**2*Xc/compact_denominator**4),
            KU(change_m**3*Yc/compact_denominator**6),
        )
        exact_polynomial_points.append((node_type, int(record["signed_index"]), point))


norm_six_raw = matrix(ZZ, pari(child_frame).qfminim(6)[2]).transpose().rows()
lattice_norm_six = []
seen_norm_six = set()
for raw in norm_six_raw:
    for short_vector in (vector(ZZ, raw), -vector(ZZ, raw)):
        key = tuple(short_vector)
        if key in seen_norm_six or short_vector*child_frame*short_vector != 6:
            continue
        seen_norm_six.add(key)
        section = vector(ZZ, [2, 1] + list(short_vector))
        hits = [
            [int(section*gram*component) for component in cycle]
            for cycle in cycles
        ]
        if not all(sorted(row) == [0, 0, 0, 1] for row in hits):
            continue
        lattice_norm_six.append((section, tuple(row.index(1) for row in hits)))
assert len(lattice_norm_six) == 618


def map_exact_branches(solution):
    orientation, anchor_assignment = solution
    anchor_map = {
        anchors[index]: lattice_anchors[anchor_assignment[index]][0]
        for index in range(16)
    }
    records = []
    for branch in exact_lifts["exact_QQ_horizontal_sections"]:
        for sign in (1, -1):
            point = reduced_exact_branch(branch, sign)
            profile = component_profile(point, orientation)
            matches = [
                section for section, candidate_profile in lattice_norm_six
                if candidate_profile == profile
                and all(
                    section_intersection(point, anchor)
                    == int(section*gram*anchor_class)
                    for anchor, anchor_class in anchor_map.items()
                )
            ]
            assert len(matches) == 1
            records.append({
                "stored_branch_index": int(branch["branch_index"]),
                "Y_sign_relative_to_stored": sign,
                "component_profile": list(profile),
                "NS_coordinates": [int(value) for value in matches[0]],
            })
    return records


exact_branch_maps = [map_exact_branches(solution) for solution in anchor_solutions]
log("EXACT_BRANCH_MAP", signed=22, norm_six=618, orientations=2)


def map_shell(points, pair, orientation, anchor_map):
    candidates = [
        (section, profile) for section, profile in lattice_sections
        if all(profile[index] != identity_indices[index] for index in pair)
    ]
    answer = {}
    for point in points:
        profile = component_profile(point, orientation)
        matches = [
            section for section, candidate_profile in candidates
            if candidate_profile == profile
            and all(
                section_intersection(point, anchor) == int(section * gram * anchor_class)
                for anchor, anchor_class in anchor_map.items()
            )
        ]
        assert len(matches) == 1
        answer[point] = matches[0]
    return answer


old_fibre = vector(ZZ, [1, 0] + [0] * 17)
old_zero = vector(ZZ, [-1, 1] + [0] * 17)
root_components = []
for cycle, identity_index in zip(cycles, identity_indices):
    root_components.extend(
        component for index, component in enumerate(cycle) if index != identity_index
    )
trivial_rows = [old_fibre, old_zero] + root_components
target_class = vector(ZZ, suffix_audit["wall_correction"]["physical_horizontal"])
expected_profile = tuple(
    [int(target_class * gram * component) for component in cycle].index(1)
    for cycle in cycles
)


def recover_target(solution):
    orientation, anchor_assignment = solution
    anchor_map = {
        anchors[index]: lattice_anchors[anchor_assignment[index]][0]
        for index in range(16)
    }
    mapped = map_shell(shell_01, (0, 1), orientation, anchor_map)
    mapped.update(map_shell(shell_02, (0, 2), orientation, anchor_map))
    mapped.update(map_shell(shell_12, (1, 2), orientation, anchor_map))
    exact_polynomial_maps = []
    for node_type, signed_index, point in exact_polynomial_points:
        profile = component_profile(point, orientation)
        matches = [
            section for section, candidate_profile in lattice_sections
            if candidate_profile == profile
            and all(
                section_intersection(point, anchor)
                == int(section*gram*anchor_class)
                for anchor, anchor_class in anchor_map.items()
            )
        ]
        assert len(matches) == 1
        exact_polynomial_maps.append({
            "node_type": node_type,
            "signed_index": signed_index,
            "component_profile": list(profile),
            "NS_coordinates": [int(value) for value in matches[0]],
        })
    points = list(mapped)
    span_matrix = matrix(
        ZZ, trivial_rows + [mapped[point] for point in points]
    ).transpose()
    span_rank = span_matrix.rank()
    augmented_rank = span_matrix.augment(
        matrix(ZZ, [target_class]).transpose()
    ).rank()
    assert span_rank == 18 and augmented_rank == 19
    return {
        "orientation": orientation, "mapped_count": len(mapped),
        "mapped_union_NS_coordinates": sorted({
            tuple(int(value) for value in section)
            for section in mapped.values()
        }),
        "exact_QQ_polynomial_section_maps": exact_polynomial_maps,
        "span_rank_with_trivial_lattice": span_rank,
        "augmented_rank_with_target": augmented_rank,
    }


target_recoveries = [recover_target(solution) for solution in anchor_solutions]
mapped_count = target_recoveries[0]["mapped_count"]
assert target_recoveries[1]["mapped_count"] == mapped_count
target_hits_by_orientation = [
    [
        record for record in records
        if vector(ZZ, record["NS_coordinates"]) == target_class
    ]
    for records in exact_branch_maps
]
target_hit_counts = [len(records) for records in target_hits_by_orientation]
log("EXACT_TARGET_HITS", counts=",".join(map(str, target_hit_counts)))
log("SHELL_MATCH", mapped=mapped_count, polynomial_MW_rank=7, target_outside=1)


def finite_record(value):
    value = KU(value)
    return {
        "numerator_coefficients_low_to_high": [int(entry) for entry in value.numerator().list()],
        "denominator_coefficients_low_to_high": [int(entry) for entry in value.denominator().list()],
        "degrees_numerator_denominator": [
            int(value.numerator().degree()), int(value.denominator().degree()),
        ],
    }


def exact_record(value):
    value = KUQ(value)
    return {
        "numerator_coefficients_low_to_high": [str(entry) for entry in value.numerator().list()],
        "denominator_coefficients_low_to_high": [str(entry) for entry in value.denominator().list()],
        "degrees_numerator_denominator": [
            int(value.numerator().degree()), int(value.denominator().degree()),
        ],
    }


payload = {
    "schema": "elkies-k3.h3-q4o208-physical-q4o323-horizontal-mod131.v1",
    "status": "PASS_EXACT_Q4O323_POLYNOMIAL_SECTION_SUBGROUP_OBSTRUCTION",
    "prime": 131,
    "exact_inherited_sections": {
        "unordered_C7_second_affine_pair": [
            {"x": exact_record(point[0]), "y": exact_record(point[1])}
            for point in gamma_points
        ],
        "exact_group_relation": "gamma_point_0 + gamma_point_1 = first_affine",
        "remaining_label_ambiguity": "C7 and second affine are exchanged",
    },
    "search": {
        "pair_01_tests": 131**3,
        "pair_01_sections": len(shell_01),
        "pair_02_tests": 131**3,
        "pair_02_sections": len(shell_02),
        "pair_12_tests": 131**3,
        "pair_12_sections": len(shell_12),
        "anchor_sections": len(anchors),
        "lattice_norm4_sections": len(lattice_sections),
        "lattice_norm6_sections": len(lattice_norm_six),
        "mapped_union_sections": mapped_count,
        "component_orientations": [
            list(recovery["orientation"]) for recovery in target_recoveries
        ],
        "large_Groebner_required": False,
    },
    "target": {
        "NS_coordinates": [int(value) for value in target_class],
        "P_dot_O": 1,
        "component_profile": list(expected_profile),
        "in_polynomial_section_subgroup": False,
        "span_rank_with_trivial_lattice": 18,
        "augmented_rank_with_target": 19,
    },
    "exact_simple_pole_mapping": {
        "maps_by_global_component_orientation": exact_branch_maps,
        "target_hits_by_global_component_orientation": target_hits_by_orientation,
        "target_hit_counts_by_global_component_orientation": target_hit_counts,
        "interpretation": (
            "The two entries differ by the single global inversion ambiguity of the "
            "mod-131 anchor marking. The exact signed characteristic-zero branches are "
            "reported without assuming that the target is represented by one of them."
        ),
    },
    "mapped_polynomial_section_subgroup": [
        {
            "global_component_orientation": list(recovery["orientation"]),
            "union_NS_coordinates": recovery["mapped_union_NS_coordinates"],
            "exact_QQ_polynomial_section_maps": recovery[
                "exact_QQ_polynomial_section_maps"
            ],
        }
        for recovery in target_recoveries
    ],
    "proof_boundary": (
        "The unordered pair of inherited constant-old-base sections is exact over QQ; labeling "
        "the pair as C7 versus second affine remains a resolved-component gate. The complete "
        "two-node polynomial-section shells over GF(131) map to an MW-rank-seven subgroup. "
        "The q4/o323 simple-pole horizontal class is outside it, so a direct simple-pole seed "
        "or a different presentation is still required."
    ),
    "runtime_seconds": time.monotonic() - started,
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in INPUTS},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
log(
    "DONE", pole=1, profile="".join(map(str, expected_profile)),
    polynomial_MW_rank=7,
    status=payload["status"], output=OUTPUT,
)
