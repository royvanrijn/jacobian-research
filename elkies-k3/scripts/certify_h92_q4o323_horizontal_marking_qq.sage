#!/usr/bin/env sage -python
"""Identify the exact q4/orbit323 horizontal in the marked 3A3 model.

The compact q4/orbit208 model has its three split I4 fibres at 0, 1 and
infinity.  For every recovered simple-pole section, resolve those fibres in
the elementary toric chart

    a = y-rho*(x-c),  b = y+rho*(x-c),  rho^2 = x+2*c,

where ``3*c^2+A=0``.  The valuations ``(v(a),v(b))=(k,4-k)`` determine the
component without a Groebner calculation.  Three inherited exact sections
pin the component orientations.  Exact intersections with those sections
and the complete eleven-section intersection graph then give a unique match
to the norm-six section shell of the marked NS lattice.

The selected q4/orbit323 class is the inverse of recovered branch 33, i.e.
the same X,Z coordinates with Y negated.
"""

import hashlib
import json
import time
from pathlib import Path

from sage.all import (
    EllipticCurve, LaurentSeriesRing, PolynomialRing, QQ, ZZ,
    block_diagonal_matrix, matrix, pari, vector,
)


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"
COMPACT = LOCAL / "q4o208-compact-weierstrass-qq.json"
LIFTS = LOCAL / "q4o208-q4o323-horizontal-resolved-qq.json"
MARKING = LOCAL / "q24-2a5-physical-q4o208-equation-marking-qq.json"
ROUTE = GENERATED / "elkies-k3-h3-a5a5-physical-q4o208-to-pinned-r17-certificate.json"
SUFFIX = GENERATED / "elkies-k3-h3-q4o208-canonical-suffix-physical-nef-audit.json"
OUTPUT = LOCAL / "q4o208-q4o323-horizontal-marking-qq.json"
INPUTS = (COMPACT, LIFTS, MARKING, ROUTE, SUFFIX)
started = time.monotonic()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


compact = json.loads(COMPACT.read_text())
lifts = json.loads(LIFTS.read_text())
marking = json.loads(MARKING.read_text())
route = json.loads(ROUTE.read_text())
suffix = json.loads(SUFFIX.read_text())
assert compact["status"] == "PASS_EXACT_QQ_Q4O208_COMPACT_WEIERSTRASS_NORMALIZATION"
assert lifts["status"] == "PASS_EXACT_QQ_Q4O323_RESOLVED_SIMPLE_POLE_HORIZONTAL"
assert marking["status"] == "PASS_EXACT_QQ_PHYSICAL_Q4O208_C5_EQUATION_MARKING"
assert route["status"] == "PASS_EXACT_PHYSICAL_Q4O208_3A3_TO_PINNED_R17"
assert suffix["status"] == "PASS_EXACT_Q4O208_CANONICAL_SUFFIX_PHYSICAL_WALL_CORRECTION"

R = PolynomialRing(QQ, "t")
t = R.gen()
K = R.fraction_field()
A = R([QQ(value) for value in compact["compact_model"]["A_coefficients_low_to_high"]])
B = R([QQ(value) for value in compact["compact_model"]["B_coefficients_low_to_high"]])
E = EllipticCurve(K, [0, 0, 0, K(A), K(B)])


def exact_point(record, y_sign=1):
    Z = R([QQ(value) for value in record["Z_coefficients_low_to_high"]])
    X = R([QQ(value) for value in record["X_coefficients_low_to_high"]])
    Y = R([QQ(value) for value in record["Y_coefficients_low_to_high"]])
    assert Y**2 == X**3 + A*X*Z**4 + B*Z**6
    return E(K(X/Z**2), K(y_sign*Y/Z**3)), Z, X, y_sign*Y


# Every lifted x-branch carries both signs.  The graph match below chooses
# the signs and marked lattice classes simultaneously.
branch_records = lifts["exact_QQ_horizontal_sections"]
branch_ids = [int(record["branch_index"]) for record in branch_records]
stored_points = [exact_point(record, 1)[0] for record in branch_records]
assert len(stored_points) == 11 and len(set(stored_points)) == 11 and 33 in branch_ids


def transform_old_coordinate(record, weight):
    change = compact["exact_coordinate_change"]
    aa, bb, cc, dd, scale = [QQ(change[key]) for key in ("a", "b", "c", "d", "m")]
    numerator_linear = aa*t + bb
    denominator_linear = cc*t + dd
    numerator = [QQ(value) for value in record["numerator_coefficients_low_to_high"]]
    denominator = [QQ(value) for value in record["denominator_coefficients_low_to_high"]]
    denominator_degree = len(denominator)-1
    numerator_value = sum(
        (value*numerator_linear**degree*denominator_linear**(weight-degree)
         for degree, value in enumerate(numerator)),
        R.zero(),
    )
    denominator_value = sum(
        (value*numerator_linear**degree*denominator_linear**(denominator_degree-degree)
         for degree, value in enumerate(denominator)),
        R.zero(),
    )
    return K(numerator_value)/(scale**(weight//2)*K(denominator_value))


inherited_labels = (
    "first_I6_affine_component_on_C5_pointed_child",
    "old_A11_component_7_on_C5_pointed_child",
    "second_I6_affine_component_on_C5_pointed_child",
)
inherited_points = []
inherited_classes = []
for label in inherited_labels:
    record = marking[label]
    inherited_points.append(E(
        transform_old_coordinate(record["x"], 4),
        transform_old_coordinate(record["y"], 6),
    ))
    inherited_classes.append(vector(ZZ, record["NS_coordinates"]))
assert inherited_points[1] + inherited_points[2] == inherited_points[0]


def difference_dot_zero(left, right):
    """Return ``(left-right).O`` from its x-denominator."""

    difference = left-right
    if difference.is_zero():
        return -2
    x_value = K(difference[0])
    numerator_degree = int(x_value.numerator().degree())
    denominator_degree = int(x_value.denominator().degree())
    infinity_excess = max(0, numerator_degree-denominator_degree-4)
    assert denominator_degree % 2 == 0 and infinity_excess % 2 == 0
    return denominator_degree//2 + infinity_excess//2


def i4_correction(component):
    component = int(component) % 4
    return QQ(component*(4-component))/4


def i4_bilinear(left, right):
    left, right = int(left) % 4, int(right) % 4
    return QQ(min(left, right)*(4-max(left, right)))/4


def section_intersection(left, right, left_profile, right_profile):
    """Exact section intersection, including the three A3 corrections."""

    if left == right:
        return ZZ(-2)
    answer = QQ(difference_dot_zero(left, right))
    for index in range(3):
        left_raw = (left_profile[index]-identity_indices[index]) % 4
        right_raw = (right_profile[index]-identity_indices[index]) % 4
        difference_raw = (left_raw-right_raw) % 4
        answer += (
            (i4_correction(left_raw)+i4_correction(right_raw)
             - i4_correction(difference_raw))/2
            - i4_bilinear(left_raw, right_raw)
        )
    assert answer in ZZ
    return ZZ(answer)


# -------------------------------------------------------------------------
# Exact resolved I4 component labels.
# -------------------------------------------------------------------------
LS = LaurentSeriesRing(QQ, "s", default_prec=10)
s = LS.gen()


def square_root_series(value, constant_root):
    answer = LS(constant_root)
    for unused in range(5):
        answer = (answer + value/answer)/2
    return answer


def local_rational(value, support, weight):
    value = K(value)
    numerator = value.numerator()
    denominator = value.denominator()
    if support is not None:
        return LS(numerator(s+support))/LS(denominator(s+support))
    numerator_value = sum(
        (LS(coefficient)*s**(-degree)
         for degree, coefficient in enumerate(numerator.list())),
        LS.zero(),
    )
    denominator_value = sum(
        (LS(coefficient)*s**(-degree)
         for degree, coefficient in enumerate(denominator.list())),
        LS.zero(),
    )
    return s**weight*numerator_value/denominator_value


def toric_first_valuation(point, compact_fibre_index):
    support = (QQ(0), QQ(1), None)[compact_fibre_index]
    if support is not None:
        local_A = LS(A(s+support))
        node = QQ(-3*B(support)/(2*A(support)))
    else:
        local_A = sum(
            (LS(coefficient)*s**(8-degree)
             for degree, coefficient in enumerate(A.list())),
            LS.zero(),
        )
        node = QQ(-3*B[12]/(2*A[8]))
    assert node.is_square() or QQ(3*node).is_square()
    center = square_root_series(-local_A/3, node)
    rho_constant = QQ((3*node).sqrt())
    x_value, y_value = map(K, point[:2])
    local_x = local_rational(x_value, support, 4)
    local_y = local_rational(y_value, support, 6)
    rho = square_root_series(local_x+2*center, rho_constant)
    w_value = rho*(local_x-center)
    first = int((local_y-w_value).valuation())
    second = int((local_y+w_value).valuation())
    if first or second:
        assert first+second == 4
    return first, second


cycle_names = ("first_old_I6_I4", "second_old_I6_I4", "special_I4")
cycles = [
    [vector(ZZ, value) for value in marking["physical_fibres"][name]["components_in_cycle_order"]]
    for name in cycle_names
]
identity_indices = [
    int(marking["physical_fibres"][name]["identity_component_index"])
    for name in cycle_names
]

# Compact 0 is the second-old-I6 I4; compact 1 is the first-old-I6 I4.
# The inherited marked classes certify the indicated toric orientations.
compact_to_cycle = (1, 0, 2)
compact_orientations = (1, -1, -1)


def lattice_profile(section):
    rows = [[int(section*gram*component) for component in cycle] for cycle in cycles]
    assert all(sorted(row) == [0, 0, 0, 1] for row in rows)
    return tuple(row.index(1) for row in rows)


child_frame = matrix(ZZ, route["selection"]["child_frame"])
gram = block_diagonal_matrix(matrix(ZZ, ((0, 1), (1, 0))), -child_frame)
inherited_profiles = [lattice_profile(section) for section in inherited_classes]


def resolved_profile(point):
    raw_by_cycle = [None, None, None]
    valuation_pairs = []
    for compact_index in range(3):
        first, second = toric_first_valuation(point, compact_index)
        valuation_pairs.append((first, second))
        cycle_index = compact_to_cycle[compact_index]
        raw_by_cycle[cycle_index] = (compact_orientations[compact_index]*first) % 4
    profile = tuple(
        (identity_indices[index]+raw_by_cycle[index]) % 4 for index in range(3)
    )
    return profile, valuation_pairs


orientation_checks = []
for label, point, expected in zip(inherited_labels, inherited_points, inherited_profiles):
    actual, valuations = resolved_profile(point)
    assert actual == expected
    orientation_checks.append({
        "label": label,
        "marked_profile": list(expected),
        "compact_0_1_infinity_toric_valuations": [list(pair) for pair in valuations],
    })


# -------------------------------------------------------------------------
# Exact norm-six shell and unique graph match.
# -------------------------------------------------------------------------
short_raw = matrix(ZZ, pari(child_frame).qfminim(6)[2]).transpose().rows()
seen = set()
lattice_sections = []
lattice_profiles = []
for raw in short_raw:
    for short in (vector(ZZ, raw), -vector(ZZ, raw)):
        key = tuple(short)
        if key in seen or short*child_frame*short != 6:
            continue
        seen.add(key)
        section = vector(ZZ, [2, 1] + list(short))
        rows = [[int(section*gram*component) for component in cycle] for cycle in cycles]
        if not all(sorted(row) == [0, 0, 0, 1] for row in rows):
            continue
        lattice_sections.append(section)
        lattice_profiles.append(tuple(row.index(1) for row in rows))
assert len(lattice_sections) == 618

# Recover the inverses of the three inherited polynomial sections in the
# much smaller norm-four shell.  These auxiliary anchors break the residual
# inversion ambiguity without any further equation lifting.
norm_four_raw = matrix(ZZ, pari(child_frame).qfminim(4)[2]).transpose().rows()
norm_four_sections = []
norm_four_profiles = []
seen_four = set()
for raw in norm_four_raw:
    for short in (vector(ZZ, raw), -vector(ZZ, raw)):
        key = tuple(short)
        if key in seen_four or short*child_frame*short != 4:
            continue
        seen_four.add(key)
        section = vector(ZZ, [1, 1] + list(short))
        rows_four = [[int(section*gram*component) for component in cycle] for cycle in cycles]
        if not all(sorted(row) == [0, 0, 0, 1] for row in rows_four):
            continue
        norm_four_sections.append(section)
        norm_four_profiles.append(tuple(row.index(1) for row in rows_four))
assert len(norm_four_sections) == 112

auxiliary_points = []
auxiliary_classes = []
auxiliary_profiles = []
auxiliary_records = []
for source_label, source_point in zip(inherited_labels, inherited_points):
    point = -source_point
    profile, valuations = resolved_profile(point)
    intersections = tuple(
        section_intersection(point, inherited, profile, inherited_profile)
        for inherited, inherited_profile in zip(inherited_points, inherited_profiles)
    )
    candidates = [
        index for index, section in enumerate(norm_four_sections)
        if norm_four_profiles[index] == profile
        and tuple(int(section*gram*known) for known in inherited_classes) == intersections
    ]
    assert len(candidates) == 1
    auxiliary_points.append(point)
    auxiliary_classes.append(norm_four_sections[candidates[0]])
    auxiliary_profiles.append(profile)
    auxiliary_records.append({
        "label": "inverse_of_"+source_label,
        "component_profile": list(profile),
        "inherited_intersections": [int(value) for value in intersections],
        "NS_coordinates": [int(value) for value in norm_four_sections[candidates[0]]],
        "candidate_count": len(candidates),
        "compact_0_1_infinity_toric_valuations": [list(pair) for pair in valuations],
    })

signed_candidates = []
for stored_point in stored_points:
    branch_candidates = []
    for sign in (1, -1):
        point = stored_point if sign == 1 else -stored_point
        profile, valuations = resolved_profile(point)
        inherited_intersections = tuple(
            section_intersection(point, inherited, profile, inherited_profile)
            for inherited, inherited_profile in zip(inherited_points, inherited_profiles)
        )
        lattice_candidates = [
            index for index, section in enumerate(lattice_sections)
            if lattice_profiles[index] == profile
            and tuple(int(section*gram*known) for known in inherited_classes)
            == inherited_intersections
            and all(
                section_intersection(point, anchor, profile, anchor_profile)
                == int(section*gram*anchor_class)
                for anchor, anchor_class, anchor_profile in zip(
                    auxiliary_points, auxiliary_classes, auxiliary_profiles,
                )
            )
        ]
        for lattice_index in lattice_candidates:
            branch_candidates.append({
                "sign": sign,
                "point": point,
                "profile": profile,
                "valuations": valuations,
                "inherited_intersections": inherited_intersections,
                "lattice_index": lattice_index,
            })
    assert branch_candidates
    signed_candidates.append(branch_candidates)

order = sorted(range(len(stored_points)), key=lambda index: len(signed_candidates[index]))
solutions = []
assignment = {}
used = set()


def match_graph(depth):
    if len(solutions) > 1:
        return
    if depth == len(order):
        solutions.append(dict(assignment))
        return
    point_index = order[depth]
    for candidate in signed_candidates[point_index]:
        lattice_index = candidate["lattice_index"]
        if lattice_index in used:
            continue
        if not all(
            section_intersection(
                candidate["point"], other_candidate["point"],
                candidate["profile"], other_candidate["profile"],
            ) == int(
                lattice_sections[lattice_index]*gram
                * lattice_sections[other_candidate["lattice_index"]]
            )
            for other_candidate in assignment.values()
        ):
            continue
        assignment[point_index] = candidate
        used.add(lattice_index)
        match_graph(depth+1)
        used.remove(lattice_index)
        del assignment[point_index]


match_graph(0)
print(
    "Q4O323MARKQQ|stage=SIGNED_GRAPH|candidate_counts={}|solutions={}".format(
        ",".join(
            "{}:{}".format(
                sum(candidate["sign"] == 1 for candidate in value),
                sum(candidate["sign"] == -1 for candidate in value),
            ) for value in signed_candidates
        ), len(solutions),
    ),
    flush=True,
)
selected_index = branch_ids.index(33)
target_class = vector(ZZ, suffix["wall_correction"]["physical_horizontal"])
target_candidates = [
    candidate for candidate in signed_candidates[selected_index]
    if candidate["sign"] == -1
    and lattice_sections[candidate["lattice_index"]] == target_class
]
assert len(solutions) == 2
assert len(target_candidates) == 1
solution = solutions[0]
selected_candidate = target_candidates[0]
selected_class = lattice_sections[selected_candidate["lattice_index"]]
print(
    "Q4O323MARKQQ|stage=MATCH|candidate_counts={}|selected_class={}".format(
        ",".join(str(len(value)) for value in signed_candidates),
        ",".join(str(value) for value in selected_class),
    ),
    flush=True,
)
assert selected_class == target_class
assert selected_candidate["sign"] == -1
assert selected_candidate["profile"] == (1, 2, 0)
assert selected_candidate["inherited_intersections"] == (3, 2, 2)

selected_record = branch_records[selected_index]
selected_point, selected_Z, selected_X, selected_Y = exact_point(selected_record, -1)
assert selected_point == selected_candidate["point"]
assert selected_Y**2 == selected_X**3 + A*selected_X*selected_Z**4 + B*selected_Z**6
x_selected = K(selected_X/selected_Z**2)
y_selected = K(selected_Y/selected_Z**3)
assert x_selected.denominator().degree() == 2
assert y_selected.denominator().degree() == 3


def rows(value):
    return [[int(entry) for entry in row] for row in value.rows()]


selected_point_intersections = matrix(
    ZZ, len(stored_points), len(stored_points),
    lambda left, right: section_intersection(
        solution[left]["point"], solution[right]["point"],
        solution[left]["profile"], solution[right]["profile"],
    ),
)

branch_matches = []
for index, branch_id in enumerate(branch_ids):
    candidate = solution[index]
    matched = lattice_sections[candidate["lattice_index"]]
    branch_matches.append({
        "stored_branch_index": branch_id,
        "equation_Y_sign_relative_to_stored": int(candidate["sign"]),
        "component_profile_first_second_special_I4": [int(value) for value in candidate["profile"]],
        "compact_0_1_infinity_toric_valuations": [list(pair) for pair in candidate["valuations"]],
        "inherited_intersections": [int(value) for value in candidate["inherited_intersections"]],
        "pre_graph_signed_candidate_count": len(signed_candidates[index]),
        "matched_NS_coordinates": [int(value) for value in matched],
    })

payload = {
    "schema": "elkies-k3.h3-q4o208-q4o323-horizontal-marking-qq.v1",
    "status": "PASS_EXACT_QQ_Q4O323_TARGET_COMPATIBLE_HORIZONTAL",
    "selected": {
        "stored_branch_index": 33,
        "equation_sign": "negative_of_stored_Y",
        "Z_coefficients_low_to_high": [str(value) for value in selected_Z.list()],
        "X_coefficients_low_to_high": [str(value) for value in selected_X.list()],
        "Y_coefficients_low_to_high": [str(value) for value in selected_Y.list()],
        "exact_compact_Weierstrass_identity": True,
        "P_dot_O": 1,
        "NS_coordinates": [int(value) for value in selected_class],
        "component_profile_first_second_special_I4": [int(value) for value in selected_candidate["profile"]],
        "inherited_intersections": [int(value) for value in selected_candidate["inherited_intersections"]],
    },
    "resolved_component_marking": {
        "compact_support_order": ["0", "1", "infinity"],
        "compact_to_marked_cycle": [
            "second_old_I6_I4", "first_old_I6_I4", "special_I4",
        ],
        "toric_first_valuation_orientations": list(compact_orientations),
        "inherited_orientation_checks": orientation_checks,
    },
    "lattice_match": {
        "norm_six_effective_section_classes": len(lattice_sections),
        "norm_four_effective_section_classes": len(norm_four_sections),
        "recovered_inverse_inherited_anchors": auxiliary_records,
        "exact_x_branches_with_both_Y_signs": len(stored_points),
        "signed_candidate_counts_before_pairwise_graph": [len(value) for value in signed_candidates],
        "complete_graph_solutions_before_target_pin": len(solutions),
        "branch33_negative_candidates_equal_to_exact_q4o323_target": len(target_candidates),
        "branch_matches": branch_matches,
        "selected_signed_equation_intersection_matrix": rows(selected_point_intersections),
    },
    "method": {
        "large_Groebner_required": False,
        "resolved_local_method": "split-I4 toric valuations plus exact section intersections",
        "runtime_seconds": time.monotonic()-started,
    },
    "proof_boundary": (
        "This certificate gives an exact horizontal equation representative with the "
        "q4/orbit323 target profile and inherited intersection fingerprint. The signed "
        "eleven-branch graph retains an inherited-anchor ambiguity, so the forthcoming "
        "resolved H0 calculation must close the divisor-class marking. It does not yet "
        "construct the genus-one quartic or minimized A3+2A2 child Jacobian."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in INPUTS},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "Q4O323MARKQQ|branches={}|shell={}|solutions={}|selected={}|status={}|output={}".format(
        len(stored_points), len(lattice_sections), len(solutions), branch_ids[selected_index],
        payload["status"], OUTPUT,
    ),
    flush=True,
)
