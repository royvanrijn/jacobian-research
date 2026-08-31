#!/usr/bin/env sage -python
"""Anchor the fibre-safe q4/o323 polynomial shell to the exact child lattice.

status: ACTIVE_SEARCH
claim: exact component and intersection domains for every regular p=61 shell point
inputs: q323 child equation, component-2 pointing, reflected child marking, p=61 shell
outputs: artifacts/local/elkies-k3/q4o323-p0-shell-anchor-domains-mod61.json

The exact opposite t=0 branch is used as a P.O=1 anchor.  Component labels
are computed by specialization in the elliptic group, including the I4 fibre
at infinity.  The output retains every remaining finite-component orientation
instead of choosing one heuristically.  No Groebner basis is used.
"""

import hashlib
import json
import time
from collections import Counter
from itertools import product
from pathlib import Path

from sage.all import (
    EllipticCurve, GF, PolynomialRing, QQ, ZZ, block_diagonal_matrix,
    matrix, pari, vector,
)


ROOT = Path(__file__).resolve().parents[2]
if not (ROOT / "MATH_STATUS.json").exists():
    ROOT = Path.cwd()
LOCAL = ROOT / "artifacts/local/elkies-k3"
SHELL = LOCAL / "q4o323-p0-shell-mod61.json"
POINTING = LOCAL / "q4o323-component2-pointing-qq.json"
MARKING = LOCAL / "q4o323-reflected-fixed-suffix-component2-marking.json"
OUTPUT = LOCAL / "q4o323-p0-shell-anchor-domains-mod61.json"
INPUTS = (SHELL, POINTING, MARKING)
started = time.monotonic()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


shell = json.loads(SHELL.read_text())
pointing = json.loads(POINTING.read_text())
marking = json.loads(MARKING.read_text())
assert shell["status"] == "PASS_MODP_Q4O323_COMPLETE_POLYNOMIAL_P0_SHELL"
assert shell["prime"] == 61
assert pointing["status"] == "PASS_EXACT_QQ_Q4O323_OLD_A11_COMPONENT2_POINTING"
assert marking["status"] == "PASS_EXACT_Q4O323_REFLECTED_FIXED_SUFFIX_MARKING"

prime = ZZ(61)
F = GF(prime)
R = PolynomialRing(F, "u")
u = R.gen()
K = R.fraction_field()


def reduce_qq(value):
    value = QQ(value)
    if ZZ(value.denominator()) % prime == 0:
        raise ZeroDivisionError("p=61 is bad for an exact anchor coefficient")
    return F(value.numerator())/F(value.denominator())


def polynomial(values):
    return R([reduce_qq(value) for value in values])


def record_function(record):
    return K(polynomial(record["numerator_coefficients_low_to_high"])) / K(
        polynomial(record["denominator_coefficients_low_to_high"])
    )


A = polynomial(pointing["global_short_model"]["A_coefficients_low_to_high"])
B = polynomial(pointing["global_short_model"]["B_coefficients_low_to_high"])
E = EllipticCurve(K, [0, 0, 0, A, B])
records = shell["shell"]["records"]
points = [
    E(
        K(R(record["x_coefficients_low_to_high"])),
        K(R(record["y_coefficients_low_to_high"])),
    )
    for record in records
]
anchor_record = pointing["opposite_t0_branch_section"]
anchor = E(record_function(anchor_record["x"]), record_function(anchor_record["y"]))
assert not anchor.is_zero()

supports = [F(value) for value in shell["model"]["finite_I3_supports"]]
nodes = [F(value) for value in shell["model"]["finite_I3_nodes"]]
infinity_node = F(shell["model"]["infinity_I4_node"])


def scaled_infinity_value(function, weight):
    function = K(function)
    numerator = R(function.numerator())
    denominator = R(function.denominator())
    excess = numerator.degree()-denominator.degree()
    if excess < weight:
        return F.zero()
    if excess == weight:
        return numerator[numerator.degree()]/denominator[denominator.degree()]
    return None


def hits_node(point, fibre_index):
    if point.is_zero():
        return False
    x_coordinate, y_coordinate = point.xy()
    if fibre_index < 2:
        support = supports[fibre_index]
        if R(x_coordinate.denominator())(support) == 0:
            return False
        if R(y_coordinate.denominator())(support) == 0:
            return False
        return x_coordinate(support) == nodes[fibre_index] and y_coordinate(support) == 0
    return (
        scaled_infinity_value(x_coordinate, 4) == infinity_node
        and scaled_infinity_value(y_coordinate, 6) == 0
    )


# Any nonidentity component generates Z/3.  For I4 require an odd component,
# detected by the fact that its double remains on the node.
reference_indices = []
references = []
for fibre_index, order in enumerate((3, 3, 4)):
    candidates = [
        index for index, point in enumerate(points)
        if hits_node(point, fibre_index)
        and (order != 4 or hits_node(2*point, fibre_index))
    ]
    assert candidates
    reference_indices.append(candidates[0])
    references.append(points[candidates[0]])
    assert not hits_node(order*references[-1], fibre_index)


def raw_component_profile(point):
    answer = []
    for fibre_index, order in enumerate((3, 3, 4)):
        labels = [
            multiplier for multiplier in range(order)
            if not hits_node(point-multiplier*references[fibre_index], fibre_index)
        ]
        assert len(labels) == 1
        answer.append(labels[0])
    return tuple(answer)  # finite I3, finite I3, infinity I4


def section_intersection(left, right):
    difference = left-right
    if difference.is_zero():
        return -2
    x_coordinate = K(difference[0])
    numerator_degree = R(x_coordinate.numerator()).degree()
    denominator_degree = R(x_coordinate.denominator()).degree()
    infinity_excess = max(0, numerator_degree-denominator_degree-4)
    assert denominator_degree % 2 == 0 and infinity_excess % 2 == 0
    return int(denominator_degree//2+infinity_excess//2)


frame = matrix(ZZ, [
    list(map(int, line.split()))
    for line in (ROOT / marking["frame_output"]).read_text().splitlines()
    if line.strip() and not line.lstrip().startswith("#")
])
gram = block_diagonal_matrix(matrix(ZZ, ((0, 1), (1, 0))), -frame)
fibre = vector(ZZ, [1, 0] + [0]*17)
simple = [
    vector(ZZ, [0, 0] + [-ZZ(index == other) for other in range(17)])
    for index in range(7)
]
highest = (
    (1, 1, 1, 0, 0, 0, 0),
    (0, 0, 0, 1, 1, 0, 0),
    (0, 0, 0, 0, 0, 1, 1),
)
affine = [
    fibre-sum(
        (highest[cycle][index]*simple[index] for index in range(7)),
        vector(ZZ, 19),
    )
    for cycle in range(3)
]
cycles = (
    (affine[0], simple[0], simple[2], simple[1]),
    (affine[1], simple[3], simple[4]),
    (affine[2], simple[5], simple[6]),
)

half_shell = matrix(ZZ, pari(frame).qfminim(4)[2]).transpose().rows()
lattice_sections = []
for raw in half_shell + [-row for row in half_shell]:
    raw = vector(ZZ, raw)
    if raw*frame*raw != 4:
        continue
    section = vector(ZZ, [1, 1] + list(raw))
    pairings = [
        [int(section*gram*component) for component in cycle]
        for cycle in cycles
    ]
    if not all(sorted(row) == [0]*(len(row)-1)+[1] for row in pairings):
        continue
    lattice_sections.append({
        "NS_coordinates": section,
        "component_profile_I4_I3_I3": tuple(row.index(1) for row in pairings),
    })
assert len(lattice_sections) == 258

anchor_class = vector(
    ZZ, marking["equation_explicit_curves_in_child"]["old_A11_component_1"]
)
assert anchor_class*gram*anchor_class == -2
assert anchor_class[0]-anchor_class[1] == 1
anchor_profile = tuple(
    next(index for index, component in enumerate(cycle) if anchor_class*gram*component == 1)
    for cycle in cycles
)
anchor_raw_profile = raw_component_profile(anchor)

# Recover -anchor intrinsically from the marked Mordell--Weil quotient.  A
# P.O=1 section has [a,b]=[2,1] and positive-frame norm six.
zero = vector(ZZ, [-1, 1] + [0]*17)
trivial = matrix(QQ, [fibre, zero] + simple)
assert trivial.rank() == 9
quotient_columns = trivial.right_kernel().basis_matrix().transpose()


def mw_signature(section):
    return vector(QQ, section)*quotient_columns


anchor_signature = mw_signature(anchor_class)
half_norm_six = matrix(ZZ, pari(frame).qfminim(6)[2]).transpose().rows()
inverse_anchor_candidates = []
for raw in half_norm_six + [-row for row in half_norm_six]:
    raw = vector(ZZ, raw)
    if raw*frame*raw != 6:
        continue
    section = vector(ZZ, [2, 1] + list(raw))
    pairings = [
        [int(section*gram*component) for component in cycle]
        for cycle in cycles
    ]
    if not all(sorted(row) == [0]*(len(row)-1)+[1] for row in pairings):
        continue
    if mw_signature(section) == -anchor_signature:
        inverse_anchor_candidates.append(section)
assert len(inverse_anchor_candidates) == 1
inverse_anchor_class = inverse_anchor_candidates[0]
inverse_anchor = -anchor
anchor_inverse_difference_pole_count = section_intersection(anchor, inverse_anchor)
anchor_inverse_lattice_intersection = int(anchor_class*gram*inverse_anchor_class)
inverse_anchor_profile = tuple(
    next(
        index for index, component in enumerate(cycle)
        if inverse_anchor_class*gram*component == 1
    )
    for cycle in cycles
)
inverse_anchor_raw_profile = raw_component_profile(inverse_anchor)


def mapped_profile(raw, finite_swap, signs):
    finite = (raw[1], raw[0]) if finite_swap else (raw[0], raw[1])
    return (
        signs[0]*raw[2] % 4,
        signs[1]*finite[0] % 3,
        signs[2]*finite[1] % 3,
    )


regular_indices = [
    index for index, record in enumerate(records)
    if record["ordinary_coefficient_jacobian_rank"] == 12
]
assert len(regular_indices) == 120
shell_indices = list(range(len(records)))


def direct_affine_intersection_with_anchor(point):
    """Count common affine smooth-fibre points with the exact anchor."""
    x_coordinate, y_coordinate = point.xy()
    anchor_x, anchor_y = anchor.xy()
    x_difference = K(x_coordinate-anchor_x)
    y_difference = K(y_coordinate-anchor_y)
    common = R(x_difference.numerator()).gcd(R(y_difference.numerator()))
    for support in supports:
        while common.degree() >= 1 and common(support) == 0:
            common //= u-support
    return int(common.degree())


def direct_affine_intersection(left, right):
    x_left, y_left = left.xy()
    x_right, y_right = right.xy()
    x_difference = K(x_left-x_right)
    y_difference = K(y_left-y_right)
    common = R(x_difference.numerator()).gcd(R(y_difference.numerator()))
    for support in supports:
        while common.degree() >= 1 and common(support) == 0:
            common //= u-support
    return int(common.degree())


def is_anchor_local_safe(raw_profile):
    # Different components at every reducible fibre imply that the two strict
    # transforms have no omitted local intersection there.  The P.O=0 shell
    # also has no common zero-section intersection with the P.O=1 anchor.
    return all(left != right for left, right in zip(raw_profile, anchor_raw_profile))


def is_dual_anchor_local_safe(raw_profile):
    return (
        is_anchor_local_safe(raw_profile)
        and all(
            left != right
            for left, right in zip(raw_profile, inverse_anchor_raw_profile)
        )
    )


modular_rows = {
    index: {
        "raw_component_profile_I3_I3_I4": raw_component_profile(points[index]),
    }
    for index in shell_indices
}
for index in shell_indices:
    raw_profile = modular_rows[index]["raw_component_profile_I3_I3_I4"]
    modular_rows[index]["safe_for_direct_anchor_intersection"] = is_anchor_local_safe(
        raw_profile
    )
    modular_rows[index]["safe_for_both_direct_anchor_intersections"] = (
        is_dual_anchor_local_safe(raw_profile)
    )
    modular_rows[index]["direct_affine_intersection_with_anchor"] = (
        direct_affine_intersection_with_anchor(points[index])
        if modular_rows[index]["safe_for_direct_anchor_intersection"] else None
    )
    modular_rows[index]["direct_affine_intersection_with_inverse_anchor"] = (
        direct_affine_intersection(points[index], inverse_anchor)
        if modular_rows[index]["safe_for_both_direct_anchor_intersections"] else None
    )
anchored_shell_indices = [
    index for index in shell_indices
    if modular_rows[index]["safe_for_both_direct_anchor_intersections"]
]

orientation_records = []
for finite_swap in (False, True):
    for signs in product((1, -1), repeat=3):
        if mapped_profile(anchor_raw_profile, finite_swap, signs) != anchor_profile:
            continue
        if mapped_profile(
            inverse_anchor_raw_profile, finite_swap, signs
        ) != inverse_anchor_profile:
            continue
        domains = {}
        for index in anchored_shell_indices:
            row = modular_rows[index]
            profile = mapped_profile(
                row["raw_component_profile_I3_I3_I4"], finite_swap, signs
            )
            intersection = row["direct_affine_intersection_with_anchor"]
            inverse_intersection = row[
                "direct_affine_intersection_with_inverse_anchor"
            ]
            domains[index] = [
                lattice_index
                for lattice_index, candidate in enumerate(lattice_sections)
                if candidate["component_profile_I4_I3_I3"] == profile
                and int(candidate["NS_coordinates"]*gram*anchor_class) == intersection
                and int(candidate["NS_coordinates"]*gram*inverse_anchor_class)
                == inverse_intersection
            ]
        histogram = Counter(len(domain) for domain in domains.values())
        singleton_assignments = {
            index: domain[0] for index, domain in domains.items() if len(domain) == 1
        }
        singleton_safe_pair_checks = []
        singleton_safe_pair_mismatches = []
        singleton_indices = sorted(singleton_assignments)
        for left_position, left in enumerate(singleton_indices):
            for right in singleton_indices[left_position+1:]:
                left_profile = modular_rows[left]["raw_component_profile_I3_I3_I4"]
                right_profile = modular_rows[right]["raw_component_profile_I3_I3_I4"]
                if not all(a != b for a, b in zip(left_profile, right_profile)):
                    continue
                modular_intersection = direct_affine_intersection(
                    points[left], points[right]
                )
                lattice_intersection = int(
                    lattice_sections[singleton_assignments[left]]["NS_coordinates"]
                    * gram
                    * lattice_sections[singleton_assignments[right]]["NS_coordinates"]
                )
                check = {
                    "left_shell_index": left,
                    "right_shell_index": right,
                    "modular_intersection": modular_intersection,
                    "lattice_intersection": lattice_intersection,
                }
                singleton_safe_pair_checks.append(check)
                if modular_intersection != lattice_intersection:
                    singleton_safe_pair_mismatches.append(check)
        orientation_records.append({
            "finite_I3_swap": finite_swap,
            "signs_I4_I3_I3": list(signs),
            "domain_size_histogram": {str(key): value for key, value in sorted(histogram.items())},
            "empty_domain_count": histogram.get(0, 0),
            "singleton_domain_count": histogram.get(1, 0),
            "maximum_domain_size": max(map(len, domains.values()), default=0),
            "singleton_assignments": {
                str(index): lattice_index
                for index, lattice_index in singleton_assignments.items()
            },
            "singleton_safe_pair_check_count": len(singleton_safe_pair_checks),
            "singleton_safe_pair_mismatch_count": len(singleton_safe_pair_mismatches),
            "singleton_safe_pair_mismatches": singleton_safe_pair_mismatches,
            "domains_by_shell_index": {
                str(index): domain for index, domain in domains.items()
            },
        })
assert orientation_records

payload = {
    "schema": "elkies-k3.h92-q4o323-p0-shell-anchor-domains-mod61.v1",
    "status": "PASS_MOD61_Q4O323_REGULAR_P0_SHELL_ANCHOR_DOMAINS",
    "prime": 61,
    "exact_anchor": {
        "label": "old_A11_component_1_opposite_to_component2_zero",
        "NS_coordinates": list(map(int, anchor_class)),
        "P_dot_O": 1,
        "raw_component_profile_I3_I3_I4": list(anchor_raw_profile),
        "marked_component_profile_I4_I3_I3": list(anchor_profile),
        "inverse": {
            "NS_coordinates": list(map(int, inverse_anchor_class)),
            "raw_component_profile_I3_I3_I4": list(inverse_anchor_raw_profile),
            "marked_component_profile_I4_I3_I3": list(inverse_anchor_profile),
            "intersection_with_anchor": int(anchor_class*gram*inverse_anchor_class),
            "difference_section_naive_pole_count": anchor_inverse_difference_pole_count,
            "note": (
                "The naive pole count of P-(-P) is not P.(-P) on the resolved "
                "surface; the missing contribution is at the common zero-section pole."
            ),
            "exact_MW_signature_is_negative": True,
        },
    },
    "component_references": {
        "shell_indices_finite_I3_finite_I3_infinity_I4": reference_indices,
        "orders": [3, 3, 4],
    },
    "shell": {
        "signed_count": len(records),
        "ordinary_rank12_indices": regular_indices,
        "ordinary_rank12_count": len(regular_indices),
        "locally_safe_shell_indices": anchored_shell_indices,
        "locally_safe_shell_count": len(anchored_shell_indices),
        "locally_safe_rank12_indices": [
            index for index in anchored_shell_indices if index in set(regular_indices)
        ],
        "locally_safe_rank12_count": sum(
            index in set(regular_indices) for index in anchored_shell_indices
        ),
        "regular_point_anchor_data": {
            str(index): {
                "raw_component_profile_I3_I3_I4": list(
                    modular_rows[index]["raw_component_profile_I3_I3_I4"]
                ),
                "safe_for_direct_anchor_intersection": modular_rows[index][
                    "safe_for_direct_anchor_intersection"
                ],
                "safe_for_both_direct_anchor_intersections": modular_rows[index][
                    "safe_for_both_direct_anchor_intersections"
                ],
                "direct_affine_intersection_with_anchor": modular_rows[index][
                    "direct_affine_intersection_with_anchor"
                ],
                "direct_affine_intersection_with_inverse_anchor": modular_rows[index][
                    "direct_affine_intersection_with_inverse_anchor"
                ],
            }
            for index in shell_indices
        },
    },
    "lattice": {
        "P_dot_O_zero_section_count": len(lattice_sections),
        "sections": [
            {
                "NS_coordinates": list(map(int, row["NS_coordinates"])),
                "component_profile_I4_I3_I3": list(row["component_profile_I4_I3_I3"]),
            }
            for row in lattice_sections
        ],
    },
    "orientation_domains": orientation_records,
    "method": {
        "large_Groebner_required": False,
        "elimination_required": False,
        "runtime_seconds": time.monotonic()-started,
    },
    "proof_boundary": (
        "Every component and anchor-intersection datum is exact over GF(61), and every "
        "candidate lattice class is an exact nef norm-four section. Pairwise graph matching, "
        "target-coset naming, characteristic-zero lifting, and the q12 equation remain open."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in INPUTS},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "Q4O323P0ANCHOR|regular={}|orientations={}|empty={}|singletons={}|max_domains={}|"
    "status={}|output={}".format(
        len(regular_indices), len(orientation_records),
        [row["empty_domain_count"] for row in orientation_records],
        [row["singleton_domain_count"] for row in orientation_records],
        [row["maximum_domain_size"] for row in orientation_records],
        payload["status"], OUTPUT,
    ), flush=True,
)
