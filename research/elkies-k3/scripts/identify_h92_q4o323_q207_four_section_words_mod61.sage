#!/usr/bin/env sage -python
"""Identify the reflected q207 horizontal from four q323 P.O=0 sections.

status: ACTIVE_SEARCH
claim: exhaustive equation-side test of every minimal four-section lattice word
inputs: q323 p=61 shell, component-2 pointing, reflected fixed-suffix marking
outputs: artifacts/local/elkies-k3/q4o323-q207-four-section-words-mod61.json

The marked lattice has no one-, two-, or three-section P.O=0 expression for
the q207 MW class and exactly three unordered four-section expressions.  For
each component orientation, this script matches their component profiles and
directional inverse-pencil degrees against the complete p=61 shell.  Degree
drops under reduction are allowed, then the full sum is required to have the
exact parent degree 16 and target component profile.  No Groebner basis or
surface elimination is used.
"""

import hashlib
import json
import time
from collections import defaultdict
from itertools import product
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
if not (ROOT / "MATH_STATUS.json").exists():
    ROOT = Path.cwd()
LOCAL = ROOT / "artifacts/local/elkies-k3"
SHELL = LOCAL / "q4o323-p0-shell-mod61.json"
ANCHOR = LOCAL / "q4o323-p0-shell-anchor-domains-mod61.json"
POINTING = LOCAL / "q4o323-component2-pointing-qq.json"
MARKING = LOCAL / "q4o323-reflected-fixed-suffix-component2-marking.json"
OUTPUT = LOCAL / "q4o323-q207-four-section-words-mod61.json"
INPUTS = (SHELL, ANCHOR, POINTING, MARKING)
started = time.monotonic()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


shell = json.loads(SHELL.read_text())
anchor_data = json.loads(ANCHOR.read_text())
pointing = json.loads(POINTING.read_text())
marking = json.loads(MARKING.read_text())
assert shell["status"] == "PASS_MODP_Q4O323_COMPLETE_POLYNOMIAL_P0_SHELL"
assert anchor_data["status"] == "PASS_MOD61_Q4O323_REGULAR_P0_SHELL_ANCHOR_DOMAINS"
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
        raise ZeroDivisionError("p=61 is bad for a pointing coefficient")
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
shell_records = shell["shell"]["records"]
points = [
    E(
        K(R(record["x_coefficients_low_to_high"])),
        K(R(record["y_coefficients_low_to_high"])),
    )
    for record in shell_records
]
raw_profiles = [
    tuple(
        anchor_data["shell"]["regular_point_anchor_data"][str(index)][
            "raw_component_profile_I3_I3_I4"
        ]
    )
    for index in range(len(points))
]

a1, a2, a3, unused_a4, unused_a6 = [
    record_function(pointing["pointed_generalized_weierstrass"][name])
    for name in ("a1", "a2", "a3", "a4", "a6")
]
b2 = a1**2+4*a2
w0 = K(polynomial(pointing["quartic_square_at_t0"]["L0_coefficients_low_to_high"]))


def inverse_parent_base(point):
    if point.is_zero():
        return None
    X, Y = point.xy()
    x_general = K(X)/9-b2/12
    y_general = K(Y)/27-(a1*x_general+a3)/2
    if not y_general:
        return None
    return K(2*w0*(x_general+a2)/y_general)


def rational_degree(value):
    if value is None:
        return None
    return int(max(R(value.numerator()).degree(), R(value.denominator()).degree()))


modular_directional_degrees = [
    (rational_degree(inverse_parent_base(point)), rational_degree(inverse_parent_base(-point)))
    for point in points
]

lattice_rows = anchor_data["lattice"]["sections"]
lattice_sections = [vector(ZZ, row["NS_coordinates"]) for row in lattice_rows]
lattice_profiles = [tuple(row["component_profile_I4_I3_I3"]) for row in lattice_rows]
tail_to_index = {
    tuple(section[-10:]): index for index, section in enumerate(lattice_sections)
}
assert len(tail_to_index) == len(lattice_sections)
inverse_indices = [
    tail_to_index[tuple(-value for value in section[-10:])]
    for section in lattice_sections
]
child_in_source = matrix(ZZ, marking["basis_in_source"])
lattice_parent_degrees = [
    int((section*child_in_source)[1]) for section in lattice_sections
]

target = vector(
    ZZ,
    marking["fixed_suffix_transport"]["q207_component_reduction"][
        "equation_preflight"
    ]["horizontal_section"],
)
target_tail = tuple(target[-10:])
target_parent_degree = int(
    marking["fixed_suffix_transport"]["q207_component_reduction"][
        "equation_preflight"
    ]["q4o208_parent_degree"]
)
assert target_parent_degree == 16

# Recover the I4+A2+A2 target component profile directly from the standard
# child basis ordering used by the reflected-suffix certificate.
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
frame = matrix(ZZ, [
    list(map(int, line.split()))
    for line in (ROOT / marking["frame_output"]).read_text().splitlines()
    if line.strip() and not line.lstrip().startswith("#")
])
from sage.all import block_diagonal_matrix
gram = block_diagonal_matrix(matrix(ZZ, ((0, 1), (1, 0))), -frame)
target_profile = tuple(
    next(
        index for index, component in enumerate(cycle)
        if target*gram*component == 1
    )
    for cycle in cycles
)

# Exhaust all unordered one/two/three/four-section MW-tail words.
tail_indices = {tuple(section[-10:]): index for index, section in enumerate(lattice_sections)}
assert target_tail not in tail_indices
assert not any(
    tuple(target_tail[k]-left[-10+k] for k in range(10)) in tail_indices
    for left in lattice_sections
)
pair_sums = defaultdict(list)
for left in range(len(lattice_sections)):
    for right in range(left, len(lattice_sections)):
        pair_sums[tuple(
            lattice_sections[left][-10+k]+lattice_sections[right][-10+k]
            for k in range(10)
        )].append((left, right))
assert not any(
    tuple(target_tail[k]-section[-10+k] for k in range(10)) in pair_sums
    for section in lattice_sections
)
four_words = set()
for pair_tail, left_pairs in pair_sums.items():
    complement = tuple(target_tail[k]-pair_tail[k] for k in range(10))
    for left_pair in left_pairs:
        for right_pair in pair_sums.get(complement, ()):
            four_words.add(tuple(sorted(left_pair+right_pair)))
four_words = sorted(four_words)
assert four_words == [
    (3, 25, 63, 224),
    (4, 25, 60, 224),
    (4, 25, 208, 247),
]


def mapped_profile(raw, finite_swap, signs):
    finite = (raw[1], raw[0]) if finite_swap else (raw[0], raw[1])
    return (
        signs[0]*raw[2] % 4,
        signs[1]*finite[0] % 3,
        signs[2]*finite[1] % 3,
    )


def point_key(point):
    X, Y = point.xy()
    return (
        tuple(R(X.numerator()).list()), tuple(R(X.denominator()).list()),
        tuple(R(Y.numerator()).list()), tuple(R(Y.denominator()).list()),
    )


def point_shape(point):
    X, Y = point.xy()
    return [
        int(R(X.numerator()).degree()), int(R(X.denominator()).degree()),
        int(R(Y.numerator()).degree()), int(R(Y.denominator()).degree()),
    ]


orientation_records = []
hits_by_key = {}
tested_combinations = 0
for finite_swap in (False, True):
    for signs in product((1, -1), repeat=3):
        word_records = []
        for word in four_words:
            pools = []
            for lattice_index in word:
                forward_bound = lattice_parent_degrees[lattice_index]
                reverse_bound = lattice_parent_degrees[inverse_indices[lattice_index]]
                pools.append([
                    index for index in range(len(points))
                    if mapped_profile(raw_profiles[index], finite_swap, signs)
                    == lattice_profiles[lattice_index]
                    and modular_directional_degrees[index][0] is not None
                    and modular_directional_degrees[index][1] is not None
                    and modular_directional_degrees[index][0] <= forward_bound
                    and modular_directional_degrees[index][1] <= reverse_bound
                ])
            tested_here = 0
            hit_keys = []
            if all(pools):
                for indices in product(*pools):
                    if len(set(indices)) != len(indices):
                        continue
                    tested_here += 1
                    tested_combinations += 1
                    point = sum((points[index] for index in indices), E(0))
                    parent_degree = rational_degree(inverse_parent_base(point))
                    if parent_degree != target_parent_degree:
                        continue
                    raw_sum_profile = tuple(
                        sum(raw_profiles[index][component] for index in indices)
                        % (3, 3, 4)[component]
                        for component in range(3)
                    )
                    if mapped_profile(raw_sum_profile, finite_swap, signs) != target_profile:
                        continue
                    key = point_key(point)
                    if key not in hits_by_key:
                        X, Y = point.xy()
                        hits_by_key[key] = {
                            "point": point,
                            "records": [],
                            "shape_Xnum_Xden_Ynum_Yden": point_shape(point),
                            "inverse_parent_degree": parent_degree,
                            "raw_component_profile_I3_I3_I4": list(raw_sum_profile),
                            "mapped_component_profile_I4_I3_I3": list(target_profile),
                            "exact_mod61_weierstrass_identity": bool(
                                Y**2 == X**3+K(A)*X+K(B)
                            ),
                        }
                    hits_by_key[key]["records"].append({
                        "finite_I3_swap": finite_swap,
                        "signs_I4_I3_I3": list(signs),
                        "lattice_word": list(word),
                        "shell_indices": list(indices),
                    })
                    hit_keys.append(key)
            word_records.append({
                "lattice_word": list(word),
                "pool_sizes": list(map(len, pools)),
                "tested_distinct_combinations": tested_here,
                "hit_count_with_multiplicity": len(hit_keys),
            })
        orientation_records.append({
            "finite_I3_swap": finite_swap,
            "signs_I4_I3_I3": list(signs),
            "words": word_records,
        })

hits = []
for row in hits_by_key.values():
    point = row.pop("point")
    X, Y = point.xy()
    row["x"] = {
        "numerator_coefficients_low_to_high": list(map(int, R(X.numerator()).list())),
        "denominator_coefficients_low_to_high": list(map(int, R(X.denominator()).list())),
    }
    row["y"] = {
        "numerator_coefficients_low_to_high": list(map(int, R(Y.numerator()).list())),
        "denominator_coefficients_low_to_high": list(map(int, R(Y.denominator()).list())),
    }
    hits.append(row)

payload = {
    "schema": "elkies-k3.h92-q4o323-q207-four-section-words-mod61.v1",
    "status": (
        "PASS_MOD61_Q4O323_Q207_FOUR_SECTION_CANDIDATES"
        if hits else "REJECTED_MOD61_Q4O323_Q207_NO_FOUR_SECTION_CANDIDATE"
    ),
    "prime": 61,
    "target": {
        "MW_tail": list(map(int, target_tail)),
        "component_profile_I4_I3_I3": list(target_profile),
        "inverse_parent_degree": target_parent_degree,
        "minimal_new_P0_section_count": 4,
        "all_unordered_minimal_lattice_words": [list(word) for word in four_words],
    },
    "search": {
        "degree_drop_under_reduction_allowed": True,
        "orientation_records": orientation_records,
        "tested_distinct_combinations": tested_combinations,
        "unique_candidate_points": len(hits),
        "candidates": hits,
    },
    "method": {
        "large_Groebner_required": False,
        "elimination_required": False,
        "runtime_seconds": time.monotonic()-started,
    },
    "proof_boundary": (
        "This exhausts the three minimal marked four-section words inside the complete "
        "p=61 polynomial shell, with component and inverse-parent-degree gates. A modular "
        "candidate still requires marking uniqueness and a characteristic-zero lift before "
        "it proves the physical q207 horizontal or its q12 equation."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in INPUTS},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "Q4O323Q207WORD61|lattice_words={}|tests={}|candidates={}|status={}|output={}".format(
        len(four_words), tested_combinations, len(hits), payload["status"], OUTPUT
    ), flush=True,
)
