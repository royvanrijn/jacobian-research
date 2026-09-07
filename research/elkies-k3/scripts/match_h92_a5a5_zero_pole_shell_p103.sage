#!/usr/bin/env sage -python
"""Match the exact 54-class PO=0 shell inside the p=103 section pool.

status: EXPERIMENT
claim: exhaustive marked weighted-graph embedding over the bounded pools
inputs: exact 2A5 marking/frame and modular polynomial-section enumeration
output: artifacts/local/elkies-k3/q24-2a5-zero-pole-shell-match-p103.json

The exact shell consists of every nef class [1,1,w] with w.Frame.w=4.
Vertices are marked by both oriented I6 components and their intersection
with the known old-A11 affine section.  Distinct section intersections give
the edge colours.  The search embeds this complete marked graph into the 130
modular polynomial sections.  It is a finite marking calculation, not a QQ
lift.
"""

import hashlib
import json
import time
from collections import Counter
from pathlib import Path

from sage.all import (
    ZZ, QQ, EllipticCurve, GF, PolynomialRing, block_matrix, matrix, pari,
    vector, zero_matrix,
)


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
LOCAL = ROOT / "artifacts/local/elkies-k3"
FRAME_PATH = GENERATED / "elkies-k3-h3-a11-q8-orbit12-explicit-zero-frame.txt"
MARKING_PATH = LOCAL / "q24-a11-to-2a5-q8-equation-marking-qq.json"
MODULAR_PATH = LOCAL / "q24-2a5-zero-pole-sections-p103.json"
OUTPUT = LOCAL / "q24-2a5-zero-pole-shell-match-p103.json"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


frame = matrix(ZZ, [
    list(map(int, line.split()))
    for line in FRAME_PATH.read_text().splitlines()
    if line.strip() and not line.lstrip().startswith("#")
])
assert frame.nrows() == frame.ncols() == 17
gram = block_matrix(ZZ, [
    [matrix(ZZ, [[0, 1], [1, 0]]), zero_matrix(ZZ, 2, 17)],
    [zero_matrix(ZZ, 17, 2), -frame],
])

marking = json.loads(MARKING_PATH.read_text())
modular = json.loads(MODULAR_PATH.read_text())
assert marking["status"] == "PASS_EXACT_Q24_A11_Q8_2A5_EQUATION_MARKING"
assert modular["status"] == "PASS_BOUNDED_MOD103_ZERO_POLE_SECTION_ENUMERATION"

old_fibre = vector(ZZ, [1, 0] + [0] * 17)
old_components = [
    vector(ZZ, marking["physical_2A5"]["child_coordinates"][f"old_A11_component_{index}"])
    for index in range(11)
]
chains = marking["physical_2A5"]["chains"]
root_components = [old_components[index] for chain in chains for index in chain]
affine_components = [
    old_fibre - sum((old_components[index] for index in chain), vector(ZZ, 19))
    for chain in chains
]
fibre_components = root_components + affine_components
known_affine = vector(
    ZZ,
    marking["old_A11_affine_section_on_component9_pointed_child"][
        "NS_coordinates_in_selected_child_basis"
    ],
)
target = vector(ZZ, [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 1])

# qfminim returns one representative from each +/- pair by default.
half_shell = matrix(ZZ, pari(frame).qfminim(4)[2]).transpose().rows()
norm_four = half_shell + [-row for row in half_shell]
exact_sections = []
for row in norm_four:
    if row * frame * row != 4:
        continue
    section = vector(ZZ, [1, 1] + list(row))
    pairings = tuple(int(section * gram * component) for component in fibre_components)
    if min(pairings) < 0:
        continue
    assert sum(pairings[:5]) + pairings[10] == 1
    assert sum(pairings[5:10]) + pairings[11] == 1
    # Oriented cycle labels: T=68 is the second chain, T=89 the first.
    component_68 = next(
        component for component, position in ((4, 5), (3, 6), (2, 7), (1, 8), (0, 9), (5, 11))
        if pairings[position] == 1
    )
    component_89 = next(
        component for component, position in ((4, 0), (3, 1), (2, 2), (1, 3), (0, 4), (5, 10))
        if pairings[position] == 1
    )
    exact_sections.append({
        "NS": section,
        "profile": (component_68, component_89),
        "affine_intersection": int(section * gram * known_affine),
    })
assert len(exact_sections) == 54
exact_affine = next(index for index, row in enumerate(exact_sections) if row["NS"] == known_affine)
exact_target = next(index for index, row in enumerate(exact_sections) if row["NS"] == target)
assert exact_sections[exact_target]["profile"] == (0, 0)
assert exact_sections[exact_target]["affine_intersection"] == 3
assert sum(
    row["profile"] == (0, 0) and row["affine_intersection"] == 3
    for row in exact_sections
) == 1

# Quotient signatures modulo the rational trivial lattice.  There is no
# torsion ambiguity in this pointed 2A5 frame, so signature addition is the
# exact Mordell--Weil group law on these classes.
trivial = matrix(QQ, [old_fibre, old_components[9]] + root_components)
assert trivial.rank() == 12
quotient_columns = trivial.right_kernel().basis_matrix().transpose()


def mw_signature(section):
    return tuple(vector(QQ, section) * quotient_columns)


exact_signatures = [mw_signature(row["NS"]) for row in exact_sections]
assert len(set(exact_signatures)) == len(exact_signatures)
signature_to_exact = {signature: index for index, signature in enumerate(exact_signatures)}
exact_sum = [[None] * len(exact_sections) for unused in exact_sections]
sum_predecessors = [[] for unused in exact_sections]
for left in range(len(exact_sections)):
    for right in range(len(exact_sections)):
        signature = tuple(
            exact_signatures[left][index] + exact_signatures[right][index]
            for index in range(len(exact_signatures[left]))
        )
        result = signature_to_exact.get(signature)
        exact_sum[left][right] = result
        if result is not None:
            sum_predecessors[result].append((left, right))
exact_addition_relations = sum(result is not None for row in exact_sum for result in row)
assert exact_addition_relations > 0

F = GF(103)
RT = PolynomialRing(F, "T")


modular_sections = []
for row in modular["sections"]:
    point = (
        RT(row["X_coefficients_low_to_high"]),
        RT(row["Y_coefficients_low_to_high"]),
    )
    profile = tuple(
        int(incidence["oriented_component_in_I6_cycle"])
        for incidence in row["I6_incidence"]
    )
    modular_sections.append({
        "point": point,
        "profile": profile,
        "affine_intersection": row["intersection_with_known_old_A11_affine_section"],
    })
assert len(modular_sections) == 130
modular_affine = int(modular["known_affine_section_index"])


def inverse_pair_intersection(profile):
    """Intersection P.(-P) from Shioda's formula for the two I6 fibres."""
    height = QQ(4)
    cross_correction = QQ(0)
    for component in profile:
        if component == 0:
            continue
        height -= QQ(component * (6 - component)) / 6
        inverse_component = 6 - component
        cross_correction += QQ(
            min(component, inverse_component) * (6 - max(component, inverse_component))
        ) / 6
    answer = QQ(2) + height - cross_correction
    assert answer in ZZ
    return int(answer)

K = RT.fraction_field()
E = EllipticCurve(K, [
    0, 0, 0,
    K(RT(modular["surface_mod_103"]["A_coefficients_low_to_high"])),
    K(RT(modular["surface_mod_103"]["B_coefficients_low_to_high"])),
])
modular_points = [E(K(row["point"][0]), K(row["point"][1])) for row in modular_sections]


def padded(poly, length):
    return tuple(int(poly[index]) if index <= poly.degree() else 0 for index in range(length))


modular_key_to_index = {
    (padded(row["point"][0], 5), padded(row["point"][1], 7)): index
    for index, row in enumerate(modular_sections)
}
assert len(modular_key_to_index) == len(modular_sections)
modular_sum_cache = {}


def modular_sum(left, right):
    key = (min(left, right), max(left, right))
    if key in modular_sum_cache:
        return modular_sum_cache[key]
    result = modular_points[left] + modular_points[right]
    answer = None
    if not result.is_zero():
        x_value, y_value = result.xy()
        if RT(x_value.denominator()).degree() <= 0 and RT(y_value.denominator()).degree() <= 0:
            x_poly = RT(x_value.numerator()) / F(RT(x_value.denominator())[0])
            y_poly = RT(y_value.numerator()) / F(RT(y_value.denominator())[0])
            if x_poly.degree() <= 4 and y_poly.degree() <= 6:
                answer = modular_key_to_index.get((padded(x_poly, 5), padded(y_poly, 7)))
    modular_sum_cache[key] = answer
    return answer

domains = {}
for exact_index, row in enumerate(exact_sections):
    if exact_index == exact_affine:
        domains[exact_index] = [modular_affine]
        continue
    domains[exact_index] = [
        modular_index
        for modular_index, candidate in enumerate(modular_sections)
        if candidate["profile"] == row["profile"]
        and (
            candidate["affine_intersection"]
            if candidate["affine_intersection"] is not None
            else inverse_pair_intersection(candidate["profile"])
        ) == row["affine_intersection"]
    ]
    if not domains[exact_index]:
        available = Counter(
            candidate["affine_intersection"]
            for candidate in modular_sections
            if candidate["profile"] == row["profile"]
        )
        raise ArithmeticError(
            f"empty domain exact={exact_index} profile={row['profile']} "
            f"affine_intersection={row['affine_intersection']} available={available}"
        )
assert set(domains[exact_target]) == set(modular["P1229_mod103_candidates"])

started = time.monotonic()
mapping = {exact_affine: modular_affine}
used = {modular_affine}
solution_count = 0
target_images = Counter()
complete_mappings = []


def compatible(exact_index, modular_index):
    # Relations exact_index + assigned = result.
    for exact_other, modular_other in mapping.items():
        exact_result = exact_sum[exact_index][exact_other]
        if exact_result is None:
            continue
        observed = modular_sum(modular_index, modular_other)
        if observed is None:
            return False
        if exact_result in mapping:
            if observed != mapping[exact_result]:
                return False
        elif observed not in domains[exact_result] or observed in used:
            return False

    # Relations assigned + assigned = exact_index.
    for exact_left, exact_right in sum_predecessors[exact_index]:
        if exact_left in mapping and exact_right in mapping:
            if modular_sum(mapping[exact_left], mapping[exact_right]) != modular_index:
                return False
    return True


def search():
    global solution_count
    if len(mapping) == len(exact_sections):
        solution_count += 1
        target_images[mapping[exact_target]] += 1
        complete_mappings.append(tuple(mapping[index] for index in range(len(exact_sections))))
        return

    best_exact = None
    best_options = None
    for exact_index in range(len(exact_sections)):
        if exact_index in mapping:
            continue
        options = [
            modular_index for modular_index in domains[exact_index]
            if modular_index not in used and compatible(exact_index, modular_index)
        ]
        if not options:
            return
        if best_options is None or len(options) < len(best_options):
            best_exact = exact_index
            best_options = options
            if len(options) == 1:
                break

    for modular_index in best_options:
        mapping[best_exact] = modular_index
        used.add(modular_index)
        search()
        used.remove(modular_index)
        del mapping[best_exact]


search()
assert solution_count > 0
print(
    f"A5A5SHELLMATCH_DIAG|embeddings={solution_count}|target_images={dict(target_images)}",
    flush=True,
)
assert solution_count == len(complete_mappings) == 36
assert set(target_images) == set(modular["P1229_mod103_candidates"])
assert set(target_images.values()) == {6}
canonical_mapping = min(complete_mappings)
selected_target = canonical_mapping[exact_target]
assert selected_target in modular["P1229_mod103_candidates"]

payload = {
    "schema": "elkies-k3.q24-2a5-zero-pole-shell-match-p103.v1",
    "status": "PASS_EXHAUSTIVE_MOD103_ZERO_POLE_SHELL_EMBEDDINGS_CANONICAL_MARKING",
    "software": "SageMath 10.9 (conda-forge pinned repository environment)",
    "exact_shell": {
        "norm_four_signed_vectors": len(norm_four),
        "nef_zero_pole_classes": len(exact_sections),
        "indexed_component_nef_classes": [
            {
                "exact_index": index,
                "NS_coordinates": list(map(int, row["NS"])),
                "oriented_I6_components_at_T68_T89": list(row["profile"]),
                "intersection_with_known_affine": row["affine_intersection"],
            }
            for index, row in enumerate(exact_sections)
        ],
        "known_affine_exact_index": exact_affine,
        "P1229_exact_index": exact_target,
        "P1229_unique_exact_fingerprint": {
            "oriented_I6_components_at_T68_T89": [0, 0],
            "intersection_with_known_affine": 3,
        },
    },
    "modular_pool": {
        "sections": len(modular_sections),
        "known_affine_index": modular_affine,
        "coarse_P1229_candidates": modular["P1229_mod103_candidates"],
    },
    "embedding_search": {
        "complete_marked_shell_embeddings": solution_count,
        "P1229_image_multiplicities": dict(sorted(target_images.items())),
        "exact_MW_addition_relations": exact_addition_relations,
        "edge_invariant": "Mordell-Weil addition commuting with reduction",
        "residual_marking_stabilizer_order": 36,
        "intrinsic_P1229_image_is_unique": False,
        "canonical_choice": "lexicographically least complete exact-index-to-modular-index mapping",
        "elapsed_seconds": round(time.monotonic() - started, 6),
    },
    "selected_P1229_mod103": {
        "section_index": selected_target,
        "X_coefficients_low_to_high": modular["sections"][selected_target]["X_coefficients_low_to_high"],
        "Y_coefficients_low_to_high": modular["sections"][selected_target]["Y_coefficients_low_to_high"],
    },
    "canonical_complete_mapping_exact_index_to_modular_index": {
        str(index): canonical_mapping[index] for index in range(len(exact_sections))
    },
    "all_complete_mappings_exact_index_to_modular_index": [
        {str(index): complete[index] for index in range(len(exact_sections))}
        for complete in sorted(complete_mappings)
    ],
    "proof_boundary": (
        "The exhaustive finite marked-addition search leaves 36 shell embeddings and six "
        "equally possible images of P1229. The displayed image is a deterministic compiler "
        "marking convention, not an intrinsic identification from the currently explicit "
        "equation anchors. It does not lift the selected modular section to QQ."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in (FRAME_PATH, MARKING_PATH, MODULAR_PATH)],
        "sha256": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in (FRAME_PATH, MARKING_PATH, MODULAR_PATH)
        },
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "A5A5SHELLMATCH|"
    f"exact=54|modular=130|embeddings={solution_count}|"
    f"P1229_index={selected_target}|X={tuple(payload['selected_P1229_mod103']['X_coefficients_low_to_high'])}|"
    f"Y={tuple(payload['selected_P1229_mod103']['Y_coefficients_low_to_high'])}|"
    f"status={payload['status']}|output={OUTPUT}",
    flush=True,
)
