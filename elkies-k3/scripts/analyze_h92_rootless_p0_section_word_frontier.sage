#!/usr/bin/env sage -python
"""Compare short P.O=0 section words for every retained rootless exit.

status: ACTIVE_COMPILER
claim: exhaustive length-at-most-four P.O=0 word search modulo explicit sections
inputs: physical q8/o376 4A1 marking and rootless equation-cost artifacts
outputs: generated rootless P.O=0 section-word frontier
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
U2 = matrix(ZZ, ((0, 1), (1, 0)))
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--marking", type=Path, required=True)
parser.add_argument("--cost", type=Path, action="append", required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()


def load_matrix(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


def entries(value):
    return [int(item) for item in vector(ZZ, value)]


marking_path = args.marking.resolve()
cost_paths = [path.resolve() for path in args.cost]
output = args.output.resolve()
marking = json.loads(marking_path.read_text())
frame_path = ROOT / marking["frame_output"]
frame = load_matrix(frame_path)
root_rank = int(marking["root_data"][0])
assert root_rank == 4
gram = block_diagonal_matrix(U2, -frame)
basis_in_parent = matrix(ZZ, marking["basis_in_source"])
parent_in_basis = matrix(ZZ, marking["source_in_basis"])
parent_fibre = vector(ZZ, [1, 0] + [0] * 17) * parent_in_basis

explicit = {
    name: vector(ZZ, value)
    for name, value in marking["equation_explicit_curves_in_child"].items()
}
known_sections = [
    (name, vector(ZZ, curve[-13:]))
    for name, curve in explicit.items()
    if curve[1] == 1 and any(curve[-13:])
]
known_matrix = matrix(ZZ, [list(value) for _, value in known_sections])
smith, _, right = known_matrix.smith_form()
known_rank = int(known_matrix.rank())
diagonal = [abs(ZZ(smith[index, index])) for index in range(known_rank)]


def quotient_key(value):
    transformed = vector(ZZ, value) * right
    return tuple(
        [int(transformed[index] % diagonal[index]) for index in range(known_rank)]
        + [int(transformed[index]) for index in range(known_rank, 13)]
    )


def quotient_subtract(left, right_key):
    return (
        tuple(
            (left[index] - right_key[index]) % diagonal[index]
            for index in range(known_rank)
        )
        + tuple(left[index] - right_key[index] for index in range(known_rank, 13))
    )


shell = pari(frame).qfminim(4)
half = [vector(ZZ, column) for column in matrix(ZZ, shell[2]).columns()]
best_by_mw = {}
best_by_mw_parent_pole = {}
physical_p0_class_count = 0
for tail in half + [-value for value in half]:
    if tail * frame * tail != 4:
        continue
    component_pairings = vector(ZZ, tail) * frame[:, :root_rank]
    affine_pairings = vector(ZZ, [1] * root_rank) - component_pairings
    if min(tuple(component_pairings) + tuple(affine_pairings)) < 0:
        continue
    physical_p0_class_count += 1
    section = vector(ZZ, [1, 1] + list(tail))
    parent_degree = int(section * gram * parent_fibre)
    assert parent_degree >= 0
    mw = vector(ZZ, tail[-13:])
    ordering = (
        parent_degree,
        max(abs(int(item)) for item in tail),
        tuple(tail),
    )
    key = tuple(mw)
    parent_curve = section * basis_in_parent
    parent_pole_proxy = int(parent_curve[0] - parent_curve[1])
    candidate = {
        "current_4A1_section": section,
        "current_4A1_P_dot_O": 0,
        "mw": mw,
        "q4o164_parent_curve": parent_curve,
        "q4o164_parent_degree": parent_degree,
        "q4o164_parent_a_minus_b": parent_pole_proxy,
        "current_4A1_component_pairings": component_pairings,
        "current_4A1_affine_pairings": affine_pairings,
    }
    if key not in best_by_mw or ordering < best_by_mw[key][0]:
        best_by_mw[key] = (ordering, candidate)
    pole_ordering = (
        parent_pole_proxy,
        parent_degree,
        max(abs(int(item)) for item in tail),
        tuple(tail),
    )
    if key not in best_by_mw_parent_pole or pole_ordering < best_by_mw_parent_pole[key][0]:
        best_by_mw_parent_pole[key] = (pole_ordering, candidate)

candidates = [item[1] for item in best_by_mw.values()]
candidate_keys = [quotient_key(item["mw"]) for item in candidates]
parent_pole_candidates = [item[1] for item in best_by_mw_parent_pole.values()]

# Also retain the complete P.O=1 shell.  It is small enough to test the two
# mixed patterns most likely to beat four polynomial branches: two sections
# with pole split (0,1) or (1,1), and three sections with split (0,0,1).
shell6 = pari(frame).qfminim(6)
half6 = [vector(ZZ, column) for column in matrix(ZZ, shell6[2]).columns()]
best_p1_by_mw = {}
physical_p1_class_count = 0
for tail in half6 + [-value for value in half6]:
    if tail * frame * tail != 6:
        continue
    component_pairings = vector(ZZ, tail) * frame[:, :root_rank]
    affine_pairings = vector(ZZ, [1] * root_rank) - component_pairings
    if min(tuple(component_pairings) + tuple(affine_pairings)) < 0:
        continue
    physical_p1_class_count += 1
    section = vector(ZZ, [2, 1] + list(tail))
    parent_degree = int(section * gram * parent_fibre)
    assert parent_degree >= 0
    mw = vector(ZZ, tail[-13:])
    ordering = (
        parent_degree,
        max(abs(int(item)) for item in tail),
        tuple(tail),
    )
    key = tuple(mw)
    candidate = {
        "current_4A1_section": section,
        "current_4A1_P_dot_O": 1,
        "mw": mw,
        "q4o164_parent_curve": section * basis_in_parent,
        "q4o164_parent_degree": parent_degree,
        "q4o164_parent_a_minus_b": int((section * basis_in_parent)[0] - parent_degree),
        "current_4A1_component_pairings": component_pairings,
        "current_4A1_affine_pairings": affine_pairings,
    }
    if key not in best_p1_by_mw or ordering < best_p1_by_mw[key][0]:
        best_p1_by_mw[key] = (ordering, candidate)
p1_candidates = [item[1] for item in best_p1_by_mw.values()]
p1_by_quotient = {}
for index, item in enumerate(p1_candidates):
    key = quotient_key(item["mw"])
    ordering = (
        item["q4o164_parent_degree"],
        max(abs(int(value)) for value in item["mw"]),
        index,
    )
    if key not in p1_by_quotient or ordering < p1_by_quotient[key][0]:
        p1_by_quotient[key] = (ordering, item)

# For each quotient class retain the pair with minimum parent-degree sum, then
# minimum maximum degree.  This is sufficient for both the length-three and
# meet-in-the-middle length-four optima under the same lexicographic cost.
pairs = {}
for left_index, left in enumerate(candidates):
    for right_index in range(left_index, len(candidates)):
        right_item = candidates[right_index]
        key = quotient_key(left["mw"] + right_item["mw"])
        ordering = (
            left["q4o164_parent_degree"] + right_item["q4o164_parent_degree"],
            max(left["q4o164_parent_degree"], right_item["q4o164_parent_degree"]),
            left_index,
            right_index,
        )
        if key not in pairs or ordering < pairs[key]:
            pairs[key] = ordering

parent_pole_pairs = {}
for left_index, left in enumerate(parent_pole_candidates):
    for right_index in range(left_index, len(parent_pole_candidates)):
        right_item = parent_pole_candidates[right_index]
        key = quotient_key(left["mw"] + right_item["mw"])
        po_values = [left["q4o164_parent_a_minus_b"], right_item["q4o164_parent_a_minus_b"]]
        degrees = [left["q4o164_parent_degree"], right_item["q4o164_parent_degree"]]
        ordering = (
            sum(po_values), sum(degrees), max(po_values), max(degrees),
            left_index, right_index,
        )
        if key not in parent_pole_pairs or ordering < parent_pole_pairs[key]:
            parent_pole_pairs[key] = ordering

# Extra branches are relevant only if they strictly improve the best q12/o5867
# parent a-b total of six.  Build all optimal triple quotient states under the
# strict cutoff five, then use pair+triple and triple+triple meet-in-the-middle
# searches for exact lengths five and six.
strict_parent_pole_cutoff = 5
parent_pole_triples_under_cutoff = {}
for pair_key, pair in parent_pole_pairs.items():
    for index, item in enumerate(parent_pole_candidates):
        total_po = pair[0] + item["q4o164_parent_a_minus_b"]
        if total_po > strict_parent_pole_cutoff:
            continue
        key = quotient_key(
            parent_pole_candidates[pair[-2]]["mw"]
            + parent_pole_candidates[pair[-1]]["mw"]
            + item["mw"]
        )
        ordering = (
            total_po,
            pair[1] + item["q4o164_parent_degree"],
            max(pair[2], item["q4o164_parent_a_minus_b"]),
            max(pair[3], item["q4o164_parent_degree"]),
            pair[-2], pair[-1], index,
        )
        if (
            key not in parent_pole_triples_under_cutoff
            or ordering < parent_pole_triples_under_cutoff[key]
        ):
            parent_pole_triples_under_cutoff[key] = ordering


def word_order(indices):
    degrees = [candidates[index]["q4o164_parent_degree"] for index in indices]
    return (sum(degrees), max(degrees), tuple(indices))


def best_word(target_mw, target_key, length):
    if length == 1:
        words = [(index,) for index, key in enumerate(candidate_keys) if key == target_key]
    elif length == 2:
        words = [pairs[target_key][2:]] if target_key in pairs else []
    elif length == 3:
        words = []
        for index, candidate in enumerate(candidates):
            needed = quotient_key(target_mw - candidate["mw"])
            if needed in pairs:
                words.append((index,) + pairs[needed][2:])
    elif length == 4:
        words = []
        for left_key, left_pair in pairs.items():
            needed = quotient_subtract(target_key, left_key)
            if needed in pairs and left_key <= needed:
                words.append(left_pair[2:] + pairs[needed][2:])
    else:
        raise ValueError(length)
    return min(words, key=word_order) if words else None


def word_payload(indices, target_mw):
    word_sections = [candidates[index] for index in indices]
    residual = target_mw - sum(
        (item["mw"] for item in word_sections), vector(ZZ, 13)
    )
    coefficients = vector(
        QQ, known_matrix.transpose().solve_right(residual.column()).column(0)
    )
    assert all(value in ZZ for value in coefficients)
    coefficients = vector(ZZ, coefficients)
    assert sum(
        (item["mw"] for item in word_sections), vector(ZZ, 13)
    ) + coefficients * known_matrix == target_mw
    degrees = [item["q4o164_parent_degree"] for item in word_sections]
    return {
        "new_sections": [
            {
                "current_4A1_section": entries(item["current_4A1_section"]),
                "current_4A1_P_dot_O": int(item["current_4A1_P_dot_O"]),
                "current_4A1_mw": entries(item["mw"]),
                "q4o164_parent_curve": entries(item["q4o164_parent_curve"]),
                "q4o164_parent_degree": int(item["q4o164_parent_degree"]),
                "q4o164_parent_a_minus_b": int(item["q4o164_parent_a_minus_b"]),
                "current_4A1_component_pairings": entries(item["current_4A1_component_pairings"]),
                "current_4A1_affine_pairings": entries(item["current_4A1_affine_pairings"]),
                "physical_component_and_affine_nef_gate": True,
            }
            for item in word_sections
        ],
        "known_section_correction": [
            {"name": known_sections[index][0], "coefficient": int(value)}
            for index, value in enumerate(coefficients)
            if value
        ],
        "q4o164_parent_degree_sum": int(sum(degrees)),
        "q4o164_parent_degree_max": int(max(degrees)),
        "q4o164_parent_a_minus_b_sum": int(sum(
            item["q4o164_parent_a_minus_b"] for item in word_sections
        )),
        "q4o164_parent_a_minus_b_max": int(max(
            item["q4o164_parent_a_minus_b"] for item in word_sections
        )),
        "exact_mw_identity_pass": True,
    }


def item_word_order(items):
    degrees = [item["q4o164_parent_degree"] for item in items]
    return (
        sum(item["current_4A1_P_dot_O"] for item in items),
        sum(degrees),
        max(degrees),
        tuple(tuple(item["mw"]) for item in items),
    )


def item_word_payload(items, target_mw):
    residual = target_mw - sum((item["mw"] for item in items), vector(ZZ, 13))
    coefficients = vector(
        QQ, known_matrix.transpose().solve_right(residual.column()).column(0)
    )
    assert all(value in ZZ for value in coefficients)
    coefficients = vector(ZZ, coefficients)
    assert sum((item["mw"] for item in items), vector(ZZ, 13)) + coefficients * known_matrix == target_mw
    degrees = [item["q4o164_parent_degree"] for item in items]
    return {
        "new_sections": [
            {
                "current_4A1_section": entries(item["current_4A1_section"]),
                "current_4A1_P_dot_O": int(item["current_4A1_P_dot_O"]),
                "current_4A1_mw": entries(item["mw"]),
                "q4o164_parent_curve": entries(item["q4o164_parent_curve"]),
                "q4o164_parent_degree": int(item["q4o164_parent_degree"]),
                "q4o164_parent_a_minus_b": int(item["q4o164_parent_a_minus_b"]),
                "current_4A1_component_pairings": entries(item["current_4A1_component_pairings"]),
                "current_4A1_affine_pairings": entries(item["current_4A1_affine_pairings"]),
                "physical_component_and_affine_nef_gate": True,
            }
            for item in items
        ],
        "known_section_correction": [
            {"name": known_sections[index][0], "coefficient": int(value)}
            for index, value in enumerate(coefficients)
            if value
        ],
        "new_section_count": len(items),
        "total_P_dot_O": sum(item["current_4A1_P_dot_O"] for item in items),
        "q4o164_parent_degree_sum": int(sum(degrees)),
        "q4o164_parent_degree_max": int(max(degrees)),
        "q4o164_parent_a_minus_b_sum": int(sum(
            item["q4o164_parent_a_minus_b"] for item in items
        )),
        "q4o164_parent_a_minus_b_max": int(max(
            item["q4o164_parent_a_minus_b"] for item in items
        )),
        "exact_mw_identity_pass": True,
    }


def best_parent_pole_four_word(target_mw, target_key):
    pool = parent_pole_candidates
    words = []
    for left_key, left_pair in parent_pole_pairs.items():
        needed = quotient_subtract(target_key, left_key)
        if needed not in parent_pole_pairs or left_key > needed:
            continue
        right_pair = parent_pole_pairs[needed]
        indices = left_pair[-2:] + right_pair[-2:]
        items = [pool[index] for index in indices]
        po_values = [item["q4o164_parent_a_minus_b"] for item in items]
        degrees = [item["q4o164_parent_degree"] for item in items]
        words.append((
            sum(po_values), sum(degrees), max(po_values), max(degrees),
            tuple(indices), items,
        ))
    if not words:
        return None
    selected = min(words)
    return item_word_payload(selected[-1], target_mw)


def best_strict_parent_pole_long_word(target_mw, target_key, length):
    words = []
    if length == 5:
        left_states = parent_pole_pairs
        right_states = parent_pole_triples_under_cutoff
    elif length == 6:
        left_states = parent_pole_triples_under_cutoff
        right_states = parent_pole_triples_under_cutoff
    else:
        raise ValueError(length)
    for left_key, left in left_states.items():
        needed = quotient_subtract(target_key, left_key)
        if needed not in right_states or (length == 6 and left_key > needed):
            continue
        right = right_states[needed]
        total_po = left[0] + right[0]
        if total_po > strict_parent_pole_cutoff:
            continue
        left_count = 2 if length == 5 else 3
        left_indices = left[-left_count:]
        right_indices = right[-3:]
        indices = left_indices + right_indices
        items = [parent_pole_candidates[index] for index in indices]
        degrees = [item["q4o164_parent_degree"] for item in items]
        po_values = [item["q4o164_parent_a_minus_b"] for item in items]
        words.append((
            sum(po_values), sum(degrees), max(po_values), max(degrees),
            tuple(indices), items,
        ))
    return None if not words else item_word_payload(min(words)[-1], target_mw)


def best_mixed_words(target_mw, target_key):
    patterns = {"two_P0_P1": [], "two_P1_P1": [], "three_P0_P0_P1": []}
    for left in candidates:
        needed = quotient_key(target_mw - left["mw"])
        if needed in p1_by_quotient:
            patterns["two_P0_P1"].append([left, p1_by_quotient[needed][1]])
    for left in p1_candidates:
        needed = quotient_key(target_mw - left["mw"])
        if needed in p1_by_quotient:
            patterns["two_P1_P1"].append([left, p1_by_quotient[needed][1]])
    for p1 in p1_candidates:
        needed = quotient_key(target_mw - p1["mw"])
        if needed in pairs:
            pair = pairs[needed]
            patterns["three_P0_P0_P1"].append([
                candidates[pair[2]], candidates[pair[3]], p1,
            ])
    result = {}
    for name, words in patterns.items():
        result[name] = {
            "decomposition_count_before_word_symmetry": len(words),
            "best_word": None if not words else item_word_payload(
                min(words, key=item_word_order), target_mw
            ),
        }
    return result


targets = []
for cost_path in cost_paths:
    cost = json.loads(cost_path.read_text())
    for item in cost["retained_candidates"]:
        candidate_id = item["candidate_id"]
        # These cost files are already filtered to rootless candidates.  Retain
        # the endpoint field when present, but do not infer rootlessness from ADE.
        section = vector(ZZ, item["horizontal"]["section"])
        target_mw = vector(ZZ, section[-13:])
        target_key = quotient_key(target_mw)
        words = {}
        for length in range(1, 5):
            indices = best_word(target_mw, target_key, length)
            words[str(length)] = None if indices is None else word_payload(indices, target_mw)
        first_length = next((length for length in range(1, 5) if words[str(length)]), None)
        targets.append({
            "candidate_id": candidate_id,
            "source_cost_artifact": str(cost_path.relative_to(ROOT)),
            "horizontal_section": entries(section),
            "direct_P_dot_O": int(item["horizontal"]["P_dot_O"]),
            "direct_expected_RR_ambient": int(item["expected_RR_ambient"]),
            "mw": entries(target_mw),
            "minimum_word_length_at_most_four": first_length,
            "best_words_by_exact_length": words,
            "best_four_P0_word_by_parent_a_minus_b": best_parent_pole_four_word(
                target_mw, target_key
            ),
            "strict_parent_a_minus_b_improvements_with_more_P0_sections": {
                "cutoff": strict_parent_pole_cutoff,
                "five_sections": best_strict_parent_pole_long_word(
                    target_mw, target_key, 5
                ),
                "six_sections": best_strict_parent_pole_long_word(
                    target_mw, target_key, 6
                ),
            },
            "mixed_low_pole_words": best_mixed_words(target_mw, target_key),
        })

inputs = [marking_path, frame_path] + cost_paths
payload = {
    "schema": "elkies-k3.h3-rootless-p0-section-word-frontier.v1",
    "status": "PASS_EXACT_ROOTLESS_P0_SECTION_WORD_FRONTIER",
    "known_section_subgroup": {
        "rank": known_rank,
        "smith_diagonal": [int(value) for value in diagonal],
        "named_generators": [name for name, _ in known_sections],
    },
    "complete_P_dot_O_zero_shell": {
        "signed_norm_at_most_4_vector_count": int(shell[0]),
        "unique_mw_vector_count": len(candidates),
        "physical_P_dot_O_zero_class_count": physical_p0_class_count,
        "unique_pair_quotient_count": len(pairs),
        "signed_norm_at_most_6_vector_count": int(shell6[0]),
        "unique_P_dot_O_one_mw_vector_count": len(p1_candidates),
        "physical_P_dot_O_one_class_count": physical_p1_class_count,
    },
    "targets": targets,
    "compiler_interpretation": (
        "For each already retained rootless lattice exit, compare the direct high-pole "
        "horizontal with exact group-law words in polynomial P.O=0 sections. Parent degrees "
        "measure the expected difficulty of recovering those sections from q4/o164 data."
    ),
    "proof_boundary": (
        "The norm<=4 shell, physical simple/affine component gate, and word searches through "
        "length four are exhaustive in the exact marked MW quotient modulo the displayed "
        "explicit subgroup. Every emitted branch is an effective chamber section rather than "
        "a vertical-root pseudo-section. The length-five and length-six meet-in-the-middle "
        "searches are exhaustive under the displayed parent-a-minus-b cutoff. These are equation compiler "
        "plans, not QQ section constructions. Rootless and pinned-R17 status remains delegated "
        "to the source route certificates and must be checked separately before promotion."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in inputs],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in inputs
        },
    },
}
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
for target in targets:
    first = target["minimum_word_length_at_most_four"]
    best = target["best_words_by_exact_length"].get(str(first)) if first else None
    print(
        "ROOTLESSP0|q={}|orbit={}|direct_po={}|word_length={}|parent_sum={}|parent_max={}".format(
            target["candidate_id"]["q"], target["candidate_id"]["orbit_index"],
            target["direct_P_dot_O"], first,
            None if best is None else best["q4o164_parent_degree_sum"],
            None if best is None else best["q4o164_parent_degree_max"],
        )
    )
print("ROOTLESSP0|status={}|output={}".format(payload["status"], output))
