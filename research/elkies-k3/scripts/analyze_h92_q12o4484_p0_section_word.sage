#!/usr/bin/env sage -python
"""Find a minimum low-parent-degree P.O=0 word for the q12/o4484 section.

status: ACTIVE_COMPILER
claim: exhaustive length-at-most-four P.O=0 word search modulo explicit sections
inputs: physical q8/o376 4A1 marking and q12/o4484 equation-cost artifact
outputs: generated q12/o4484 four-P.O=0-section word certificate
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
parser.add_argument("--q12-cost", type=Path, required=True)
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
cost_path = args.q12_cost.resolve()
output = args.output.resolve()
marking = json.loads(marking_path.read_text())
cost = json.loads(cost_path.read_text())
frame_path = ROOT / marking["frame_output"]
frame = load_matrix(frame_path)
root_rank = int(marking["root_data"][0])
assert root_rank == 4
gram = block_diagonal_matrix(U2, -frame)
basis_in_parent = matrix(ZZ, marking["basis_in_source"])
parent_in_basis = matrix(ZZ, marking["source_in_basis"])
parent_fibre = vector(ZZ, [1, 0] + [0] * 17) * parent_in_basis

selected = next(
    item for item in cost["retained_candidates"]
    if item["candidate_id"] == {"q": 12, "old_fibre_degree": 2, "orbit_index": 4484}
)
target_section = vector(ZZ, selected["horizontal"]["section"])
target_mw = vector(ZZ, target_section[-13:])

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
physical_class_count = 0
for tail in half + [-value for value in half]:
    if tail * frame * tail != 4:
        continue
    component_pairings = vector(ZZ, tail) * frame[:, :root_rank]
    affine_pairings = vector(ZZ, [1] * root_rank) - component_pairings
    if min(tuple(component_pairings) + tuple(affine_pairings)) < 0:
        continue
    physical_class_count += 1
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
    candidate = {
        "current_4A1_section": section,
        "current_4A1_tail": tail,
        "mw": mw,
        "q4o164_parent_curve": section * basis_in_parent,
        "q4o164_parent_degree": parent_degree,
        "current_4A1_component_pairings": component_pairings,
        "current_4A1_affine_pairings": affine_pairings,
    }
    if key not in best_by_mw or ordering < best_by_mw[key][0]:
        best_by_mw[key] = (ordering, candidate)

candidates = [item[1] for item in best_by_mw.values()]
target_key = quotient_key(target_mw)
one_section_count = sum(quotient_key(item["mw"]) == target_key for item in candidates)

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

two_section_exists = target_key in pairs
three_section_count = 0
for candidate in candidates:
    needed = quotient_key(target_mw - candidate["mw"])
    three_section_count += int(needed in pairs)

four_section_words = []
for left_key, left_pair in pairs.items():
    needed = quotient_subtract(target_key, left_key)
    if needed not in pairs or left_key > needed:
        continue
    right_pair = pairs[needed]
    indices = left_pair[2:] + right_pair[2:]
    four_section_words.append((
        left_pair[0] + right_pair[0],
        max(left_pair[1], right_pair[1]),
        indices,
    ))
four_section_words.sort()
assert not one_section_count and not two_section_exists and not three_section_count
assert four_section_words

minimum_parent_degree_sum, minimum_parent_degree_max, selected_indices = four_section_words[0]
word_sections = [candidates[index] for index in selected_indices]
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

known_word = [
    {"name": known_sections[index][0], "coefficient": int(value)}
    for index, value in enumerate(coefficients)
    if value
]
section_payload = []
for index, item in enumerate(word_sections, 1):
    current = item["current_4A1_section"]
    parent = item["q4o164_parent_curve"]
    assert current * gram * current == -2 and current[1] == 1 and current[0] - current[1] == 0
    section_payload.append({
        "name": f"Q{index}",
        "current_4A1_section": entries(current),
        "current_4A1_P_dot_O": 0,
        "current_4A1_mw": entries(item["mw"]),
        "q4o164_parent_curve": entries(parent),
        "q4o164_parent_degree": int(item["q4o164_parent_degree"]),
        "current_4A1_component_pairings": entries(item["current_4A1_component_pairings"]),
        "current_4A1_affine_pairings": entries(item["current_4A1_affine_pairings"]),
        "physical_component_and_affine_nef_gate": True,
    })

inputs = (marking_path, cost_path, frame_path)
payload = {
    "schema": "elkies-k3.h3-q12o4484-p0-section-word.v1",
    "status": "PASS_EXACT_Q12O4484_FOUR_P0_SECTION_WORD",
    "target": {
        "candidate_id": selected["candidate_id"],
        "horizontal_section": entries(target_section),
        "P_dot_O": int(selected["horizontal"]["P_dot_O"]),
        "expected_RR_ambient": int(selected["expected_RR_ambient"]),
        "mw": entries(target_mw),
    },
    "known_section_subgroup": {
        "rank": known_rank,
        "smith_diagonal": [int(value) for value in diagonal],
        "named_generators": [name for name, _ in known_sections],
    },
    "complete_P_dot_O_zero_shell": {
        "signed_norm_at_most_4_vector_count": int(shell[0]),
        "physical_P_dot_O_zero_class_count": physical_class_count,
        "unique_mw_vector_count": len(candidates),
        "one_section_decomposition_count": one_section_count,
        "two_section_decomposition_exists": two_section_exists,
        "three_section_decomposition_count": three_section_count,
        "four_section_decomposition_count": len(four_section_words),
    },
    "selected_word": {
        "new_sections": section_payload,
        "known_section_correction": known_word,
        "minimum_q4o164_parent_degree_sum": int(minimum_parent_degree_sum),
        "minimum_q4o164_parent_degree_max": int(minimum_parent_degree_max),
        "mw_identity": "P4484 = Q1+Q2+Q3+Q4 + named known-section correction",
        "exact_mw_identity_pass": True,
    },
    "compiler_interpretation": (
        "Replace one direct P.O=10 section reconstruction by four P.O=0 sections on the "
        "current 4A1 equation, of q4/o164 parent degrees 3,2,3,2, followed by exact group law "
        "and the displayed already-explicit section correction."
    ),
    "proof_boundary": (
        "The norm<=4 shell, physical simple/affine component gate, and word-length exclusions "
        "are exhaustive in the exact marked MW quotient. Every emitted branch is an effective "
        "chamber section, not a vertical-root pseudo-section. This is an exact lattice/group-law "
        "compilation plan, not yet four QQ section "
        "equations or a measured equation-lift timing. The existing q12/o4484 lattice route and "
        "pinned-R17 endpoint certificate remain unchanged."
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
print(
    "Q12P0WORD|shell={}|mw={}|length_lt4=0|length4={}|parent_sum={}|parent_max={}|"
    "status={}|output={}".format(
        payload["complete_P_dot_O_zero_shell"]["signed_norm_at_most_4_vector_count"],
        len(candidates), len(four_section_words), minimum_parent_degree_sum,
        minimum_parent_degree_max, payload["status"], output,
    )
)
