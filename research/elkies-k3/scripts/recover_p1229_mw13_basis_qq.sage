#!/usr/bin/env sage -python
"""Fail-closed recovery gate for a P1229-pointed equation-level MW13 basis.

status: ACTIVE_SEARCH
claim: exact obstruction to the inverse-parent-degree-one MW13 shortcut
inputs: pinned P1229 marking, three modular shell classifications, 42 controls
outputs: artifacts/generated-results/elkies-k3-p1229-mw13-degree1-recovery-gate.json

The proposed shortcut starts with polynomial sections on the P1229-pointed
4A1 model whose inverse q4-parent map has degree one.  This script checks that
proposal against the complete marked lattice shell and against independently
enumerated shells at p=83, 89, and 137.  It also reconstructs the exact
Shioda height lattice of the parent from the root-adapted frame.

The script deliberately stops before Hensel lifting when the filtered classes
cannot span rank 13.  In that case lifting, saturation, control-point words,
and a point-factory complexity verdict would be mathematically unjustified.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import QQ, ZZ, matrix


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
FRAME = GENERATED / "elkies-k3-h3-q4o164-c8-q8o376-4a1-p1229-frame.txt"
MARKING = GENERATED / "elkies-k3-h3-q4o164-c8-q8o376-4a1-p1229-marking.json"
CONTROLS = GENERATED / "elkies-k3-q12o5867-genus-one-point-factory-controls.json"
LOCAL = ROOT / "artifacts/local/elkies-k3"
PRIMES = (83, 89, 137)
CLASSIFICATIONS = {
    prime: LOCAL / f"q12o5867-p0-shell-lattice-classification-mod{prime}.json"
    for prime in PRIMES
}
SHELLS = {
    prime: LOCAL / f"q12o5867-p0-shell-all-records-mod{prime}.json"
    for prime in PRIMES
}
DEFAULT_OUTPUT = GENERATED / "elkies-k3-p1229-mw13-degree1-recovery-gate.json"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_matrix(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


def rank_of_classes(classes, indices):
    indices = sorted(set(indices))
    if not indices:
        return 0
    return int(matrix(QQ, [
        classes[index]["current_4A1_mw"] for index in indices
    ]).rank())


def refined_regular_realizations(data):
    """Return class -> shell records after exact smooth-intersection refinement."""
    records = data["polynomial_shell"]["records"]
    answer = {}
    for record in records:
        alternatives = record["profile_compatible_lattice_class_indices"]
        if record["ordinary_coefficient_jacobian_rank"] == 12 and len(alternatives) == 1:
            answer.setdefault(alternatives[0], []).append(record)
    for result in data["polynomial_shell"].get(
        "complete_pairwise_intersection_disambiguation", []
    ):
        if not result["resolved_uniquely"]:
            continue
        record = records[result["shell_index"]]
        if record["ordinary_coefficient_jacobian_rank"] != 12:
            continue
        answer.setdefault(result["surviving_class_alternatives"][0], []).append(record)
    return answer


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
output = args.output if args.output.is_absolute() else ROOT / args.output

marking = json.loads(MARKING.read_text())
controls = json.loads(CONTROLS.read_text())
classifications = {
    prime: json.loads(path.read_text()) for prime, path in CLASSIFICATIONS.items()
}
shells = {prime: json.loads(path.read_text()) for prime, path in SHELLS.items()}

assert marking["status"] == "PASS_EXACT_PHYSICAL_AN_EFFECTIVE_ZERO_MARKING"
assert controls["status"] == "PASS_EXACT_Q12O5867_BIRATIONAL_POINT_MAP_AND_CONTROL_ROUNDTRIPS"
assert sum(len(control["points"]) for control in controls["controls"]) == 42
assert all(
    point["exact_forward_inverse_roundtrip"]
    and point["exact_parent_and_child_equation_checks"]
    for control in controls["controls"] for point in control["points"]
)

for prime in PRIMES:
    assert classifications[prime]["prime"] == prime
    assert shells[prime]["prime"] == prime
    assert len(shells[prime]["all_records"]) == len(
        classifications[prime]["polynomial_shell"]["records"]
    )

# The lattice payload is prime-independent and is protected by the classifier's
# input hashes.  Assert literal equality before using the p=89 copy below.
lattice_payloads = [
    classifications[prime]["lattice_shell"]["classes"] for prime in PRIMES
]
prime_dependent_class_fields = {
    "realizing_shell_indices", "possible_realizing_shell_indices"
}
def stable_lattice_payload(items):
    return [
        {key: value for key, value in item.items()
         if key not in prime_dependent_class_fields}
        for item in items
    ]
assert (
    stable_lattice_payload(lattice_payloads[0])
    == stable_lattice_payload(lattice_payloads[1])
    == stable_lattice_payload(lattice_payloads[2])
)
classes = lattice_payloads[1]
assert len(classes) == 938

# Reconstruct the exact saturated MW13 height Gram as the Schur complement of
# the four A1 roots in the positive root-adapted frame.
frame = load_matrix(FRAME)
assert frame.nrows() == frame.ncols() == 17
root_gram = frame[:4, :4]
root_mw = frame[:4, 4:]
mw_block = frame[4:, 4:]
mw_gram = mw_block-root_mw.transpose()*root_gram.inverse()*root_mw
assert root_gram.det() == 16
assert frame.det() == 948
assert mw_gram.det() == QQ(237)/4

# This is the decisive exhaustive gate.  It uses all physical P.O=0 lattice
# classes, not merely the classes that happened to appear at one prime.
all_degree_one = [
    item["class_index"] for item in classes
    if item["q4o164_parent_degree"] == 1
]
all_degree_one_rank = rank_of_classes(classes, all_degree_one)
assert len(all_degree_one) == 23
assert all_degree_one_rank == 12

prime_records = {}
for prime in PRIMES:
    records = classifications[prime]["polynomial_shell"]["records"]
    degree_one_records = [record for record in records if record["inverse_parent_degree"] == 1]
    compatible = {
        index for record in degree_one_records
        for index in record["profile_compatible_lattice_class_indices"]
    }
    unique = {
        record["profile_compatible_lattice_class_indices"][0]
        for record in degree_one_records
        if len(record["profile_compatible_lattice_class_indices"]) == 1
    }
    refined = refined_regular_realizations(classifications[prime])
    refined_degree_one = {
        index for index, rows in refined.items()
        if classes[index]["q4o164_parent_degree"] == 1 and rows
    }
    prime_records[str(prime)] = {
        "signed_shell_record_count": len(records),
        "inverse_degree_one_record_count": len(degree_one_records),
        "profile_compatible_class_count": len(compatible),
        "profile_compatible_class_rank": rank_of_classes(classes, compatible),
        "unique_class_count": len(unique),
        "unique_class_rank": rank_of_classes(classes, unique),
        "regular_pairwise_refined_class_count": len(refined_degree_one),
        "regular_pairwise_refined_class_rank": rank_of_classes(
            classes, refined_degree_one
        ),
        "regular_pairwise_refined_class_indices": sorted(refined_degree_one),
    }

refined_maps = {
    prime: refined_regular_realizations(classifications[prime]) for prime in PRIMES
}
common_refined = set.intersection(*(set(mapping) for mapping in refined_maps.values()))
common_refined_degree_one = {
    index for index in common_refined
    if classes[index]["q4o164_parent_degree"] == 1
}
assert rank_of_classes(classes, common_refined_degree_one) <= 12

# For comparison, the complete physical P.O=0 lattice shell is saturated.  The
# presently enumerated two-profile modular shell is not: its regular refined
# classes have an index-eight defect.  Thus a viable follow-up must widen the
# component-profile enumeration as well as admit higher inverse degree.
all_p0_mw = matrix(ZZ, [item["current_4A1_mw"] for item in classes])
assert all_p0_mw.rank() == 13
p0_elementary_divisors = [
    int(value) for value in all_p0_mw.elementary_divisors() if value
]
assert p0_elementary_divisors == [1]*13
p89_refined = refined_maps[89]
p89_refined_mw = matrix(ZZ, [
    classes[index]["current_4A1_mw"] for index in sorted(p89_refined)
])
p89_refined_elementary_divisors = [
    int(value) for value in p89_refined_mw.elementary_divisors() if value
]
assert p89_refined_mw.rank() == 13
assert p89_refined_elementary_divisors == [1]*10+[2, 2, 2]

input_paths = [FRAME, MARKING, CONTROLS]
for prime in PRIMES:
    input_paths.extend((CLASSIFICATIONS[prime], SHELLS[prime]))

payload = {
    "schema": "elkies-k3.p1229-mw13-degree1-recovery-gate.v1",
    "status": "REJECTED_EXACT_DEGREE1_SHELL_CANNOT_SPAN_MW13",
    "parent_lattice": {
        "root_type": "4A1",
        "mw_rank": 13,
        "ns_frame_determinant": int(frame.det()),
        "root_gram_determinant": int(root_gram.det()),
        "saturated_shioda_height_gram": [
            [str(value) for value in row] for row in mw_gram.rows()
        ],
        "saturated_shioda_height_determinant": str(mw_gram.det()),
    },
    "complete_lattice_degree_one_gate": {
        "physical_P_dot_O_zero_class_count": len(classes),
        "inverse_parent_degree_one_class_count": len(all_degree_one),
        "inverse_parent_degree_one_class_indices": all_degree_one,
        "inverse_parent_degree_one_span_rank": all_degree_one_rank,
        "required_rank": 13,
        "pass": False,
    },
    "three_good_prime_shells": prime_records,
    "three_prime_common_pairwise_refined_degree_one": {
        "class_indices": sorted(common_refined_degree_one),
        "class_count": len(common_refined_degree_one),
        "span_rank": rank_of_classes(classes, common_refined_degree_one),
    },
    "component_profile_coverage_and_saturation": {
        "complete_physical_P_dot_O_zero_shell_rank": int(all_p0_mw.rank()),
        "complete_physical_P_dot_O_zero_shell_smith_nonzero_diagonal": p0_elementary_divisors,
        "complete_physical_P_dot_O_zero_shell_index": 1,
        "p89_regular_pairwise_refined_two_profile_shell_rank": int(p89_refined_mw.rank()),
        "p89_regular_pairwise_refined_two_profile_shell_smith_nonzero_diagonal": p89_refined_elementary_divisors,
        "p89_regular_pairwise_refined_two_profile_shell_index": 8,
        "required_extra_certificate": (
            "At least one inverse-degree-greater-than-one class is needed for rank. "
            "The modular enumeration must also include the missing I2 component profiles "
            "or supply three exact 2-divisions to remove its current index-8 defect."
        ),
    },
    "control_point_gate": {
        "exact_roundtrip_control_count": 42,
        "mw13_words_computed": 0,
        "reason": (
            "No saturated equation-level MW13 basis was produced; assigning words in an "
            "index-defective or rank-12 set would be invalid."
        ),
    },
    "point_factory_decision": {
        "decision": "DO_NOT_ENUMERATE_AND_DO_NOT_CLOSE_FROM_THIS_GATE",
        "reason": (
            "The proposed basis-recovery premise fails before the statistical complexity "
            "test. The ten invisible rank-28 directions therefore remain inconclusive."
        ),
    },
    "proof_boundary": (
        "This is an exact exhaustive obstruction to the inverse-parent-degree-one recovery "
        "workflow inside the complete displayed P.O=0 physical lattice shell. It is not a "
        "nonexistence theorem for a saturated equation-level MW13 basis obtained from higher "
        "inverse degree, positive P.O, or exact division."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in input_paths],
        "sha256": {
            str(path.relative_to(ROOT)): sha256(path) for path in input_paths
        },
    },
}

output.parent.mkdir(parents=True, exist_ok=True)
output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "P1229MW13GATE|degree1_classes={}|degree1_rank={}|required_rank=13|"
    "p0_smith={}|controls=42|status={}|output={}".format(
        len(all_degree_one), all_degree_one_rank, p0_elementary_divisors,
        payload["status"], output,
    ),
    flush=True,
)
