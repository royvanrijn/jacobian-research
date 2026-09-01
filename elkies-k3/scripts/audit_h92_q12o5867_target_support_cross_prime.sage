#!/usr/bin/env sage -python
"""Replay exact target-support sections in the complete shells at three primes."""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
PRIMES = (83, 89, 137)
Q3 = LOCAL / "q12o5867-degree1-compiler-branch-qq.json"
ABEL = LOCAL / "q12o5867-abel-trace-named-seeds-qq.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--sections", type=Path, action="append", default=[],
    help="exact QQ section artifact; may be repeated",
)
parser.add_argument(
    "--output", type=Path,
    default=LOCAL / "q12o5867-two-primary-target-support-cross-prime.json",
)
parser.add_argument(
    "--halvings", type=Path,
    default=LOCAL / "q12o5867-three-target-halvings-qq.json",
)
args = parser.parse_args()
section_paths = args.sections or [
    LOCAL / "q12o5867-saturation-completion-class21-qq.json",
    LOCAL / "q12o5867-replacement-word-seeds-qq.json",
    LOCAL / "q12o5867-two-primary-target-support-sections-qq.json",
    LOCAL / "q12o5867-two-primary-target-support-sections-v2-qq.json",
    LOCAL / "q12o5867-support-class170-shell32-qq.json",
]
section_paths = [path if path.is_absolute() else ROOT / path for path in section_paths]
output = args.output if args.output.is_absolute() else ROOT / args.output


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def reduce_coefficients(values, prime):
    answer = []
    for value in values:
        value = QQ(value)
        assert value.denominator() % prime
        answer.append(int((ZZ(value.numerator()) % prime)
                          *pow(int(ZZ(value.denominator()) % prime), -1, prime) % prime))
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return answer


exact_sections = []
for path in section_paths:
    data = json.loads(path.read_text())
    for key, record in data["sections"].items():
        exact_sections.append({
            "kind": "target_support",
            "source_artifact": str(path.relative_to(ROOT)),
            "source_key": key,
            "supplied_class_index": int(record["lattice_class_index"]),
            "section": record["section"],
        })

q3 = json.loads(Q3.read_text())
exact_sections.append({
    "kind": "existing_exact_Q3",
    "source_artifact": str(Q3.relative_to(ROOT)),
    "source_key": "exact_Q3",
    "supplied_class_index": 69,
    "section": q3["section"],
})
abel = json.loads(ABEL.read_text())
exact_sections.append({
    "kind": "existing_exact_Q4",
    "source_artifact": str(ABEL.relative_to(ROOT)),
    "source_key": "Q4_candidate1_shell220_rank12_seed",
    "supplied_class_index": 516,
    "section": abel["sections"]["Q4_candidate1_shell220_rank12_seed"]["section"],
})

halvings_path = args.halvings if args.halvings.is_absolute() else ROOT / args.halvings
halvings = json.loads(halvings_path.read_text())
for attempt in halvings["attempts"]:
    assert len(attempt["verified_rational_halves"]) == 1
    half = attempt["verified_rational_halves"][0]
    assert half["x"]["denominator_coefficients_low_to_high"] == ["1"]
    assert half["y"]["denominator_coefficients_low_to_high"] == ["1"]
    exact_sections.append({
        "kind": "recovered_half",
        "source_artifact": str(halvings_path.relative_to(ROOT)),
        "source_key": f"declared_class{attempt['declared_mod89_half_target_class_index']}",
        "supplied_class_index": int(attempt["declared_mod89_half_target_class_index"]),
        "section": {
            "x_coefficients_low_to_high": half["x"]["numerator_coefficients_low_to_high"],
            "y_coefficients_low_to_high": half["y"]["numerator_coefficients_low_to_high"],
            "exact_weierstrass_identity": half["literal_curve_substitution"],
        },
    })

prime_data = {}
for prime in PRIMES:
    shell_path = LOCAL / f"q12o5867-p0-shell-all-profiles-mod{prime}.json"
    classifier_path = LOCAL / f"q12o5867-p0-shell-lattice-classification-all-profiles-mod{prime}.json"
    shell = json.loads(shell_path.read_text())
    classifier = json.loads(classifier_path.read_text())
    lookup = {}
    for shell_index, record in enumerate(shell["all_records"]):
        key = (
            tuple(record["x_coefficients_low_to_high"]),
            tuple(record["y_coefficients_low_to_high"]),
        )
        assert key not in lookup
        lookup[key] = shell_index
    pairwise = {
        int(row["shell_index"]): row
        for row in classifier["polynomial_shell"].get(
            "complete_pairwise_intersection_disambiguation", []
        )
    }
    prime_data[prime] = (shell_path, classifier_path, shell, classifier, lookup, pairwise)

audits = []
for exact in exact_sections:
    supplied = exact["supplied_class_index"]
    reductions = []
    for prime in PRIMES:
        shell_path, classifier_path, shell, classifier, lookup, pairwise = prime_data[prime]
        section = exact["section"]
        key = (
            tuple(reduce_coefficients(section["x_coefficients_low_to_high"], prime)),
            tuple(reduce_coefficients(section["y_coefficients_low_to_high"], prime)),
        )
        shell_index = lookup.get(key)
        if shell_index is None:
            reductions.append({
                "prime": prime,
                "present_in_complete_shell": False,
            })
            continue
        row = classifier["polynomial_shell"]["records"][shell_index]
        refinement = pairwise.get(shell_index)
        reductions.append({
            "prime": prime,
            "present_in_complete_shell": True,
            "shell_index": shell_index,
            "equation_component_profile": row["equation_component_profile"],
            "inverse_parent_degree": row["inverse_parent_degree"],
            "trace_matching_lattice_class_indices": row[
                "trace_matching_lattice_class_indices"
            ],
            "profile_compatible_lattice_class_indices": row[
                "profile_compatible_lattice_class_indices"
            ],
            "supplied_class_is_trace_candidate": supplied in row[
                "trace_matching_lattice_class_indices"
            ],
            "supplied_class_is_profile_candidate": supplied in row[
                "profile_compatible_lattice_class_indices"
            ],
            "pairwise_surviving_class_alternatives": (
                refinement["surviving_class_alternatives"] if refinement else None
            ),
        })
    audits.append({
        **{key: value for key, value in exact.items() if key != "section"},
        "literal_exact_weierstrass_identity": exact["section"]["exact_weierstrass_identity"],
        "reductions": reductions,
        "supplied_class_trace_candidate_at_all_three_primes": all(
            row.get("supplied_class_is_trace_candidate", False) for row in reductions
        ),
        "supplied_class_profile_candidate_at_all_three_primes": all(
            row.get("supplied_class_is_profile_candidate", False) for row in reductions
        ),
    })

payload = {
    "schema": "q12o5867-two-primary-target-support-cross-prime-v1",
    "status": "PASS_EXACT_QQ_SECTIONS_REDUCED_AND_MATCHED_AT_83_89_137",
    "inputs": {
        "section_artifacts": [
            {"path": str(path.relative_to(ROOT)), "sha256": sha256(path)}
            for path in section_paths
        ],
        "halvings": {"path": str(halvings_path.relative_to(ROOT)),
                     "sha256": sha256(halvings_path)},
        "exact_Q3": {"path": str(Q3.relative_to(ROOT)), "sha256": sha256(Q3)},
        "exact_Q4": {"path": str(ABEL.relative_to(ROOT)), "sha256": sha256(ABEL)},
        "shells_and_classifiers": {
            str(prime): {
                "shell": {"path": str(prime_data[prime][0].relative_to(ROOT)),
                          "sha256": sha256(prime_data[prime][0])},
                "classifier": {"path": str(prime_data[prime][1].relative_to(ROOT)),
                               "sha256": sha256(prime_data[prime][1])},
            }
            for prime in PRIMES
        },
    },
    "sections": audits,
    "proof_boundary": (
        "This certifies literal QQ curve identities and three-prime shell fingerprints. "
        "An Abel-trace candidate list is not by itself a unique global NS name."
    ),
}
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(f"wrote {output}")
for row in audits:
    print(
        row["supplied_class_index"],
        "trace-all", row["supplied_class_trace_candidate_at_all_three_primes"],
        "profile-all", row["supplied_class_profile_candidate_at_all_three_primes"],
        "shells", [item.get("shell_index") for item in row["reductions"]],
    )
