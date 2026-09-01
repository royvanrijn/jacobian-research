#!/usr/bin/env python3
"""Bind the exact artifacts for the declared bounded wgxli rebasing rejection.

This checker does no new search.  It verifies the hashes and status fields of
the sign/permutation audit, the one-shear audit, and their complete projective
first-jet eliminations.  Its conclusion is only about the explicitly recorded
rebasing bound; it is not a rejection of arbitrary GL(17,Z) changes of basis.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ARTIFACT_DIR = ROOT / "artifacts/generated-results/elliptic-curves"
DEFAULT_OUTPUT = ARTIFACT_DIR / "icarm_wgxli_rank17_bounded_rejection_v1.json"
FILES = {
    "lineage": ARTIFACT_DIR / "icarm_wgxli_rank17_lineage_v1.json",
    "rebasing": ARTIFACT_DIR / "icarm_wgxli_rank17_signed_permutation_rebasing_v1.json",
    "mutation": ARTIFACT_DIR / "icarm_wgxli_rank17_mutation_p4_minus_p1_v1.json",
    "literal_mod17": ARTIFACT_DIR / "icarm_wgxli_rank17_first_jet_mod17_v2.json",
    "literal_mod53": ARTIFACT_DIR / "icarm_wgxli_rank17_first_jet_mod53_v2.json",
    "literal_mod67": ARTIFACT_DIR / "icarm_wgxli_rank17_first_jet_mod67_v1.json",
    "mutation_mod17": ARTIFACT_DIR
    / "icarm_wgxli_rank17_mutation_p4_minus_p1_first_jet_mod17_v1.json",
    "mutation_mod53": ARTIFACT_DIR
    / "icarm_wgxli_rank17_mutation_p4_minus_p1_first_jet_mod53_v1.json",
}
TARGETS = (351, 356, 376, 377, 385)
NEGATIVE_CONTROLS = ((363, 364, 378), (389, 390, 391))


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    return parser.parse_args()


def sha256(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def load_inputs():
    raw = {name: path.read_bytes() for name, path in FILES.items()}
    return raw, {name: json.loads(value) for name, value in raw.items()}


def assert_complete_elimination(record, prime, input_hash):
    assert record["schema"] == "icarm.wgxli-rank17-first-jet-elimination.v1"
    assert record["status"] == "PASS_COMPLETE_DISTINCT_CHART_EMPTY"
    assert record["prime"] == prime
    assert record["input"]["sha256"] == input_hash
    chart = record["chart"]
    assert chart["complete"] is True
    assert chart["parameter_space"] == (
        "ordered pairs in P1(F_p) minus the three fixed nodes"
    )
    assert chart["ordered_pair_count"] == (prime - 2) * (prime - 3)
    assert chart["solution_pair_count"] == 0
    assert chart["timeout_pair_count"] == 0
    assert chart["normalized_chart_types"] == [
        "both residual nodes finite",
        "t_376 at infinity",
        "t_377 at infinity",
    ]


def main():
    arguments = parse_args()
    raw, data = load_inputs()
    lineage_hash = sha256(raw["lineage"])
    rebasing_hash = sha256(raw["rebasing"])
    mutation_hash = sha256(raw["mutation"])

    rebasing = data["rebasing"]
    assert rebasing["status"] == "PASS_BOUNDED_SIGN_AND_PERMUTATION_ALIGNMENT"
    retained = rebasing["retained_signed_permutation_candidates"]
    assert len(retained) == 1
    for curve_id in TARGETS:
        row = retained[0][str(curve_id)]
        assert row["signs"] == "+" * 17
        assert row["permutation_new_label_to_old_point"] == list(range(1, 18))
    assert rebasing["exact_signed_permutation_group_law_replay"][
        "all_85_transformed_points_verified"
    ] is True
    components = tuple(
        tuple(component) for component in rebasing["negative_controls"]["components"]
    )
    assert all(component in components for component in NEGATIVE_CONTROLS)
    assert rebasing["negative_controls"]["status"] == "PASS_CONTROLS_REMAIN_SEPARATE"

    mutation = data["mutation"]
    assert mutation["status"] == "PASS_EXACT_POINT_CONSTRUCTION_FOR_ONE_BOUNDED_MUTATION"
    assert mutation["inputs"][str(FILES["lineage"].relative_to(ROOT))] == lineage_hash
    assert mutation["inputs"][str(FILES["rebasing"].relative_to(ROOT))] == rebasing_hash
    bound = mutation["declared_search_bound"]
    assert bound["enumerated_transform_count"] == 352
    assert bound["retained_transform_count"] == 1
    assert bound["elementary_mutation_count"] == 1
    assert bound["maximum_column_l1_norm"] == 2
    assert bound["unimodular_determinant"] == 1
    assert mutation["retained_proposal"]["new_point_word"] == [[4, 1], [1, -1]]
    assert mutation["exact_group_law"] == {
        "all_five_mutated_points_verified_on_their_short_curves": True,
        "transformation_determinant": 1,
        "word": "P4-P1",
    }

    assert_complete_elimination(data["literal_mod17"], 17, lineage_hash)
    assert_complete_elimination(data["literal_mod53"], 53, lineage_hash)
    assert_complete_elimination(data["literal_mod67"], 67, lineage_hash)
    assert_complete_elimination(data["mutation_mod17"], 17, mutation_hash)
    assert_complete_elimination(data["mutation_mod53"], 53, mutation_hash)
    for name in ("literal_mod53", "literal_mod67"):
        control = data[name]["positive_control"]
        assert control["available"] is True
        assert control["verified_section_count"] == 17
        assert control["status"] == "PASS_EXACT_PROJECTIVE_FIRST_JET_SYSTEM"
        assert len(control["normalized_chart_controls"]) == 3
        assert all(
            chart["status"] == "PASS_ELIMINATED_SYSTEM_WITNESS"
            for chart in control["normalized_chart_controls"]
        )

    payload = {
        "schema": "icarm.wgxli-rank17-bounded-rejection.v1",
        "status": "PASS_FORMAL_BOUNDED_REJECTION_AT_TWO_CLEAN_PRIMES",
        "inputs": {
            str(FILES[name].relative_to(ROOT)): sha256(raw[name])
            for name in FILES
        },
        "target_curves": list(TARGETS),
        "literal_projective_first_jet_rejections": {
            "17": data["literal_mod17"]["chart"]["ordered_pair_count"],
            "53": data["literal_mod53"]["chart"]["ordered_pair_count"],
            "67": data["literal_mod67"]["chart"]["ordered_pair_count"],
        },
        "retained_signed_permutation_count": 1,
        "retained_signed_permutation": "displayed signs and displayed order",
        "bounded_mutation_enumerated_count": 352,
        "bounded_mutation_retained_count": 1,
        "retained_bounded_mutation": "P4-P1",
        "retained_mutation_projective_first_jet_rejections": {
            "17": data["mutation_mod17"]["chart"]["ordered_pair_count"],
            "53": data["mutation_mod53"]["chart"]["ordered_pair_count"],
        },
        "negative_controls": [list(component) for component in NEGATIVE_CONTROLS],
        "conclusion": (
            "The five records share the pinned numerical/private-pipeline fingerprint, "
            "but have no rootless-K3 (8,12;4,6) realization under the declared signed-"
            "permutation and one-anchor-preserving-elementary-mutation rebasing."
        ),
        "proof_boundary": (
            "This is a finite exact rejection inside the recorded bounds. It is not an "
            "unrestricted GL(17,Z) search, does not exclude larger or repeated basis "
            "mutations, and does not exclude a different family shape or reductions bad "
            "at every tested prime. Numerical Gram alignment is used only to define the "
            "finite candidates; every retained candidate is constructed exactly and "
            "rejected by algebraically closed modular elimination."
        ),
    }
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if arguments.output.read_text() != rendered:
            raise SystemExit("stale bounded-rejection artifact")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered)
        print(
            f"WGXLIREJECT|output={arguments.output}|"
            f"sha256={sha256(rendered.encode())}"
        )
    print("WGXLIREJECT|status=PASS_FORMAL_BOUNDED_REJECTION_AT_TWO_CLEAN_PRIMES")


if __name__ == "__main__":
    main()
