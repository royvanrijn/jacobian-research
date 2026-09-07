#!/usr/bin/env python3
"""Build the durable expanded-q12 compiler and pinned-route comparison.

status: ACTIVE_PROOF
claim: exact compiler Pareto audit for the cap-50000 q12 rootless sample
inputs: rootless search/cost/frontier artifacts and exact marked route certificates
outputs: generated expanded-q12 compiler Pareto artifact
"""

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GEN = ROOT / "artifacts/generated-results"
FILES = {
    "neighbors": GEN / "elkies-k3-h3-q4o164-q8o376-4a1-p1229-q12d2-cap50000-mw17-neighbors.json",
    "cost": GEN / "elkies-k3-h3-q4o164-q8o376-4a1-p1229-q12d2-cap50000-rootless-equation-cost.json",
    "frontier": GEN / "elkies-k3-h3-q4o164-q8o376-rootless-p0-section-word-frontier.json",
    "cap_stability": GEN / "elkies-k3-h3-q4o164-q8o376-q12-cap50000-cap100000-rootless-stability-audit.json",
    "old_edge": GEN / "elkies-k3-h3-q4o164-q8o376-q12o5867-rootless-certificate.json",
    "old_route": GEN / "elkies-k3-h3-q4o208-q323-free-q4o1584-q4o164-q8o376-q12o5867-pinned-r17-route-certificate.json",
    "winner_edge": GEN / "elkies-k3-h3-q4o164-q8o376-q12-fp-6-2-m13-10-rootless-certificate.json",
    "winner_route": GEN / "elkies-k3-h3-q4o208-q323-free-q4o1584-q4o164-q8o376-q12-fp-m13-10-pinned-r17-route-certificate.json",
    "tie_edge": GEN / "elkies-k3-h3-q4o164-q8o376-q12-fp-6-2-m12-10-rootless-certificate.json",
    "tie_route": GEN / "elkies-k3-h3-q4o208-q323-free-q4o1584-q4o164-q8o376-q12-fp-m12-10-pinned-r17-route-certificate.json",
}


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path):
    return str(path.relative_to(ROOT))


data = {name: json.loads(path.read_text()) for name, path in FILES.items()}
summary = data["neighbors"]["summaries"][0]
assert summary["q"] == 12
assert summary["mw_vector_cap"] == 50000
assert summary["mw_projection_representatives"] == 50000
assert summary["mw_pari_vector_count"] == 114347416
assert summary["primitive_neighbors"] == 61352
assert sum(
    item["orbit_count"]
    for item in summary["root_histogram"]
    if item["root_rank"] == 0 and item["mw_rank"] == 17
) == 28
assert not summary["mw_enumeration_complete"]
assert data["cost"]["status"] == "PASS_EXACT_MARKED_FRONTIER_EQUATION_COST_SCORING"
assert data["cost"]["retained_count"] == 28
assert data["frontier"]["status"] == "PASS_EXACT_ROOTLESS_P0_SECTION_WORD_FRONTIER"
assert data["cap_stability"]["status"] == "PASS_EXACT_BOUNDED_Q12_ROOTLESS_SET_STABLE_ACROSS_DOUBLED_CAP"
assert data["cap_stability"]["rootless_fibre_sets_equal"]

OLD_FIBRE = [6, 2, -16, 18, 2, 5, -11, 2, 13, 11, -6, -41, 6, 5, -14, 0, -3, -1, -2]
WINNER_FIBRE = [6, 2, -13, 10, 2, 5, -12, 2, 9, 7, -6, -22, 6, 2, -10, 0, -1, -1, -2]
TIE_FIBRE = [6, 2, -12, 10, 3, 5, -12, 2, 9, 5, -4, -22, 6, 2, -10, 0, -1, -1, -2]


def target(fibre):
    matches = [item for item in data["frontier"]["targets"] if item.get("fibre") == fibre]
    assert len(matches) == 1
    return matches[0]


old = target(OLD_FIBRE)
winner = target(WINNER_FIBRE)
tie = target(TIE_FIBRE)
old_four = old["best_four_P0_word_by_parent_a_minus_b"]
winner_mixed = winner["mixed_low_pole_words"]["three_P0_P0_P1"]["best_word"]
tie_mixed = tie["mixed_low_pole_words"]["three_P0_P0_P1"]["best_word"]
assert old_four["total_P_dot_O"] == 0
assert old_four["new_section_count"] == 4
assert old_four["q4o164_parent_a_minus_b_sum"] == 6
assert old_four["q4o164_parent_degree_sum"] == 8
for mixed in (winner_mixed, tie_mixed):
    assert mixed["total_P_dot_O"] == 1
    assert mixed["new_section_count"] == 3
    assert mixed["q4o164_parent_a_minus_b_sum"] == 3
    assert mixed["q4o164_parent_degree_sum"] == 8
    assert mixed["q4o164_parent_a_minus_b_max"] == 2
    assert mixed["q4o164_parent_degree_max"] == 3
assert winner_mixed["new_sections"] == tie_mixed["new_sections"]


def mixed_metrics(item):
    word = item["mixed_low_pole_words"]["three_P0_P0_P1"]["best_word"]
    return (
        word["q4o164_parent_a_minus_b_sum"],
        word["q4o164_parent_degree_sum"],
        word["q4o164_parent_a_minus_b_max"],
        word["q4o164_parent_degree_max"],
        len(word["known_section_correction"]),
    )


q12_targets = [item for item in data["frontier"]["targets"] if item["candidate_id"]["q"] == 12]
pareto = []
for item in q12_targets:
    metrics = mixed_metrics(item)
    dominated = any(
        all(left <= right for left, right in zip(mixed_metrics(other), metrics))
        and any(left < right for left, right in zip(mixed_metrics(other), metrics))
        for other in q12_targets
        if other is not item
    )
    if not dominated:
        pareto.append({
            "candidate_id_sample_local": item["candidate_id"],
            "fibre": item["fibre"],
            "fibre_fingerprint_sha256": item["fibre_fingerprint_sha256"],
            "metrics": list(metrics),
        })

for key in ("old_edge", "winner_edge", "tie_edge"):
    assert data[key]["status"] == "PASS_EXACT_MARKED_DEGREE_TWO_CANDIDATE_CERTIFICATE"
    assert data[key]["child"]["root_data"] == [0, 0, 1]
    assert data[key]["child"]["mw_rank"] == 17
for key in ("old_route", "winner_route", "tie_route"):
    assert data[key]["route"][-1]["child_root_data"] == [0, 0, 1]
    assert abs(data[key]["isometry_determinant"]) == 1

payload = {
    "schema": "elkies-k3.h3-q12-expanded-compiler-pareto.v1",
    "status": "PASS_EXACT_Q12_EXPANDED_COMPILER_PARETO_AND_PINNED_R17_ALTERNATIVE",
    "search_scope": {
        "q": 12,
        "old_fibre_degree": 2,
        "mw_vector_cap": 50000,
        "full_mw_vector_count": 114347416,
        "mw_enumeration_complete": False,
        "primitive_candidates_in_sample": 61352,
        "rootless_mw17_candidates_in_sample": 28,
        "claim": (
            "The compiler comparison is exact on the 28 rootless fibres found at cap 50,000; "
            "an independently stored cap-100,000 audit finds the identical rootless set. This "
            "is not a global q12 optimum theorem."
        ),
        "doubled_cap_stability": {
            "larger_cap": data["cap_stability"]["larger_cap"],
            "larger_primitive_candidate_count": data["cap_stability"]["larger_primitive_candidate_count"],
            "new_rootless_fibres": data["cap_stability"]["new_rootless_fibres_in_larger_sample"],
            "rootless_fibre_sets_equal": data["cap_stability"]["rootless_fibre_sets_equal"],
        },
    },
    "stable_identity_rule": (
        "Bounded-search orbit indices change when the MW-vector cap changes. Identify a route "
        "by its exact fibre vector and SHA-256 fingerprint, not by the sample-local orbit index."
    ),
    "promoted_preferred": {
        "sample_local_candidate_id": old["candidate_id"],
        "fibre": old["fibre"],
        "fibre_fingerprint_sha256": old["fibre_fingerprint_sha256"],
        "compiler": old_four,
        "edge_certificate": relative(FILES["old_edge"]),
        "route_certificate": relative(FILES["old_route"]),
        "promotion_scope": (
            "Preferred optional final equation-compiler target. All four new branches have "
            "current P.O=0; their q4/o164 parent degrees are (3,2,1,2) and parent "
            "a-b values are (2,2,1,1)."
        ),
    },
    "low_pole_alternative": {
        "sample_local_candidate_id": winner["candidate_id"],
        "fibre": winner["fibre"],
        "fibre_fingerprint_sha256": winner["fibre_fingerprint_sha256"],
        "compiler": winner_mixed,
        "compiler_word_profile": {
            "current_4A1_P_dot_O": [
                item["current_4A1_P_dot_O"] for item in winner_mixed["new_sections"]
            ],
            "q4o164_parent_a_minus_b": [
                item["q4o164_parent_a_minus_b"] for item in winner_mixed["new_sections"]
            ],
            "q4o164_parent_degrees": [
                item["q4o164_parent_degree"] for item in winner_mixed["new_sections"]
            ],
        },
        "edge_certificate": relative(FILES["winner_edge"]),
        "route_certificate": relative(FILES["winner_route"]),
        "endpoint": "rootless MW17 with exact determinant-one isometry to pinned R17",
        "retention_scope": (
            "Certified alternative from the retained sample. It uses one fewer new branch "
            "and a smaller parent pole proxy, but introduces a P.O=1 branch and four named "
            "corrections; q12/o5867 remains preferred for its all-polynomial P.O=0 compiler."
        ),
    },
    "exact_tie": {
        "sample_local_candidate_id": tie["candidate_id"],
        "fibre": tie["fibre"],
        "fibre_fingerprint_sha256": tie["fibre_fingerprint_sha256"],
        "same_three_new_sections": True,
        "compiler": tie_mixed,
        "edge_certificate": relative(FILES["tie_edge"]),
        "route_certificate": relative(FILES["tie_route"]),
    },
    "alternative_comparison": {
        "new_section_count": {"before": 4, "after": 3},
        "total_current_P_dot_O": {"before": 0, "after": 1},
        "q4o164_parent_a_minus_b_sum": {"before": 6, "after": 3},
        "q4o164_parent_degree_sum": {"before": 8, "after": 8},
        "q4o164_parent_a_minus_b_max": {"before": 2, "after": 2},
        "q4o164_parent_degree_max": {"before": 3, "after": 3},
        "known_section_correction_count": {"before": 2, "after": 4},
        "interpretation": (
            "The alternative exchanges one polynomial branch for a P.O=1 branch and halves the "
            "inherited-parent pole proxy at unchanged total parent degree, but it requires four "
            "named corrections. This does not supersede the all-P.O=0 q12/o5867 compiler."
        ),
    },
    "mixed_word_pareto_frontier": pareto,
    "proof_boundary": (
        "The section shells, physical component/affine gates, MW identities, marked U, nef and "
        "finite horizontal-wall audits, bidirectional unimodular NS transports, rootless data, "
        "and pinned-R17 endpoint identifications are exact. The rootless set is unchanged between "
        "the 50,000- and 100,000-vector q12 samples, but the full shell is not exhausted. The "
        "comparison retains q12/o5867 as the preferred optional "
        "all-P.O=0 compiler; characteristic-zero section equations and final RR pencils are not claimed."
    ),
    "inputs": {
        "paths": [relative(path) for path in FILES.values()],
        "sha256": {relative(path): sha256(path) for path in FILES.values()},
    },
}

output = GEN / "elkies-k3-h3-q4o164-q8o376-q12-expanded-compiler-pareto.json"
output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(output)
