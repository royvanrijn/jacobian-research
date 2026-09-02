#!/usr/bin/env python3
"""Build the compact result for the prospective frozen-R17 shell search."""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results/elliptic-curves"
LOCAL = ROOT / "artifacts/local/elliptic-curves"
RANKING = GENERATED / "r17_frozen_nagao_shell_h10001_30000_v1.json"
COHORT = LOCAL / "r17-frozen-shell-h10001-30000-cohort.jsonl"
BASE_LABELS = LOCAL / "r17-frozen-shell-h10001-30000-bisection-labels.jsonl"
BASE_SUMMARY = LOCAL / "r17-frozen-shell-h10001-30000-bisection-labels-summary.json"
MEDIUM_LABELS = LOCAL / "r17-frozen-shell-h10001-30000-medium-labels.jsonl"
MEDIUM_SUMMARY = LOCAL / "r17-frozen-shell-h10001-30000-medium-summary.json"
DEEPEST_LABELS = LOCAL / "r17-frozen-shell-h10001-30000-deepest-labels.jsonl"
DEEPEST_SUMMARY = LOCAL / "r17-frozen-shell-h10001-30000-deepest-summary.json"
OUTPUT = GENERATED / "r17_frozen_nagao_shell_search_v1.json"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def read_jsonl(path: Path) -> list[dict[str, object]]:
    return [json.loads(line) for line in path.read_text().splitlines() if line]


def main() -> None:
    ranking = json.loads(RANKING.read_text())
    if ranking.get("status") != "PASS_COMPLETE_FROZEN_RULE_REGION_RANKING":
        raise SystemExit("the prospective population ranking is missing")
    cohort = read_jsonl(COHORT)
    base = read_jsonl(BASE_LABELS)
    medium = {row["parameter"]: row for row in read_jsonl(MEDIUM_LABELS)}
    deepest = {row["parameter"]: row for row in read_jsonl(DEEPEST_LABELS)}
    if len(cohort) != 384 or len(base) != len(cohort):
        raise SystemExit("the matched search cohort is incomplete")
    if [row["parameter"] for row in cohort] != [row["parameter"] for row in base]:
        raise SystemExit("cohort and base outcomes are not in the same order")

    lane_counts = Counter(row["selection_lanes"][0] for row in cohort)
    if lane_counts != {
        "frozen_weakest_block": 128,
        "ordinary_nagao_control": 128,
        "random_control": 128,
    }:
        raise SystemExit("the matched lanes changed")

    successes = []
    censored = []
    lane_gains = Counter()
    lane_successes = Counter()
    for selection, label in zip(cohort, base):
        lane = selection["selection_lanes"][0]
        rank = int(selection["lane_rank"])
        gain = label["outcomes"]["finite_quotient_gain_lower_bound"]
        strongest_bound = 199
        stronger = None
        if lane == "frozen_weakest_block" and rank <= 16:
            strongest_bound = 997
            stronger = deepest.get(label["parameter"])
        elif lane == "frozen_weakest_block" and rank <= 64:
            strongest_bound = 499
            stronger = medium.get(label["parameter"])
        if stronger is not None:
            stronger_gain = stronger["outcomes"]["finite_quotient_gain_lower_bound"]
            if stronger_gain != gain:
                raise AssertionError("a stronger frozen-rank replay changed the gain")
        if gain is None:
            censored.append(
                {
                    "parameter": label["parameter"],
                    "lane": lane,
                    "lane_rank": rank,
                    "frozen_population_rank": selection["frozen_population_rank"],
                    "split_bisection_count": label["bisection_census"][
                        "nonzero_split_bisection_count"
                    ],
                    "blocker": label["outcomes"]["censoring_or_blocker"],
                }
            )
            continue
        if int(gain) <= 0:
            continue
        lane_gains[lane] += int(gain)
        lane_successes[lane] += 1
        hits = label["bisection_census"]["hits"]
        successes.append(
            {
                "parameter": label["parameter"],
                "projective_pair": label["projective_pair"],
                "lane": lane,
                "lane_rank": rank,
                "frozen_population_rank": selection["frozen_population_rank"],
                "quotient_rank_lower_bound_beyond_generic_17": int(gain),
                "specialized_rank_lower_bound": 17 + int(gain),
                "split_bisection_count": label["bisection_census"][
                    "nonzero_split_bisection_count"
                ],
                "independent_bisection_labels": label["finite_quotient_audit"][
                    "independent_hit_labels"
                ],
                "minimum_exceptional_x_naive_height": label["outcomes"][
                    "minimum_exceptional_naive_height"
                ],
                "strongest_frozen_depth_reduction_prime_bound": strongest_bound,
                "all_displayed_points_verified_exactly": all(
                    hit["exact_verification"]["both_branches_on_source_fibre"]
                    and hit["exact_verification"]["sum_of_branches_equals_stored_trace"]
                    and hit["exact_verification"][
                        "stored_trace_equals_published_basis_word"
                    ]
                    for hit in hits
                ),
            }
        )

    if len(successes) != 7 or sum(lane_gains.values()) != 7 or len(censored) != 1:
        raise SystemExit("the prospective exact outcome counts changed")
    base_summary = json.loads(BASE_SUMMARY.read_text())
    medium_summary = json.loads(MEDIUM_SUMMARY.read_text())
    deepest_summary = json.loads(DEEPEST_SUMMARY.read_text())
    output = {
        "schema": "elliptic-curves.r17-frozen-nagao-shell-search.v1",
        "status": "PASS_PROSPECTIVE_R17_SHELL_ATLAS_GAINS_NO_DEEP_POINT_SEARCH",
        "population": {
            "definition": ranking["search"],
            "complete_primitive_parameter_count": ranking["population_count"],
            "all_parameters_scored_without_presieve": True,
            "ranking_artifact": relative(RANKING),
            "ranking_artifact_sha256": digest(RANKING),
        },
        "frozen_rule": {
            "reference": ranking["frozen_reference"],
            "ranking": ranking["frozen_ranking"],
            "changed_after_prospective_shell_opened": False,
        },
        "matched_lane_outcomes": {
            lane: {
                "rows": lane_counts[lane],
                "rows_with_certified_gain": lane_successes[lane],
                "certified_quotient_gain_sum": lane_gains[lane],
                "observed_success_fraction": lane_successes[lane] / lane_counts[lane],
            }
            for lane in (
                "frozen_weakest_block",
                "ordinary_nagao_control",
                "random_control",
            )
        },
        "certified_new_directions": successes,
        "censored_split_rows": censored,
        "progressive_depth": {
            "base": {
                "rows": 384,
                "complete_preexisting_bisection_atlas_size": 39_120,
                "finite_reduction_prime_bound": 199,
                "summary": relative(BASE_SUMMARY),
                "summary_sha256": digest(BASE_SUMMARY),
            },
            "medium": {
                "selection": "frozen target ranks 17 through 64 only",
                "rows": 48,
                "finite_reduction_prime_bound": 499,
                "summary": relative(MEDIUM_SUMMARY),
                "summary_sha256": digest(MEDIUM_SUMMARY),
            },
            "deepest": {
                "selection": "frozen target ranks 1 through 16 only",
                "rows": 16,
                "finite_reduction_prime_bound": 997,
                "summary": relative(DEEPEST_SUMMARY),
                "summary_sha256": digest(DEEPEST_SUMMARY),
            },
            "allocation_depended_only_on_frozen_rank": True,
            "point_mechanism": (
                "complete evaluation of the preexisting certified 39,120-bisection atlas"
            ),
        },
        "unrestricted_point_search": {
            "executed": False,
            "reason": (
                "No candidate has a completed same-fibre residual 2-Selmer quotient "
                "of dimension at least 15; all retained ratpoints, eclib, slope-chart, "
                "and two-cover entry points therefore fail closed. Magma is absent on "
                "the current host."
            ),
            "interpretation": (
                "This is a policy/backend boundary, not a rank upper bound or a "
                "negative point-search outcome."
            ),
        },
        "inputs": {
            relative(path): digest(path)
            for path in (
                COHORT,
                BASE_LABELS,
                BASE_SUMMARY,
                MEDIUM_LABELS,
                MEDIUM_SUMMARY,
                DEEPEST_LABELS,
                DEEPEST_SUMMARY,
            )
        },
        "summary_replays": {
            "base_outcomes": base_summary["outcomes"],
            "medium_outcomes": medium_summary["outcomes"],
            "deepest_outcomes": deepest_summary["outcomes"],
        },
        "generation": {
            "script": relative(Path(__file__)),
            "script_sha256": digest(Path(__file__)),
            "command": (
                ".venv/bin/python "
                "elliptic-curves/scripts/summarize_r17_frozen_nagao_shell.py"
            ),
        },
        "proof_boundary": [
            "Each positive row has an exact specialized point and an exact finite-reduction quotient-rank lower bound.",
            "The seven gains occur on seven different fibres; they do not combine into one rank-24 fibre.",
            "A zero atlas outcome is not a rank-17 claim, and the one censored row has no promoted gain.",
            "Nagao scores and population ranks remain heuristics.",
            "No unrestricted rational-point search or residual-Selmer computation was completed.",
        ],
    }
    OUTPUT.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(
        f"R17FROZENSHELLRESULT|population={ranking['population_count']}|"
        f"gains={len(successes)}|censored={len(censored)}|output={OUTPUT}"
    )


if __name__ == "__main__":
    main()
