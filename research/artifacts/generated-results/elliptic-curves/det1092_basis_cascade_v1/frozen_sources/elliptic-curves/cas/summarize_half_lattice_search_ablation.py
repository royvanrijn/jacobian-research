#!/usr/bin/env python3
"""Build the compact certificate for the equal-budget chart-selection ablation."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
from statistics import mean


ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / "artifacts/generated-results/elliptic-curves"
DEVELOPMENT_BLIND = ART / "half_lattice_search_ablation_r17_development_blind_v1.json"
DEVELOPMENT_VERIFIED = ART / "half_lattice_search_ablation_r17_development_verification_v1.json"
HOLDOUT_BLIND = ART / "half_lattice_search_ablation_rank29_holdout_blind_v1.json"
HOLDOUT_VERIFIED = ART / "half_lattice_search_ablation_rank29_holdout_verification_v1.json"
OUTPUT = ART / "half_lattice_search_ablation_summary_v1.json"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def arm_compact(arm: dict) -> dict:
    return {
        "id": arm["id"],
        "class_count": arm["class_count"],
        "exact_quotient_rank_over_Q": arm["exact_quotient_rank_over_Q"],
        "exact_quotient_rank_mod2": arm["exact_quotient_rank_mod2"],
        "rank_normalized_per_43_covers": arm["rank_normalized_per_43_covers"],
        "cover_cpu_seconds": arm["cover_cpu_seconds"],
        "quotient_rank_per_cpu_second": arm["quotient_rank_per_cpu_second"],
        "generic_deepest_overlap": arm["generic_deepest_overlap"],
        "specialized_deepest_overlap": arm["specialized_deepest_overlap"],
    }


def compact_case(case: dict) -> dict:
    arms = [arm_compact(arm) for arm in case["arms"]]
    random_arms = [arm for arm in arms if arm["id"].startswith("random43-")]
    generic = next(arm for arm in arms if arm["id"] == "generic-deepest43")
    specialized = next(arm for arm in arms if arm["id"] == "specialized-deepest43")
    return {
        "label": case["label"],
        "parameter": case["parameter"],
        "target_quotient_dimension": case["target_quotient_dimension"],
        "generic_specialized_intersection_count": case["ranking"][
            "generic_specialized_intersection_count"
        ],
        "maximum_generic_norm": case["ranking"]["maximum_generic_norm"],
        "maximum_generic_norm_stratum_count": case["ranking"].get(
            "maximum_generic_norm_stratum_count", 43
        ),
        "specialized_top43_stable_at_scales_1e5_and_1e6": case["ranking"][
            "specialized"
        ]["top43_set_stable"],
        "arms": arms,
        "random43_summary": {
            "ranks": [arm["exact_quotient_rank_over_Q"] for arm in random_arms],
            "mean_rank": mean(
                arm["exact_quotient_rank_over_Q"] for arm in random_arms
            ),
            "best_rank": max(
                arm["exact_quotient_rank_over_Q"] for arm in random_arms
            ),
            "generic_deepest_minus_random_mean": generic[
                "exact_quotient_rank_over_Q"
            ]
            - mean(arm["exact_quotient_rank_over_Q"] for arm in random_arms),
            "specialized_deepest_minus_random_mean": specialized[
                "exact_quotient_rank_over_Q"
            ]
            - mean(arm["exact_quotient_rank_over_Q"] for arm in random_arms),
        },
    }


def phase_summary(cases: list[dict]) -> dict:
    ids = [arm["id"] for arm in cases[0]["arms"]]
    arm_rows = []
    for arm_id in ids:
        arms = [
            next(arm for arm in case["arms"] if arm["id"] == arm_id)
            for case in cases
        ]
        total_rank = sum(arm["exact_quotient_rank_over_Q"] for arm in arms)
        total_cpu = sum(arm["cover_cpu_seconds"] for arm in arms)
        arm_rows.append(
            {
                "id": arm_id,
                "case_count": len(cases),
                "total_exact_quotient_rank_over_Q": total_rank,
                "mean_exact_quotient_rank_over_Q": total_rank / len(cases),
                "mean_rank_normalized_per_43_covers": mean(
                    arm["rank_normalized_per_43_covers"] for arm in arms
                ),
                "pooled_quotient_rank_per_cpu_second": total_rank / total_cpu,
            }
        )
    random_rows = [row for row in arm_rows if row["id"].startswith("random43-")]
    random_arms = [
        arm
        for case in cases
        for arm in case["arms"]
        if arm["id"].startswith("random43-")
    ]
    random_total_rank = sum(
        arm["exact_quotient_rank_over_Q"] for arm in random_arms
    )
    random_total_cpu = sum(arm["cover_cpu_seconds"] for arm in random_arms)
    return {
        "arms": arm_rows,
        "random43_across_all_five_arms": {
            "mean_of_arm_mean_ranks": mean(
                row["mean_exact_quotient_rank_over_Q"] for row in random_rows
            ),
            "best_arm_mean_rank": max(
                row["mean_exact_quotient_rank_over_Q"] for row in random_rows
            ),
            "pooled_quotient_rank_per_cpu_second": random_total_rank
            / random_total_cpu,
        },
    }


def main() -> None:
    paths = (
        DEVELOPMENT_BLIND,
        DEVELOPMENT_VERIFIED,
        HOLDOUT_BLIND,
        HOLDOUT_VERIFIED,
    )
    documents = [json.loads(path.read_text()) for path in paths]
    development_blind, development, holdout_blind, holdout = documents
    if development_blind["status"] != "PASS_BLIND_ABLATION_SEARCH":
        raise ValueError("development blind search is not complete")
    if holdout_blind["status"] != "PASS_BLIND_ABLATION_SEARCH":
        raise ValueError("holdout blind search is not complete")
    if any(
        document["status"] != "PASS_EXACT_PUBLIC_QUOTIENT_ABLATION"
        for document in (development, holdout)
    ):
        raise ValueError("an exact fixture verification is not passing")
    if development_blind["declared_budget"] != holdout_blind["declared_budget"]:
        raise ArithmeticError("development and holdout budgets differ")
    if development_blind["arm_definition"] != holdout_blind["arm_definition"]:
        raise ArithmeticError("development and holdout arm definitions differ")
    if development["phase_boundary"]["blind_artifact_sha256_before_fixture_import"] != digest(
        DEVELOPMENT_BLIND
    ):
        raise ArithmeticError("development verifier did not hash these blind bytes")
    if holdout["phase_boundary"]["blind_artifact_sha256_before_fixture_import"] != digest(
        HOLDOUT_BLIND
    ):
        raise ArithmeticError("holdout verifier did not hash these blind bytes")
    expected_ids = [arm["id"] for arm in development["results"][0]["arms"]]
    for phase in (development, holdout):
        for case in phase["results"]:
            if [arm["id"] for arm in case["arms"]] != expected_ids:
                raise ArithmeticError("arm order changed")
            for arm in case["arms"]:
                if arm["exact_quotient_rank_over_Q"] != arm["exact_quotient_rank_mod2"]:
                    raise ArithmeticError("Q and F2 quotient ranks differ")
                if arm["id"] != "deep-union" and arm["class_count"] != 43:
                    raise ArithmeticError("a fixed-size arm stopped having 43 covers")
    blind_cases = development_blind["results"] + holdout_blind["results"]
    random_reference = {
        arm["id"]: arm["masks"]
        for arm in blind_cases[0]["arms"]
        if arm["id"].startswith("random43-")
    }
    if len({mask for masks in random_reference.values() for mask in masks}) != 5 * 43:
        raise ArithmeticError("the five random arms stopped being disjoint")
    for case in blind_cases:
        observed = {
            arm["id"]: arm["masks"]
            for arm in case["arms"]
            if arm["id"].startswith("random43-")
        }
        if observed != random_reference:
            raise ArithmeticError("a deterministic random arm changed between fibres")

    development_cases = [compact_case(case) for case in development["results"]]
    holdout_cases = [compact_case(case) for case in holdout["results"]]
    payload = {
        "schema": "elliptic-curves.half-lattice-search-ablation-summary.v1",
        "status": "PASS_EQUAL_BUDGET_BLIND_ABLATION",
        "evidence_labels": {
            "exact": (
                "All reported quotient ranks and displayed relations are exact in the "
                "published point subgroups; Q-rank equals the retained mod-2 rank."
            ),
            "numerical": (
                "Specialized CVP order is numerical, with its top-43 set reproduced at "
                "height-matrix rounding scales 10^5 and 10^6."
            ),
            "bounded": (
                "Each cover received one height-100000, 15-second search; every miss is "
                "only a bounded-search miss."
            ),
            "heuristic": (
                "The comparisons measure enrichment and search efficiency, not a theorem "
                "that depth determines solubility or rank."
            ),
        },
        "blind_boundary": {
            "development_blind_sha256": digest(DEVELOPMENT_BLIND),
            "holdout_blind_sha256_before_fixture_import": holdout["phase_boundary"][
                "blind_artifact_sha256_before_fixture_import"
            ],
            "holdout_search_source_sha256": holdout_blind["input_hashes"][
                "elliptic-curves/cas/replay_half_lattice_search_ablation.sage"
            ],
            "source_history": holdout_blind["source_history"],
            "public_holdout_points_loaded_only_after_blind_artifact_frozen": holdout[
                "phase_boundary"
            ]["public_points_loaded_only_after_blind_artifact_frozen"],
        },
        "identical_budget": holdout_blind["declared_budget"],
        "arm_definition": holdout_blind["arm_definition"],
        "development": {
            "cases": development_cases,
            "aggregate": phase_summary(development_cases),
        },
        "holdout": {
            "cases": holdout_cases,
            "aggregate": phase_summary(holdout_cases),
        },
        "interpretation": [
            (
                "The generic deepest 43 are genuinely enriched: on the three sealed +12 "
                "holdouts their mean exact quotient gain is 8.33, versus 3.20 averaged "
                "over the five deterministic random arms."
            ),
            (
                "The generic deepest arm also has the best aggregate holdout CPU "
                "efficiency. The specialized deepest arm is useful but is weaker than "
                "the generic arm on these holdouts."
            ),
            (
                "The deep union maximizes raw recovery only when complementary directions "
                "occur, but its 64--86 covers make it less efficient per 43 covers and per "
                "CPU second than the generic deepest arm."
            ),
            (
                "Every shallowest arm recovers quotient dimension zero across all seven "
                "fibres. This rejects the explanation that merely trying 43 reduced "
                "quartic charts accounts for the result."
            ),
            (
                "Median and random arms can recover many directions, especially on the "
                "easier development fibres, so minimization/reduction and multi-chart "
                "search remain substantial parts of the mechanism."
            ),
            (
                "The data favor generic half-lattice geometry as the primary selector, "
                "not specialization-specific CVP. This is a calibrated bounded-search "
                "finding from five fixed random arms, not a distributional theorem."
            ),
        ],
        "input_hashes": {
            str(path.relative_to(ROOT)): digest(path) for path in paths
        }
        | {str(Path(__file__).resolve().relative_to(ROOT)): digest(Path(__file__).resolve())},
        "reproducing_command": "python3 elliptic-curves/cas/summarize_half_lattice_search_ablation.py",
    }
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(f"HALFABLATESUMMARY|status={payload['status']}|output={OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
