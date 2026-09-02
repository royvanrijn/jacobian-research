"""Evaluate frozen cheap-score baselines on the recovered Fermigier corpus."""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
import gzip
import json
from math import ceil
from pathlib import Path
from statistics import median
from typing import Any, Iterable


CONFIG_SCHEMA = "elliptic-curves.fermigier-baseline-rankers.v1"
RESULT_SCHEMA = "elliptic-curves.fermigier-baseline-evaluation.v1"


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_config(config: dict[str, Any]) -> None:
    if config.get("schema") != CONFIG_SCHEMA:
        raise ValueError("unsupported Fermigier baseline configuration")
    if config.get("evaluation_role") != "retrospective_development_not_prospective_holdout":
        raise ValueError("baseline evaluation role must remain retrospective development")
    positives = config.get("admitted_positive_ids")
    if not isinstance(positives, list) or len(positives) != 2 or len(set(positives)) != 2:
        raise ValueError("exactly two distinct admitted positives are required")
    budgets = config.get("budgets")
    if not isinstance(budgets, list) or budgets != sorted(set(budgets)) or any(
        not isinstance(value, int) or isinstance(value, bool) or value <= 0
        for value in budgets
    ):
        raise ValueError("budgets must be distinct increasing positive integers")
    rankers = config.get("rankers")
    if not isinstance(rankers, list) or not rankers:
        raise ValueError("rankers must be a nonempty list")
    identifiers: set[str] = set()
    for ranker in rankers:
        identifier = ranker.get("id")
        if not isinstance(identifier, str) or not identifier or identifier in identifiers:
            raise ValueError("ranker ids must be nonempty and unique")
        identifiers.add(identifier)
        if ranker.get("direction") not in {"higher", "lower"}:
            raise ValueError(f"{identifier}: unknown direction")
        for field in ("feature_source", "field", "role", "selection", "leakage_boundary"):
            if not isinstance(ranker.get(field), str) or not ranker[field]:
                raise ValueError(f"{identifier}: missing {field}")
        initial = ranker.get("initial_population_count")
        materialized = ranker.get("materialized_population_count")
        if not isinstance(initial, int) or not isinstance(materialized, int):
            raise ValueError(f"{identifier}: invalid population counts")
        if not 0 < materialized <= initial:
            raise ValueError(f"{identifier}: inconsistent population counts")


def _candidate_sort_key(candidate: dict[str, Any], direction: str) -> tuple[Any, ...]:
    value = float(candidate["value"])
    primary = -value if direction == "higher" else value
    numerator, denominator = candidate["projective_pair_T"]
    return primary, max(abs(numerator), denominator), denominator, numerator


def _budget_metrics(
    ranked: list[dict[str, Any]],
    admitted_positive_ids: list[str],
    budgets: Iterable[int],
) -> list[dict[str, Any]]:
    present = {
        row["positive_id"]
        for row in ranked
        if row["positive_id"] is not None
    }
    answer = []
    for requested in budgets:
        evaluated = min(requested, len(ranked))
        hits = [
            row["positive_id"]
            for row in ranked[:evaluated]
            if row["positive_id"] is not None
        ]
        expected_present = evaluated * len(present) / len(ranked)
        answer.append(
            {
                "requested_budget": requested,
                "evaluated_budget": evaluated,
                "known_positive_hits": hits,
                "recall_of_all_admitted_positives": len(hits) / len(admitted_positive_ids),
                "recall_of_present_positives": None if not present else len(hits) / len(present),
                "random_order_expected_present_positive_hits": expected_present,
                "enrichment_over_random_for_present_positives": (
                    None if expected_present == 0 else len(hits) / expected_present
                ),
            }
        )
    return answer


def _stratum_metrics(
    ranked: list[dict[str, Any]],
    predicate,
    budgets: Iterable[int],
) -> dict[str, Any]:
    positions = [index for index, row in enumerate(ranked, 1) if predicate(row)]
    return {
        "count": len(positions),
        "best_position": min(positions, default=None),
        "median_position": None if not positions else median(positions),
        "budget_hits": {
            str(budget): sum(position <= min(budget, len(ranked)) for position in positions)
            for budget in budgets
        },
    }


def evaluate_ranker(
    ranker: dict[str, Any],
    candidates: list[dict[str, Any]],
    admitted_positive_ids: list[str],
    budgets: list[int],
) -> dict[str, Any]:
    ranked = sorted(candidates, key=lambda row: _candidate_sort_key(row, ranker["direction"]))
    positive_positions: dict[str, dict[str, Any]] = {}
    for position, row in enumerate(ranked, 1):
        positive_id = row["positive_id"]
        if positive_id is None:
            continue
        same_value = sum(other["value"] == row["value"] for other in ranked)
        positive_positions[positive_id] = {
            "position_one_based": position,
            "value": row["value"],
            "top_fraction_of_materialized_cohort": position / len(ranked),
            "materialized_search_space_reduction_factor": len(ranked) / position,
            "operational_full_point_search_fraction_of_initial_population": (
                position / ranker["initial_population_count"]
            ),
            "operational_full_point_search_reduction_factor": (
                ranker["initial_population_count"] / position
            ),
            "equal_primary_value_count": same_value,
        }
    missing = [value for value in admitted_positive_ids if value not in positive_positions]
    outcome_strata = {
        "legacy_rank_at_least_13_uncertified": _stratum_metrics(
            ranked,
            lambda row: row["positive_id"] is None and (row["legacy_rank"] or -1) >= 13,
            budgets,
        ),
        "legacy_rank_at_least_15_uncertified": _stratum_metrics(
            ranked,
            lambda row: row["positive_id"] is None and (row["legacy_rank"] or -1) >= 15,
            budgets,
        ),
        "legacy_rank_at_least_17_uncertified": _stratum_metrics(
            ranked,
            lambda row: row["positive_id"] is None and (row["legacy_rank"] or -1) >= 17,
            budgets,
        ),
        "positive_quartic_point_count": _stratum_metrics(
            ranked,
            lambda row: row["positive_id"] is None and (row["quartic_point_count"] or 0) > 0,
            budgets,
        ),
        "structured_numerical_rank_at_least_13": _stratum_metrics(
            ranked,
            lambda row: row["positive_id"] is None and (row["numerical_rank"] or -1) >= 13,
            budgets,
        ),
    }
    return {
        "id": ranker["id"],
        "role": ranker["role"],
        "feature_source": ranker["feature_source"],
        "field": ranker["field"],
        "direction": ranker["direction"],
        "initial_population_count": ranker["initial_population_count"],
        "materialized_population_count": len(ranked),
        "materialized_fraction_of_initial_population": (
            len(ranked) / ranker["initial_population_count"]
        ),
        "population_count_matches_declaration": (
            len(ranked) == ranker["materialized_population_count"]
        ),
        "selection": ranker["selection"],
        "leakage_boundary": ranker["leakage_boundary"],
        "positive_coverage": {
            "present_count": len(positive_positions),
            "admitted_count": len(admitted_positive_ids),
            "missing_positive_ids": missing,
        },
        "positive_positions": positive_positions,
        "budget_metrics": _budget_metrics(ranked, admitted_positive_ids, budgets),
        "censored_outcome_overlap": {
            "interpretation": "descriptive only; searches were selected adaptively and missing outcomes are not negatives",
            "rows_with_any_recorded_outcome": sum(
                row["positive_id"] is None
                and (
                    row["legacy_rank"] is not None
                    or row["quartic_point_count"] is not None
                    or row["numerical_rank"] is not None
                )
                for row in ranked
            ),
            "strata": outcome_strata,
        },
    }


def load_candidates(
    corpus_path: Path,
    rankers: list[dict[str, Any]],
    admitted_positive_ids: list[str],
) -> tuple[dict[str, list[dict[str, Any]]], dict[str, int]]:
    needed = {ranker["feature_source"] for ranker in rankers}
    answer = {source: [] for source in needed}
    label_counts: Counter[str] = Counter()
    positive_ids: set[str] = set()
    with gzip.open(corpus_path, "rt", encoding="utf-8") as stream:
        for line in stream:
            row = json.loads(line)
            label = row["label"]
            label_counts[label["state"]] += 1
            positive_id = None
            if label["state"] == "certified_positive":
                if row["split"] != "positive_holdout":
                    raise ValueError("a certified positive escaped the positive holdout")
                positive_id = label["certificate"]["id"]
                positive_ids.add(positive_id)
            outcome = row["outcome_summary"]
            for source, features in row["cheap_features"].items():
                if source not in needed:
                    continue
                answer[source].append(
                    {
                        "parameter": row["parameter"]["normalized_T"],
                        "projective_pair_T": row["parameter"]["projective_pair_T"],
                        "positive_id": positive_id,
                        "features": features,
                        "legacy_rank": outcome[
                            "maximum_legacy_reported_rank_floor_uncertified"
                        ],
                        "quartic_point_count": outcome[
                            "maximum_reported_quartic_point_count"
                        ],
                        "numerical_rank": outcome["maximum_stable_numerical_rank"],
                    }
                )
    if positive_ids != set(admitted_positive_ids):
        raise ValueError("the corpus positive set does not match the frozen configuration")
    return answer, dict(label_counts)


def evaluate(config: dict[str, Any], repository_root: Path) -> dict[str, Any]:
    validate_config(config)
    summary_path = repository_root / config["corpus"]["summary"]
    corpus_path = repository_root / config["corpus"]["path"]
    summary = json.loads(summary_path.read_text())
    expected_hash = summary["output"]["sha256"]
    if file_sha256(corpus_path) != expected_hash:
        raise ValueError("the corpus hash does not match its summary")
    by_source, label_counts = load_candidates(
        corpus_path, config["rankers"], config["admitted_positive_ids"]
    )
    results = []
    for ranker in config["rankers"]:
        candidates = []
        for row in by_source[ranker["feature_source"]]:
            value = row["features"].get(ranker["field"])
            if not isinstance(value, (int, float)) or isinstance(value, bool):
                raise ValueError(f"{ranker['id']}: missing numeric field {ranker['field']}")
            candidates.append({**row, "value": value})
        result = evaluate_ranker(
            ranker, candidates, config["admitted_positive_ids"], config["budgets"]
        )
        if not result["population_count_matches_declaration"]:
            raise ValueError(f"{ranker['id']}: materialized population count changed")
        results.append(result)
    return {
        "schema": RESULT_SCHEMA,
        "status": "RETROSPECTIVE_BASELINES_EVALUATED_NO_MODEL_FIT",
        "evaluation_role": config["evaluation_role"],
        "corpus": {
            "path": config["corpus"]["path"],
            "sha256": expected_hash,
            "summary": config["corpus"]["summary"],
            "summary_sha256": file_sha256(summary_path),
            "row_count": summary["output"]["row_count"],
            "label_counts": label_counts,
        },
        "budgets": config["budgets"],
        "rankers": results,
        "interpretation": config["interpretation"],
        "proof_boundary": [
            "No baseline is fitted to a rank or point-search label.",
            "All results are retrospective development diagnostics, not prospective holdouts.",
            "Ranks are exact only inside the materialized finalist cohort named by each baseline.",
            "Censored-outcome overlap is selection-biased and is not classification accuracy.",
            "Only the two attached certificates establish exceptional quotient-rank lower bounds.",
        ],
    }
