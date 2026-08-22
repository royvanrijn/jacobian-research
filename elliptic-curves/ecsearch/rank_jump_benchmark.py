#!/usr/bin/env python3
"""Leakage-resistant ranking evaluation for exceptional quotient-rank data.

The input is a labelled *finite* manifest.  Each candidate supplies a lower
bound for its exceptional quotient rank, its four structural identities, and
the same predeclared feature vector on disjoint discovery and held-forward
prime bands.  The evaluator never reports classification accuracy.  Instead,
for each leave-one-group-out fold it reports the top-k enrichment of known
quotient-rank lower bounds +6, +8, and +10.

Labels must be backed by an independence/quotient certificate.  A numerical
rank or a raw trace score is not an admissible label.  This program does not
manufacture a labelled corpus and writes no pinned result by default.
"""

from __future__ import annotations

import argparse
import json
from math import isfinite, sqrt
from pathlib import Path
from statistics import fmean
from typing import Any, Iterable, Sequence


SCHEMA = "elliptic-curves.rank-jump-benchmark.v1"
GROUP_FIELDS = (
    "family",
    "root_shape",
    "parametrization_component",
    "quadratic_twist_class",
)
TARGET_JUMPS = (6, 8, 10)
REQUIRED_FEATURE_KINDS = (
    "family_residual_s0",
    "family_residual_s5",
    "conductor_scaled",
    "root_number",
    "predicted_local_conductor",
    "quotient_escape",
)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _numeric_feature_vector(
    candidate: dict[str, Any], key: str, feature_names: Sequence[str]
) -> None:
    features = candidate.get(key)
    _require(isinstance(features, dict), f"{candidate['id']}: missing {key}")
    _require(
        set(features) == set(feature_names),
        f"{candidate['id']}: {key} does not match declared feature_names",
    )
    for name in feature_names:
        value = features[name]
        _require(
            isinstance(value, (int, float)) and not isinstance(value, bool) and isfinite(value),
            f"{candidate['id']}: non-finite numeric feature {key}.{name}",
        )


def validate_manifest(manifest: dict[str, Any]) -> tuple[str, ...]:
    """Validate the small, explicit manifest format before fitting any score."""

    _require(manifest.get("schema") == SCHEMA, "unsupported benchmark schema")
    _require(
        manifest.get("discovery_primes_disjoint_from_held_forward") is True,
        "the manifest must attest disjoint discovery and held-forward prime bands",
    )
    feature_names = manifest.get("feature_names")
    _require(
        isinstance(feature_names, list)
        and feature_names
        and all(isinstance(name, str) and name for name in feature_names)
        and len(set(feature_names)) == len(feature_names),
        "feature_names must be a nonempty duplicate-free string list",
    )
    for kind in REQUIRED_FEATURE_KINDS:
        _require(
            any(kind in name for name in feature_names),
            f"missing required {kind} feature",
        )
    candidates = manifest.get("candidates")
    _require(isinstance(candidates, list) and candidates, "candidates must be nonempty")
    identifiers: set[str] = set()
    for candidate in candidates:
        _require(isinstance(candidate, dict), "candidate must be an object")
        identifier = candidate.get("id")
        _require(
            isinstance(identifier, str) and identifier and identifier not in identifiers,
            "candidate ids must be nonempty and unique",
        )
        identifiers.add(identifier)
        for field in GROUP_FIELDS:
            _require(
                isinstance(candidate.get(field), str) and candidate[field],
                f"{identifier}: missing structural group {field}",
            )
        jump = candidate.get("exceptional_quotient_rank_lower_bound")
        _require(
            isinstance(jump, int) and not isinstance(jump, bool) and jump >= 0,
            f"{identifier}: invalid exceptional quotient-rank lower bound",
        )
        _numeric_feature_vector(candidate, "discovery_features", feature_names)
        _numeric_feature_vector(candidate, "held_forward_features", feature_names)
    return tuple(feature_names)


def _feature_statistics(
    training: Sequence[dict[str, Any]], feature_names: Sequence[str]
) -> tuple[dict[str, float], dict[str, float]]:
    means = {
        name: fmean(float(row["discovery_features"][name]) for row in training)
        for name in feature_names
    }
    scales = {}
    for name in feature_names:
        variance = fmean(
            (float(row["discovery_features"][name]) - means[name]) ** 2
            for row in training
        )
        scales[name] = max(sqrt(variance), 1e-12)
    return means, scales


def fit_fixed_ranker(
    training: Sequence[dict[str, Any]],
    feature_names: Sequence[str],
    target_jump: int,
) -> dict[str, Any]:
    """Fit one fixed diagonal contrast score using discovery bands only.

    This is a ranker, not a calibrated probability model: each coefficient is
    the standardized positive-minus-control mean.  There is no feature or
    hyperparameter selection, so held-forward bands remain strictly unused by
    fitting.
    """

    positives = [
        row
        for row in training
        if row["exceptional_quotient_rank_lower_bound"] >= target_jump
    ]
    controls = [
        row
        for row in training
        if row["exceptional_quotient_rank_lower_bound"] < target_jump
    ]
    if not positives or not controls:
        raise ValueError("training fold has no positive or no control")
    means, scales = _feature_statistics(training, feature_names)
    weights = {}
    for name in feature_names:
        positive_mean = fmean(float(row["discovery_features"][name]) for row in positives)
        control_mean = fmean(float(row["discovery_features"][name]) for row in controls)
        weights[name] = (positive_mean - control_mean) / scales[name]
    return {
        "fit_band": "discovery",
        "target_exceptional_quotient_rank_lower_bound": target_jump,
        "training_positive_count": len(positives),
        "training_control_count": len(controls),
        "means": means,
        "scales": scales,
        "weights": weights,
    }


def score_candidate(
    model: dict[str, Any], candidate: dict[str, Any], feature_band: str
) -> float:
    """Score a candidate without accessing its label or any held-out row."""

    return sum(
        model["weights"][name]
        * (float(candidate[feature_band][name]) - model["means"][name])
        / model["scales"][name]
        for name in model["weights"]
    )


def top_k_enrichment(
    ranked: Sequence[dict[str, Any]], target_jump: int, k: int
) -> dict[str, Any]:
    """Report the predeclared rare-event metric, never accuracy."""

    selected = ranked[: min(k, len(ranked))]
    positives = sum(
        row["exceptional_quotient_rank_lower_bound"] >= target_jump
        for row in ranked
    )
    hits = [
        row["id"]
        for row in selected
        if row["exceptional_quotient_rank_lower_bound"] >= target_jump
    ]
    expected = len(selected) * positives / len(ranked)
    return {
        "requested_k": k,
        "evaluated_k": len(selected),
        "population_count": len(ranked),
        "population_positive_count": positives,
        "top_k_positive_count": len(hits),
        "top_k_positive_ids": hits,
        "random_ranking_expected_positive_count": expected,
        "enrichment_over_random": None if expected == 0 else len(hits) / expected,
    }


def _rank(
    test: Sequence[dict[str, Any]], model: dict[str, Any], feature_band: str
) -> list[dict[str, Any]]:
    return sorted(
        test,
        key=lambda row: (-score_candidate(model, row, feature_band), row["id"]),
    )


def evaluate_protocol(
    candidates: Sequence[dict[str, Any]],
    feature_names: Sequence[str],
    group_field: str,
    top_ks: Sequence[int],
) -> dict[str, Any]:
    """Run one leave-one-structural-group-out protocol without leakage."""

    _require(group_field in GROUP_FIELDS, "unknown structural group")
    groups = sorted({str(row[group_field]) for row in candidates})
    folds = []
    for group in groups:
        test = [row for row in candidates if row[group_field] == group]
        training = [row for row in candidates if row[group_field] != group]
        _require(test and training, f"degenerate {group_field} fold {group}")
        fold: dict[str, Any] = {
            "held_out_group": group,
            "held_out_candidate_ids": [row["id"] for row in test],
            "training_candidate_ids": [row["id"] for row in training],
            "leakage_check": {
                "training_excludes_entire_held_out_group": all(
                    row[group_field] != group for row in training
                ),
                "held_forward_features_used_for_fitting": False,
            },
            "targets": {},
        }
        for target in TARGET_JUMPS:
            try:
                model = fit_fixed_ranker(training, feature_names, target)
            except ValueError as error:
                fold["targets"][str(target)] = {
                    "status": "unavailable",
                    "reason": str(error),
                }
                continue
            target_record: dict[str, Any] = {
                "status": "evaluated",
                "model": model,
                "rankings": {},
            }
            for feature_band in ("discovery_features", "held_forward_features"):
                ranked = _rank(test, model, feature_band)
                target_record["rankings"][feature_band] = {
                    "ordered_candidate_ids": [row["id"] for row in ranked],
                    "top_k_enrichment": [
                        top_k_enrichment(ranked, target, k) for k in top_ks
                    ],
                }
            fold["targets"][str(target)] = target_record
        folds.append(fold)
    return {"group_field": group_field, "folds": folds}


def evaluate_manifest(manifest: dict[str, Any], top_ks: Iterable[int]) -> dict[str, Any]:
    """Evaluate all four structural holdout protocols."""

    feature_names = validate_manifest(manifest)
    normalized_top_ks = tuple(sorted(set(top_ks)))
    _require(
        normalized_top_ks and all(isinstance(k, int) and k > 0 for k in normalized_top_ks),
        "top-k values must be positive integers",
    )
    candidates = manifest["candidates"]
    return {
        "schema": "elliptic-curves.rank-jump-benchmark-result.v1",
        "objective": (
            "rank exceptional quotient jumps, not total rank; report top-k "
            "enrichment for certified quotient-rank lower bounds"
        ),
        "feature_names": list(feature_names),
        "top_k_values": list(normalized_top_ks),
        "prime_band_protocol": {
            "discovery_primes_disjoint_from_held_forward": True,
            "model_fit_band": "discovery_features",
            "evaluation_bands": ["discovery_features", "held_forward_features"],
        },
        "protocols": [
            evaluate_protocol(candidates, feature_names, field, normalized_top_ks)
            for field in GROUP_FIELDS
        ],
        "interpretation": {
            "proved": "the declared finite manifest and its leakage exclusions were evaluated exactly",
            "not_proved": "a rank upper bound, a probabilistic rank law, or a new exceptional quotient rank",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--top-k",
        default="1,5,10",
        help="comma-separated positive top-k values (default: 1,5,10)",
    )
    args = parser.parse_args()
    try:
        top_ks = tuple(int(value) for value in args.top_k.split(","))
    except ValueError as error:
        raise SystemExit("--top-k must be comma-separated integers") from error
    manifest = json.loads(args.input.read_text(encoding="utf-8"))
    result = evaluate_manifest(manifest, top_ks)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
