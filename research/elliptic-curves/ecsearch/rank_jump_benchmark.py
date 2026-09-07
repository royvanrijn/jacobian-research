#!/usr/bin/env python3
"""Leakage-resistant ranking evaluation for exceptional quotient-rank data.

Two input protocols are supported.  The original finite candidate manifest
fits a fixed contrast score on disjoint discovery and held-forward prime
bands.  The laboratory registry inventories certified within-family labels
and evaluates complete pre-point-search ranking artifacts by search budget.
Neither protocol reports classification accuracy.

Labels must be backed by an independence/quotient certificate.  A numerical
rank or a raw trace score is not an admissible label.  This program does not
manufacture a labelled corpus and writes no pinned result by default.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from math import isfinite, sqrt
from pathlib import Path
from statistics import fmean
from typing import Any, Iterable, Sequence


SCHEMA = "elliptic-curves.rank-jump-benchmark.v1"
LAB_SCHEMA = "elliptic-curves.rank-jump-laboratory.v1"
LAB_RESULT_SCHEMA = "elliptic-curves.rank-jump-laboratory-result.v1"
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
ALLOWED_LAB_READINESS = {
    "ranked_development_corpus",
    "labels_only",
    "pending_quotient_audit",
    "excluded_pending_lineage",
}
ALLOWED_EVALUATION_ROLES = {"development", "held_out"}
ALLOWED_FEATURE_STAGES = {
    "parameter_only",
    "local_arithmetic",
    "cover_incidence",
}


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


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _validate_source(
    source: dict[str, Any], repository_root: Path, context: str
) -> Path:
    _require(isinstance(source, dict), f"{context}: source must be an object")
    relative = source.get("path")
    digest = source.get("sha256")
    _require(
        isinstance(relative, str) and relative and not Path(relative).is_absolute(),
        f"{context}: source path must be nonempty and repository-relative",
    )
    _require(
        isinstance(digest, str)
        and len(digest) == 64
        and all(character in "0123456789abcdef" for character in digest),
        f"{context}: source sha256 must be 64 lowercase hexadecimal characters",
    )
    relative_path = Path(relative)
    _require(".." not in relative_path.parts, f"{context}: source path escapes repository")
    path = repository_root / relative_path
    _require(path.is_file(), f"{context}: missing source {relative}")
    _require(_sha256(path) == digest, f"{context}: stale source hash for {relative}")
    return path


def validate_lab_manifest(
    manifest: dict[str, Any], repository_root: Path
) -> None:
    """Validate the laboratory inventory and every proof/search dependency.

    Failed point searches are deliberately not binary negative labels.  They
    may appear only as censored control sources in the family inventory.
    """

    _require(manifest.get("schema") == LAB_SCHEMA, "unsupported laboratory schema")
    budget_fractions = manifest.get("budget_fractions")
    _require(
        isinstance(budget_fractions, list) and budget_fractions,
        "budget_fractions must be a nonempty list",
    )
    parsed_fractions: list[Fraction] = []
    for value in budget_fractions:
        _require(isinstance(value, str), "budget fractions must be exact strings")
        try:
            fraction = Fraction(value)
        except (ValueError, ZeroDivisionError) as error:
            raise ValueError(f"invalid budget fraction {value}") from error
        _require(0 < fraction <= 1, f"budget fraction outside (0,1]: {value}")
        parsed_fractions.append(fraction)
    _require(
        parsed_fractions == sorted(set(parsed_fractions)),
        "budget_fractions must be strictly increasing",
    )

    families = manifest.get("families")
    _require(isinstance(families, list) and families, "families must be nonempty")
    family_ids: set[str] = set()
    positive_ids: set[str] = set()
    for family in families:
        _require(isinstance(family, dict), "family must be an object")
        family_id = family.get("id")
        _require(
            isinstance(family_id, str) and family_id and family_id not in family_ids,
            "family ids must be nonempty and unique",
        )
        family_ids.add(family_id)
        _require(
            family.get("readiness") in ALLOWED_LAB_READINESS,
            f"{family_id}: unknown readiness",
        )
        generic_rank = family.get("proved_arithmetic_generic_rank")
        if generic_rank is not None:
            _require(
                isinstance(generic_rank, int)
                and not isinstance(generic_rank, bool)
                and generic_rank >= 0,
                f"{family_id}: invalid proved arithmetic generic rank",
            )

        positives = family.get("known_positives", [])
        _require(isinstance(positives, list), f"{family_id}: known_positives must be a list")
        local_positive_ids: set[str] = set()
        for positive in positives:
            _require(isinstance(positive, dict), f"{family_id}: positive must be an object")
            positive_id = positive.get("id")
            _require(
                isinstance(positive_id, str)
                and positive_id
                and positive_id not in positive_ids,
                "positive ids must be nonempty and globally unique",
            )
            positive_ids.add(positive_id)
            local_positive_ids.add(positive_id)
            _require(
                isinstance(positive.get("parameter"), str) and positive["parameter"],
                f"{positive_id}: missing canonical parameter",
            )
            rank = positive.get("certified_rank_lower_bound")
            gain = positive.get("exceptional_quotient_rank_lower_bound")
            _require(
                isinstance(rank, int) and not isinstance(rank, bool) and rank >= 0,
                f"{positive_id}: invalid certified rank lower bound",
            )
            _require(
                isinstance(gain, int) and not isinstance(gain, bool) and gain > 0,
                f"{positive_id}: invalid exceptional quotient-rank lower bound",
            )
            _require(
                generic_rank is not None and gain <= rank - generic_rank,
                f"{positive_id}: quotient gain exceeds rank lower bound minus exact generic rank",
            )
            sources = positive.get("certificate_sources")
            _require(
                isinstance(sources, list) and sources,
                f"{positive_id}: certificate_sources must be nonempty",
            )
            for index, source in enumerate(sources):
                _validate_source(source, repository_root, f"{positive_id}.certificate_sources[{index}]")

        control_sources = family.get("censored_control_sources", [])
        _require(
            isinstance(control_sources, list),
            f"{family_id}: censored_control_sources must be a list",
        )
        for index, source_record in enumerate(control_sources):
            _require(
                isinstance(source_record, dict)
                and source_record.get("outcome_type")
                in {
                    "unlabelled_population_background",
                    "bounded_search_miss",
                    "mixed_bounded_search_outcomes",
                    "bounded_search_with_certified_lower_bounds",
                },
                f"{family_id}: controls must be explicitly censored",
            )
            _validate_source(
                source_record.get("source"),
                repository_root,
                f"{family_id}.censored_control_sources[{index}]",
            )

        runs = family.get("ranking_runs", [])
        _require(isinstance(runs, list), f"{family_id}: ranking_runs must be a list")
        run_ids: set[str] = set()
        for run in runs:
            _require(isinstance(run, dict), f"{family_id}: ranking run must be an object")
            run_id = run.get("id")
            _require(
                isinstance(run_id, str) and run_id and run_id not in run_ids,
                f"{family_id}: ranking run ids must be nonempty and unique",
            )
            run_ids.add(run_id)
            _require(
                run.get("evaluation_role") in ALLOWED_EVALUATION_ROLES,
                f"{family_id}/{run_id}: unknown evaluation role",
            )
            _require(
                run.get("latest_feature_stage") in ALLOWED_FEATURE_STAGES,
                f"{family_id}/{run_id}: feature stage is not pre-point-search",
            )
            _require(
                run.get("uses_point_search_features") is False,
                f"{family_id}/{run_id}: point-search features are forbidden",
            )
            _require(
                run.get("uses_rank_labels_for_scoring") is False,
                f"{family_id}/{run_id}: rank-label leakage is forbidden",
            )
            extractor = run.get("extractor")
            _require(isinstance(extractor, dict), f"{family_id}/{run_id}: missing extractor")
            _require(
                extractor.get("kind")
                in {
                    "elkies-positive-control-worst-block-v1",
                    "r17-bisection-gain-quarantined-replay-v1",
                    "fermigier-rank-jump-replay-v1",
                    "nagao-section7-rank-jump-replay-v1",
                },
                f"{family_id}/{run_id}: unsupported extractor",
            )
            if extractor.get("kind") == "fermigier-rank-jump-replay-v1":
                _require(
                    extractor.get("metric")
                    in {
                        "discovery_rank",
                        "discovery_composite",
                        "held_rank",
                        "held_composite",
                    },
                    f"{family_id}/{run_id}: unsupported Fermigier replay metric",
                )
            if extractor.get("kind") == "nagao-section7-rank-jump-replay-v1":
                _require(
                    extractor.get("metric") in {"training", "validation"},
                    f"{family_id}/{run_id}: unsupported Nagao replay metric",
                )
            if extractor.get("kind") == "r17-bisection-gain-quarantined-replay-v1":
                _require(
                    extractor.get("metric")
                    in {"learned_contrast", "weakest_block_nagao"},
                    f"{family_id}/{run_id}: unsupported R17 training replay metric",
                )
                if extractor.get("metric") == "learned_contrast":
                    group_gate_source = extractor.get("arithmetic_group_gate")
                    _validate_source(
                        group_gate_source,
                        repository_root,
                        f"{family_id}/{run_id}.arithmetic_group_gate",
                    )
                    group_gate = json.loads(
                        (repository_root / group_gate_source["path"]).read_text(
                            encoding="utf-8"
                        )
                    )
                    _require(
                        group_gate.get("schema")
                        == "elliptic-curves.r17-training-arithmetic-group-audit.v1"
                        and group_gate.get("status")
                        == "PASS_EXACT_GROUPING_BEFORE_FURTHER_LEARNED_SCORE_REUSE",
                        f"{family_id}/{run_id}: exact arithmetic grouping is not closed",
                    )
                    reuse_gate = group_gate.get("definition", {}).get("gate", {})
                    _require(
                        reuse_gate.get("status")
                        == "PASS_EXACT_ISOMORPHISM_TWIST_GROUPING"
                        and reuse_gate.get("learned_score_reuse_authorized") is True,
                        f"{family_id}/{run_id}: learned-score reuse is not authorized",
                    )
                    _require(
                        reuse_gate.get("authorized_score_artifact")
                        == extractor.get("source", {}).get("path")
                        and reuse_gate.get("authorized_score_artifact_sha256")
                        == extractor.get("source", {}).get("sha256"),
                        f"{family_id}/{run_id}: grouping gate authorizes another score",
                    )
            _validate_source(
                extractor.get("source"), repository_root, f"{family_id}/{run_id}.extractor"
            )
            declared_ids = run.get("positive_ids")
            _require(
                isinstance(declared_ids, list)
                and len(declared_ids) == len(set(declared_ids))
                and set(declared_ids) == local_positive_ids,
                f"{family_id}/{run_id}: run must cover every admitted family positive",
            )
        if family["readiness"] == "ranked_development_corpus":
            _require(
                bool(positives) and bool(runs),
                f"{family_id}: ranked corpus must have positives and ranking runs",
            )
        else:
            _require(
                not runs,
                f"{family_id}: only a ranked corpus may contain ranking runs",
            )


def _load_lab_run(
    family: dict[str, Any], run: dict[str, Any], repository_root: Path
) -> tuple[int, dict[str, int], dict[str, Any]]:
    extractor = run["extractor"]
    path = repository_root / extractor["source"]["path"]
    artifact = json.loads(path.read_text(encoding="utf-8"))
    kind = extractor["kind"]
    expected_schema = {
        "elkies-positive-control-worst-block-v1": (
            "elkies-2026-positive-control-worst-block-nagao-v1"
        ),
        "r17-bisection-gain-quarantined-replay-v1": (
            "elliptic-curves.r17-bisection-gain-ranker-quarantined-replay.v1"
        ),
        "fermigier-rank-jump-replay-v1": "elliptic-curves.fermigier-rank-jump-replay.v1",
        "nagao-section7-rank-jump-replay-v1": (
            "elliptic-curves.nagao-section7-rank-jump-replay.v1"
        ),
    }[kind]
    _require(
        artifact.get("schema") == expected_schema,
        f"{family['id']}/{run['id']}: unexpected ranking artifact schema",
    )
    if kind == "elkies-positive-control-worst-block-v1":
        population_count = artifact.get("population_count")
    elif kind == "r17-bisection-gain-quarantined-replay-v1":
        population_count = artifact.get("evaluation", {}).get("population_count")
    else:
        population_count = artifact.get("population", {}).get("primitive_parameter_count")
    _require(
        isinstance(population_count, int)
        and not isinstance(population_count, bool)
        and population_count > 0,
        f"{family['id']}/{run['id']}: invalid population count",
    )
    parameter_to_id = {
        positive["parameter"]: positive["id"] for positive in family["known_positives"]
    }
    ranks: dict[str, int] = {}
    if kind == "elkies-positive-control-worst-block-v1":
        records = artifact.get("positive_controls", [])
    elif kind == "r17-bisection-gain-quarantined-replay-v1":
        records = artifact.get("evaluation", {}).get("quarantined_controls", [])
    else:
        records = artifact.get("anchors", [])
    metric = extractor.get("metric")
    for record in records:
        if kind == "elkies-positive-control-worst-block-v1":
            parameter = record.get("score", {}).get("parameter")
        elif kind == "r17-bisection-gain-quarantined-replay-v1":
            parameter = record.get("parameter")
        elif kind == "nagao-section7-rank-jump-replay-v1":
            parameter = record.get("constructor_parameter_T")
        else:
            parameter = record.get("canonical_parameter_u")
        if parameter not in parameter_to_id:
            continue
        if kind == "elkies-positive-control-worst-block-v1":
            rank = record.get("population_rank")
        elif kind == "r17-bisection-gain-quarantined-replay-v1":
            rank = (
                record.get("methods", {})
                .get(metric, {})
                .get("empirical_rank_among_100000_development_rows")
            )
        else:
            rank = (
                record.get("scores", {})
                .get(metric, {})
                .get("rank_position_one_based")
            )
        _require(
            isinstance(rank, int)
            and not isinstance(rank, bool)
            and 1 <= rank <= population_count,
            f"{family['id']}/{run['id']}: invalid positive rank for {parameter}",
        )
        positive_id = parameter_to_id[parameter]
        _require(positive_id not in ranks, f"{family['id']}/{run['id']}: duplicate positive")
        ranks[positive_id] = rank
    _require(
        set(ranks) == set(run["positive_ids"]),
        f"{family['id']}/{run['id']}: ranking artifact does not cover declared positives",
    )
    if kind == "elkies-positive-control-worst-block-v1":
        metadata = {
            "ranking_key": artifact.get("scoring", {}).get("primary_ranking_key"),
            "prime_ensembles": artifact.get("scoring", {}).get("prime_ensembles"),
            "population_definition": artifact.get("search"),
            "proof_boundary": artifact.get("proof_boundary"),
        }
    elif kind == "r17-bisection-gain-quarantined-replay-v1":
        metadata = {
            "ranking_key": metric,
            "frozen_model_sha256": artifact.get("frozen_model_sha256_before_quarantine_open"),
            "population_definition": {
                "sampled_parameter_count": population_count,
                "feature_names": artifact.get("feature_names"),
            },
            "proof_boundary": artifact.get("proof_boundary"),
        }
    else:
        metadata = {
            "ranking_key": metric,
            "prime_ensembles": artifact.get("frozen_feature_specification"),
            "population_definition": artifact.get("population"),
            "proof_boundary": artifact.get("interpretation"),
        }
    return population_count, ranks, metadata


def _fraction_record(numerator: int, denominator: int) -> dict[str, Any]:
    fraction = Fraction(numerator, denominator)
    return {
        "exact": f"{fraction.numerator}/{fraction.denominator}",
        "decimal": float(fraction),
    }


def _retrieval_metrics(
    positive_ids: Sequence[str],
    ranks: dict[str, int],
    population_count: int,
    budget_fractions: Sequence[Fraction],
) -> dict[str, Any]:
    ordered = sorted(((ranks[positive_id], positive_id) for positive_id in positive_ids))
    count = len(ordered)
    _require(count > 0, "retrieval metrics require at least one positive")
    average_precision = fmean(index / rank for index, (rank, _) in enumerate(ordered, 1))
    budgets = []
    for fraction in budget_fractions:
        budget = (population_count * fraction.numerator + fraction.denominator - 1) // fraction.denominator
        hits = [positive_id for rank, positive_id in ordered if rank <= budget]
        expected = Fraction(budget * count, population_count)
        budgets.append(
            {
                "population_fraction": f"{fraction.numerator}/{fraction.denominator}",
                "candidate_budget": budget,
                "known_positive_hits": len(hits),
                "known_positive_ids": hits,
                "known_positive_recall": len(hits) / count,
                "random_ranking_expected_known_hits": float(expected),
                "known_positive_enrichment_over_random": float(
                    Fraction(len(hits), 1) / expected
                ),
            }
        )
    return {
        "known_positive_count": count,
        "best_known_positive_rank": ordered[0][0],
        "worst_known_positive_rank": ordered[-1][0],
        "worst_known_positive_population_fraction": _fraction_record(
            ordered[-1][0], population_count
        ),
        "average_precision_against_known_controls_only": average_precision,
        "budget_metrics": budgets,
    }


def evaluate_lab_manifest(
    manifest: dict[str, Any], repository_root: Path
) -> dict[str, Any]:
    """Evaluate every available within-family ranking run in the lab registry."""

    validate_lab_manifest(manifest, repository_root)
    budget_fractions = [Fraction(value) for value in manifest["budget_fractions"]]
    family_results = []
    for family in manifest["families"]:
        positives = {positive["id"]: positive for positive in family["known_positives"]}
        run_results = []
        for run in family["ranking_runs"]:
            population_count, ranks, metadata = _load_lab_run(family, run, repository_root)
            observations = []
            for positive_id in run["positive_ids"]:
                rank = ranks[positive_id]
                observations.append(
                    {
                        "positive_id": positive_id,
                        "parameter": positives[positive_id]["parameter"],
                        "certified_rank_lower_bound": positives[positive_id][
                            "certified_rank_lower_bound"
                        ],
                        "exceptional_quotient_rank_lower_bound": positives[positive_id][
                            "exceptional_quotient_rank_lower_bound"
                        ],
                        "population_rank": rank,
                        "population_fraction": _fraction_record(rank, population_count),
                        "search_space_reduction_factor": population_count / rank,
                    }
                )
            distinct_thresholds = sorted(
                {positive["exceptional_quotient_rank_lower_bound"] for positive in positives.values()}
            )
            threshold_metrics = {}
            for threshold in distinct_thresholds:
                threshold_ids = [
                    positive_id
                    for positive_id, positive in positives.items()
                    if positive["exceptional_quotient_rank_lower_bound"] >= threshold
                ]
                threshold_metrics[str(threshold)] = _retrieval_metrics(
                    threshold_ids, ranks, population_count, budget_fractions
                )
            run_results.append(
                {
                    "id": run["id"],
                    "evaluation_role": run["evaluation_role"],
                    "feature_boundary": {
                        "latest_feature_stage": run["latest_feature_stage"],
                        "uses_point_search_features": False,
                        "uses_rank_labels_for_scoring": False,
                    },
                    "population_count": population_count,
                    "source": run["extractor"]["source"],
                    "metadata": metadata,
                    "positive_observations": observations,
                    "all_known_positives": _retrieval_metrics(
                        list(positives), ranks, population_count, budget_fractions
                    ),
                    "minimum_quotient_gain_metrics": threshold_metrics,
                }
            )
        family_results.append(
            {
                "id": family["id"],
                "readiness": family["readiness"],
                "proved_arithmetic_generic_rank": family.get(
                    "proved_arithmetic_generic_rank"
                ),
                "known_positive_count": len(positives),
                "censored_control_source_count": len(
                    family.get("censored_control_sources", [])
                ),
                "ranking_runs": run_results,
                "next_action": family.get("next_action"),
            }
        )
    return {
        "schema": LAB_RESULT_SCHEMA,
        "objective": manifest["objective"],
        "label_policy": manifest["label_policy"],
        "families": family_results,
        "summary": {
            "family_count": len(family_results),
            "ranked_family_count": sum(bool(row["ranking_runs"]) for row in family_results),
            "certified_positive_count": sum(row["known_positive_count"] for row in family_results),
            "ranked_positive_count": sum(
                len(run["positive_observations"])
                for row in family_results
                for run in row["ranking_runs"]
            ),
        },
        "claim_boundary": (
            "This evaluates retrieval of already-known certified positive controls. "
            "A bounded search miss is censored, not a rank upper bound; development "
            "runs do not estimate out-of-family generalization."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--repository-root",
        type=Path,
        default=Path(__file__).resolve().parents[2],
        help="repository root for provenance checks",
    )
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
    if manifest.get("schema") == LAB_SCHEMA:
        result = evaluate_lab_manifest(manifest, args.repository_root.resolve())
    else:
        result = evaluate_manifest(manifest, top_ks)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
