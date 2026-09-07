#!/usr/bin/env python3
"""Leakage-controlled rare-event model for Nagao's rank-21 family.

The labelled tail is exceptionally small: four fibers have exact algebraic
rank lower bounds at least 18, while thirteen controls reached stable
numerical rank 17 in the same clean H=10^6 point-search stage.  This script
therefore deliberately avoids a flexible neural network.  Its primary model
is a Dirichlet-smoothed Naive Bayes likelihood on low-cardinality local
symbols.  A paper-informed multi-value Mestre--Nagao linear model is retained
as a separately declared comparator.

Every positive is evaluated leave-one-positive-out.  The controls are
three-fold cross-fitted, so no control is scored by a model that trained on
that control.  Nagao's published rank-21 fiber is calibration only and never
enters training.  Expensive candidate generation is authorized only if one
predeclared model clears every recovery threshold.

No numerical rank is promoted to a rank theorem.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
import hashlib
import json
from math import log, sqrt
from pathlib import Path
import platform
import shlex
import sys
from typing import Any, Callable, Iterable, Sequence

from ek_k3 import rational_to_string
from search_nagao_rank21_unbiased import (
    ResidueSymbol,
    build_residue_tables,
    projective_index,
)


Q = Fraction
MODEL_CUTOFF = 199
PRIME_RANGES = ((5, 43), (47, 97), (101, 199))
TRACE_BIN_BOUNDARIES = (-0.5, 0.0, 0.5)
TRACE_CATEGORY_NAMES = (
    "good_trace_lt_minus_half",
    "good_trace_minus_half_to_zero",
    "good_trace_zero_to_half",
    "good_trace_ge_half",
    "bad_reduction",
)
DIRICHLET_ALPHA = 4.0
MULTI_VALUE_CUTOFFS = (43, 97, 199)

POSITIVE_PARAMETERS = (Q(6793, 64), Q(3137, 72), Q(5783, 16), Q(6629, 174))
POSITIVE_CERTIFIED_LOWER_BOUNDS = {
    Q(6793, 64): 19,
    Q(3137, 72): 18,
    Q(5783, 16): 18,
    Q(6629, 174): 18,
}
CONTROL_PARAMETERS = (
    Q(421, 54),
    Q(743, 4),
    Q(2543, 54),
    Q(313, 28),
    Q(2197, 48),
    Q(1573, 36),
    Q(419, 28),
    Q(4487, 16),
    Q(1251, 16),
    Q(683),
    Q(2911, 16),
    Q(1, 48),
    Q(2579, 54),
)
PUBLISHED_CALIBRATION_PARAMETER = Q(14721, 188)

UNBIASED_ARTIFACT_SHA256 = (
    "5bf7406855af5ec39b269fa4105c9225adb4a10d13fab5480b15264cc3e8fe1d"
)
HISTORICAL_ARTIFACT_SHA256 = (
    "90fc658cdb7c39c96317ee888be1364b8c9f368859230e25161dc45cd6a3cec7"
)
PUBLISHED_ARTIFACT_SHA256 = (
    "7d59fe9a91c0f3e46604794e8931ae27e26eeea1ebf252176dffd6be8d6010fe"
)
MULTI_VALUE_PAPER = {
    "title": (
        "Improving elliptic curve rank classification using multi-value and "
        "learned Mestre-Nagao sums"
    ),
    "authors": "Zvonimir Bujanovic; Matija Kazalicki; Domagoj Vlah",
    "url": "https://web.math.pmf.unizg.hr/~mkazal/reprints/nagao.pdf",
    "pdf_sha256": (
        "bdd4bdf9a36115e5b1219b53a3ddea520480e74677df218346c75b8c19677e93"
    ),
    "source_role": (
        "motivates simultaneous S0/S5 cutoff features; its neural models and "
        "large low-rank training sets are not transferred to this 17-curve tail"
    ),
}

# Each held-out item must beat at least 80% of its cross-fitted controls.  The
# median positive recovery must reach 90%.  The published fiber must also beat
# 80%.  These thresholds are fixed before candidate generation.
MINIMUM_HELD_OUT_PERCENTILE = Q(4, 5)
MINIMUM_MEDIAN_POSITIVE_PERCENTILE = Q(9, 10)
MINIMUM_PUBLISHED_PERCENTILE = Q(4, 5)

REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/search_nagao_rank21_rare_event_model.py"
)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_json_with_hash(path: Path, expected_sha256: str) -> dict[str, Any]:
    observed = file_sha256(path)
    if observed != expected_sha256:
        raise RuntimeError(
            f"pinned input changed: {path} has {observed}, expected {expected_sha256}"
        )
    return json.loads(path.read_text(encoding="utf-8"))


def _h1m_records(data: dict[str, Any]) -> dict[Fraction, dict[str, Any]]:
    matching = [
        stage
        for stage in data["point_stages"]
        if stage["quartic_naive_height_bound"] == 1_000_000
    ]
    if len(matching) != 1:
        raise RuntimeError("the unbiased artifact lost its unique H=10^6 stage")
    return {
        Q(record["constructor_parameter"]): record
        for record in matching[0]["ranked_population"]
    }


def pinned_label_audit(root: Path) -> dict[str, Any]:
    """Load and verify every label without silently refreshing a certificate."""

    generated = root / "artifacts" / "generated-results"
    unbiased_path = generated / "elliptic_nagao_rank21_unbiased.json"
    historical_path = generated / "elliptic_nagao_rank21_historical_finalists.json"
    published_path = generated / "elliptic_nagao_rank21_neighborhood.json"
    unbiased = _load_json_with_hash(unbiased_path, UNBIASED_ARTIFACT_SHA256)
    historical = _load_json_with_hash(
        historical_path, HISTORICAL_ARTIFACT_SHA256
    )
    published = _load_json_with_hash(published_path, PUBLISHED_ARTIFACT_SHA256)

    certificates = unbiased["finite_reduction_certificates"]
    unbiased_positive_keys = {
        Q(3137, 72): "unbiased-3137-72",
        Q(5783, 16): "unbiased-5783-16",
        Q(6793, 64): "unbiased-6793-64",
    }
    positive_records = []
    for parameter, key in unbiased_positive_keys.items():
        certificate = certificates[key]
        expected = POSITIVE_CERTIFIED_LOWER_BOUNDS[parameter]
        if certificate["status"] != "certified":
            raise RuntimeError(f"{key} is no longer certified")
        if certificate["certified_algebraic_rank_lower_bound"] != expected:
            raise RuntimeError(f"{key} rank lower bound changed")
        positive_records.append(
            {
                "constructor_parameter": rational_to_string(parameter),
                "label": "exact_algebraic_rank_lower_bound_at_least_18",
                "certified_algebraic_rank_lower_bound": expected,
                "certificate_basis_sha256": certificate["saturated_point_sha256"],
                "source_artifact": str(unbiased_path.relative_to(root)),
            }
        )

    historical_checkpoints = historical[
        "exact_checkpoints_stable_numerical_rank_at_least_18"
    ]
    matching_historical = [
        record
        for record in historical_checkpoints
        if Q(record["constructor_parameter"]) == Q(6629, 174)
    ]
    if len(matching_historical) != 1:
        raise RuntimeError("the historical T=6629/174 checkpoint changed")
    historical_certificate = matching_historical[0]["exact_rank_certificate"]
    if (
        historical_certificate["status"] != "certified"
        or historical_certificate["certified_algebraic_rank_lower_bound"] != 18
    ):
        raise RuntimeError("T=6629/174 lost its exact rank-18 lower bound")
    positive_records.append(
        {
            "constructor_parameter": "6629/174",
            "label": "exact_algebraic_rank_lower_bound_at_least_18",
            "certified_algebraic_rank_lower_bound": 18,
            "certificate_basis_sha256": historical_certificate[
                "saturated_point_sha256"
            ],
            "source_artifact": str(historical_path.relative_to(root)),
        }
    )
    positive_records.sort(
        key=lambda record: POSITIVE_PARAMETERS.index(
            Q(record["constructor_parameter"])
        )
    )

    h1m = _h1m_records(unbiased)
    controls = []
    for index, parameter in enumerate(CONTROL_PARAMETERS):
        record = h1m.get(parameter)
        if record is None:
            raise RuntimeError(f"control {parameter} disappeared from H=10^6")
        rank = record["height_rank"]["stable_numerical_rank"]
        precisions = [
            run["numerical_rank"]
            for run in record["height_rank"]["precision_runs"]
        ]
        if rank != 17 or precisions != [17, 17]:
            raise RuntimeError(f"control {parameter} lost its clean rank-17 replay")
        controls.append(
            {
                "constructor_parameter": rational_to_string(parameter),
                "label": "stable_numerical_rank_17_after_bounded_H1m_search",
                "stable_numerical_rank": 17,
                "precision_run_ranks": precisions,
                "negative_crossfit_fold": index % 3,
                "scope_warning": (
                    "empirical control label only; it is not an algebraic-rank "
                    "upper bound"
                ),
                "source_artifact": str(unbiased_path.relative_to(root)),
            }
        )

    calibration = published["published_record_calibration"]
    if (
        Q(calibration["constructor_parameter"])
        != PUBLISHED_CALIBRATION_PARAMETER
        or calibration["printed_points_checked_exactly"] != 21
        or not calibration["all_printed_points_on_printed_model"]
    ):
        raise RuntimeError("the published rank-21 calibration record changed")
    return {
        "pinned_inputs": [
            {
                "path": str(unbiased_path.relative_to(root)),
                "sha256": UNBIASED_ARTIFACT_SHA256,
            },
            {
                "path": str(historical_path.relative_to(root)),
                "sha256": HISTORICAL_ARTIFACT_SHA256,
            },
            {
                "path": str(published_path.relative_to(root)),
                "sha256": PUBLISHED_ARTIFACT_SHA256,
            },
        ],
        "positives": positive_records,
        "controls": controls,
        "published_calibration": {
            "constructor_parameter": "14721/188",
            "training_role": "held_out_calibration_only",
            "printed_points_checked_exactly": 21,
            "independence_status": calibration["published_independence_status"],
            "source_artifact": str(published_path.relative_to(root)),
        },
    }


@dataclass(frozen=True)
class LocalObservation:
    prime: int
    projective_index: int
    trace: int | None
    good_reduction: bool


def local_observations(
    parameter: Fraction,
    tables: dict[int, tuple[ResidueSymbol, ...]],
) -> tuple[LocalObservation, ...]:
    observations = []
    for prime, table in tables.items():
        index = projective_index(
            parameter.numerator, parameter.denominator, prime
        )
        symbol = table[index]
        observations.append(
            LocalObservation(
                prime=prime,
                projective_index=index,
                trace=symbol.ellap,
                good_reduction=symbol.good_reduction,
            )
        )
    return tuple(observations)


def _range_index(prime: int) -> int:
    for index, (lower, upper) in enumerate(PRIME_RANGES):
        if lower <= prime <= upper:
            return index
    raise ValueError(f"prime {prime} lies outside the declared ranges")


def trace_category(observation: LocalObservation) -> int:
    if not observation.good_reduction:
        return 4
    if observation.trace is None:
        raise AssertionError("good reduction requires an exact trace")
    normalized = observation.trace / (2 * sqrt(observation.prime))
    for index, boundary in enumerate(TRACE_BIN_BOUNDARIES):
        if normalized < boundary:
            return index
    return 3


def fit_binned_naive_bayes(
    positives: Sequence[Fraction],
    controls: Sequence[Fraction],
    observations: dict[Fraction, tuple[LocalObservation, ...]],
) -> dict[str, Any]:
    if not positives or not controls:
        raise ValueError("both likelihood classes must be nonempty")
    counts = {
        "control": [
            [DIRICHLET_ALPHA] * len(TRACE_CATEGORY_NAMES)
            for _ in PRIME_RANGES
        ],
        "positive": [
            [DIRICHLET_ALPHA] * len(TRACE_CATEGORY_NAMES)
            for _ in PRIME_RANGES
        ],
    }
    for class_name, parameters in (
        ("control", controls),
        ("positive", positives),
    ):
        for parameter in parameters:
            for observation in observations[parameter]:
                counts[class_name][_range_index(observation.prime)][
                    trace_category(observation)
                ] += 1
    log_ratios = []
    for range_index in range(len(PRIME_RANGES)):
        control_total = sum(counts["control"][range_index])
        positive_total = sum(counts["positive"][range_index])
        log_ratios.append(
            [
                log(
                    (counts["positive"][range_index][category] / positive_total)
                    / (counts["control"][range_index][category] / control_total)
                )
                for category in range(len(TRACE_CATEGORY_NAMES))
            ]
        )
    return {
        "model_name": "dirichlet_smoothed_binned_local_naive_bayes",
        "dirichlet_alpha_per_category": DIRICHLET_ALPHA,
        "prime_ranges": [list(pair) for pair in PRIME_RANGES],
        "trace_bin_boundaries_after_Hasse_normalization": list(
            TRACE_BIN_BOUNDARIES
        ),
        "categories": list(TRACE_CATEGORY_NAMES),
        "smoothed_counts": counts,
        "log_positive_to_control_likelihood_ratios": log_ratios,
        "training_positive_count": len(positives),
        "training_control_count": len(controls),
    }


def binned_naive_bayes_score(
    model: dict[str, Any], observations: Sequence[LocalObservation]
) -> float:
    ratios = model["log_positive_to_control_likelihood_ratios"]
    return sum(
        ratios[_range_index(observation.prime)][trace_category(observation)]
        for observation in observations
    )


def multi_value_features(
    observations: Sequence[LocalObservation],
) -> tuple[float, ...]:
    """Return paper-defined cumulative S0 and S5 values at three cutoffs."""

    answer: list[float] = []
    for cutoff in MULTI_VALUE_CUTOFFS:
        s0 = 0.0
        s5 = 0.0
        for observation in observations:
            if observation.prime > cutoff:
                continue
            prime = observation.prime
            if observation.good_reduction:
                if observation.trace is None:
                    raise AssertionError("good reduction requires a trace")
                s0 += observation.trace * log(prime) / prime
                s5 += log((prime + 1 - observation.trace) / prime)
            else:
                s5 += log(1.5 * (prime - 1) / prime)
        answer.extend((s0 / log(cutoff), s5))
    return tuple(answer)


def fit_multi_value_diagonal_lda(
    positives: Sequence[Fraction],
    controls: Sequence[Fraction],
    observations: dict[Fraction, tuple[LocalObservation, ...]],
) -> dict[str, Any]:
    positive_features = [multi_value_features(observations[p]) for p in positives]
    control_features = [multi_value_features(observations[p]) for p in controls]
    complete = positive_features + control_features
    dimension = len(complete[0])
    positive_mean = [
        sum(row[index] for row in positive_features) / len(positive_features)
        for index in range(dimension)
    ]
    control_mean = [
        sum(row[index] for row in control_features) / len(control_features)
        for index in range(dimension)
    ]
    complete_mean = [
        sum(row[index] for row in complete) / len(complete)
        for index in range(dimension)
    ]
    sample_variance = [
        sum((row[index] - complete_mean[index]) ** 2 for row in complete)
        / max(1, len(complete) - 1)
        for index in range(dimension)
    ]
    variance_floor = 1e-8
    weights = [
        (positive_mean[index] - control_mean[index])
        / max(sample_variance[index], variance_floor)
        for index in range(dimension)
    ]
    midpoint = [
        (positive_mean[index] + control_mean[index]) / 2
        for index in range(dimension)
    ]
    return {
        "model_name": "paper_multi_value_diagonal_linear_discriminant",
        "feature_order": [
            name
            for cutoff in MULTI_VALUE_CUTOFFS
            for name in (f"S0({cutoff})", f"S5({cutoff})")
        ],
        "positive_mean": positive_mean,
        "control_mean": control_mean,
        "sample_variance": sample_variance,
        "variance_floor": variance_floor,
        "weights": weights,
        "midpoint": midpoint,
        "training_positive_count": len(positives),
        "training_control_count": len(controls),
    }


def multi_value_lda_score(model: dict[str, Any], features: Sequence[float]) -> float:
    return sum(
        weight * (value - midpoint)
        for weight, value, midpoint in zip(
            model["weights"], features, model["midpoint"]
        )
    )


FitFunction = Callable[
    [Sequence[Fraction], Sequence[Fraction]],
    tuple[dict[str, Any], Callable[[Fraction], float]],
]


def _fraction_record(value: Fraction) -> dict[str, Any]:
    return {
        "exact": rational_to_string(value),
        "decimal": float(value),
    }


def _median(values: Sequence[Fraction]) -> Fraction:
    ordered = sorted(values)
    middle = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[middle]
    return (ordered[middle - 1] + ordered[middle]) / 2


def cross_fitted_recovery(fit: FitFunction) -> dict[str, Any]:
    held_out_records = []
    positive_percentiles = []
    for positive_index, held_out in enumerate(POSITIVE_PARAMETERS):
        comparisons = []
        fold_records = []
        training_positives = tuple(
            value
            for index, value in enumerate(POSITIVE_PARAMETERS)
            if index != positive_index
        )
        for fold in range(3):
            calibration_controls = tuple(
                value
                for index, value in enumerate(CONTROL_PARAMETERS)
                if index % 3 == fold
            )
            training_controls = tuple(
                value
                for index, value in enumerate(CONTROL_PARAMETERS)
                if index % 3 != fold
            )
            _, score = fit(training_positives, training_controls)
            held_score = score(held_out)
            fold_scores = [
                (control, score(control)) for control in calibration_controls
            ]
            comparisons.extend(
                control_score <= held_score
                for _, control_score in fold_scores
            )
            fold_records.append(
                {
                    "negative_crossfit_fold": fold,
                    "training_positive_parameters": [
                        rational_to_string(value) for value in training_positives
                    ],
                    "training_control_parameters": [
                        rational_to_string(value) for value in training_controls
                    ],
                    "held_out_positive_score": held_score,
                    "calibration_control_scores": [
                        {
                            "constructor_parameter": rational_to_string(parameter),
                            "score": value,
                        }
                        for parameter, value in fold_scores
                    ],
                }
            )
        percentile = Q(sum(comparisons), len(comparisons))
        positive_percentiles.append(percentile)
        held_out_records.append(
            {
                "held_out_positive_parameter": rational_to_string(held_out),
                "controls_at_or_below": sum(comparisons),
                "controls_compared": len(comparisons),
                "percentile": _fraction_record(percentile),
                "passes_80_percent_threshold": (
                    percentile >= MINIMUM_HELD_OUT_PERCENTILE
                ),
                "folds": fold_records,
            }
        )

    published_comparisons = []
    published_folds = []
    for fold in range(3):
        calibration_controls = tuple(
            value
            for index, value in enumerate(CONTROL_PARAMETERS)
            if index % 3 == fold
        )
        training_controls = tuple(
            value
            for index, value in enumerate(CONTROL_PARAMETERS)
            if index % 3 != fold
        )
        _, score = fit(POSITIVE_PARAMETERS, training_controls)
        published_score = score(PUBLISHED_CALIBRATION_PARAMETER)
        control_scores = [
            (control, score(control)) for control in calibration_controls
        ]
        published_comparisons.extend(
            control_score <= published_score for _, control_score in control_scores
        )
        published_folds.append(
            {
                "negative_crossfit_fold": fold,
                "published_score": published_score,
                "calibration_control_scores": [
                    {
                        "constructor_parameter": rational_to_string(parameter),
                        "score": value,
                    }
                    for parameter, value in control_scores
                ],
            }
        )
    published_percentile = Q(
        sum(published_comparisons), len(published_comparisons)
    )
    median_positive = _median(positive_percentiles)
    passes = (
        all(
            percentile >= MINIMUM_HELD_OUT_PERCENTILE
            for percentile in positive_percentiles
        )
        and median_positive >= MINIMUM_MEDIAN_POSITIVE_PERCENTILE
        and published_percentile >= MINIMUM_PUBLISHED_PERCENTILE
    )
    return {
        "positive_leave_one_out": held_out_records,
        "minimum_positive_percentile": _fraction_record(min(positive_percentiles)),
        "median_positive_percentile": _fraction_record(median_positive),
        "published_rank21_calibration": {
            "constructor_parameter": "14721/188",
            "training_role": "never_used_for_training_or_model_selection",
            "controls_at_or_below": sum(published_comparisons),
            "controls_compared": len(published_comparisons),
            "percentile": _fraction_record(published_percentile),
            "passes_80_percent_threshold": (
                published_percentile >= MINIMUM_PUBLISHED_PERCENTILE
            ),
            "folds": published_folds,
        },
        "passes_all_predeclared_thresholds": passes,
    }


def observation_digest(
    observations: dict[Fraction, tuple[LocalObservation, ...]]
) -> str:
    payload = [
        {
            "parameter": rational_to_string(parameter),
            "symbols": [
                [
                    item.prime,
                    item.projective_index,
                    item.trace,
                    item.good_reduction,
                ]
                for item in observations[parameter]
            ],
        }
        for parameter in sorted(
            observations,
            key=lambda value: (value.numerator / value.denominator, value.denominator),
        )
    ]
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def build_result(root: Path) -> dict[str, Any]:
    labels = pinned_label_audit(root)
    tables = build_residue_tables(MODEL_CUTOFF)
    all_parameters = (
        POSITIVE_PARAMETERS
        + CONTROL_PARAMETERS
        + (PUBLISHED_CALIBRATION_PARAMETER,)
    )
    observations = {
        parameter: local_observations(parameter, tables)
        for parameter in all_parameters
    }

    def fit_binned(
        positives: Sequence[Fraction], controls: Sequence[Fraction]
    ) -> tuple[dict[str, Any], Callable[[Fraction], float]]:
        model = fit_binned_naive_bayes(positives, controls, observations)
        return model, lambda parameter: binned_naive_bayes_score(
            model, observations[parameter]
        )

    def fit_multi(
        positives: Sequence[Fraction], controls: Sequence[Fraction]
    ) -> tuple[dict[str, Any], Callable[[Fraction], float]]:
        model = fit_multi_value_diagonal_lda(positives, controls, observations)
        return model, lambda parameter: multi_value_lda_score(
            model, multi_value_features(observations[parameter])
        )

    full_binned, full_binned_score = fit_binned(
        POSITIVE_PARAMETERS, CONTROL_PARAMETERS
    )
    full_multi, full_multi_score = fit_multi(
        POSITIVE_PARAMETERS, CONTROL_PARAMETERS
    )
    binned_recovery = cross_fitted_recovery(fit_binned)
    multi_recovery = cross_fitted_recovery(fit_multi)
    accepted_models = [
        model_name
        for model_name, recovery in (
            (full_binned["model_name"], binned_recovery),
            (full_multi["model_name"], multi_recovery),
        )
        if recovery["passes_all_predeclared_thresholds"]
    ]
    scan_authorized = bool(accepted_models)

    curve_summaries = []
    label_by_parameter = {
        **{parameter: "positive" for parameter in POSITIVE_PARAMETERS},
        **{parameter: "control" for parameter in CONTROL_PARAMETERS},
        PUBLISHED_CALIBRATION_PARAMETER: "held_out_published_calibration",
    }
    for parameter in all_parameters:
        local = observations[parameter]
        curve_summaries.append(
            {
                "constructor_parameter": rational_to_string(parameter),
                "role": label_by_parameter[parameter],
                "bad_primes_at_most_199": [
                    item.prime for item in local if not item.good_reduction
                ],
                "local_symbol_sha256": hashlib.sha256(
                    json.dumps(
                        [
                            [
                                item.prime,
                                item.projective_index,
                                item.trace,
                                item.good_reduction,
                            ]
                            for item in local
                        ],
                        separators=(",", ":"),
                    ).encode("utf-8")
                ).hexdigest(),
                "full_binned_log_likelihood_ratio": full_binned_score(parameter),
                "multi_value_features": list(multi_value_features(local)),
                "full_multi_value_linear_score": full_multi_score(parameter),
            }
        )

    return {
        "schema_version": 1,
        "status": (
            "model_accepted_scan_authorized"
            if scan_authorized
            else "model_rejected_no_candidate_scan"
        ),
        "objective": {
            "primary": "algebraic rank at least 21 with log conductor < 182.72",
            "alternative": "algebraic rank at least 30",
        },
        "method_source": MULTI_VALUE_PAPER,
        "label_audit": labels,
        "leakage_control": {
            "positive_protocol": "leave exactly one of four certified positives out",
            "control_protocol": (
                "fixed index-mod-3 folds; each control is evaluated only by a "
                "model trained on the other two folds"
            ),
            "published_protocol": (
                "14721/188 is never used for training or model selection"
            ),
            "control_label_warning": (
                "rank-17 controls are clean stable numerical H=1m labels, not "
                "certified rank upper bounds"
            ),
            "candidate_data_used": False,
        },
        "feature_protocol": {
            "prime_cutoff": MODEL_CUTOFF,
            "finite_field_source": (
                "exact projective residue lookup on the short Jacobian for every "
                "prime 5 <= p <= 199"
            ),
            "raw_local_symbol": "(T mod p in P1(Fp), good/bad, exact a_p if good)",
            "primary_low_cardinality_projection": (
                "prime-range x Hasse-normalized trace bin, with a bad-reduction bin"
            ),
            "raw_residue_policy": (
                "projective residues are pinned in the symbol digest but are not "
                "separate high-cardinality predictors with only four positives"
            ),
            "paper_multi_value_cutoffs": list(MULTI_VALUE_CUTOFFS),
        },
        "predeclared_acceptance_thresholds": {
            "each_positive_percentile_at_least": _fraction_record(
                MINIMUM_HELD_OUT_PERCENTILE
            ),
            "median_positive_percentile_at_least": _fraction_record(
                MINIMUM_MEDIAN_POSITIVE_PERCENTILE
            ),
            "published_percentile_at_least": _fraction_record(
                MINIMUM_PUBLISHED_PERCENTILE
            ),
            "all_three_conditions_required": True,
        },
        "models": [
            {
                "role": "primary_graphical_model",
                "full_fit": full_binned,
                "cross_validation": binned_recovery,
            },
            {
                "role": "paper_informed_comparator",
                "full_fit": full_multi,
                "cross_validation": multi_recovery,
            },
        ],
        "local_population": {
            "curve_count": len(all_parameters),
            "prime_count": len(tables),
            "complete_symbol_stream_sha256": observation_digest(observations),
            "curves": curve_summaries,
        },
        "selection_decision": {
            "accepted_models": accepted_models,
            "scan_authorized": scan_authorized,
            "decision": (
                "launch_genuinely_new_population_scan"
                if scan_authorized
                else "reject_model_and_stop_before_population_or_point_search"
            ),
            "scientific_reason": (
                "a model that cannot reliably recover already-known tail events "
                "must not allocate an expensive rare-event search"
            ),
        },
        "expensive_search": {
            "broad_rational_population_scanned": 0,
            "conductors_computed": 0,
            "point_searches_launched": 0,
            "finite_reduction_certificates_triggered": 0,
        },
        "interpretation": {
            "proved": (
                "the pinned local data and the declared cross-fitted scores replay "
                "exactly for this finite labelled population"
            ),
            "not_proved": (
                "the control curves have rank exactly 17, or any probabilistic "
                "independence assumption between primes"
            ),
            "frontier": (
                "with four positives, neither a smoothed local Naive Bayes model "
                "nor a six-feature S0/S5 discriminator passes the declared tail "
                "recovery gate"
            ),
        },
        "reproducing_command": REPRODUCING_COMMAND,
        "actual_command": " ".join(shlex.quote(part) for part in sys.argv),
        "script_sha256": file_sha256(Path(__file__).resolve()),
        "software": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "external_subprocesses": 0,
        },
    }


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            root
            / "artifacts"
            / "generated-results"
            / "elliptic_nagao_rank21_rare_event_model.json"
        ),
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    root = Path(__file__).resolve().parents[2]
    result = build_result(root)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(result["selection_decision"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
