"""Interpretable ranker for exact R17 known-bisection gain labels.

The model is deliberately small: an overlap-adjusted, standardized difference
of train-fold class means.  It has no hyperparameter search and consumes no
validation, internal-test, or quarantined-control labels during fitting.
"""

from __future__ import annotations

from hashlib import sha256
import json
from math import log, sqrt
from typing import Any, Iterable, Sequence


FEATURE_NAMES = (
    "log_projective_height",
    "nagao_block_0",
    "nagao_block_1",
    "nagao_block_2",
    "nagao_worst_block",
    "nagao_mean_block",
    "nagao_bad_prime_count",
    "conductor_quality_proxy",
    "conductor_known_prime_log_saving",
    "conductor_log_discriminant_after_known_scaling",
    "quotient_code_rarity",
    "quotient_p19_mod2_dimension",
    "quotient_p19_mod3_dimension",
    "quotient_p23_mod2_dimension",
    "quotient_p23_mod3_dimension",
    "quotient_p29_mod2_dimension",
    "quotient_p29_mod3_dimension",
    "quotient_p31_mod2_dimension",
    "quotient_p31_mod3_dimension",
    "quotient_p37_mod2_dimension",
    "quotient_p37_mod3_dimension",
    "cover_distinct_character_patterns",
    "cover_character_pattern_entropy",
    "cover_branch_zero_count",
)


def semantic_label_sha256(rows: Iterable[dict[str, Any]]) -> str:
    """Hash arithmetic labels while excluding deliberately variable timings."""

    digest = sha256()
    for row in rows:
        semantic = dict(row)
        outcomes = dict(row["outcomes"])
        for name in ("cpu_seconds", "cpu_seconds_per_new_direction", "peak_rss_bytes"):
            outcomes.pop(name, None)
        semantic["outcomes"] = outcomes
        digest.update(
            (json.dumps(semantic, sort_keys=True, separators=(",", ":")) + "\n").encode()
        )
    return digest.hexdigest()


def feature_vector(row: dict[str, Any]) -> dict[str, float]:
    features = row["features"]
    nagao = features["level1_nagao"]
    conductor = features["level2_conductor_proxy"]
    quotient = features["level2_quotient_code"]
    cover = features["level2_cover_diversity"]
    blocks = nagao["standardized_block_signals"]
    code = quotient["local_E_mod_2_and_mod_3_dimensions"]
    if len(blocks) != 3 or len(code) != 5:
        raise ValueError("the frozen R17 feature shape changed")
    values = {
        "log_projective_height": log(float(row["height"])),
        "nagao_block_0": float(blocks[0]),
        "nagao_block_1": float(blocks[1]),
        "nagao_block_2": float(blocks[2]),
        "nagao_worst_block": float(nagao["worst_block_signal"]),
        "nagao_mean_block": float(nagao["mean_block_signal"]),
        "nagao_bad_prime_count": float(nagao["bad_prime_count"]),
        "conductor_quality_proxy": float(conductor["quality_proxy"]),
        "conductor_known_prime_log_saving": float(
            conductor["known_prime_log_discriminant_saving"]
        ),
        "conductor_log_discriminant_after_known_scaling": float(
            conductor["log_discriminant_after_known_scaling"]
        ),
        "quotient_code_rarity": float(quotient["rarity"]),
        "cover_distinct_character_patterns": float(cover["distinct_character_patterns"]),
        "cover_character_pattern_entropy": float(cover["character_pattern_entropy"]),
        "cover_branch_zero_count": float(cover["branch_zero_count"]),
    }
    for index, prime in enumerate((19, 23, 29, 31, 37)):
        values[f"quotient_p{prime}_mod2_dimension"] = float(code[index][0])
        values[f"quotient_p{prime}_mod3_dimension"] = float(code[index][1])
    if set(values) != set(FEATURE_NAMES):
        raise AssertionError("feature extraction disagrees with the frozen feature list")
    return values


def fit_weighted_contrast(
    rows: Sequence[dict[str, Any]], labels: Sequence[int], weights: Sequence[float]
) -> dict[str, Any]:
    if not (len(rows) == len(labels) == len(weights)) or not rows:
        raise ValueError("training rows, labels, and weights must have one nonempty length")
    if any(label not in (0, 1) for label in labels):
        raise ValueError("the contrast target must be binary")
    if any(weight <= 0 for weight in weights):
        raise ValueError("training weights must be positive")
    vectors = [feature_vector(row) for row in rows]

    def weighted_mean(name: str, selected: Iterable[int]) -> float:
        indices = list(selected)
        denominator = sum(weights[index] for index in indices)
        if denominator == 0:
            raise ValueError("a training class has zero total weight")
        return sum(weights[index] * vectors[index][name] for index in indices) / denominator

    all_indices = range(len(rows))
    positive_indices = [index for index, label in enumerate(labels) if label]
    negative_indices = [index for index, label in enumerate(labels) if not label]
    if not positive_indices or not negative_indices:
        raise ValueError("training requires positive and negative labels")
    means = {name: weighted_mean(name, all_indices) for name in FEATURE_NAMES}
    scales = {}
    total_weight = sum(weights)
    for name in FEATURE_NAMES:
        variance = sum(
            weights[index] * (vectors[index][name] - means[name]) ** 2
            for index in range(len(rows))
        ) / total_weight
        scales[name] = max(sqrt(variance), 1e-12)
    positive_means = {
        name: weighted_mean(name, positive_indices) for name in FEATURE_NAMES
    }
    negative_means = {
        name: weighted_mean(name, negative_indices) for name in FEATURE_NAMES
    }
    coefficients = {
        name: (positive_means[name] - negative_means[name]) / scales[name]
        for name in FEATURE_NAMES
    }
    return {
        "kind": "overlap-adjusted standardized train-class-mean contrast",
        "feature_names": list(FEATURE_NAMES),
        "means": means,
        "scales": scales,
        "positive_means": positive_means,
        "negative_means": negative_means,
        "coefficients": coefficients,
        "training_row_count": len(rows),
        "training_positive_count": len(positive_indices),
        "training_negative_count": len(negative_indices),
        "training_weight_sum": total_weight,
    }


def score(model: dict[str, Any], row: dict[str, Any]) -> float:
    vector = feature_vector(row)
    return sum(
        float(model["coefficients"][name])
        * (vector[name] - float(model["means"][name]))
        / float(model["scales"][name])
        for name in model["feature_names"]
    )


def ranking_metrics(scores: Sequence[float], labels: Sequence[int], budgets: Sequence[int]) -> dict[str, Any]:
    if len(scores) != len(labels) or not scores:
        raise ValueError("scores and labels must have one nonempty length")
    order = sorted(range(len(scores)), key=lambda index: (-scores[index], index))
    positives = sum(labels)
    metrics = {
        "row_count": len(labels),
        "positive_count": positives,
        "positive_rate": positives / len(labels),
        "budgets": {},
    }
    for requested in budgets:
        count = min(requested, len(labels))
        hits = sum(labels[index] for index in order[:count])
        expected = count * positives / len(labels)
        metrics["budgets"][str(requested)] = {
            "evaluated_count": count,
            "positive_count": hits,
            "recall": None if positives == 0 else hits / positives,
            "enrichment_over_random": None if expected == 0 else hits / expected,
        }
    return metrics
