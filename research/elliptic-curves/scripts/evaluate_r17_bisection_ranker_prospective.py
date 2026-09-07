#!/usr/bin/env python3
"""Evaluate the frozen R17 ranker on the label-after-freeze random holdout."""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import log
from pathlib import Path
import sys
from typing import Any, Callable


ROOT = Path(__file__).resolve().parents[2]
ECSEARCH = ROOT / "elliptic-curves/ecsearch"
sys.path.insert(0, str(ECSEARCH))

from r17_bisection_ranker import score, semantic_label_sha256  # noqa: E402


DEFAULT_FEATURES = ROOT / "artifacts/local/elliptic-curves/r17-prospective-holdout.jsonl"
DEFAULT_COMMITMENT = (
    ROOT / "artifacts/local/elliptic-curves/r17-prospective-holdout-summary.json"
)
DEFAULT_LABELS = (
    ROOT / "artifacts/local/elliptic-curves/r17-prospective-holdout-bisection-labels.jsonl"
)
DEFAULT_MODEL = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "r17_bisection_gain_ranker_quarantined_replay_v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "r17_bisection_gain_ranker_prospective_holdout_v1.json"
)
BUDGETS = (Fraction(1, 100), Fraction(1, 20), Fraction(1, 10))


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_sha256(value: Any) -> str:
    return sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def average_precision(scores: list[float], labels: list[int]) -> float | None:
    positives = sum(labels)
    if positives == 0:
        return None
    order = sorted(range(len(scores)), key=lambda index: (-scores[index], index))
    hits = 0
    total = 0.0
    for rank, index in enumerate(order, 1):
        if labels[index]:
            hits += 1
            total += hits / rank
    return total / positives


def roc_auc(scores: list[float], labels: list[int]) -> float | None:
    positives = sum(labels)
    negatives = len(labels) - positives
    if positives == 0 or negatives == 0:
        return None
    ordered = sorted(zip(scores, labels), key=lambda item: item[0])
    rank_sum = 0.0
    index = 0
    while index < len(ordered):
        end = index + 1
        while end < len(ordered) and ordered[end][0] == ordered[index][0]:
            end += 1
        average_rank = ((index + 1) + end) / 2
        rank_sum += average_rank * sum(label for _score, label in ordered[index:end])
        index = end
    return (rank_sum - positives * (positives + 1) / 2) / (positives * negatives)


def method_metrics(
    values: list[float], labels: list[int], gains: list[int], parameters: list[str]
) -> dict[str, Any]:
    order = sorted(range(len(values)), key=lambda index: (-values[index], parameters[index]))
    positives = sum(labels)
    budgets = []
    for fraction in BUDGETS:
        count = (len(values) * fraction.numerator + fraction.denominator - 1) // fraction.denominator
        chosen = order[:count]
        hits = sum(labels[index] for index in chosen)
        expected = Fraction(count * positives, len(values))
        budgets.append(
            {
                "population_fraction": f"{fraction.numerator}/{fraction.denominator}",
                "candidate_count": count,
                "positive_count": hits,
                "certified_gain_sum": sum(gains[index] for index in chosen),
                "positive_recall": None if positives == 0 else hits / positives,
                "enrichment_over_random": None if expected == 0 else float(Fraction(hits, 1) / expected),
            }
        )
    return {
        "average_precision": average_precision(values, labels),
        "roc_auc": roc_auc(values, labels),
        "budgets": budgets,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", type=Path, default=DEFAULT_FEATURES)
    parser.add_argument("--commitment", type=Path, default=DEFAULT_COMMITMENT)
    parser.add_argument("--labels", type=Path, default=DEFAULT_LABELS)
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    features = read_jsonl(args.features)
    label_rows = read_jsonl(args.labels)
    labels_by_parameter = {row["parameter"]: row for row in label_rows}
    if set(labels_by_parameter) != {row["parameter"] for row in features}:
        raise AssertionError("prospective features and labels have different parameters")
    commitment = json.loads(args.commitment.read_text())
    if commitment["status"] != "FROZEN_BEFORE_BISECTION_LABEL_EVALUATION":
        raise AssertionError("the holdout lacks a pre-label commitment")
    if commitment["selection"]["output_sha256"] != file_sha256(args.features):
        raise AssertionError("the prospective feature cohort changed after commitment")
    model_artifact = json.loads(args.model.read_text())
    frozen_model = model_artifact["frozen_model"]
    if canonical_sha256(frozen_model) != model_artifact[
        "frozen_model_sha256_before_quarantine_open"
    ]:
        raise AssertionError("the frozen ranker model hash changed")
    if commitment["commitment"]["frozen_model_sha256"] != model_artifact[
        "frozen_model_sha256_before_quarantine_open"
    ]:
        raise AssertionError("the prospective commitment names another model")

    usable_features = []
    binary_labels = []
    gains = []
    censored = []
    for row in features:
        gain = labels_by_parameter[row["parameter"]]["outcomes"][
            "finite_quotient_gain_lower_bound"
        ]
        if gain is None:
            censored.append(row["parameter"])
            continue
        usable_features.append(row)
        gains.append(int(gain))
        binary_labels.append(int(gain > 0))
    parameters = [row["parameter"] for row in usable_features]
    methods: dict[str, Callable[[dict[str, Any]], float]] = {
        "frozen_learned_contrast": lambda row: score(frozen_model["model"], row),
        "weakest_block_nagao": lambda row: float(
            row["features"]["level1_nagao"]["worst_block_signal"]
        ),
        "partial_conductor_quality_proxy": lambda row: float(
            row["features"]["level2_conductor_proxy"]["quality_proxy"]
        ),
        "negative_log_projective_height": lambda row: -log(float(row["height"])),
    }
    metric_rows = {
        name: method_metrics(
            [method(row) for row in usable_features], binary_labels, gains, parameters
        )
        for name, method in methods.items()
    }
    payload = {
        "schema": "elliptic-curves.r17-bisection-gain-ranker-prospective-holdout.v1",
        "status": "EXPERIMENTAL_PROSPECTIVE_BISECTION_TARGET_EVALUATION",
        "population": {
            "committed_row_count": len(features),
            "usable_row_count": len(usable_features),
            "censored_row_count": len(censored),
            "censored_parameters": censored,
            "positive_count": sum(binary_labels),
            "positive_rate": sum(binary_labels) / len(binary_labels),
            "certified_gain_sum": sum(gains),
        },
        "methods": metric_rows,
        "commitment": commitment["commitment"],
        "frozen_model_sha256": model_artifact["frozen_model_sha256_before_quarantine_open"],
        "semantic_label_sha256_excluding_timings": semantic_label_sha256(label_rows),
        "inputs": {
            str(args.features.resolve()): file_sha256(args.features),
            str(args.commitment.resolve()): file_sha256(args.commitment),
            str(args.model.resolve()): file_sha256(args.model),
        },
        "proof_boundary": [
            "The cohort and four evaluation methods were committed before its bisection labels were computed.",
            "This evaluates known-bisection visibility, not total rank or residual Selmer dimension.",
            "Censored rows are omitted rather than treated as negatives.",
            "Timing-bearing raw-label hashes are excluded from the deterministic evaluation artifact.",
        ],
        "generation": {
            "command": "python3 elliptic-curves/scripts/evaluate_r17_bisection_ranker_prospective.py",
            "script_sha256": file_sha256(Path(__file__)),
        },
    }
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.is_file() or args.output.read_text() != rendered:
            raise SystemExit("stale prospective R17 bisection evaluation")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(
        "R17PROSPECTIVE"
        f"|rows={len(features)}|usable={len(usable_features)}"
        f"|positives={sum(binary_labels)}|gain={sum(gains)}"
    )
    for name, metrics in metric_rows.items():
        top = metrics["budgets"][0]
        print(
            f"METHOD|name={name}|auc={metrics['roc_auc']:.6f}"
            f"|ap={metrics['average_precision']:.6f}"
            f"|top1_hits={top['positive_count']}|top1_enrichment={top['enrichment_over_random']:.6f}"
        )


if __name__ == "__main__":
    main()
