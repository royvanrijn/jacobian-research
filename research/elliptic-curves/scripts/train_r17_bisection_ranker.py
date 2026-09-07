#!/usr/bin/env python3
"""Fit a train-only R17 bisection-gain ranker and open the quarantined replay.

The frozen model is an interpretable standardized class-mean contrast.  It is
fit only on exact, noncensored train-fold labels from the selected cohort.
Validation, internal test, and the four published controls are scored only
after the serialized model hash has been fixed.
"""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from math import log
from pathlib import Path
import sys
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC_ROOT = ROOT / "elliptic-curves"
ECSEARCH = ELLIPTIC_ROOT / "ecsearch"
K3_SCRIPTS = ROOT / "elkies-k3/scripts"
sys.path[:0] = [str(ECSEARCH), str(K3_SCRIPTS)]

from r17_bisection_ranker import (  # noqa: E402
    FEATURE_NAMES,
    fit_weighted_contrast,
    ranking_metrics,
    score,
    semantic_label_sha256,
)
from r17_training_data import (  # noqa: E402
    EMBARGOED_PARAMETERS,
    build_cheap_prime_tables,
    conductor_proxy_features,
    cover_character_tables,
    cover_diversity_features,
    nagao_features,
    parameter_text,
    quotient_code_features,
    select_cover_panel,
)
from search_h92_q12o5867_rootless_nagao import is_prime, load_family_model  # noqa: E402


MODEL = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
PANEL_SOURCE = ROOT / "artifacts/generated-results/elkies-2026-equation-bisections.json"
DEFAULT_POPULATION = ROOT / "artifacts/local/elliptic-curves/r17-training-population.jsonl"
DEFAULT_SELECTED = ROOT / "artifacts/local/elliptic-curves/r17-training-selected.jsonl"
DEFAULT_LABELS = ROOT / "artifacts/local/elliptic-curves/r17-training-bisection-labels.jsonl"
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "r17_bisection_gain_ranker_quarantined_replay_v1.json"
)


CONTROL_RANK_LOWER_BOUNDS = {
    (-2, 377): 25,
    (-308, 251): 26,
    (2456, 135): 27,
    (-9529, 5471): 28,
}


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def canonical_sha256(value: Any) -> str:
    return sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def read_jsonl(path: Path) -> Iterable[dict[str, Any]]:
    with path.open() as handle:
        for line in handle:
            if line.strip():
                yield json.loads(line)


def fold_metrics(
    feature_rows: list[dict[str, Any]],
    label_by_parameter: dict[str, dict[str, Any]],
    model: dict[str, Any],
    *,
    split: str,
    required_lane: str | None = None,
) -> dict[str, Any]:
    rows = []
    labels = []
    for row in feature_rows:
        if row["split"] != split:
            continue
        if required_lane is not None and required_lane not in row["selection_lanes"]:
            continue
        gain = label_by_parameter[row["parameter"]]["outcomes"][
            "finite_quotient_gain_lower_bound"
        ]
        if gain is None:
            continue
        rows.append(row)
        labels.append(int(gain > 0))
    budgets = [10, 50, 100]
    return {
        "learned_contrast": ranking_metrics([score(model, row) for row in rows], labels, budgets),
        "weakest_block_nagao": ranking_metrics(
            [float(row["features"]["level1_nagao"]["worst_block_signal"]) for row in rows],
            labels,
            budgets,
        ),
    }


def build_control_features(
    *,
    training_code_counts: Counter[str],
    smoothed_population: int,
) -> dict[tuple[int, int], dict[str, Any]]:
    family_model = load_family_model(MODEL)
    model_document = json.loads(MODEL.read_text())
    a_coefficients = [int(value) for value in model_document["A_coefficients_low_to_high"]]
    b_coefficients = [int(value) for value in model_document["B_coefficients_low_to_high"]]
    score_primes = [value for value in range(19, 600) if is_prime(value)]
    prime_blocks = [score_primes[index::3] for index in range(3)]
    quotient_primes = [19, 23, 29, 31, 37]
    conductor_primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43]
    prime_tables = build_cheap_prime_tables(
        family_model, sorted(set(score_primes + quotient_primes))
    )
    panel_document = json.loads(PANEL_SOURCE.read_text())
    cover_panel = select_cover_panel(panel_document["bisections"], 128)
    cover_primes = [101, 103, 107, 109, 113]
    cover_tables = cover_character_tables(cover_panel, cover_primes)
    answer = {}
    for parameter in sorted(EMBARGOED_PARAMETERS):
        quotient = quotient_code_features(parameter, prime_tables, quotient_primes)
        frequency = training_code_counts.get(str(quotient["code"]), 0)
        quotient["training_frequency"] = frequency
        quotient["rarity"] = log(smoothed_population / (frequency + 1))
        quotient["rarity_reference"] = "train split with add-one smoothing"
        height = max(abs(parameter[0]), parameter[1])
        answer[parameter] = {
            "parameter": parameter_text(parameter),
            "projective_pair": list(parameter),
            "height": height,
            "features": {
                "level0": {
                    "numerator": parameter[0],
                    "denominator": parameter[1],
                    "projective_height": height,
                },
                "level1_nagao": nagao_features(parameter, prime_tables, prime_blocks),
                "level2_conductor_proxy": conductor_proxy_features(
                    parameter, a_coefficients, b_coefficients, conductor_primes
                ),
                "level2_quotient_code": quotient,
                "level2_cover_diversity": cover_diversity_features(
                    parameter, cover_tables, cover_primes
                ),
            },
        }
    return answer


def empirical_rank(
    population_scores: list[tuple[str, float]], control_parameter: str, control_score: float
) -> int:
    return 1 + sum(
        candidate_score > control_score
        or (candidate_score == control_score and candidate_parameter < control_parameter)
        for candidate_parameter, candidate_score in population_scores
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--population", type=Path, default=DEFAULT_POPULATION)
    parser.add_argument("--selected", type=Path, default=DEFAULT_SELECTED)
    parser.add_argument("--labels", type=Path, default=DEFAULT_LABELS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    selected_rows = list(read_jsonl(args.selected))
    label_rows = list(read_jsonl(args.labels))
    label_semantic_sha256 = semantic_label_sha256(label_rows)
    label_by_parameter = {row["parameter"]: row for row in label_rows}
    if set(label_by_parameter) != {row["parameter"] for row in selected_rows}:
        raise AssertionError("selected features and exact label rows do not have identical keys")
    if any(tuple(row["projective_pair"]) in EMBARGOED_PARAMETERS for row in selected_rows):
        raise AssertionError("an embargoed control entered model development")

    training_rows = []
    training_labels = []
    training_weights = []
    for row in selected_rows:
        if row["split"] != "train":
            continue
        gain = label_by_parameter[row["parameter"]]["outcomes"][
            "finite_quotient_gain_lower_bound"
        ]
        if gain is None:
            continue
        training_rows.append(row)
        training_labels.append(int(gain > 0))
        training_weights.append(1.0 / len(row["selection_lanes"]))
    fitted = fit_weighted_contrast(training_rows, training_labels, training_weights)
    frozen_model = {
        "schema": "elliptic-curves.r17-bisection-gain-ranker.v1",
        "target": "finite_quotient_gain_lower_bound > 0 from the complete known-bisection atlas",
        "fit_split": "train",
        "censored_rows_excluded": True,
        "row_weight": "1 / selected-lane membership count",
        "validation_or_control_used_for_fit": False,
        "model": fitted,
        "development_inputs": {
            relative(args.selected): file_sha256(args.selected),
            relative(args.labels): {
                "semantic_sha256_excluding_timings": label_semantic_sha256,
                "reason": "per-row CPU timings are measured data and vary across exact replays",
            },
        },
        "training_parameter_sha256": sha256(
            "\n".join(row["parameter"] for row in training_rows).encode()
        ).hexdigest(),
    }
    frozen_model_sha256 = canonical_sha256(frozen_model)

    # Everything below this line evaluates the already-frozen score.
    training_code_counts: Counter[str] = Counter()
    training_population_count = 0
    population_count = 0
    for row in read_jsonl(args.population):
        population_count += 1
        if row["split"] == "train":
            training_population_count += 1
            training_code_counts[row["features"]["level2_quotient_code"]["code"]] += 1
    smoothed_population = training_population_count + len(training_code_counts)
    controls = build_control_features(
        training_code_counts=training_code_counts,
        smoothed_population=smoothed_population,
    )
    control_scores = {
        parameter: {
            "learned_contrast": score(fitted, row),
            "weakest_block_nagao": float(
                row["features"]["level1_nagao"]["worst_block_signal"]
            ),
        }
        for parameter, row in controls.items()
    }
    population_scores = {"learned_contrast": [], "weakest_block_nagao": []}
    for row in read_jsonl(args.population):
        population_scores["learned_contrast"].append((row["parameter"], score(fitted, row)))
        population_scores["weakest_block_nagao"].append(
            (row["parameter"], float(row["features"]["level1_nagao"]["worst_block_signal"]))
        )

    replay = []
    for parameter in sorted(controls):
        parameter_label = parameter_text(parameter)
        methods = {}
        for method, control_score in control_scores[parameter].items():
            rank = empirical_rank(population_scores[method], parameter_label, control_score)
            methods[method] = {
                "score": control_score,
                "empirical_rank_among_100000_development_rows": rank,
                "empirical_population_fraction_at_or_above": rank / population_count,
                "retrieved_at_top_0_01_percent": rank <= max(1, round(0.0001 * population_count)),
                "retrieved_at_top_0_1_percent": rank <= max(1, round(0.001 * population_count)),
                "retrieved_at_top_1_percent": rank <= max(1, round(0.01 * population_count)),
            }
        replay.append(
            {
                "parameter": parameter_label,
                "published_rank_lower_bound_labels_only": CONTROL_RANK_LOWER_BOUNDS[parameter],
                "methods": methods,
            }
        )

    fold_results = {}
    for split in ("train", "validation", "internal_test"):
        fold_results[split] = {
            "selected_union": fold_metrics(
                selected_rows, label_by_parameter, fitted, split=split
            ),
            "random_control_lane": fold_metrics(
                selected_rows,
                label_by_parameter,
                fitted,
                split=split,
                required_lane="random_controls",
            ),
        }
    primary = {row["parameter"]: row for row in replay}
    primary_success = all(
        primary[parameter]["methods"]["learned_contrast"]["retrieved_at_top_1_percent"]
        for parameter in ("2456/135", "-9529/5471")
    )
    payload = {
        "schema": "elliptic-curves.r17-bisection-gain-ranker-quarantined-replay.v1",
        "status": "EXPERIMENTAL_MECHANICALLY_QUARANTINED_REPLAY",
        "frozen_model": frozen_model,
        "frozen_model_sha256_before_quarantine_open": frozen_model_sha256,
        "evaluation": {
            "fold_metrics": fold_results,
            "population_count": population_count,
            "quarantined_controls": replay,
            "primary_rank27_rank28_top_one_percent_success": primary_success,
        },
        "quarantine_audit": {
            "all_four_controls_absent_from_feature_population_and_label_cohort": True,
            "controls_scored_only_after_model_serialization_hash": True,
            "control_rank_lower_bounds_are_labels_only_and_never_score_inputs": True,
            "human_blind": False,
            "description": (
                "Mechanically quarantined replay: the public parameters were known to the "
                "investigators, but the executable fit path had no control row or label."
            ),
        },
        "inputs": {
            relative(path): file_sha256(path)
            for path in (args.population, args.selected, MODEL, PANEL_SOURCE)
        },
        "semantic_label_input": {
            "path": relative(args.labels),
            "sha256_excluding_timings": label_semantic_sha256,
        },
        "feature_names": list(FEATURE_NAMES),
        "proof_boundary": [
            "The target is visibility through the known bisection atlas, not total Mordell--Weil rank.",
            "Empirical control ranks use the sampled 100,000-row population, not all height-at-most-10,000 parameters.",
            "The selected cohort is stratified; random-control-lane metrics are reported separately.",
            "No rank upper bound, Selmer dimension, or unrestricted point-search result is inferred.",
        ],
        "generation": {
            "command": "python3 elliptic-curves/scripts/train_r17_bisection_ranker.py",
            "script_sha256": file_sha256(Path(__file__)),
        },
    }
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.is_file() or args.output.read_text() != rendered:
            raise SystemExit("stale R17 bisection-gain ranker replay artifact")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(
        "R17BISECTIONRANKER"
        f"|model={frozen_model_sha256}"
        f"|primary_top1={str(primary_success).lower()}"
        f"|artifact={args.output}"
    )
    for row in replay:
        learned = row["methods"]["learned_contrast"]
        nagao = row["methods"]["weakest_block_nagao"]
        print(
            f"CONTROL|t={row['parameter']}|rank_lb={row['published_rank_lower_bound_labels_only']}"
            f"|learned_rank={learned['empirical_rank_among_100000_development_rows']}"
            f"|nagao_rank={nagao['empirical_rank_among_100000_development_rows']}"
        )


if __name__ == "__main__":
    main()
