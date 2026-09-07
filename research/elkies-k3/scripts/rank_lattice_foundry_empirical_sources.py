#!/usr/bin/env python3
"""Rank MW1 foundry sources by transparent equation-attempt analogues.

This is deliberately not a calibrated equation-success classifier.  The
current foundry ledger has one marked positive-dimensional modular locus and
no characteristic-zero MW1 equation success.  The script therefore records
all requested features, computes exact finite-field gate rates, and ranks
primitive rows by similarity to the six attempted MW1 sources.
"""

from __future__ import annotations

import argparse
import json
import math
import statistics
from fractions import Fraction
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RANKING = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-source-ranking-v2.json"
DEFAULT_POLES = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-rank1-section-poles-v1.json"
DEFAULT_ATTEMPTS = ROOT / "elkies-k3/data/lattice-foundry/source-equation-attempts-v1.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-empirical-source-ranking-v1.json"

MODEL_FEATURES = (
    "support_count",
    "minimum_pole_order",
    "log_height",
    "component_correction_total",
    "component_correction_max",
    "nonzero_component_corrections",
    "repeated_fibre_excess",
    "semistable_indicator",
    "log_discriminant",
)


def load(path: Path) -> Any:
    return json.loads(path.read_text())


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def source_key(source_artifact: str, source_id: str) -> str:
    return f"{source_artifact}::{source_id}"


def fraction_string(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def solve_fraction_matrix(matrix: list[list[Fraction]], vector: list[Fraction]) -> list[Fraction]:
    size = len(matrix)
    augmented = [row[:] + [vector[index]] for index, row in enumerate(matrix)]
    for column in range(size):
        pivot = next(row for row in range(column, size) if augmented[row][column])
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        scale = augmented[column][column]
        augmented[column] = [entry / scale for entry in augmented[column]]
        for row in range(size):
            if row == column or not augmented[row][column]:
                continue
            scale = augmented[row][column]
            augmented[row] = [
                left - scale * right
                for left, right in zip(augmented[row], augmented[column])
            ]
    return [augmented[index][-1] for index in range(size)]


def connected_components(root_gram: list[list[Fraction]]) -> list[list[int]]:
    unseen = set(range(len(root_gram)))
    components = []
    while unseen:
        first = min(unseen)
        unseen.remove(first)
        stack = [first]
        component = []
        while stack:
            left = stack.pop()
            component.append(left)
            for right in sorted(unseen):
                if root_gram[left][right]:
                    unseen.remove(right)
                    stack.append(right)
        components.append(sorted(component))
    return sorted(components, key=lambda indices: indices[0])


def determinant_fraction(matrix: list[list[Fraction]]) -> Fraction:
    work = [row[:] for row in matrix]
    determinant = Fraction(1)
    for column in range(len(work)):
        pivot = next((row for row in range(column, len(work)) if work[row][column]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            determinant *= -1
        value = work[column][column]
        determinant *= value
        for row in range(column + 1, len(work)):
            scale = work[row][column] / value
            for index in range(column, len(work)):
                work[row][index] -= scale * work[column][index]
    return determinant


def component_type(rank: int, determinant: int) -> str:
    if determinant == rank + 1:
        return f"A{rank}"
    if rank >= 4 and determinant == 4:
        return f"D{rank}"
    exceptional = {(6, 3): "E6", (7, 2): "E7", (8, 1): "E8"}
    return exceptional.get((rank, determinant), f"rank{rank}_det{determinant}")


def kodaira_discriminant_order(root_type: str) -> int | None:
    if root_type.startswith("A") and root_type[1:].isdigit():
        return int(root_type[1:]) + 1
    if root_type.startswith("D") and root_type[1:].isdigit():
        return int(root_type[1:]) + 2
    return {"E6": 8, "E7": 9, "E8": 10}.get(root_type)


def repeated_fibre_features(root_types: list[str]) -> dict[str, Any]:
    orders: dict[int, int] = {}
    for root_type in root_types:
        order = kodaira_discriminant_order(root_type)
        if order is not None:
            orders[order] = orders.get(order, 0) + 1
    repeated = [
        {"kodaira_discriminant_order": order, "support_count": count}
        for order, count in sorted(orders.items())
        if count > 1
    ]
    return {
        "repeated_fibre_multiplicities": repeated,
        "repeated_fibre_excess": sum(row["support_count"] - 1 for row in repeated),
    }


def root_features(source: dict[str, Any]) -> dict[str, Any]:
    gram = [[Fraction(value) for value in row] for row in source["root_adapted_gram"]]
    root_gram = [row[:16] for row in gram[:16]]
    cross = [gram[index][16] for index in range(16)]
    corrections = []
    for indices in connected_components(root_gram):
        block = [[root_gram[left][right] for right in indices] for left in indices]
        block_cross = [cross[index] for index in indices]
        solution = solve_fraction_matrix(block, block_cross)
        correction = sum(left * right for left, right in zip(block_cross, solution))
        determinant = abs(int(determinant_fraction(block)))
        root_type = component_type(len(indices), determinant)
        corrections.append(
            {
                "root_type": root_type,
                "kodaira_discriminant_order": kodaira_discriminant_order(root_type),
                "correction": fraction_string(correction),
            }
        )
    values = [Fraction(row["correction"]) for row in corrections]
    return {
        "component_corrections": corrections,
        "component_correction_total": fraction_string(sum(values, Fraction(0))),
        "component_correction_max": fraction_string(max(values, default=Fraction(0))),
        "nonzero_component_corrections": sum(value != 0 for value in values),
        **repeated_fibre_features([row["root_type"] for row in corrections]),
    }


def attempt_summary(attempt: dict[str, Any]) -> dict[str, Any]:
    eligible = [trial for trial in attempt["trials"] if trial["model_eligible"]]
    fibres = sum(bool(trial["fibre_success"]) for trial in eligible)
    marked = sum(bool(trial["marked_section_success"]) for trial in eligible)
    fibre_eligible = [trial for trial in eligible if trial["fibre_success"]]
    marked_degrees = [
        int(trial["extension_degree"])
        for trial in attempt["trials"]
        if trial["marked_section_success"]
    ]
    locus = bool(attempt["positive_dimensional_marked_locus"])
    stages = []
    for trial in eligible:
        stage = 0
        if trial["fibre_success"]:
            stage = 1
        if trial["marked_section_success"]:
            stage = 2
            if locus:
                stage = 3
        stages.append(stage)
    progress = sum(stages) / (3 * len(stages)) if stages else None
    return {
        "eligible_finite_field_trials": len(eligible),
        "finite_field_fibre_successes": fibres,
        "finite_field_fibre_success_rate": fibres / len(eligible) if eligible else None,
        "finite_field_marked_section_successes": marked,
        "finite_field_marked_section_success_rate": marked / len(eligible) if eligible else None,
        "finite_field_marked_section_success_rate_given_fibre": (
            marked / len(fibre_eligible) if fibre_eligible else None
        ),
        "minimum_successful_extension_degree": min(marked_degrees) if marked_degrees else None,
        "positive_dimensional_marked_locus": locus,
        "equation_success": bool(attempt["equation_success"]),
        "observed_precursor_progress": progress,
    }


def median_scale(values: list[float]) -> float:
    median = statistics.median(values)
    deviations = [abs(value - median) for value in values]
    mad = statistics.median(deviations)
    if mad > 0:
        return mad
    spread = max(values) - min(values)
    return spread / 4 if spread else 1.0


def model_vector(row: dict[str, Any]) -> dict[str, float]:
    return {
        "support_count": float(row["support_count"]),
        "minimum_pole_order": float(row["minimum_pole_order"]),
        "log_height": math.log1p(float(Fraction(row["height"]))),
        "component_correction_total": float(Fraction(row["component_correction_total"])),
        "component_correction_max": float(Fraction(row["component_correction_max"])),
        "nonzero_component_corrections": float(row["nonzero_component_corrections"]),
        "repeated_fibre_excess": float(row["repeated_fibre_excess"]),
        "semistable_indicator": float(bool(row["semistable_compatible"])),
        "log_discriminant": math.log1p(float(row["discriminant"])),
    }


def distance(left: dict[str, float], right: dict[str, float], scales: dict[str, float]) -> float:
    return sum(
        abs(left[name] - right[name]) / scales[name]
        for name in MODEL_FEATURES
    ) / len(MODEL_FEATURES)


def analogue_score(
    vector: dict[str, float],
    training: list[tuple[str, dict[str, float], float]],
    scales: dict[str, float],
) -> tuple[float, list[dict[str, Any]]]:
    neighbours = []
    for key, other, outcome in training:
        separation = distance(vector, other, scales)
        weight = 1.0 / (0.25 + separation) ** 2
        neighbours.append(
            {
                "source_key": key,
                "distance": separation,
                "weight": weight,
                "observed_precursor_progress": outcome,
            }
        )
    global_mean = sum(outcome for _, _, outcome in training) / len(training)
    numerator = 2.0 * global_mean + sum(
        row["weight"] * row["observed_precursor_progress"] for row in neighbours
    )
    denominator = 2.0 + sum(row["weight"] for row in neighbours)
    return numerator / denominator, sorted(neighbours, key=lambda row: row["distance"])


def build_payload(ranking_path: Path, poles_path: Path, attempts_path: Path) -> dict[str, Any]:
    ranking = load(ranking_path)
    poles = load(poles_path)
    attempt_ledger = load(attempts_path)

    pole_map = {
        source_key(row["source_artifact"], row["source_id"]): row
        for row in poles["sources"]
    }
    source_cache: dict[str, dict[str, Any]] = {}

    def source_entry(candidate: dict[str, Any]) -> dict[str, Any]:
        artifact = candidate["source_artifact"]
        if artifact not in source_cache:
            payload = load(ROOT / artifact)
            source_cache[artifact] = {row["source_id"]: row for row in payload["sources"]}
        return source_cache[artifact][candidate["source_id"]]

    attempt_map = {}
    for attempt in attempt_ledger["attempts"]:
        key = source_key(attempt["source_artifact"], attempt["source_id"])
        attempt_map[key] = {**attempt_summary(attempt), "trial_records": attempt["trials"]}

    feature_rows = []
    for candidate in ranking["candidates"]:
        if int(candidate["source_mw_rank"]) != 1:
            continue
        key = source_key(candidate["source_artifact"], candidate["source_id"])
        entry = source_entry(candidate)
        source = entry["source"]
        pole = pole_map[key]
        primitive = bool(source["root_lattice_primitive"])
        if primitive:
            root = root_features(source)
        else:
            root = {
                "component_corrections": None,
                "component_correction_total": None,
                "component_correction_max": None,
                "nonzero_component_corrections": None,
                **repeated_fibre_features(
                    [component["type"] for component in source["root_components"]]
                ),
            }
        minimum_pole = pole["minimum_section_pole_order"]
        height = (
            str(source["mw_height_gram"][0][0])
            if source["mw_height_gram"] is not None
            else None
        )
        row = {
            "source_key": key,
            "source_artifact": candidate["source_artifact"],
            "source_id": candidate["source_id"],
            "source_gram_sha256": source["gram_sha256"],
            "ns_id": candidate["ns_id"],
            "mw_rank": 1,
            "root_type": candidate["source_root_type"],
            "support_count": int(candidate["reducible_fibre_support_count"]),
            "minimum_pole_order": minimum_pole,
            "minimum_pole_multiple": pole.get("minimizing_multiple"),
            "height": height,
            "root_lattice_primitive": primitive,
            "semistable_compatible": bool(candidate["semistable_configuration_compatible"]),
            "discriminant": int(candidate["determinant"]),
            **root,
            "attempt_evidence": attempt_map.get(key),
        }
        feature_rows.append(row)

    eligible = [
        row for row in feature_rows
        if row["root_lattice_primitive"] and row["minimum_pole_order"] is not None
    ]
    vectors = {row["source_key"]: model_vector(row) for row in eligible}
    scales = {
        name: median_scale([vector[name] for vector in vectors.values()])
        for name in MODEL_FEATURES
    }
    training = []
    for row in eligible:
        evidence = row["attempt_evidence"]
        if evidence is None or evidence["observed_precursor_progress"] is None:
            continue
        training.append(
            (
                row["source_key"],
                vectors[row["source_key"]],
                float(evidence["observed_precursor_progress"]),
            )
        )
    if len(training) != 6:
        raise ValueError(f"expected six MW1 training sources, found {len(training)}")

    ranked = []
    for row in eligible:
        score, neighbours = analogue_score(vectors[row["source_key"]], training, scales)
        ranked.append(
            {
                **row,
                "empirical_equation_precursor_score": score,
                "nearest_attempted_analogues": neighbours[:3],
            }
        )
    ranked.sort(
        key=lambda row: (
            -row["empirical_equation_precursor_score"],
            row["minimum_pole_order"],
            row["support_count"],
            float(Fraction(row["height"])),
            row["source_key"],
        )
    )
    for index, row in enumerate(ranked, 1):
        row["empirical_rank"] = index

    diversified = []
    seen_ns = set()
    for row in ranked:
        if row["ns_id"] in seen_ns:
            continue
        diversified.append(row)
        seen_ns.add(row["ns_id"])
        if len(diversified) == 10:
            break

    recommended_next = []
    seen_ns = set()
    for row in ranked:
        evidence = row["attempt_evidence"]
        if evidence is not None and not evidence["positive_dimensional_marked_locus"]:
            continue
        if row["ns_id"] in seen_ns:
            continue
        recommended_next.append(row)
        seen_ns.add(row["ns_id"])
        if len(recommended_next) == 10:
            break

    leave_one_out = []
    for key, vector, outcome in training:
        reduced = [item for item in training if item[0] != key]
        prediction, _ = analogue_score(vector, reduced, scales)
        leave_one_out.append(
            {
                "source_key": key,
                "observed_precursor_progress": outcome,
                "predicted_precursor_score": prediction,
                "absolute_error": abs(outcome - prediction),
            }
        )

    unknown = [row for row in feature_rows if row not in eligible]
    return {
        "schema": "elkies-k3-lattice-foundry-empirical-source-ranking-v1",
        "status": "PASS_REPRODUCIBLE_DESCRIPTIVE_PRECURSOR_RANKING",
        "inputs": {
            "source_ranking": relative(ranking_path),
            "rank1_section_poles": relative(poles_path),
            "attempt_ledger": relative(attempts_path),
        },
        "feature_coverage": {
            "mw1_candidates": len(feature_rows),
            "primitive_candidates_with_complete_model_features": len(eligible),
            "nonprimitive_candidates_left_unranked": len(unknown),
            "attempted_mw1_sources": len(training),
            "characteristic_zero_equation_successes_in_training": 0,
            "positive_dimensional_marked_loci_in_training": sum(
                bool(attempt["positive_dimensional_marked_locus"])
                for attempt in attempt_ledger["attempts"]
            ),
        },
        "model": {
            "name": "shrunk_inverse_distance_equation_precursor_analogue",
            "target": "Mean exact-chart progress through no fibre=0, fibre=1/3, marked section=2/3, positive-dimensional marked locus=1.",
            "features": list(MODEL_FEATURES),
            "feature_scales_median_absolute_deviation": scales,
            "distance": "Equal-weight mean standardized absolute difference.",
            "analogue_weight": "1/(0.25+distance)^2",
            "shrinkage": "Two pseudo-observations at the six-source global mean.",
            "interpretation": "Relative triage score only; it is not a calibrated probability that a characteristic-zero equation exists or will be found.",
        },
        "training_sources": [
            next(row for row in ranked if row["source_key"] == key)
            for key, _, _ in training
        ],
        "leave_one_source_out_diagnostic": {
            "rows": leave_one_out,
            "mean_absolute_error": sum(row["absolute_error"] for row in leave_one_out) / len(leave_one_out),
            "warning": "The only positive-dimensional marked source cannot be validated from the remaining all-semistable attempts; this is evidence scarcity, not model performance.",
        },
        "top_ten": ranked[:10],
        "diversified_top_ten_one_per_ns": diversified,
        "recommended_next_ten": recommended_next,
        "ranked_primitive_candidates": ranked,
        "unranked_nonprimitive_candidates": unknown,
        "proof_boundary": "All lattice features are exact consequences of the stored candidate and pole artifacts. Finite-field rates use only model_eligible charts from the curated attempt ledger. The empirical score ranks similarity to observed equation precursors; with six selected attempts, one marked-locus success, and no MW1 equation success, it cannot estimate absolute equation probabilities or support theorem claims.",
        "reproduce": "python3 elkies-k3/scripts/rank_lattice_foundry_empirical_sources.py --check",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ranking", type=Path, default=DEFAULT_RANKING)
    parser.add_argument("--poles", type=Path, default=DEFAULT_POLES)
    parser.add_argument("--attempts", type=Path, default=DEFAULT_ATTEMPTS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    payload = build_payload(args.ranking, args.poles, args.attempts)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.exists():
            raise SystemExit(f"missing output for --check: {args.output}")
        if args.output.read_text() != rendered:
            raise SystemExit(f"stale output: {args.output}")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(
        "EMPIRICAL_SOURCE_RANKING|"
        f"candidates={payload['feature_coverage']['mw1_candidates']}|"
        f"ranked={payload['feature_coverage']['primitive_candidates_with_complete_model_features']}|"
        f"attempts={payload['feature_coverage']['attempted_mw1_sources']}|"
        f"leader={payload['top_ten'][0]['source_key']}|"
        f"status={payload['status']}"
    )


if __name__ == "__main__":
    main()
