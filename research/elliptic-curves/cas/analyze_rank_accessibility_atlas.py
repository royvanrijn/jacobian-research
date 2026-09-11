#!/usr/bin/env sage -python
"""Descriptive, point-level analysis of the frozen Curve302/11952 atlas.

The unit of presentation is a known exceptional direction (or an adjacent
known-subgroup transition), never an individual anchor/point pair.  Thus the
many exact anchor evaluations are used to obtain minima but are not treated as
independent observations or fed into p-values.
"""

from __future__ import annotations

import argparse
import gzip
from hashlib import sha256
import json
from math import log2
from pathlib import Path
from statistics import mean, median
from typing import Iterable

from sympy import Matrix


ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / "artifacts" / "generated-results" / "elliptic-curves"
ATLAS = ART / "rank_accessibility_atlas_curve302_11952_v1.json.gz"
ATLAS_SUMMARY = ART / "rank_accessibility_atlas_curve302_11952_v1.summary.json"
OUTPUT = ART / "rank_accessibility_atlas_analysis_v1.json"
REPORT = ART / "rank_accessibility_atlas_analysis_v1.md"


def canonical(value) -> bytes:
    return (json.dumps(value, sort_keys=True, indent=2, ensure_ascii=True) + "\n").encode()


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def summary(values: Iterable[float]) -> dict:
    values = list(values)
    if not values:
        return {"count": 0}
    return {
        "count": len(values), "minimum": min(values), "median": median(values),
        "mean": mean(values), "maximum": max(values),
    }


def survival(values: Iterable[float]) -> list[dict]:
    """The full empirical A(c)=#{value >= c}, evaluated at every observed c."""
    values = sorted(values)
    return [{"c": value, "A": sum(candidate >= value for candidate in values)} for value in sorted(set(values))]


def residual_r2(observations: list[dict]) -> float | None:
    """R^2 for kappa explained by half the residual canonical height, with intercept."""
    x = [0.5 * entry["residual_canonical_height"] for entry in observations]
    y = [entry["kappa_bits"] for entry in observations]
    if len(x) < 2 or len(set(x)) < 2 or len(set(y)) < 2:
        return None
    xbar, ybar = mean(x), mean(y)
    slope = sum((left - xbar) * (right - ybar) for left, right in zip(x, y)) / sum((left - xbar) ** 2 for left in x)
    fitted = [ybar + slope * (left - xbar) for left in x]
    sse = sum((right - estimate) ** 2 for right, estimate in zip(y, fitted))
    sst = sum((right - ybar) ** 2 for right in y)
    return 1.0 - sse / sst


def finite_x_height(row: dict) -> float:
    numerator = abs(int(row["raw_x_map"]["N_homogeneous"]))
    denominator = abs(int(row["raw_x_map"]["F_homogeneous"]))
    cancellation = int(row["finite_cancellation_gcd"])
    return log2(max(numerator // cancellation, denominator // cancellation))


def kappa_observation(panel: dict, bank: str, class_of_word) -> dict:
    winner = panel[bank]["measurement"]["minima"]["kappa_parameter_height_bits"]
    kappa = float(winner["parameter"]["height_log2"])
    residual = float(winner["residual_canonical_height"])
    x_quarter = 0.25 * finite_x_height(winner)
    finite = float(winner["lambda_finite_bits"])
    arch = float(winner["real_distortion_bits"])
    # Exactly: kappa = 1/4 h_x + arch + finite.  The model term compares
    # h_x/4 with the leading 1/2 hhat residual proxy.
    model = x_quarter - 0.5 * residual
    reconstruction_error = kappa - (0.5 * residual + model + arch + finite)
    if abs(reconstruction_error) > 1e-8:
        raise ArithmeticError("cover-height decomposition no longer reconstructs kappa")
    return {
        "kappa_bits": kappa,
        "residual_canonical_height": residual,
        "quarter_x_height": x_quarter,
        "residual_half_height_component": 0.5 * residual,
        "model_height_difference_component": model,
        "finite_cancellation_component": finite,
        "archimedean_distortion_component": arch,
        "reconstruction_error": reconstruction_error,
        "winning_anchor_word_in_D": winner["anchor_word_in_D"],
        "winning_anchor_mod2_quotient_class": class_of_word(winner["anchor_word_in_D"]),
        "normalization_maximum_coefficient_bits": winner["normalization_maximum_coefficient_bits"],
    }


def class_map(curve: dict):
    generic = [entry["word_in_D"] for entry in curve["M17_specialization"]]
    exceptional = [entry["word_in_D"] for entry in curve["exceptional_directions"]]
    basis = Matrix.hstack(*(Matrix(word) for word in [*generic, *exceptional]))
    if abs(int(basis.det())) != 1:
        raise ArithmeticError("displayed M17/exceptional basis ceased to be unimodular")
    inverse = basis.inv()
    rank = len(exceptional)
    def convert(word):
        coordinates = inverse * Matrix(word)
        if any(value.q != 1 for value in coordinates):
            raise ArithmeticError("unexpected nonintegral displayed quotient coordinates")
        return [int(coordinates[17 + index]) % 2 for index in range(rank)]
    return convert


def class_key(vector: list[int]) -> str:
    return "".join(map(str, vector))


def curve_analysis(curve: dict) -> dict:
    class_of_word = class_map(curve)
    panels = sorted(curve["nested_target_blind_panels"], key=lambda row: (row["target_exceptional_index_zero_based"], row["stage"]))
    observations = []
    by_target: dict[int, list[dict]] = {}
    for panel in panels:
        native = kappa_observation(panel, "native_bank", class_of_word)
        generic = kappa_observation(panel, "generic_only_control", class_of_word)
        target_class = class_of_word(panel["target_word_in_D"])
        record = {
            "stage": panel["stage"],
            "target_exceptional_index_zero_based": panel["target_exceptional_index_zero_based"],
            "target_mod2_quotient_class": target_class,
            "native": native,
            "generic_control": generic,
            "exceptional_bank_advantage_G_bits": generic["kappa_bits"] - native["kappa_bits"],
        }
        observations.append(record)
        by_target.setdefault(record["target_exceptional_index_zero_based"], []).append(record)

    initial = [entry for entry in observations if entry["stage"] == 0]
    initial_differences = [entry["exceptional_bank_advantage_G_bits"] for entry in initial]
    # The protocol makes M17 native and generic-only banks byte-for-byte the
    # same.  Making that null comparison explicit is more informative than a
    # fabricated 'cheap' threshold.
    if any(abs(value) > 1e-10 for value in initial_differences):
        raise ArithmeticError("initial equal-bank control unexpectedly differs")

    transitions = []
    for target, path in sorted(by_target.items()):
        path.sort(key=lambda row: row["stage"])
        for before, after in zip(path, path[1:]):
            added = before["stage"]
            if after["stage"] != added + 1:
                raise ArithmeticError("nested panel stage is missing")
            native_delta = before["native"]["kappa_bits"] - after["native"]["kappa_bits"]
            generic_delta = before["generic_control"]["kappa_bits"] - after["generic_control"]["kappa_bits"]
            transitions.append({
                "target_exceptional_index_zero_based": target,
                "added_exceptional_index_zero_based": added,
                "target_mod2_quotient_class": before["target_mod2_quotient_class"],
                "added_mod2_quotient_class": class_of_word(curve["exceptional_directions"][added]["word_in_D"]),
                "native_kappa_collapse_delta_bits": native_delta,
                "generic_control_kappa_collapse_delta_bits": generic_delta,
                "exceptional_excess_collapse_bits": native_delta - generic_delta,
                "G_before_bits": before["exceptional_bank_advantage_G_bits"],
                "G_after_bits": after["exceptional_bank_advantage_G_bits"],
            })

    native_deltas = [entry["native_kappa_collapse_delta_bits"] for entry in transitions]
    positive = [value for value in native_deltas if value > 0]
    advantage = [entry["exceptional_bank_advantage_G_bits"] for entry in observations if entry["stage"] > 0]
    target_strata = {}
    for key in sorted({class_key(entry["target_mod2_quotient_class"]) for entry in observations}):
        selected = [entry for entry in observations if class_key(entry["target_mod2_quotient_class"]) == key]
        target_strata[key] = {
            "observation_count": len(selected),
            "kappa_bits": summary(entry["native"]["kappa_bits"] for entry in selected),
            "G_bits": summary(entry["exceptional_bank_advantage_G_bits"] for entry in selected),
            "residual_component_bits": summary(entry["native"]["residual_half_height_component"] for entry in selected),
            "finite_cancellation_bits": summary(entry["native"]["finite_cancellation_component"] for entry in selected),
        }
    transition_strata = {}
    for key in sorted({class_key(entry["added_mod2_quotient_class"]) for entry in transitions}):
        selected = [entry for entry in transitions if class_key(entry["added_mod2_quotient_class"]) == key]
        transition_strata[key] = {
            "transition_count": len(selected),
            "native_delta_bits": summary(entry["native_kappa_collapse_delta_bits"] for entry in selected),
            "excess_delta_bits": summary(entry["exceptional_excess_collapse_bits"] for entry in selected),
        }
    native_observations = [entry["native"] for entry in observations]
    return {
        "label": curve["label"],
        "known_exceptional_dimensions": curve["known_quotient"]["known_exceptional_quotient_rank"],
        "analysis_unit": "known exceptional target or adjacent known-subgroup transition",
        "initial_accessibility": {
            "cheap_threshold": None,
            "reason": "No threshold was selected; retain the complete empirical kappa distribution.",
            "native_kappa_bits": summary(entry["native"]["kappa_bits"] for entry in initial),
            "native_residual_canonical_height": summary(entry["native"]["residual_canonical_height"] for entry in initial),
            "equal_generic_control_difference_bits": summary(initial_differences),
            "kappa_survival_curve": survival(entry["native"]["kappa_bits"] for entry in initial),
        },
        "winning_chart_decomposition": {
            "observation_count": len(native_observations),
            "residual_half_height_component_bits": summary(entry["residual_half_height_component"] for entry in native_observations),
            "model_height_difference_component_bits": summary(entry["model_height_difference_component"] for entry in native_observations),
            "finite_cancellation_component_bits": summary(entry["finite_cancellation_component"] for entry in native_observations),
            "archimedean_distortion_component_bits": summary(entry["archimedean_distortion_component"] for entry in native_observations),
            "residual_height_explained_fraction_R2": residual_r2(native_observations),
        },
        "avalanche": {
            "transition_count": len(transitions),
            "native_kappa_delta_bits": summary(native_deltas),
            "largest_kappa_collapse_bits": max(native_deltas) if native_deltas else None,
            "median_positive_kappa_collapse_bits": median(positive) if positive else None,
            "positive_collapse_count": len(positive),
            "collapse_survival_A_of_c": survival(native_deltas),
            "collapse_AUC_bit_directions": sum(positive),
            "generic_control_delta_bits": summary(entry["generic_control_kappa_collapse_delta_bits"] for entry in transitions),
            "exceptional_excess_collapse_bits": summary(entry["exceptional_excess_collapse_bits"] for entry in transitions),
        },
        "exceptional_vs_generic_information": {
            "G_definition": "kappa_generic_control - kappa_exceptional_native",
            "enlarged_bank_G_bits": summary(advantage),
            "positive_G_count": sum(value > 0 for value in advantage),
            "G_survival_curve": survival(advantage),
        },
        "mod2_quotient_stratification": {
            "target_class": target_strata,
            "added_direction_class": transition_strata,
            "boundary": "Classes label displayed D/M17 mod-2 coordinates. They are stratification labels, not claimed predictors or independent samples.",
        },
        "point_level_records": observations,
        "adjacent_transition_records": transitions,
    }


def report(analysis: dict) -> str:
    rows = analysis["curves"]
    def cell(curve, path):
        value = curve
        for key in path:
            value = value[key]
        return "—" if value is None else f"{value:.3f}" if isinstance(value, float) else str(value)
    lines = [
        "# Rank-accessibility atlas: descriptive comparison\n",
        "The unit is a known exceptional direction or an adjacent subgroup transition; no anchor-level p-values are computed. `11952` remains a rank-at-least-25 fibre.\n",
        "| Metric | Curve302 | 11952 |",
        "|---|---:|---:|",
        "| Known exceptional dimensions | " + " | ".join(str(row["known_exceptional_dimensions"]) for row in rows) + " |",
        "| Initially cheap | not thresholded | not thresholded |",
        "| Initial median κ (bits) | " + " | ".join(cell(row, ["initial_accessibility", "native_kappa_bits", "median"]) for row in rows) + " |",
        "| Initial generic-control difference (bits) | " + " | ".join(cell(row, ["initial_accessibility", "equal_generic_control_difference_bits", "maximum"]) for row in rows) + " |",
        "| Median enlarged exceptional-bank advantage G (bits) | " + " | ".join(cell(row, ["exceptional_vs_generic_information", "enlarged_bank_G_bits", "median"]) for row in rows) + " |",
        "| Positive G panels / enlarged panels | " + " | ".join(str(row["exceptional_vs_generic_information"]["positive_G_count"]) + "/" + str(row["exceptional_vs_generic_information"]["enlarged_bank_G_bits"]["count"]) for row in rows) + " |",
        "| Maximum exceptional-bank advantage G (bits) | " + " | ".join(cell(row, ["exceptional_vs_generic_information", "enlarged_bank_G_bits", "maximum"]) for row in rows) + " |",
        "| Largest κ collapse (bits) | " + " | ".join(cell(row, ["avalanche", "largest_kappa_collapse_bits"]) for row in rows) + " |",
        "| Median positive κ collapse (bits) | " + " | ".join(cell(row, ["avalanche", "median_positive_kappa_collapse_bits"]) for row in rows) + " |",
        "| Collapse AUC (bit-directions) | " + " | ".join(cell(row, ["avalanche", "collapse_AUC_bit_directions"]) for row in rows) + " |",
        "| Residual-height R² for κ winners | " + " | ".join(cell(row, ["winning_chart_decomposition", "residual_height_explained_fraction_R2"]) for row in rows) + " |",
        "| Median finite-cancellation term (bits) | " + " | ".join(cell(row, ["winning_chart_decomposition", "finite_cancellation_component_bits", "median"]) for row in rows) + " |",
        "\nThe initial equal-bank difference is exactly zero by design: both use the same M17 signed unit/pair bank. The JSON contains the full survival curves, per-winning-chart decomposition, and mod-2 class strata.\n",
    ]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true", help="recompute and compare deterministic outputs")
    args = parser.parse_args()
    payload = json.loads(gzip.decompress(ATLAS.read_bytes()))
    if payload["status"] != "COMPLETE_FINITE_ATLAS_MEASUREMENT":
        raise ArithmeticError("atlas is not a complete frozen measurement")
    output = {
        "schema": "elliptic-curves.rank-accessibility-atlas-analysis.v1",
        "status": "DESCRIPTIVE_COMPLETE",
        "claim_boundary": "Descriptive analysis of known displayed lattices only. No point search, rank upper bound, saturation assertion in E(Q), causal mechanism, or statistical significance claim.",
        "input_hashes": {
            str(ATLAS.relative_to(ROOT)): digest(ATLAS),
            str(ATLAS_SUMMARY.relative_to(ROOT)): digest(ATLAS_SUMMARY),
            str(Path(__file__).relative_to(ROOT)): digest(Path(__file__)),
        },
        "atlas_payload_sha256": read_summary_payload_hash(),
        "curves": [curve_analysis(curve) for curve in payload["curves"]],
    }
    data = canonical(output)
    markdown = report(output)
    if args.check:
        if not OUTPUT.is_file() or OUTPUT.read_bytes() != data or not REPORT.is_file() or REPORT.read_text() != markdown:
            raise SystemExit("RANK_ACCESSIBILITY_ANALYSIS_CHECK|status=FAIL")
        print("RANK_ACCESSIBILITY_ANALYSIS_CHECK|status=PASS")
        return
    OUTPUT.write_bytes(data)
    REPORT.write_text(markdown)
    print(f"RANK_ACCESSIBILITY_ANALYSIS|status=PASS|output={OUTPUT}")


def read_summary_payload_hash() -> str:
    return json.loads(ATLAS_SUMMARY.read_text())["payload_sha256"]


if __name__ == "__main__":
    main()
