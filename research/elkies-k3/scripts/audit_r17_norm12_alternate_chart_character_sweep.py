#!/usr/bin/env python3
"""Audit exact character searches across all alternate-Q80 norm-12 charts."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
OUTPUT = GENERATED / "elkies-k3-r17-norm12-alternate-chart-character-sweep-v1.json"

INHERITED_COUNTS = {
    "11952": 121,
    "08ab4": 131,
    "08f72": 86,
    "091e4": 155,
    "135b7": 127,
    "10f72": 118,
    "1183a": 120,
    "09952": 125,
    "098fc": 95,
    "0ae21": 120,
}
FULL_LABELS = tuple(INHERITED_COUNTS)
HALVED_LABELS = ("135b7", "10f72", "09952", "0ae21")
NEGATIVE_STATUS = "PASS_EXACT_NO_EQUAL_COVER_OR_THREE_CHARACTER_CLOSURE"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def load(path: Path):
    return json.loads(path.read_text())


def inherited_cover_path(label: str) -> Path:
    suffix = "covers-v1.json" if label == "11952" else "covers-only-v1.json"
    return GENERATED / f"elkies-k3-r17-norm12-{label}-inherited-bisection-{suffix}"


def direct_path(label: str) -> Path:
    suffix = "-saturated" if label in HALVED_LABELS else ""
    return GENERATED / (
        f"elkies-k3-r17-norm12-orbit{label}-direct-fibration{suffix}-v1.json"
    )


def closure_path(label: str, inherited: bool) -> Path:
    layer = "inherited" if inherited else "streaming"
    return GENERATED / f"elkies-k3-r17-norm12-{label}-{layer}-character-closure-v1.json"


def require_negative_closure(payload, expected_count: int, label: str) -> None:
    if payload["status"] != NEGATIVE_STATUS:
        raise ArithmeticError(f"{label}: character search is not an exact negative")
    if (
        payload["character_count"] != expected_count
        or payload["distinct_polynomial_support_count"] != expected_count
        or payload["equal_cover_collision_count"] != 0
        or payload["three_character_relation_count"] != 0
        or payload["committed_character_catalog"]["comparison_count"] != 12
        or payload["committed_character_catalog"][
            "formal_variable_rename_product_match_count"
        ]
        != 0
    ):
        raise ArithmeticError(f"{label}: negative character counts changed")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    inputs = []
    inherited_records = []
    for label, expected_count in INHERITED_COUNTS.items():
        cover_file = inherited_cover_path(label)
        closure_file = closure_path(label, inherited=True)
        cover = load(cover_file)
        closure = load(closure_file)
        if len(cover["bisections"]) != expected_count:
            raise ArithmeticError(f"{label}: inherited cover count changed")
        require_negative_closure(closure, expected_count, f"{label} inherited")
        inherited_records.append(
            {
                "label": f"norm12-orbit-{label}",
                "cover_count": expected_count,
                "distinct_character_count": expected_count,
                "collision_count": 0,
                "internal_relation_count": 0,
                "committed_catalog_match_count": 0,
            }
        )
        inputs.extend((cover_file, closure_file))

    full_records = []
    for label in FULL_LABELS:
        direct_file = direct_path(label)
        closure_file = closure_path(label, inherited=False)
        direct = load(direct_file)
        closure = load(closure_file)
        if (
            direct["status"]
            != "PASS_EXACT_DIRECT_TWO_NEIGHBOR_EQUATION_FRAME_AND_SECTIONS"
            or direct["weierstrass_model"]["fibre_configuration"] != "24 I1"
            or direct["sections"]["status"] != "PASS_EXACT_SATURATED_RANK17_BASIS"
            or direct["sections"].get("index_in_saturated_mw_lattice", 1) != 1
        ):
            raise ArithmeticError(f"{label}: full chart is not saturated 24I1 rank 17")
        require_negative_closure(closure, 39147, f"{label} full")
        full_records.append(
            {
                "label": f"norm12-orbit-{label}",
                "cover_count": 39147,
                "distinct_character_count": 39147,
                "collision_count": 0,
                "internal_relation_count": 0,
                "committed_catalog_match_count": 0,
            }
        )
        inputs.extend((direct_file, closure_file))

    saturation_records = []
    for label in HALVED_LABELS:
        finite_index_file = GENERATED / (
            f"elkies-k3-r17-norm12-orbit{label}-direct-fibration-v1.json"
        )
        saturated_file = direct_path(label)
        plan_file = GENERATED / (
            f"elkies-k3-r17-norm12-{label}-direct-section-basis-obstruction-v1.json"
        )
        finite_index = load(finite_index_file)
        saturated = load(saturated_file)
        plan = load(plan_file)
        if (
            finite_index["status"]
            != "PASS_EXACT_DIRECT_TWO_NEIGHBOR_EQUATION_FRAME_AND_RANK17_SUBLATTICE"
            or finite_index["weierstrass_model"]["fibre_configuration"] != "24 I1"
            or finite_index["sections"]["status"]
            != "PASS_EXACT_RANK17_FINITE_INDEX_SUBLATTICE"
            or finite_index["sections"]["index_in_saturated_mw_lattice"] != 2
            or plan["status"] != "NO_UNIMODULAR_MARKING_IN_SEARCHED_CURVES"
            or plan["height_bound"] != 12
            or plan["old_degree_one_span_rank"] != 17
            or plan["old_degree_one_lattice_saturation_index"] != 2
            or abs(plan["selected_child_coordinate_determinant"]) != 2
            or plan["rational_bisection_glue_candidate_count"] != 0
            or saturated["status"]
            != "PASS_EXACT_DIRECT_TWO_NEIGHBOR_EQUATION_FRAME_AND_SECTIONS"
            or saturated["sections"]["status"]
            != "PASS_EXACT_SATURATED_RANK17_BASIS"
            or abs(saturated["sections"]["coordinate_matrix_determinant"]) != 1
            or saturated["sections"]["height_gram_determinant"] != 948
            or saturated["sections"]["saturation"]["input_index"] != 2
            or saturated["sections"]["records"][-1]["new_height"] != 4
            or saturated["sections"]["records"][-1]["doubling_verified"] is not True
        ):
            raise ArithmeticError(f"{label}: exact halving saturation changed")
        saturation_records.append(
            {
                "label": f"norm12-orbit-{label}",
                "input_section_sublattice_index": 2,
                "complete_old_degree_one_span_rank": 17,
                "complete_old_degree_one_lattice_saturation_index": 2,
                "committed_rational_bisection_glue_count": 0,
                "exact_rational_half_height": 4,
                "output_section_lattice_index": 1,
            }
        )
        inputs.extend((finite_index_file, plan_file))

    result = {
        "schema": "elkies-k3.r17-norm12-alternate-chart-character-sweep.v1",
        "status": "PASS_EXACT_ALTERNATE_CHART_SWEEP_NO_RANK19_OR_RANK20_CHARACTER_INPUT",
        "inherited_layer": {
            "chart_count": len(inherited_records),
            "cover_count_sum_across_distinct_base_charts": sum(INHERITED_COUNTS.values()),
            "records": inherited_records,
        },
        "complete_smooth_layer": {
            "chart_count": len(full_records),
            "cover_count_sum_across_distinct_base_charts": 39147 * len(full_records),
            "records": full_records,
        },
        "exact_halving_saturations": saturation_records,
        "inputs": {relative(path): digest(path) for path in inputs},
        "proof_boundary": (
            "Equal-cover, internal three-character, and twelve-character catalog "
            "comparisons are exact within each base chart; characters on different "
            "base coordinates are not compared as if they belonged to one function "
            "field. The inherited layer is complete for all ten alternate-Q80 norm-12 "
            "charts. The full 39147-class layer is complete for all ten charts with "
            "saturated explicit equation markings. On four charts the old degree-one "
            "section span first reduces to index two; an exact rational height-four half, "
            "verified by doubling and a determinant-one coordinate matrix, supplies the "
            "missing saturation section."
        ),
        "reproducing_command": (
            ".venv/bin/python elkies-k3/scripts/"
            "audit_r17_norm12_alternate_chart_character_sweep.py --check"
        ),
    }
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.exists() or args.output.read_text() != serialized:
            raise ArithmeticError("stored alternate-chart sweep differs from replay")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized)
    print(
        "R17ALTCHARTSWEEP"
        f"|inherited_charts={len(inherited_records)}"
        f"|inherited_covers={sum(INHERITED_COUNTS.values())}"
        f"|full_charts={len(full_records)}"
        f"|full_covers={39147 * len(full_records)}"
        f"|halved_full_charts={len(saturation_records)}"
        f"|status={result['status']}"
    )


if __name__ == "__main__":
    main()
