#!/usr/bin/env sage-python
"""Retrospectively evaluate the sealed V2 active-subgroup beam design.

This file is intentionally separate from ``design_curve302_v2_active_subgroup_beam``.
It reads the fourteen fixed exceptional directions only after the geometry-only
selector has been sealed, and has no centre-selection or point-search code.

The evaluation transports the *same* entrywise-rounded 31-dimensional metric
used by the immutable 2^14 diagnostic to the off-grid V2 subgroups.  Thus its
values extend the frozen diagnostic's local-CVP definition; they are not new
pointed-quartic searches, exact CVP claims, or prospective policy inputs.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.machinery
import json
import math
from fractions import Fraction as F
from pathlib import Path

import numpy as np
from sage.all import RealField, ZZ, matrix


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
ART = ROOT / "artifacts/generated-results/elliptic-curves"
V2 = ROOT / "artifacts/local/elliptic-curves/adaptive-visibility-cascade-v2"
BEAM = ART / "curve302_v2_active_subgroup_beam_v1.json"
LANDSCAPE = ART / "curve302_exceptional_subgroup_landscape_v1.json"
VISIBILITY = ART / "curve302_residual_visibility_geometry_v1.json"
BASE_SELECTION = V2 / "calibration302/epoch-07/selection.json"
BATCH_AUDITS = (
    V2 / "calibration302/epoch-07/mod2-007.json",
    V2 / "calibration302/epoch-08/mod2-040.json",
    V2 / "calibration302/epoch-09/mod2-029.json",
)
OUTPUT = ART / "curve302_v2_active_subgroup_beam_retrospective_v1.json"
GENERIC_RANK = 17
DIMENSION = 14
FULL_RANK = GENERIC_RANK + DIMENSION


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def read(path: Path) -> dict:
    return json.loads(path.read_text())


def load(name: str):
    return importlib.machinery.SourceFileLoader("v2_beam_evaluation_" + name.replace(".", "_"), str(CAS / name)).load_module()


def qpoint(row) -> tuple[F, F]:
    return tuple(F(value) for value in row)


def metric_sha(gram) -> str:
    return hashlib.sha256(json.dumps([list(map(int, row)) for row in gram.rows()], separators=(",", ":")).encode()).hexdigest()


def primary_target(entry: dict, base17, model, group) -> tuple[F, F]:
    """Undo a chart's M17 translation to recover its fixed diagnostic root."""

    roots = set()
    for row in entry["exact_pointed_quartic_rows"]:
        translated = qpoint(row["exact_target_point_short_model"])
        inverse = group.linear_combination(model, base17, [-int(value) for value in row["target_translation_m17_word"]])
        roots.add(group.short_add(model, translated, inverse))
    if len(roots) != 1:
        raise ArithmeticError("diagnostic chart rows do not encode one target")
    return roots.pop()


def vector_word_for_point(geometry, group, model, ambient_points, high_gram, point) -> list[int]:
    """Recognize and exact-group-law verify a word in the fixed diagnostic basis."""

    pairing = geometry.canonical_height_gram(model, (*ambient_points, point))[0]
    field = high_gram.base_ring()
    column = matrix(field, FULL_RANK, 1, [field(str(pairing[index][FULL_RANK]).strip("()")) for index in range(FULL_RANK)])
    coefficients = high_gram.solve_right(column)
    word = [int(round(float(coefficients[index, 0]))) for index in range(FULL_RANK)]
    error = max(abs(coefficients[index, 0] - word[index]) for index in range(FULL_RANK))
    if error > field("1e-30") or group.linear_combination(model, ambient_points, word) != point:
        raise ArithmeticError("V2 certified point does not have an exact diagnostic-basis word")
    return word


def exact_integer_member(rows, coordinate) -> bool:
    """Whether a diagnostic basis vector is already in the active Z-lattice."""

    row_matrix = matrix(ZZ, rows)
    try:
        solution = row_matrix.transpose().solve_right(matrix(ZZ, FULL_RANK, 1, coordinate))
    except ValueError:
        return False
    return all(value.denominator() == 1 for value in solution)


def decimal(value: float) -> str:
    if value < 0 and value > -1e-7:
        value = 0.0
    if value < 0 or not math.isfinite(value):
        raise ArithmeticError("invalid nonnegative retrospective value")
    return format(value, ".12g")


def percentile(values: list[int], value: int) -> str:
    """Empirical lower-tail percentile in the frozen 16,384-state table."""

    return format(sum(other <= value for other in values) / len(values), ".12g")


def summary(costs: list[int], retained: list[int], active_rows) -> dict:
    remaining = []
    for target in range(DIMENSION):
        coordinate = [0] * FULL_RANK
        coordinate[GENERIC_RANK + target] = 1
        if not exact_integer_member(active_rows, coordinate):
            remaining.append(target)
    if not remaining:
        raise ArithmeticError("evaluation unexpectedly absorbed every diagnostic direction")
    fresh_remaining = [costs[index] for index in remaining]
    retained_remaining = [retained[index] for index in remaining]
    return {
        "remaining_direction_indices": remaining,
        "remaining_direction_count": len(remaining),
        "fresh_max_numerator": max(fresh_remaining),
        "fresh_mean_numerator": decimal(sum(fresh_remaining) / len(fresh_remaining)),
        "retained_from_v2_M24_max_numerator": max(retained_remaining),
        "retained_from_v2_M24_mean_numerator": decimal(sum(retained_remaining) / len(retained_remaining)),
    }


def build() -> dict:
    beam, landscape, visibility = read(BEAM), read(LANDSCAPE), read(VISIBILITY)
    if beam["status"] != "PASS_GEOMETRY_ONLY_BASIS_INVARIANT_BEAM_DESIGN":
        raise ArithmeticError("the geometry-only beam artefact is not sealed")
    if landscape["status"] != "PASS_RETROSPECTIVE_FINITE_SUBGROUP_LANDSCAPE":
        raise ArithmeticError("the frozen 2^14 diagnostic is unavailable")
    if visibility["status"] != "PASS_RETROSPECTIVE_VETTED_VISIBILITY_DIAGNOSTIC":
        raise ArithmeticError("the frozen diagnostic target source is unavailable")
    if landscape["inputs"].get(rel(VISIBILITY)) != sha(VISIBILITY):
        raise ArithmeticError("visibility source no longer matches the frozen 2^14 diagnostic input hash")
    if beam["point_searches_run"] != 0:
        raise ArithmeticError("beam artefact boundary was violated")

    core = load("audit_curve302_exceptional_subgroup_landscape.sage")
    geometry = load("prospective_half_lattice_v3.sage")
    group = load("half_lattice_pointed_sieve.py")
    base_selection = read(BASE_SELECTION)
    base = tuple(qpoint(point) for point in base_selection["basis"])
    model = tuple(F(value) for value in read(BATCH_AUDITS[0])["curve"])
    base17 = base[:GENERIC_RANK]
    target_names = [row["id"] for row in visibility["directions"]]
    targets = tuple(primary_target(row, base17, model, group) for row in visibility["directions"])
    ambient_points = (*base17, *targets)
    ambient, asymmetry = core.rounded_metric(geometry, model, ambient_points)
    if metric_sha(ambient) != landscape["fixed_basis"]["rounded_metric_sha256"]:
        raise ArithmeticError("reconstructed metric differs from the frozen 2^14 diagnostic")
    field = RealField(192)
    height = geometry.canonical_height_gram(model, ambient_points)[0]
    high_gram = matrix(field, [[field(str(value).strip("()")) for value in row] for row in height])

    batch = []
    current = base
    for path in BATCH_AUDITS:
        audit = read(path)
        points = tuple(qpoint(point) for point in audit["independent_points"])
        if audit.get("status") != "COMPLETE_DECLARED_FINITE_AUDIT" or points[:len(current)] != current or len(points) != len(current) + 1:
            raise ArithmeticError("V2 batch certificate does not form the frozen rank-24 prefix chain")
        batch.append(points[-1])
        current = points
    base_words = [vector_word_for_point(geometry, group, model, ambient_points, high_gram, point) for point in base]
    batch_words = [vector_word_for_point(geometry, group, model, ambient_points, high_gram, point) for point in batch]
    if matrix(ZZ, base_words).rank() != 24:
        raise ArithmeticError("V2 M24 basis did not retain rank 24 in the fixed diagnostic lattice")

    ambient_array = np.array([[int(value) for value in row] for row in ambient.rows()], dtype=np.int64)
    base_prepared = core.prepare_subgroup(ambient, ambient_array, base_words)
    base_costs = [core.subgroup_score(ambient, ambient_array, base_prepared, target)["fresh_numerator"] for target in range(DIMENSION)]
    base_summary = summary(base_costs, base_costs, base_words)
    # The only nested V2 atlas available to a candidate is M24 plus that
    # candidate; taking its minimum records the finite cumulative value C^*.
    # The full 2^14 table is used separately for its fixed metric and target
    # distribution, never for selection.
    frozen_values = [[] for _ in range(DIMENSION)]
    for state in landscape["subset_states"]:
        for target, value in enumerate(state["stage_local_numerators"]):
            if value is not None:
                frozen_values[target].append(int(value))
    rows = []
    selector_lookup = {entry["branch_mask"]: index + 1 for index, entry in enumerate(beam["all_branch_scores_in_selection_order"])}
    for candidate in beam["all_branch_scores_in_selection_order"]:
        mask = int(candidate["branch_mask"])
        active_words = [*base_words, *(batch_words[index] for index in range(len(batch_words)) if (mask >> index) & 1)]
        prepared = core.prepare_subgroup(ambient, ambient_array, active_words)
        costs = [core.subgroup_score(ambient, ambient_array, prepared, target)["fresh_numerator"] for target in range(DIMENSION)]
        retained = [min(before, after) for before, after in zip(base_costs, costs)]
        unlock = sum(max(0.0, math.log1p(before) - math.log1p(after)) for before, after in zip(base_costs, costs))
        improved = sum(after < before for before, after in zip(base_costs, costs))
        item = {
            "branch_mask": mask,
            "increment_rank": candidate["increment_rank"],
            "selector_order": selector_lookup[mask],
            "in_selector_beam": any(mask == selected["branch_mask"] for selected in beam["beam"]),
            "stage_local_numerators_by_fixed_direction": {name: cost for name, cost in zip(target_names, costs)},
            "retained_from_v2_M24_numerators_by_fixed_direction": {name: cost for name, cost in zip(target_names, retained)},
            "frozen_16384_lower_tail_percentile_by_fixed_direction": {name: percentile(frozen_values[index], costs[index]) for index, name in enumerate(target_names)},
            "summary": summary(costs, retained, active_words),
            "retrospective_unlock_from_v2_M24_log1p_sum": decimal(unlock),
            "fixed_directions_improved_vs_v2_M24": improved,
        }
        rows.append(item)
    selected = [row for row in rows if row["in_selector_beam"]]
    inputs = {rel(path): sha(path) for path in (BEAM, LANDSCAPE, VISIBILITY, BASE_SELECTION, *BATCH_AUDITS, Path(__file__))}
    return {
        "schema": "elliptic-curves.curve302-v2-active-subgroup-beam-retrospective.v1",
        "status": "PASS_SEPARATE_RETROSPECTIVE_EVALUATION",
        "inputs": inputs,
        "boundary": {
            "selector": "The sealed beam selector did not read this evaluator, exceptional-direction inputs, the 2^14 state table, or any target cost.",
            "evaluator": "This post-sealing evaluator uses the known fourteen directions only to assess proposed subgroups. It performs no centre selection, chart scheduling, point search, or rank certification.",
            "metric": "All costs use the exact entrywise-rounded metric whose SHA-256 is fixed in the immutable 2^14 diagnostic; values at V2 subgroups are an off-grid extension of its documented nearest-plane plus exact coordinate-descent upper bound.",
            "not_claimed": "These finite local-CVP upper bounds are not exact closest-vector values, global visibility minima, a proof of a future search cost, or prospective knowledge.",
        },
        "fixed_diagnostic_reference": {
            "state_count": 1 << DIMENSION,
            "direction_ids": target_names,
            "rounded_metric_sha256": metric_sha(ambient),
            "height_pairing_asymmetry": asymmetry,
            "frozen_landscape_sha256": sha(LANDSCAPE),
        },
        "v2_M24_baseline": {
            "stage_local_numerators_by_fixed_direction": {name: cost for name, cost in zip(target_names, base_costs)},
            "summary": base_summary,
        },
        "diagnostic_coordinate_words": {
            "v2_M24_rows": base_words,
            "coalesced_batch_rows": batch_words,
            "meaning": "High-precision pairing recognition followed by exact elliptic-group-law verification; supplied for replay only, never to the selector.",
        },
        "branches_in_frozen_selector_order": rows,
        "selected_beam_retrospective_rows": selected,
        "reproducing_command": "sage -python elliptic-curves/cas/evaluate_curve302_v2_active_subgroup_beam.sage --check",
        "point_searches_run": 0,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.build == args.check:
        parser.error("choose exactly one of --build or --check")
    payload = build()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.build:
        if OUTPUT.exists():
            raise FileExistsError("preserve the immutable retrospective evaluation artifact")
        OUTPUT.write_text(rendered)
    elif not OUTPUT.is_file() or OUTPUT.read_text() != rendered:
        raise ArithmeticError("retrospective beam evaluation replay differs")
    print("V2_ACTIVE_SUBGROUP_BEAM_RETROSPECTIVE|states=7|point_searches=0|status=PASS", flush=True)


if __name__ == "__main__":
    main()
