#!/usr/bin/env sage-python
"""Enumerate the finite 302 exceptional-direction subgroup landscape.

This is deliberately a *retrospective diagnostic*.  The fourteen known
directions provide a fixed complement to the certified generic rank-17
subgroup.  For every one of the 2^14 displayed subgroups it measures the
same half-lattice problem in one rounded canonical-height metric, then solves
finite order-optimization problems on those measured scores.

The score is a reproducible nearest-plane-plus-coordinate-descent CVP upper bound
for the half-lattice distance.  It is not a global pointed-quartic-coordinate
minimum and it is not a target-free centre selector.  The continuous
orthogonal-projection lower bound is recorded alongside it.  Thus this file
is suitable for asking how favourable an ordering of the *known* directions
would have been, but cannot be used to claim that such an ordering was
available prospectively.
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
from sage.all import ZZ, matrix


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
ART = ROOT / "artifacts/generated-results/elliptic-curves"
OUTPUT = ART / "curve302_exceptional_subgroup_landscape_v1.json"
VISIBILITY = ART / "curve302_residual_visibility_geometry_v1.json"
CHAIN = ART / "curve302_residual_strict_bootstrap_chain_v1.json"
AUTONOMOUS = ART / "adaptive_visibility_cascade_v1.json"
M24 = ART / "curve302_recovered_followup_wave_03_mod2_v1.json"
SCALE = 1_000_000
GENERIC_RANK = 17
DIMENSION = 14
FULL_RANK = GENERIC_RANK + DIMENSION
PRECISION = 192


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def read(path: Path) -> dict:
    return json.loads(path.read_text())


def load_module(name: str):
    return importlib.machinery.SourceFileLoader(
        "curve302_subgroup_landscape_" + name.replace(".", "_"), str(CAS / name)
    ).load_module()


def qpoint(row) -> tuple[F, F]:
    return tuple(F(value) for value in row)


def decimal(value: float) -> str:
    """Stable human-readable nonnegative float rendering for diagnostics."""

    if value < 0 and value > -1.0e-7:
        value = 0.0
    if value < 0 or not math.isfinite(value):
        raise ArithmeticError("invalid floating diagnostic value")
    return format(value, ".12g")


def bit_indices(mask: int) -> list[int]:
    return [index for index in range(DIMENSION) if (mask >> index) & 1]


def popcount(mask: int) -> int:
    return mask.bit_count()


def primary_target(entry: dict, base17, model, group) -> tuple[F, F]:
    """Recover the direction itself, undoing its chart-specific M17 translate."""

    roots = set()
    for row in entry["exact_pointed_quartic_rows"]:
        target = qpoint(row["exact_target_point_short_model"])
        inverse_translation = group.linear_combination(
            model, base17, [-int(value) for value in row["target_translation_m17_word"]]
        )
        roots.add(group.short_add(model, target, inverse_translation))
    if len(roots) != 1:
        raise ArithmeticError("vetted target charts do not encode one fixed direction")
    return roots.pop()


def rounded_metric(geometry, model, points):
    height, asymmetry = geometry.canonical_height_gram(model, points)
    rounded = matrix(ZZ, geometry.rounded_gram(height, SCALE))
    if rounded.nrows() != FULL_RANK or rounded.ncols() != FULL_RANK:
        raise ArithmeticError("fixed diagnostic basis has the wrong dimension")
    if not rounded.is_positive_definite():
        raise ArithmeticError("rounded 31-dimensional metric is not positive definite")
    return rounded, str(asymmetry)


def nearest_plane(upper, target):
    """Cholesky nearest-plane word for ``(z-target)^T G (z-target)``."""

    answer = [0] * len(target)
    for index in range(len(target) - 1, -1, -1):
        residual = sum(upper[index, column] * (answer[column] - target[column]) for column in range(index + 1, len(target)))
        answer[index] = int(np.rint(target[index] - residual / upper[index, index]))
    return answer


def exact_coordinate_descent(gram, cross, initial, h_ii):
    """Deterministically improve a nearest-plane word in coordinate neighbours.

    This is intentionally bounded local CVP work, not a claim of global CVP
    optimality.  Every returned word is nevertheless verified exactly in the
    rounded integral Gram form by the caller.
    """

    current = list(map(int, initial))
    dimension = len(current)
    gram_times_current = [sum(int(gram[row, column]) * current[column] for column in range(dimension)) for row in range(dimension)]
    current_value = 4 * h_ii - 4 * sum(current[index] * int(cross[index]) for index in range(dimension))
    current_value += sum(current[index] * gram_times_current[index] for index in range(dimension))
    changed = True
    while changed:
        changed = False
        for index in range(dimension):
            winner_delta = 0
            winner_value = current_value
            for delta in (-1, 1):
                difference = 2 * delta * gram_times_current[index] + delta * delta * int(gram[index, index])
                difference -= 4 * delta * int(cross[index])
                value = current_value + difference
                if value < winner_value or (value == winner_value and delta < winner_delta):
                    winner_delta, winner_value = delta, value
            if winner_delta:
                current[index] += winner_delta
                current_value = winner_value
                for row in range(dimension):
                    gram_times_current[row] += winner_delta * int(gram[row, index])
                changed = True
    if current_value < 0:
        raise ArithmeticError("local CVP descent produced a negative score")
    return current, current_value


def prepare_subgroup(ambient, ambient_array, subgroup_rows) -> dict:
    """Prepare one height lattice for all remaining targets at one state."""

    rows = matrix(ZZ, subgroup_rows)
    rows_array = np.array([[int(value) for value in row] for row in rows.rows()], dtype=np.int64)
    gram_array = rows_array @ ambient_array @ rows_array.T
    gram = matrix(ZZ, [[int(value) for value in row] for row in gram_array])
    if not gram.is_positive_definite() or not np.allclose(gram_array, gram_array.T):
        raise ArithmeticError("active subgroup Gram matrix is not positive definite")
    float_gram = gram_array.astype(np.float64)
    try:
        upper = np.linalg.cholesky(float_gram).T
    except np.linalg.LinAlgError as error:
        raise ArithmeticError("float Cholesky failed on an exact positive-definite Gram matrix") from error
    return {
        "rows": rows,
        "rows_array": rows_array,
        "gram": gram,
        "gram_array": gram_array,
        "upper": upper,
    }


def subgroup_score(ambient, ambient_array, prepared, target_index: int) -> dict:
    """Score dist(R_i, 1/2 M) in the common rounded height metric.

    If rows are the coordinate rows of M in the fixed 31-point basis and e_i
    is the target coordinate, this minimizes approximately

        h(R_i - q/2), q in M.

    ``fresh_numerator`` is an actual integral-lattice candidate and therefore
    an upper bound for exact CVP in the rounded metric.  ``continuous_lower``
    is the real-span lower bound.  They share the same units.
    """

    rows, gram = prepared["rows"], prepared["gram"]
    cross = rows * ambient.column(GENERIC_RANK + target_index)
    h_ii = int(ambient[GENERIC_RANK + target_index, GENERIC_RANK + target_index])
    cross_array = prepared["rows_array"] @ ambient_array[:, GENERIC_RANK + target_index]
    target = 2.0 * np.linalg.solve(prepared["gram_array"].astype(np.float64), cross_array.astype(np.float64))
    nearest = nearest_plane(prepared["upper"], target)
    active_word, numerator = exact_coordinate_descent(gram, cross, nearest, h_ii)
    # 4 h(R-q/2) = 4h(R) - 4<R,q> + h(q), an exact integer in the
    # entrywise-rounded metric.
    if numerator < 0:
        raise ArithmeticError("negative half-lattice CVP score")
    # A real linear-algebra lower bound.  This is only a floating diagnostic;
    # the integral word and numerator above are exact in the rounded metric.
    gram_float = prepared["gram_array"].astype(np.float64)
    cross_float = cross_array.astype(np.float64)
    lower = float(h_ii) - float(cross_float @ np.linalg.solve(gram_float, cross_float))
    if lower < -1.0e-3:
        raise ArithmeticError("continuous half-lattice lower bound became materially negative")
    lower = max(0.0, lower)
    if numerator / 4.0 + 1.0e-4 < lower:
        raise ArithmeticError("integral CVP candidate fell below real-span lower bound")
    return {
        "fresh_numerator": numerator,
        "continuous_lower_bound_scaled": decimal(lower),
        "active_word": active_word,
        "active_parity": [int(value) & 1 for value in active_word],
    }


def state_rows(mask: int):
    rows = []
    for index in range(GENERIC_RANK):
        rows.append([1 if column == index else 0 for column in range(FULL_RANK)])
    for index in bit_indices(mask):
        position = GENERIC_RANK + index
        rows.append([1 if column == position else 0 for column in range(FULL_RANK)])
    return rows


def fresh_scores(ambient):
    """Compute the local CVP estimate at every fixed-basis subgroup state."""

    values = [[None for _ in range(DIMENSION)] for _ in range(1 << DIMENSION)]
    details = {}
    ambient_array = np.array([[int(value) for value in row] for row in ambient.rows()], dtype=np.int64)
    for mask in range(1 << DIMENSION):
        prepared = prepare_subgroup(ambient, ambient_array, state_rows(mask))
        for target in range(DIMENSION):
            if (mask >> target) & 1:
                continue
            detail = subgroup_score(ambient, ambient_array, prepared, target)
            values[mask][target] = detail["fresh_numerator"]
            details[mask, target] = detail
        if mask % 512 == 0:
            print("LANDSCAPE CVP", mask, "of", (1 << DIMENSION) - 1, flush=True)
    return values, details


def retain_nested_candidates(fresh):
    """Retain every exact rounded-metric candidate from every sub-subgroup.

    This produces the finite-atlas analogue of C^*: a candidate belonging to
    T is also admissible at every S containing T.  Its score therefore cannot
    increase under subgroup inclusion, independent of LLL/Babai instability.
    """

    retained = [[None for _ in range(DIMENSION)] for _ in range(1 << DIMENSION)]
    source = [[None for _ in range(DIMENSION)] for _ in range(1 << DIMENSION)]
    for mask in range(1 << DIMENSION):
        for target in range(DIMENSION):
            if (mask >> target) & 1:
                continue
            candidates = [(fresh[mask][target], mask)]
            for removed in bit_indices(mask):
                parent = mask ^ (1 << removed)
                candidates.append((retained[parent][target], source[parent][target]))
            value, origin = min(candidates, key=lambda pair: (pair[0], pair[1]))
            retained[mask][target], source[mask][target] = value, origin
    return retained, source


def reconstruct_path(final_mask: int, previous: dict[int, tuple[int, int]]):
    order = []
    mask = final_mask
    while mask:
        parent, target = previous[mask]
        order.append(target)
        mask = parent
    return list(reversed(order))


def optimize_orders(retained):
    """Dynamic programming for the exact finite 2^14 diagnostic graph."""

    size = 1 << DIMENSION
    bottleneck = [None] * size
    bottleneck_previous = {}
    total = [None] * size
    total_previous = {}
    bottleneck[0] = 0
    total[0] = 0
    for mask in range(size):
        if bottleneck[mask] is None:
            raise ArithmeticError("dynamic programming state was unreachable")
        for target in range(DIMENSION):
            if (mask >> target) & 1:
                continue
            next_mask = mask | (1 << target)
            local = retained[mask][target]
            candidate_bottleneck = max(bottleneck[mask], local)
            candidate_total = total[mask] + local
            # deterministic tie break: smaller predecessor bitmask, then target
            key = (mask, target)
            current_key = bottleneck_previous.get(next_mask, (1 << 99, 1 << 99))
            if bottleneck[next_mask] is None or (candidate_bottleneck, key) < (bottleneck[next_mask], current_key):
                bottleneck[next_mask] = candidate_bottleneck
                bottleneck_previous[next_mask] = key
            current_key = total_previous.get(next_mask, (1 << 99, 1 << 99))
            if total[next_mask] is None or (candidate_total, key) < (total[next_mask], current_key):
                total[next_mask] = candidate_total
                total_previous[next_mask] = key
    final = size - 1
    return {
        "minimax": {
            "objective_numerator": bottleneck[final],
            "order_indices": reconstruct_path(final, bottleneck_previous),
        },
        "minimum_total": {
            "objective_numerator": total[final],
            "order_indices": reconstruct_path(final, total_previous),
        },
    }


def path_rows(order, names, fresh, retained, source, details):
    mask = 0
    rows = []
    for target in order:
        if mask >> target & 1:
            raise ArithmeticError("order repeats a target")
        origin = source[mask][target]
        witness = details[origin, target]
        rows.append({
            "rank_before": GENERIC_RANK + popcount(mask),
            "state_mask": mask,
            "add_direction": names[target],
            "stage_local_numerator": fresh[mask][target],
            "retained_numerator": retained[mask][target],
            "retained_witness_source_mask": origin,
            "retained_witness_active_word": witness["active_word"],
            "retained_witness_active_parity": witness["active_parity"],
            "fresh_continuous_lower_bound_scaled": details[mask, target]["continuous_lower_bound_scaled"],
        })
        mask |= 1 << target
    return rows


def unlock_rows(retained, names):
    """Full retrospective unlock values, but only top four directions per state."""

    rows = []
    for mask in range(1 << DIMENSION):
        candidates = []
        remaining = [index for index in range(DIMENSION) if not ((mask >> index) & 1)]
        for added in remaining:
            next_mask = mask | (1 << added)
            unlock = 0.0
            exposed = 0
            for target in remaining:
                if target == added:
                    continue
                before = retained[mask][target]
                after = retained[next_mask][target]
                # A zero means the diagnostic target itself is a half-lattice
                # point of the active subgroup.  log(0) is not defined, so the
                # finite unlock statistic uses the documented log(1+C)
                # regularization rather than discarding that endpoint.
                improvement = max(0.0, math.log1p(before) - math.log1p(after))
                unlock += improvement
                if after < before:
                    exposed += 1
            candidates.append({
                "add_direction": names[added],
                "unlock_log_ratio_sum": decimal(unlock),
                "other_directions_improved": exposed,
            })
        candidates.sort(key=lambda row: (-float(row["unlock_log_ratio_sum"]), -row["other_directions_improved"], row["add_direction"]))
        rows.append({"state_mask": mask, "rank": GENERIC_RANK + popcount(mask), "top_unlock_choices": candidates[:4]})
    return rows


def autonomous_coordinates(geometry, group, ambient_points, ambient, model, base17):
    """Express autonomous-V1 points in the displayed 31-point diagnostic basis.

    Pairing coordinates propose an integral relation; the elliptic group law
    then verifies it exactly.  This gives a comparison of *subgroups*, not a
    fictional order of the fourteen diagnostic coordinate axes.
    """

    report = read(AUTONOMOUS)
    if report["status"] != "PASS_BOUNDED_AUTONOMOUS_CASCADE_AND_FINITE_ATLAS_REPLAY":
        raise ArithmeticError("autonomous V1 report is not a passed bounded replay")
    calibration = next(row for row in report["rows"] if row["id"] == "calibration302")
    # High-precision real pairings are used only to recognize integral words;
    # every such word is then checked in the exact elliptic group law.
    from sage.all import RealField
    field = RealField(192)
    high = matrix(field, [[field(str(value).strip("()")) for value in row] for row in geometry.canonical_height_gram(model, ambient_points)[0]])
    rows = []
    ambient_array = np.array([[int(value) for value in row] for row in ambient.rows()], dtype=np.int64)
    prior_count = GENERIC_RANK
    for stage in calibration["stages"]:
        cloud = read(ROOT / stage["mod2"])
        points = tuple(qpoint(point) for point in cloud["independent_points"])
        if len(points) != stage["after"] or points[:GENERIC_RANK] != base17:
            raise ArithmeticError("autonomous basis lost the common generic prefix")
        new_words = []
        for point in points[prior_count:]:
            pairing = geometry.canonical_height_gram(model, (*ambient_points, point))[0]
            column = matrix(field, FULL_RANK, 1, [field(str(pairing[index][FULL_RANK]).strip("()")) for index in range(FULL_RANK)])
            coefficients = high.solve_right(column)
            word = [int(round(float(coefficients[index, 0]))) for index in range(FULL_RANK)]
            error = max(abs(coefficients[index, 0] - word[index]) for index in range(FULL_RANK))
            if error > field("1e-30") or group.linear_combination(model, ambient_points, word) != point:
                raise ArithmeticError("autonomous point did not have an exactly verified diagnostic-basis word")
            new_words.append(word)
        prior_count = len(points)
        active = state_rows(0) + [word for earlier in rows for word in earlier["new_words"]] + new_words
        prepared = prepare_subgroup(ambient, ambient_array, active)
        costs = {}
        for target in range(DIMENSION):
            detail = subgroup_score(ambient, ambient_array, prepared, target)
            costs[target] = detail
        rows.append({
            "wave": stage["wave"],
            "rank_before": stage["before"],
            "rank_after": stage["after"],
            "new_words": new_words,
            "direction_costs": {
                str(target): {
                    "fresh_numerator": costs[target]["fresh_numerator"],
                    "continuous_lower_bound_scaled": costs[target]["continuous_lower_bound_scaled"],
                }
                for target in range(DIMENSION)
            },
        })
    return rows


def build() -> dict:
    visibility, chain, m24 = read(VISIBILITY), read(CHAIN), read(M24)
    if visibility["status"] != "PASS_RETROSPECTIVE_VETTED_VISIBILITY_DIAGNOSTIC":
        raise ArithmeticError("missing completed visibility diagnostic")
    if not visibility["target_basis_boundary"]["M17_to_M24_to_D_integral_direct_sum"]:
        raise ArithmeticError("fixed 14-direction diagnostic basis boundary changed")
    if chain["status"] != "PASS_EXACT_CHAIN_PRESERVED_RANK_AT_LEAST_31":
        raise ArithmeticError("missing completed exact rank-31 chain")
    if m24["rank_lower_bound"] != 24 or len(m24["independent_points"]) != 24:
        raise ArithmeticError("M24 prefix cloud unavailable")
    inputs = {rel(path): sha(path) for path in (VISIBILITY, CHAIN, AUTONOMOUS, M24, Path(__file__))}
    group = load_module("half_lattice_pointed_sieve.py")
    geometry = load_module("prospective_half_lattice_v3.sage")
    model = tuple(F(value) for value in m24["curve"])
    base17 = tuple(qpoint(point) for point in m24["independent_points"][:GENERIC_RANK])
    directions = visibility["directions"]
    names = [row["id"] for row in directions]
    if len(names) != DIMENSION or len(set(names)) != DIMENSION:
        raise ArithmeticError("expected fourteen distinct fixed diagnostic directions")
    targets = tuple(primary_target(row, base17, model, group) for row in directions)
    ambient_points = (*base17, *targets)
    ambient, asymmetry = rounded_metric(geometry, model, ambient_points)
    fresh, fresh_details = fresh_scores(ambient)
    retained, source = retain_nested_candidates(fresh)
    optimization = optimize_orders(retained)
    for entry in optimization.values():
        entry["order"] = [names[index] for index in entry.pop("order_indices")]
        entry["objective_scaled_height"] = decimal(entry["objective_numerator"] / (4.0 * SCALE))
    historical_tail = [arm["direction"] for arm in chain["recovery_arms"]]
    name_index = {name: index for index, name in enumerate(names)}
    recovered = [index for index, row in enumerate(directions) if row["cohort"] in {"recovered_local", "recovered_strict"}]
    residual = [name_index[name] for name in historical_tail]
    if len(recovered) != 7 or len(residual) != 7 or set(recovered) & set(residual):
        raise ArithmeticError("unexpected 7+7 diagnostic cohort split")
    # Only the M24->M31 tail has an attested single-direction historical order.
    # The first seven are an independently certified M24 block, so we do not
    # manufacture an order between them.
    tail_mask = sum(1 << index for index in recovered)
    historical_tail_rows = path_rows(residual, names, fresh, retained, source, fresh_details)
    for row in historical_tail_rows:
        row["state_mask"] |= tail_mask
        row["rank_before"] += 7
    # Recompute row values at the actual M24+tail states.  path_rows starts at
    # M17, while the retained source may still be any contained substate.
    actual_tail_rows = []
    mask = tail_mask
    for target in residual:
        origin = source[mask][target]
        actual_tail_rows.append({
            "rank_before": GENERIC_RANK + popcount(mask),
            "state_mask": mask,
            "add_direction": names[target],
            "stage_local_numerator": fresh[mask][target],
            "retained_numerator": retained[mask][target],
            "retained_witness_source_mask": origin,
            "retained_witness_active_word": fresh_details[origin, target]["active_word"],
            "retained_witness_active_parity": fresh_details[origin, target]["active_parity"],
            "fresh_continuous_lower_bound_scaled": fresh_details[mask, target]["continuous_lower_bound_scaled"],
        })
        mask |= 1 << target
    optimal_rows = {
        label: path_rows([name_index[name] for name in entry["order"]], names, fresh, retained, source, fresh_details)
        for label, entry in optimization.items()
    }
    autonomous = autonomous_coordinates(geometry, group, ambient_points, ambient, model, base17)
    subset_rows = []
    for mask in range(1 << DIMENSION):
        subset_rows.append({
            "state_mask": mask,
            "rank": GENERIC_RANK + popcount(mask),
            "stage_local_numerators": [fresh[mask][index] for index in range(DIMENSION)],
            "retained_numerators": [retained[mask][index] for index in range(DIMENSION)],
            "retained_source_masks": [source[mask][index] for index in range(DIMENSION)],
        })
    return {
        "schema": "elliptic-curves.curve302-exceptional-subgroup-landscape.v1",
        "status": "PASS_RETROSPECTIVE_FINITE_SUBGROUP_LANDSCAPE",
        "inputs": inputs,
        "fixed_basis": {
            "generic_rank": GENERIC_RANK,
            "direction_count": DIMENSION,
            "direction_ids": names,
            "cohorts": {row["id"]: row["cohort"] for row in directions},
            "basis_boundary": visibility["target_basis_boundary"],
            "rounded_metric_scale": SCALE,
            "rounded_metric_sha256": hashlib.sha256(json.dumps([list(map(int, row)) for row in ambient.rows()], separators=(",", ":")).encode()).hexdigest(),
            "height_pairing_asymmetry": asymmetry,
        },
        "cost_definition": {
            "formula": "C_i(M)=4*min_{q in M} hhat_rounded(R_i-q/2), approximated by Cholesky nearest-plane plus deterministic exact coordinate descent",
            "units": "numerator in four times the entrywise-1e6-rounded canonical-height metric",
            "stage_local": "new nearest-plane/local-CVP candidate at exactly M",
            "retained": "minimum of all stage-local candidates from subgroups T subseteq M; it is a verified finite-atlas upper bound and is monotone under inclusion",
            "continuous_lower_bound": "squared distance from R_i to span_R(M), in the same rounded-height units before multiplying by four",
            "unlock_regularization": "U uses log(1+C) so exact half-lattice endpoints with C=0 remain finite",
            "not_claimed": [
                "exact closest-vector optimality in the rounded metric",
                "a global canonical-height minimum",
                "a pointed-quartic coordinate-height minimum",
                "a target-free selection rule or a rank result",
            ],
        },
        "subset_states": subset_rows,
        "optimization": optimization,
        "historical_comparison": {
            "M17_to_M24": "The source certifies a recovered seven-direction block, not a total order on its fixed diagnostic basis; this analysis deliberately does not manufacture one.",
            "attested_M24_to_M31_order": historical_tail,
            "attested_tail_rows": actual_tail_rows,
            "autonomous_V1": {
                "interpretation": "Each autonomous point is expressed and exactly verified in the fixed 31-point basis. These are arbitrary subgroup directions, not a permutation of the fourteen diagnostic coordinate axes, so they are compared as subgroup states rather than falsely as a 14-direction order.",
                "stages": autonomous,
            },
        },
        "optimal_path_rows": optimal_rows,
        "unlock_value_top_four_per_state": unlock_rows(retained, names),
        "reproducing_command": "sage -python elliptic-curves/cas/audit_curve302_exceptional_subgroup_landscape.sage --check",
        "boundary": "All fourteen directions and the autonomous point words are known data used only retrospectively. The finite 2^14 enumeration studies the displayed direct-sum diagnostic lattice; it neither selects centres prospectively nor proves a rank, a search-cost law, or the optimality of any historical experimental policy.",
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
            raise FileExistsError("preserve the immutable subgroup-landscape artifact")
        OUTPUT.write_text(rendered)
    elif not OUTPUT.is_file() or OUTPUT.read_text() != rendered:
        raise ArithmeticError("subgroup-landscape replay differs")
    print("CURVE302_SUBGROUP_LANDSCAPE|states=16384|directions=14|status=PASS", flush=True)


if __name__ == "__main__":
    main()
