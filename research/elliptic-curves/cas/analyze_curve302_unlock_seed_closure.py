#!/usr/bin/env python3
"""Retrospective single-seed closure census for curve 302.

Given the immutable 2^14 exceptional-subgroup landscape, compute for the empty
seed and for each single exceptional direction the exact minimax bottleneck of
*that finite retained-atlas graph* needed to reach the full displayed M31
subgroup.  This is a retrospective oracle diagnostic, not a prospective point
selector and not a claim about exact CVP or search runtime.

Optionally pass --xi-direction ID once the arithmetic unlock class has been
identified with one of the fourteen displayed directions; the report then
highlights that row without changing any computation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / "artifacts/generated-results/elliptic-curves"
INPUT = ART / "curve302_exceptional_subgroup_landscape_v1.json"
OUTPUT = ART / "curve302_unlock_seed_closure_v1.json"
DIM = 14
GENERIC_RANK = 17


def read(path: Path):
    return json.loads(path.read_text())


def sha(path: Path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def reconstruct(final_mask: int, start_mask: int, prev: dict[int, tuple[int, int]]):
    order = []
    mask = final_mask
    while mask != start_mask:
        parent, target = prev[mask]
        order.append(target)
        mask = parent
    return list(reversed(order))


def minimax_from(states, start_mask: int):
    """Exact DP on the monotone subset DAG using retained finite-atlas costs."""
    full = (1 << DIM) - 1
    inf = 1 << 200
    best = [inf] * (1 << DIM)
    prev: dict[int, tuple[int, int]] = {}
    best[start_mask] = 0
    # Every edge increases popcount, so integer mask order is a valid topological
    # order for supersets of start_mask: parent = child ^ bit < child.
    for mask in range(start_mask, full + 1):
        if (mask & start_mask) != start_mask or best[mask] == inf:
            continue
        costs = states[mask]["retained_numerators"]
        for target in range(DIM):
            if mask >> target & 1:
                continue
            value = costs[target]
            if value is None:
                raise ArithmeticError(f"missing retained edge cost at {mask=} {target=}")
            child = mask | (1 << target)
            candidate = max(best[mask], int(value))
            old_key = prev.get(child, (1 << 99, 1 << 99))
            new_key = (mask, target)
            if candidate < best[child] or (candidate == best[child] and new_key < old_key):
                best[child] = candidate
                prev[child] = new_key
    if best[full] == inf:
        raise ArithmeticError("full subgroup is unreachable in monotone subset graph")
    return best[full], reconstruct(full, start_mask, prev)


def closure_at(states, start_mask: int, threshold: int):
    """Least fixed point obtained by adding every currently visible direction."""
    mask = start_mask
    waves = []
    while True:
        costs = states[mask]["retained_numerators"]
        add = [i for i in range(DIM) if not (mask >> i & 1) and costs[i] is not None and int(costs[i]) <= threshold]
        if not add:
            return mask, waves
        waves.append(add)
        for target in add:
            mask |= 1 << target


def row_for_seed(states, names, seed_index: int | None, thresholds):
    start = 0 if seed_index is None else 1 << seed_index
    bottleneck, order = minimax_from(states, start)
    named_order = [names[i] for i in order]
    result = {
        "seed_direction": None if seed_index is None else names[seed_index],
        "start_rank": GENERIC_RANK + (0 if seed_index is None else 1),
        "minimax_bottleneck_numerator": bottleneck,
        "minimax_bottleneck_scaled_height": bottleneck / 4_000_000.0,
        # Canonical field name plus a compatibility alias. Some handoff/report
        # tooling used the older completion terminology; both must remain exact.
        "optimal_followup_order": named_order,
        "optimal_completion_order": named_order,
    }
    closures = {}
    for label, value in thresholds.items():
        mask, waves = closure_at(states, start, value)
        closures[label] = {
            "threshold_numerator": value,
            "final_rank": GENERIC_RANK + mask.bit_count(),
            "full_M31_closure": mask == (1 << DIM) - 1,
            "waves": [[names[i] for i in wave] for wave in waves],
        }
    result["calibrated_closures"] = closures
    return result


def build(xi_direction: str | None):
    source = read(INPUT)
    if source.get("status") != "PASS_RETROSPECTIVE_FINITE_SUBGROUP_LANDSCAPE":
        raise ArithmeticError("immutable subgroup landscape is not a passed artifact")
    fixed = source["fixed_basis"]
    names = fixed["direction_ids"]
    if len(names) != DIM or len(set(names)) != DIM:
        raise ArithmeticError("expected fourteen distinct direction IDs")
    states = source["subset_states"]
    if len(states) != 1 << DIM or any(row["state_mask"] != i for i, row in enumerate(states)):
        raise ArithmeticError("subset-state ledger is incomplete or reordered")

    # Two empirical diagnostic thresholds already present in the immutable
    # source: the root global minimax and the attested M24->M31 tail bottleneck.
    global_minimax = int(source["optimization"]["minimax"]["objective_numerator"])
    actual = source["historical_comparison"].get("attested_tail_rows", [])
    if not actual:
        raise ArithmeticError("attested tail rows missing from landscape artifact")
    historical_tail = max(int(row["retained_numerator"]) for row in actual)
    thresholds = {"global_minimax": global_minimax, "historical_tail_bottleneck": historical_tail}

    rows = [row_for_seed(states, names, None, thresholds)]
    rows += [row_for_seed(states, names, i, thresholds) for i in range(DIM)]
    singles = rows[1:]
    singles.sort(key=lambda r: (r["minimax_bottleneck_numerator"], r["seed_direction"]))

    if xi_direction is not None and xi_direction not in names:
        raise ValueError("--xi-direction must be one of: " + ", ".join(names))
    highlighted = next((r for r in singles if r["seed_direction"] == xi_direction), None)

    return {
        "schema": "elliptic-curves.curve302-unlock-seed-closure.v1",
        "status": "PASS_RETROSPECTIVE_SINGLE_SEED_CLOSURE_CENSUS",
        "input": {str(INPUT.relative_to(ROOT)): sha(INPUT)},
        "direction_ids": names,
        "thresholds": thresholds,
        "unseeded": rows[0],
        "single_seed_ranking": singles,
        "xi_direction": xi_direction,
        "xi_result": highlighted,
        "interpretation_boundary": (
            "Exact dynamic programming on the already-published retained finite-atlas costs only. "
            "A low bottleneck says that a known seed makes the displayed half-lattice diagnostic "
            "favourable; it does not prove frozen V3 would find the next point, does not establish "
            "a point-search runtime bound, and uses the known M31 complement retrospectively."
        ),
        "reproducing_command": "python3 elliptic-curves/cas/analyze_curve302_unlock_seed_closure.py --check",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xi-direction")
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.write == args.check:
        parser.error("choose exactly one of --write or --check")
    result = build(args.xi_direction)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.write:
        if OUTPUT.exists():
            raise FileExistsError("preserve immutable unlock-seed census")
        OUTPUT.write_text(rendered)
    else:
        if not OUTPUT.exists() or OUTPUT.read_text() != rendered:
            raise ArithmeticError("unlock-seed census replay differs")
    best = result["single_seed_ranking"][0]
    print("CURVE302_UNLOCK_CENSUS|best_seed={}|bottleneck={}|scaled={:.9f}|status=PASS".format(
        best["seed_direction"], best["minimax_bottleneck_numerator"], best["minimax_bottleneck_scaled_height"]), flush=True)
    print("CURVE302_UNLOCK_ORDER|" + " -> ".join(best["optimal_followup_order"]), flush=True)
    if result["xi_result"]:
        row = result["xi_result"]
        print("CURVE302_XI_SEED|direction={}|bottleneck={}|scaled={:.9f}".format(
            row["seed_direction"], row["minimax_bottleneck_numerator"], row["minimax_bottleneck_scaled_height"]), flush=True)


if __name__ == "__main__":
    main()
