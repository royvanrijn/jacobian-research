#!/usr/bin/env python3
"""Search neutral gate perturbations around the unique rank-35 slice.

The frozen ``qb+x2s`` cleanup has slice profile ``(rank JK, excess)=(17,1)``.
This bounded search adds each existing neutral low-degree circuit atom,
replays the frozen cleanup where possible, and continues at most ten residual
shared-factor steps.  Terminals are scored by the identity-slice Hessian
objective from ``search_identity_slice_hessian_rank.py``.

The output is finite-field diagnostic evidence only.  Any rank below 35 must
be frozen and certified over characteristic zero.
"""

from __future__ import annotations

from collections import Counter
import json
from pathlib import Path

from restricted_rank_profiles import PowerRankProfile, correction_profile
from search_identity_slice_hessian_rank import (
    identity_slice_profile,
    objective,
)
from search_rank_aware_bcw import state_key
from search_restricted_bcw_circuits import (
    CircuitState,
    base_circuit,
    diverse_profiled,
    diverse_states,
    encoded_plan,
    high_terms,
    partial_key,
    replay_encoded_plan,
    structural_key,
    transitions,
)


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "identity_slice_hessian_rank_search.json"
)
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "identity_slice_local_perturbation_search.json"
)
WIDTH = 64
MAX_STEPS = 10
PREBEAM_FACTOR = 3
PARTIAL_POWER_DEPTH = 8


def prefix_states() -> tuple[
    dict[str, object],
    list[CircuitState],
    list[tuple[CircuitState, dict[str, object]]],
    list[str],
]:
    stored = json.loads(SOURCE.read_text())
    frozen_plan = stored["best_terminal"]["plan"]
    _, _, atoms = base_circuit()
    names = [
        atom.name
        for atom in atoms
        if atom.name.startswith(("npert_", "aspert_", "qqpert_"))
    ]
    frontier: list[CircuitState] = []
    terminals: list[tuple[CircuitState, dict[str, object]]] = []
    failed: list[str] = []
    for name in names:
        plan = {
            "circuit_atoms": list(frozen_plan["circuit_atoms"]) + [name],
            "monomial_plan": frozen_plan["monomial_plan"],
        }
        try:
            state = replay_encoded_plan(plan)
        except (AssertionError, ValueError):
            failed.append(name)
            continue
        if high_terms(state.expressions, state.variables)[0][0] <= 3:
            terminals.append((state, identity_slice_profile(state)))
        else:
            frontier.append(state)
    return frozen_plan, frontier, terminals, failed


def main() -> None:
    frozen_plan, frontier, terminals, failed = prefix_states()
    partial_cache: dict[tuple[object, ...], PowerRankProfile] = {}
    depth_log = []

    for depth in range(1, MAX_STEPS + 1):
        deduplicated: dict[tuple[object, ...], CircuitState] = {}
        completed: list[CircuitState] = []
        generated = 0
        for state in frontier:
            for candidate in transitions(state):
                generated += 1
                if high_terms(candidate.expressions, candidate.variables)[0][0] <= 3:
                    completed.append(candidate)
                    continue
                key = state_key(candidate.legacy())
                previous = deduplicated.get(key)
                if previous is None or candidate.plan_key < previous.plan_key:
                    deduplicated[key] = candidate

        for candidate in completed:
            terminals.append((candidate, identity_slice_profile(candidate)))

        prebeam = diverse_states(
            deduplicated.values(),
            PREBEAM_FACTOR * WIDTH,
            structural_key,
        )
        profiled = []
        for candidate in prebeam:
            key = state_key(candidate.legacy())
            profile = partial_cache.get(key)
            if profile is None:
                profile = correction_profile(
                    candidate.expressions,
                    candidate.variables,
                    max_power=PARTIAL_POWER_DEPTH,
                )
                partial_cache[key] = profile
            profiled.append((candidate, profile))
        profiled = diverse_profiled(profiled, WIDTH, mode="mixed")
        frontier = [candidate for candidate, _ in profiled]
        lead = (
            min(profiled, key=lambda pair: partial_key(pair[0], pair[1]))
            if profiled
            else None
        )
        record = {
            "depth_after_frozen_cleanup": depth,
            "generated": generated,
            "unique": len(deduplicated),
            "kept": len(frontier),
            "terminal_count": len(terminals),
            "lead_high_signature": (
                list(high_terms(lead[0].expressions, lead[0].variables)[0])
                if lead
                else None
            ),
        }
        depth_log.append(record)
        print(record, flush=True)
        if not frontier:
            break

    histogram = Counter(objective(profile) for _, profile in terminals)
    best_objective = min(histogram)
    best = min(
        (
            (state, profile)
            for state, profile in terminals
            if objective(profile) == best_objective
        ),
        key=lambda pair: pair[0].plan_key,
    )
    family_leaders: dict[
        tuple[str, ...],
        tuple[CircuitState, dict[str, object]],
    ] = {}
    for state, profile in terminals:
        previous = family_leaders.get(state.circuit_atoms)
        if previous is None or (
            objective(profile),
            state.plan_key,
        ) < (
            objective(previous[1]),
            previous[0].plan_key,
        ):
            family_leaders[state.circuit_atoms] = (state, profile)

    payload = {
        "format": "identity-slice-local-perturbation-search-v1",
        "certification_status": (
            "bounded deterministic finite-field diagnostics; not a lower "
            "bound and not an exact characteristic-zero rank certificate"
        ),
        "frozen_rank_35_plan": frozen_plan,
        "search_scope": {
            "neutral_atom_count": 64,
            "failed_frozen_cleanup_replays": failed,
            "nonterminal_prefix_count": 44,
            "immediate_terminal_count": 19,
            "width": WIDTH,
            "max_steps_after_frozen_cleanup": MAX_STEPS,
            "prebeam_factor": PREBEAM_FACTOR,
            "partial_power_depth": PARTIAL_POWER_DEPTH,
            "terminal_count": len(terminals),
        },
        "objective_histogram": [
            {
                "slice_hessian_rank": key[0],
                "slice_rank_JK": key[1],
                "slice_kernel_excess": key[2],
                "slice_dimension": key[3],
                "count": count,
            }
            for key, count in sorted(histogram.items())
        ],
        "depth_log": depth_log,
        "best_terminal": {
            "objective": list(best_objective),
            "profile": best[1],
            "plan": encoded_plan(best[0]),
        },
        "family_leaders": [
            {
                "circuit_atoms": list(family),
                "objective": list(objective(profile)),
                "profile": profile,
                "plan": encoded_plan(state),
            }
            for family, (state, profile) in sorted(family_leaders.items())
        ],
        "conclusion": (
            "strict modular improvement found"
            if best_objective[0] < 35
            else "no searched neutral perturbation beats slice Hessian rank 35"
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"BEST objective={best_objective}")
    print(f"BEST atoms={best[0].circuit_atoms}")
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
