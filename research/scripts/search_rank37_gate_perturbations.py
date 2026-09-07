#!/usr/bin/env python3
"""Continue low-degree gate perturbations from the rank-37 circuit prefix.

The certified ``qb+x2s`` discovery exposes ``a=xy`` and ``q=xy^2``.  A
target shear linear in ``q`` is killed exactly by affine-jet normalization.
The first nontrivial no-new-gate deformation is therefore ``lambda*a*q``.
This script replays the twenty-step rank-37 prefix for several rational
``lambda`` values and searches only the residual degree-five cleanup.

All ranks and power ranks in the resulting JSON are modular diagnostics.
Any strict improvement must still be frozen and certified over QQ.
"""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path

from restricted_rank_profiles import PowerRankProfile, correction_profile
from search_rank_aware_bcw import state_key
from search_restricted_bcw_circuits import (
    INCUMBENTS,
    PERTURBATION_COEFFICIENTS,
    CircuitState,
    beats_incumbent,
    diverse_profiled,
    diverse_states,
    encoded_plan,
    high_terms,
    pareto_terminals,
    partial_key,
    replay_encoded_plan,
    structural_key,
    terminal_key,
    terminal_objective,
    terminal_profile,
    transitions,
)


ROOT = Path(__file__).resolve().parents[1]
DISCOVERY = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "restricted_bcw_circuit_search_v2_w64.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "rank37_gate_perturbation_search.json"
)


def perturbation_name(kind: str, coefficient: int) -> str:
    sign = "m" if coefficient < 0 else "p"
    return f"{kind}_{sign}{abs(coefficient)}"


def rank37_plan() -> dict[str, object]:
    stored = json.loads(DISCOVERY.read_text())
    matches = [
        terminal["plan"]
        for terminal in stored["pareto_terminals"]
        if terminal["objective"][2] == 37
    ]
    if len(matches) != 1:
        raise AssertionError("rank-37 discovery plan is not unique")
    return matches[0]


def prefix_states() -> tuple[dict[str, object], list[CircuitState]]:
    plan = rank37_plan()
    baseline = replay_encoded_plan(plan)
    baseline_key = state_key(baseline.legacy())

    # Linear target shears in the exposed q-output are exactly removed by
    # the affine normalization.  Keep this as a regression, not a heuristic.
    for coefficient in PERTURBATION_COEFFICIENTS:
        name = perturbation_name("qpert", coefficient)
        perturbed_plan = {
            "circuit_atoms": list(plan["circuit_atoms"]) + [name],
            "monomial_plan": plan["monomial_plan"],
        }
        state = replay_encoded_plan(perturbed_plan)
        assert state_key(state.legacy()) == baseline_key

    states = []
    for coefficient in PERTURBATION_COEFFICIENTS:
        name = perturbation_name("aqpert", coefficient)
        perturbed_plan = {
            "circuit_atoms": list(plan["circuit_atoms"]) + [name],
            "monomial_plan": plan["monomial_plan"],
        }
        state = replay_encoded_plan(perturbed_plan)
        assert high_terms(state.expressions, state.variables)[0] == (5, 8, 2)
        states.append(state)
    return plan, states


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--width", type=int, default=16)
    parser.add_argument("--max-steps", type=int, default=10)
    parser.add_argument("--prebeam-factor", type=int, default=3)
    parser.add_argument("--partial-power-depth", type=int, default=8)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--skip-terminal-hessian-power", action="store_true")
    args = parser.parse_args()
    if not args.output.is_absolute():
        args.output = ROOT / args.output
    if min(
        args.width,
        args.max_steps,
        args.prebeam_factor,
        args.partial_power_depth,
    ) < 1:
        parser.error("all search limits must be positive")

    frozen_plan, frontier = prefix_states()
    partial_cache: dict[tuple[object, ...], PowerRankProfile] = {}
    terminals: list[tuple[CircuitState, dict[str, object]]] = []
    depth_log = []
    histogram: Counter[tuple[int, int | None, int, int]] = Counter()

    for depth in range(1, args.max_steps + 1):
        deduplicated: dict[tuple[object, ...], CircuitState] = {}
        completed: list[CircuitState] = []
        generated = 0
        for state in frontier:
            if high_terms(state.expressions, state.variables)[0][0] <= 3:
                completed.append(state)
                continue
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
            profile = terminal_profile(
                candidate,
                hessian_power=not args.skip_terminal_hessian_power,
            )
            terminals.append((candidate, profile))
            histogram[
                (
                    int(profile["cubic_rank_mod_1000003"]),
                    (
                        int(profile["cubic_sampled_index"])
                        if profile["cubic_sampled_index"] is not None
                        else None
                    ),
                    int(profile["cotangent_hessian_rank_mod_1000003"]),
                    int(profile["quotient_dimension"]),
                )
            ] += 1

        prebeam = diverse_states(
            deduplicated.values(),
            args.prebeam_factor * args.width,
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
                    max_power=args.partial_power_depth,
                )
                partial_cache[key] = profile
            profiled.append((candidate, profile))
        profiled = diverse_profiled(profiled, args.width)
        frontier = [candidate for candidate, _ in profiled]
        lead = (
            min(profiled, key=lambda pair: partial_key(pair[0], pair[1]))
            if profiled
            else None
        )
        record = {
            "depth_after_frozen_prefix": depth,
            "generated": generated,
            "unique": len(deduplicated),
            "kept": len(frontier),
            "terminal_count": len(terminals),
            "lead_high_signature": (
                list(high_terms(lead[0].expressions, lead[0].variables)[0])
                if lead
                else None
            ),
            "lead_power_ranks_mod_1000003": (
                list(lead[1].ranks) if lead else None
            ),
        }
        depth_log.append(record)
        print(record, flush=True)
        if not frontier:
            break

    archive = pareto_terminals(terminals)
    family_leaders: dict[
        tuple[str, ...], tuple[CircuitState, dict[str, object]]
    ] = {}
    for state, profile in terminals:
        previous = family_leaders.get(state.circuit_atoms)
        if previous is None or terminal_key(profile, state) < terminal_key(
            previous[1], previous[0]
        ):
            family_leaders[state.circuit_atoms] = (state, profile)

    payload = {
        "format": "rank37-gate-perturbation-search-v1",
        "certification_status": (
            "bounded modular search; strict improvements require exact QQ "
            "generation and independent replay"
        ),
        "frozen_prefix": frozen_plan,
        "coefficients": list(PERTURBATION_COEFFICIENTS),
        "exact_triviality_regression": (
            "all lambda*q target shears normalize to the unperturbed state"
        ),
        "nontrivial_family": "lambda*a*q target shear",
        "incumbents": INCUMBENTS,
        "search_scope": {
            "width": args.width,
            "max_steps_after_prefix": args.max_steps,
            "prebeam_factor": args.prebeam_factor,
            "partial_power_depth": args.partial_power_depth,
            "terminal_hessian_power_profiled": (
                not args.skip_terminal_hessian_power
            ),
        },
        "depth_log": depth_log,
        "terminal_histogram": [
            {
                "cubic_rank": key[0],
                "cubic_sampled_index": key[1],
                "quartic_hessian_rank": key[2],
                "quotient_dimension": key[3],
                "count": count,
            }
            for key, count in sorted(
                histogram.items(),
                key=lambda item: (
                    item[0][0],
                    item[0][1] if item[0][1] is not None else 10**9,
                    item[0][2:],
                ),
            )
        ],
        "pareto_terminals": [
            {
                "objective": list(terminal_objective(profile)),
                "profile": profile,
                "plan": encoded_plan(state),
                "beats_a_current_incumbent_diagnostically": beats_incumbent(
                    profile
                ),
            }
            for state, profile in archive
        ],
        "coefficient_family_leaders": [
            {
                "circuit_atoms": list(family),
                "objective": list(terminal_objective(profile)),
                "profile": profile,
                "plan": encoded_plan(state),
            }
            for family, (state, profile) in sorted(family_leaders.items())
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"wrote {args.output.relative_to(ROOT)}")
    if not archive:
        print("NO TERMINAL")
    for state, profile in archive:
        print(
            "PARETO "
            f"objective={terminal_objective(profile)} "
            f"atoms={state.circuit_atoms} "
            f"steps={len(state.monomial_plan)}"
        )


if __name__ == "__main__":
    main()
