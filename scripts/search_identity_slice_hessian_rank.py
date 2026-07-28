#!/usr/bin/env python3
"""Rescore the restricted BCW circuit census after identity slicing.

The earlier circuit search optimized homogeneous cubic and cotangent ranks.
Every retained homogeneous source has a final identity output whose collision
points lie in the slice ``s=1``.  This script repeats the frozen width-64
terminal census and scores the nonhomogeneous slice instead:

    rank(JK), rank(sum y_i Hess(K_i)), rank(Hess(y.K)), kernel excess.

All ranks are deterministic finite-field diagnostics.  A strict improvement
over Hessian rank 35 must be frozen and certified over characteristic zero
before it changes a theorem-level statement.
"""

from __future__ import annotations

from collections import Counter
import json
from pathlib import Path

import sympy as sp

from rank_compressed_bcw_homogenization import (
    extract_quadratic_cubic,
    factor_cubic_output,
    iterated_constant_kernel_quotient,
    rank_compressed_homogeneous_map,
)
from restricted_rank_profiles import (
    PowerRankProfile,
    correction_profile,
    cotangent_hessian_rank,
    polynomial_jacobian_profile,
)
from search_rank_aware_bcw import state_key
from search_restricted_bcw_circuits import (
    CORE_ATOM_NAMES,
    CircuitState,
    diverse_profiled,
    diverse_states,
    encoded_plan,
    high_terms,
    partial_key,
    seed_family,
    structural_key,
    transitions,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "identity_slice_hessian_rank_search.json"
)
WIDTH = 64
MAX_STEPS = 24
PREBEAM_FACTOR = 4
PARTIAL_POWER_DEPTH = 8
MAX_CIRCUIT_GATES = 7


def identity_slice_profile(state: CircuitState) -> dict[str, object]:
    quadratic, cubic = extract_quadratic_cubic(
        state.expressions,
        state.variables,
    )
    factorization = factor_cubic_output(cubic)
    ambient_variables, ambient_h = rank_compressed_homogeneous_map(
        state.variables,
        quadratic,
        factorization,
    )
    quotient = iterated_constant_kernel_quotient(
        ambient_variables,
        ambient_h,
    )
    variables = quotient.quotient_variables
    homogeneous_h = quotient.quotient_h
    assert homogeneous_h[-1].is_zero

    projected = []
    for point in state.collision_points:
        substitution = dict(zip(state.variables, point))
        cubic_values = [poly.eval(substitution) for poly in factorization.c]
        ambient_point = sp.Matrix(list(point) + cubic_values + [sp.Integer(1)])
        projected.append(ambient_point)
    for stage in quotient.stages:
        projected = [stage.B * point for point in projected]
    assert len({tuple(point) for point in projected}) == 3
    assert {point[-1] for point in projected} == {sp.Integer(1)}

    slice_variables = variables[:-1]
    identity_variable = variables[-1]
    sliced_h = [
        sp.Poly(
            component.as_expr().subs(identity_variable, 1),
            *slice_variables,
            domain=sp.QQ,
        )
        for component in homogeneous_h[:-1]
    ]
    sliced_points = [
        tuple(point[index] for index in range(len(slice_variables)))
        for point in projected
    ]
    images = []
    for point in sliced_points:
        substitution = dict(zip(slice_variables, point))
        images.append(
            tuple(
                coordinate + component.eval(substitution)
                for coordinate, component in zip(point, sliced_h)
            )
        )
    assert images[0] == images[1] == images[2]

    jacobian_profile = polynomial_jacobian_profile(
        sliced_h,
        slice_variables,
    )
    block_ranks = cotangent_hessian_rank(sliced_h)
    excess = block_ranks[2] - 2 * jacobian_profile.rank
    assert excess >= 0
    return {
        "circuit_atoms": list(state.circuit_atoms),
        "monomial_cleanup_steps": len(state.monomial_plan),
        "homogeneous_source_dimension": len(variables),
        "slice_dimension": len(slice_variables),
        "slice_correction_degrees": sorted(
            {
                sum(monomial)
                for component in sliced_h
                for monomial, coefficient in component.terms()
                if coefficient
            }
        ),
        "slice_JK_power_ranks_mod_1000003": list(
            jacobian_profile.ranks
        ),
        "slice_rank_JK_mod_1000003": jacobian_profile.rank,
        "slice_sampled_index": jacobian_profile.sampled_index,
        "slice_rank_A_mod_1000003": block_ranks[1],
        "slice_cotangent_hessian_rank_mod_1000003": block_ranks[2],
        "slice_cotangent_kernel_excess_mod_1000003": excess,
        "slice_collision_separated": True,
    }


def objective(profile: dict[str, object]) -> tuple[int, ...]:
    return (
        int(profile["slice_cotangent_hessian_rank_mod_1000003"]),
        int(profile["slice_rank_JK_mod_1000003"]),
        int(profile["slice_cotangent_kernel_excess_mod_1000003"]),
        int(profile["slice_dimension"]),
    )


def main() -> None:
    frontier = seed_family(
        MAX_CIRCUIT_GATES,
        frozenset(),
        CORE_ATOM_NAMES,
    )
    partial_cache: dict[tuple[object, ...], PowerRankProfile] = {}
    terminals: list[tuple[CircuitState, dict[str, object]]] = []
    depth_log: list[dict[str, object]] = []

    for depth in range(1, MAX_STEPS + 1):
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
            "depth": depth,
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
        "format": "identity-slice-hessian-rank-search-v1",
        "certification_status": (
            "bounded deterministic finite-field diagnostics; not a lower "
            "bound and not an exact characteristic-zero rank certificate"
        ),
        "search_scope": {
            "width": WIDTH,
            "max_steps": MAX_STEPS,
            "prebeam_factor": PREBEAM_FACTOR,
            "partial_power_depth": PARTIAL_POWER_DEPTH,
            "max_circuit_gates": MAX_CIRCUIT_GATES,
            "circuit_atoms": sorted(CORE_ATOM_NAMES),
            "terminal_count": len(terminals),
        },
        "incumbent_slice_hessian_rank": 35,
        "depth_log": depth_log,
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
            else "no terminal in the stored width-64 circuit census beats 35"
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"BEST objective={best_objective}")
    print(f"BEST atoms={best[0].circuit_atoms}")
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
