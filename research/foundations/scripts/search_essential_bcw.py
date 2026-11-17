#!/usr/bin/env python3
"""Search exact shared-factor BCW traces by final essential dimension.

Every completed quadratic--cubic trace is rank-compressed, homogenized,
quotiented by the constant kernel of its homogeneous Jacobian, and given a
coordinate-generated invariant-row-module diagnostic.  Output remains a
search candidate until a frozen generator and independent replay certify it.
"""

from __future__ import annotations

import argparse
from itertools import combinations
import json
from pathlib import Path

import sympy as sp

from rank_compressed_bcw_homogenization import (
    constant_kernel_quotient,
    cyclic_invariant_row_module_dimensions,
    extract_quadratic_cubic,
    factor_cubic_output,
    rank_compressed_homogeneous_map,
)
from search_rank_aware_bcw import State, initial_state, score, state_key, transitions
from verify_shared_bcw_33_route import apply_shared_step, dense_factor, high_terms


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "artifacts" / "generated-results" / "essential_bcw_candidate.json"
RESIDUAL_LABELS = ("r_gamma", "r_W", "r_C", "r_s", "r_t")


def tangent_residual_matroid_profile() -> dict[str, object]:
    """Exact collision-row matroid used to prune tangent-core terminal states."""

    x, y, z, A, B, C, W, gamma = sp.symbols("x y z A B C W gamma")
    variables = (x, y, z, A, B, C, W, gamma)
    u = 1 + x * y
    source_gamma = 1 - sp.Rational(3, 2) * x * y - sp.Rational(1, 2) * x**2 * z
    source_W = u * source_gamma
    residuals = sp.Matrix(
        [
            gamma - source_gamma,
            W - u * gamma,
            C - x * gamma,
            B * C - (2 * W - 3 * W**2 + gamma),
            A * C**2 - (W * gamma + W**2 - 2 * W**3),
        ]
    )
    lifts = []
    for point in (
        (sp.Integer(0), sp.Integer(0), -sp.Rational(1, 4)),
        (sp.Integer(1), -sp.Rational(3, 2), sp.Rational(13, 2)),
        (-sp.Integer(1), sp.Rational(3, 2), sp.Rational(13, 2)),
    ):
        source = dict(zip((x, y, z), point))
        lifts.append(
            {
                **source,
                A: -sp.Rational(1, 4),
                B: 0,
                C: 0,
                W: source_W.subs(source),
                gamma: source_gamma.subs(source),
            }
        )
    jacobian = residuals.jacobian(variables)
    admissible: dict[str, list[list[str]]] = {}
    for size in range(1, 6):
        rows = []
        for indices in combinations(range(5), size):
            ranks = [jacobian[list(indices), :].subs(lift).rank() for lift in lifts]
            if all(rank == size for rank in ranks):
                rows.append([RESIDUAL_LABELS[index] for index in indices])
        admissible[str(size)] = rows
    assert admissible["4"] == [] and admissible["5"] == []
    assert admissible["3"] == [
        ["r_gamma", "r_W", "r_C"],
        ["r_gamma", "r_C", "r_s"],
    ]
    rt_ranks = [jacobian[[4], :].subs(lift).rank() for lift in lifts]
    assert rt_ranks == [1, 0, 0]
    return {
        "residual_labels": list(RESIDUAL_LABELS),
        "admissible_terminal_subsets": admissible,
        "r_t_row_ranks_at_collision_lifts": rt_ranks,
    }


def evaluate_monomial(exponents: tuple[int, ...], point: list[sp.Expr]) -> sp.Expr:
    return sp.prod(value**exponent for value, exponent in zip(point, exponents))


def replay_collision_points(state: State) -> list[list[sp.Expr]]:
    """Replay the exact source lifts associated with a search plan."""

    replay = initial_state()
    points: list[list[sp.Expr]] = [
        [sp.Integer(0), sp.Integer(0), -sp.Rational(1, 8)],
        [sp.Integer(1), -sp.Rational(3, 4), sp.Rational(13, 4)],
        [-sp.Integer(1), sp.Rational(3, 4), sp.Rational(13, 4)],
    ]
    for component, first_support, second_support in state.plan:
        old_dimension = len(replay.variables)
        first = dense_factor(list(first_support), old_dimension)
        second = dense_factor(list(second_support), old_dimension)
        removed = tuple(a + b for a, b in zip(first, second))
        coefficient = sp.Poly(
            replay.expressions[component], *replay.variables, domain=sp.QQ
        ).coeff_monomial(removed)
        assert coefficient
        selected = (component, removed, coefficient, sum(removed))
        chosen = apply_shared_step(
            replay.expressions,
            replay.variables,
            replay.registry,
            selected,
            (first, second),
            replay.introduced,
        )
        assert chosen is not None
        metadata = chosen[4]
        for point in points:
            additions = []
            for record in metadata["new_factors"]:
                dense = [0] * old_dimension
                for index, exponent in record["factor"]:
                    dense[index] = exponent
                additions.append(-evaluate_monomial(tuple(dense), point))
            point.extend(additions)
        replay = State(chosen[0], chosen[1], chosen[2], chosen[3], replay.plan + ((component, first_support, second_support),))
    assert state_key(replay) == state_key(state)
    return points


def essential_profile(state: State) -> dict[str, object]:
    quadratic, cubic = extract_quadratic_cubic(state.expressions, state.variables)
    factorization = factor_cubic_output(cubic)
    homogeneous_variables, homogeneous_h = rank_compressed_homogeneous_map(
        state.variables, quadratic, factorization
    )
    quotient = constant_kernel_quotient(homogeneous_variables, homogeneous_h)

    source_points = replay_collision_points(state)
    source_images = [
        [
            sp.Poly(expression, *state.variables, domain=sp.QQ).eval(
                dict(zip(state.variables, point))
            )
            for expression in state.expressions
        ]
        for point in source_points
    ]
    assert source_images[0] == source_images[1] == source_images[2]
    homogeneous_points: list[sp.Matrix] = []
    for point in source_points:
        substitution = dict(zip(state.variables, point))
        y = [poly.eval(substitution) for poly in factorization.c]
        homogeneous_points.append(sp.Matrix(point + y + [sp.Integer(1)]))
    projected = [quotient.B * point for point in homogeneous_points]
    collision_separated = len({tuple(point) for point in projected}) == 3
    assert collision_separated
    ambient_images = []
    for point in homogeneous_points:
        substitution = dict(zip(homogeneous_variables, point))
        ambient_images.append(
            point
            + sp.Matrix(
                [poly.as_expr().subs(substitution, simultaneous=True) for poly in homogeneous_h]
            )
        )
    assert ambient_images[0] == ambient_images[1] == ambient_images[2]
    projected_images = [quotient.B * point for point in ambient_images]
    assert projected_images[0] == projected_images[1] == projected_images[2]
    assert constant_kernel_quotient(
        quotient.quotient_variables, quotient.quotient_h
    ).kernel.cols == 0

    module_dimensions = cyclic_invariant_row_module_dimensions(
        quotient.quotient_h, quotient.quotient_variables
    )
    cse_replacements, _ = sp.cse(state.expressions, order="canonical")
    return {
        "introduced_variables": state.introduced,
        "structural_dag_nodes": len(state.registry),
        "canonical_cse_nodes": len(cse_replacements),
        "expanded_operation_count": sum(int(sp.count_ops(value)) for value in state.expressions),
        "cubic_output_rank": len(factorization.c),
        "homogeneous_dimension": len(homogeneous_variables),
        "constant_kernel_dimension": quotient.kernel.cols,
        "final_essential_dimension": len(quotient.quotient_variables),
        "projected_collision_separated": collision_separated,
        "projected_collision_common_image_verified": True,
        "quotient_constant_kernel_dimension": 0,
        "coordinate_cyclic_invariant_row_module_dimensions_mod_1000003": list(module_dimensions),
        "invariant_module_scope": (
            "good-prime coordinate-generated cyclic row modules only; shortlisted "
            "candidates require exact QQ classification of all common invariant subspaces"
        ),
        "tangent_residual_matroid": tangent_residual_matroid_profile(),
        "search_family_scope": (
            "exact shared-factor triangular shears with exposed-factor DAG nodes; "
            "arbitrary polynomial DAG rewrites from the tangent core are not enumerated"
        ),
    }


def payload(state: State, profile: dict[str, object]) -> dict[str, object]:
    return {
        "format": "essential-bcw-search-candidate-v1",
        "certification_status": "search candidate; requires frozen generator and independent replay",
        "profile": profile,
        "plan": [
            {
                "component": component,
                "first_factor": [list(pair) for pair in first],
                "second_factor": [list(pair) for pair in second],
            }
            for component, first, second in state.plan
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--width", type=int, default=24)
    parser.add_argument("--max-steps", type=int, default=24)
    parser.add_argument("--incumbent-dimension", type=int, default=21)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if min(args.width, args.max_steps, args.incumbent_dimension) < 1:
        parser.error("width, max-steps, and incumbent dimension must be positive")

    frontier = [initial_state()]
    best_state: State | None = None
    best_profile: dict[str, object] | None = None
    profiled = 0
    for depth in range(1, args.max_steps + 1):
        deduplicated: dict[tuple[object, ...], State] = {}
        completed: list[State] = []
        generated = 0
        for state in frontier:
            for candidate in transitions(state):
                generated += 1
                signature, _ = high_terms(candidate.expressions, candidate.variables)
                if signature[0] <= 3:
                    completed.append(candidate)
                    continue
                key = state_key(candidate)
                previous = deduplicated.get(key)
                if previous is None or candidate.plan < previous.plan:
                    deduplicated[key] = candidate

        for candidate in completed:
            profile = essential_profile(candidate)
            profiled += 1
            if best_profile is None or (
                profile["final_essential_dimension"], candidate.plan
            ) < (best_profile["final_essential_dimension"], best_state.plan):
                best_state, best_profile = candidate, profile

        frontier = sorted(
            deduplicated.values(), key=lambda candidate: score(candidate, "rank-first")
        )[: args.width]
        best_dimension = (
            best_profile["final_essential_dimension"] if best_profile else args.incumbent_dimension
        )
        print(
            f"depth={depth} generated={generated} unique={len(deduplicated)} "
            f"kept={len(frontier)} completed_profiled={profiled} "
            f"best_essential_dimension={best_dimension}",
            flush=True,
        )
        if not frontier:
            break

    if best_state is not None and best_profile is not None:
        best_dimension = int(best_profile["final_essential_dimension"])
        if best_dimension < args.incumbent_dimension:
            args.output.write_text(json.dumps(payload(best_state, best_profile), indent=2) + "\n")
            print(
                f"IMPROVEMENT essential_dimension={best_dimension}; "
                f"wrote {args.output.relative_to(ROOT)}"
            )
            return
    print(
        f"NO IMPROVEMENT within width={args.width}, max_steps={args.max_steps}; "
        f"certified incumbent remains dimension {args.incumbent_dimension}"
    )


if __name__ == "__main__":
    main()
