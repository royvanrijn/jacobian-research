#!/usr/bin/env python3
"""Search shared-factor BCW traces for smaller essential dimension.

The certified trace uses s=13 exposed variables and has cubic-output rank
k=7.  The beam uses ``s+rank(C)`` as a cheap partial-state heuristic, but every
completed trace is rank-compressed, homogenized, and constant-kernel
quotiented before comparison.  It emits a JSON plan only when the final
essential dimension strictly improves the incumbent; the emitted plan must
still be frozen into a generator and independently replayed before changing
any theorem statement.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Iterable

import sympy as sp

from rank_compressed_bcw_homogenization import (
    constant_kernel_quotient,
    extract_quadratic_cubic,
    factor_cubic_output,
    homogeneous_part,
    modular_jacobian_coefficient_rank,
    rank_compressed_homogeneous_map,
)
from verify_shared_bcw_33_route import apply_shared_step, candidate_splits, high_terms


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "artifacts" / "generated-results" / "rank_aware_bcw_candidate.json"


PlanStep = tuple[int, tuple[tuple[int, int], ...], tuple[tuple[int, int], ...]]


@dataclass
class State:
    expressions: list[sp.Expr]
    variables: list[sp.Symbol]
    registry: dict[tuple[int, ...], int]
    introduced: int
    plan: tuple[PlanStep, ...]


def initial_state() -> State:
    x, y, z = sp.symbols("x y z")
    variables = [x, y, z]
    v = 1 + 2 * x * y
    long_map = [
        v**3 * z + 4 * y**2 * v * (2 + 3 * x * y),
        y + 3 * x * v**2 * z + 12 * x * y**2 * (2 + 3 * x * y),
        -x + 3 * x**2 * y + x**3 * z,
    ]
    expressions = [sp.expand(-long_map[2]), sp.expand(long_map[1]), sp.expand(long_map[0])]
    return State(expressions, variables, {}, 0, ())


def support(exponents: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    return tuple((index, exponent) for index, exponent in enumerate(exponents) if exponent)


def cubic_rank(state: State) -> int:
    cubic = [homogeneous_part(expression, state.variables, 3) for expression in state.expressions]
    monomials = sorted(
        {
            exponents
            for poly in cubic
            for exponents, coefficient in poly.terms()
            if coefficient
        }
    )
    if not monomials:
        return 0
    matrix = sp.polys.matrices.DomainMatrix.from_list_sympy(
        len(cubic),
        len(monomials),
        [[poly.coeff_monomial(exponents) for exponents in monomials] for poly in cubic],
    )
    return matrix.rank()


def state_key(state: State) -> tuple[object, ...]:
    polynomial_key = tuple(
        tuple(sp.Poly(expression, *state.variables, domain=sp.QQ).terms())
        for expression in state.expressions
    )
    registry_key = tuple(sorted(state.registry.items()))
    return polynomial_key, registry_key


def score(state: State, mode: str = "rank-first") -> tuple[object, ...]:
    signature, _ = high_terms(state.expressions, state.variables)
    rank = cubic_rank(state)
    if mode == "rank-first":
        return (
            signature[0],
            state.introduced + rank,
            state.introduced,
            rank,
            signature[1],
            signature[2],
            -len(state.registry),
            state.plan,
        )
    return (
        signature[0],
        signature[1],
        signature[2],
        state.introduced + rank,
        state.introduced,
        rank,
        -len(state.registry),
        state.plan,
    )


def transitions(state: State) -> Iterable[State]:
    signature, terms = high_terms(state.expressions, state.variables)
    maximum = signature[0]
    selected_terms = sorted(
        (term for term in terms if term[3] == maximum),
        key=lambda term: (term[0], term[1], sp.default_sort_key(term[2])),
    )
    for selected in selected_terms:
        component, _, _, _ = selected
        for first, second in candidate_splits(selected[1]):
            result = apply_shared_step(
                state.expressions,
                state.variables,
                state.registry,
                selected,
                (first, second),
                state.introduced,
            )
            if result is None:
                continue
            expressions, variables, registry, introduced, _ = result
            step: PlanStep = (component, support(first), support(second))
            yield State(expressions, variables, registry, introduced, state.plan + (step,))


def essential_dimension(state: State) -> tuple[int, int]:
    quadratic, cubic = extract_quadratic_cubic(state.expressions, state.variables)
    factorization = factor_cubic_output(cubic)
    variables, homogeneous = rank_compressed_homogeneous_map(
        state.variables, quadratic, factorization
    )
    quotient = constant_kernel_quotient(variables, homogeneous)
    return len(quotient.quotient_variables), quotient.kernel.cols


def modular_essential_lower_bound(state: State) -> int:
    quadratic, cubic = extract_quadratic_cubic(state.expressions, state.variables)
    factorization = factor_cubic_output(cubic)
    variables, homogeneous = rank_compressed_homogeneous_map(
        state.variables, quadratic, factorization
    )
    return modular_jacobian_coefficient_rank(homogeneous, variables)


def plan_json(
    state: State, rank: int, final_dimension: int, kernel_dimension: int
) -> dict[str, object]:
    return {
        "format": "rank-aware-bcw-search-candidate-v1",
        "certification_status": "search candidate; requires generator and independent replay",
        "introduced_variables": state.introduced,
        "cubic_output_rank": rank,
        "objective_s_plus_rank": state.introduced + rank,
        "homogeneous_dimension": 4 + state.introduced + rank,
        "constant_kernel_dimension": kernel_dimension,
        "final_essential_dimension": final_dimension,
        "plan": [
            {
                "component": component,
                "first_factor": [list(pair) for pair in first],
                "second_factor": [list(pair) for pair in second],
            }
            for component, first, second in state.plan
        ],
    }


def two_parameter_values(state: State) -> tuple[int, ...]:
    """Deterministic necessary tests for det(I+sJQ+tJC)=1."""
    quadratic = [homogeneous_part(expression, state.variables, 2) for expression in state.expressions]
    cubic = [homogeneous_part(expression, state.variables, 3) for expression in state.expressions]
    jq = sp.Matrix([poly.as_expr() for poly in quadratic]).jacobian(state.variables)
    jc = sp.Matrix([poly.as_expr() for poly in cubic]).jacobian(state.variables)
    sample_points = [
        [sp.Integer(1)] * len(state.variables),
        [sp.Integer(index % 2) for index in range(len(state.variables))],
    ]
    parameter_pairs = [(0, 1), (1, 0), (1, 2), (2, 1)]
    values: list[int] = []
    for point in sample_points:
        substitution = dict(zip(state.variables, point))
        jq_value = jq.subs(substitution)
        jc_value = jc.subs(substitution)
        for s_value, t_value in parameter_pairs:
            determinant = (
                sp.eye(len(state.variables)) + s_value * jq_value + t_value * jc_value
            ).det()
            assert determinant.is_Integer
            values.append(int(determinant))
    return tuple(values)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--width", type=int, default=128)
    parser.add_argument("--max-steps", type=int, default=24)
    parser.add_argument(
        "--incumbent", type=int, default=21, help="certified final essential dimension"
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--score-mode",
        choices=("rank-first", "degree-first"),
        default="rank-first",
        help="rank-first directly prioritizes s+rank after maximum degree",
    )
    parser.add_argument(
        "--check-two-parameter",
        action="store_true",
        help="test completed traces for det(I+sJQ+tJC)=1 at exact sample points",
    )
    args = parser.parse_args()
    if args.width < 1 or args.max_steps < 1:
        parser.error("width and max-steps must be positive")

    frontier = [initial_state()]
    best_complete: State | None = None
    best_objective = args.incumbent
    two_parameter_tested = 0
    two_parameter_best_failures = 9
    two_parameter_candidate: State | None = None
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
            rank = cubic_rank(candidate)
            if modular_essential_lower_bound(candidate) >= best_objective:
                continue
            objective, _ = essential_dimension(candidate)
            if objective < best_objective or (
                objective == best_objective
                and (best_complete is None or candidate.plan < best_complete.plan)
            ):
                best_objective = objective
                best_complete = candidate
            if args.check_two_parameter:
                values = two_parameter_values(candidate)
                failures = sum(value != 1 for value in values)
                two_parameter_tested += 1
                two_parameter_best_failures = min(two_parameter_best_failures, failures)
                if failures == 0 and two_parameter_candidate is None:
                    two_parameter_candidate = candidate

        ranked = sorted(deduplicated.values(), key=lambda state: score(state, args.score_mode))
        frontier = ranked[: args.width]
        if frontier:
            lead = frontier[0]
            lead_signature, _ = high_terms(lead.expressions, lead.variables)
            lead_rank = cubic_rank(lead)
            print(
                f"depth={depth} generated={generated} unique={len(deduplicated)} "
                f"kept={len(frontier)} lead_high={lead_signature} "
                f"lead_s={lead.introduced} lead_rank={lead_rank} "
                f"lead_s_plus_rank={lead.introduced + lead_rank} "
                f"best_essential_dimension={best_objective}",
                flush=True,
            )
        else:
            print(
                f"depth={depth} generated={generated} unique={len(deduplicated)} "
                f"kept=0 best_essential_dimension={best_objective}",
                flush=True,
            )
            break

    if best_complete is not None and best_objective < args.incumbent:
        rank = cubic_rank(best_complete)
        final_dimension, kernel_dimension = essential_dimension(best_complete)
        args.output.write_text(
            json.dumps(
                plan_json(best_complete, rank, final_dimension, kernel_dimension), indent=2
            )
            + "\n"
        )
        print(
            f"IMPROVEMENT essential_dimension={best_objective}; "
            f"wrote {args.output.relative_to(ROOT)}"
        )
    else:
        print(
            f"NO IMPROVEMENT within width={args.width}, max_steps={args.max_steps}; "
            f"certified incumbent remains essential_dimension={args.incumbent}"
        )
    if args.check_two_parameter:
        print(
            f"TWO_PARAMETER tested={two_parameter_tested} "
            f"best_failed_exact_samples={two_parameter_best_failures}"
        )
        if two_parameter_candidate is not None:
            rank = cubic_rank(two_parameter_candidate)
            final_dimension, kernel_dimension = essential_dimension(two_parameter_candidate)
            candidate_output = args.output.with_name("two_parameter_bcw_candidate.json")
            payload = plan_json(
                two_parameter_candidate, rank, final_dimension, kernel_dimension
            )
            payload["two_parameter_sample_status"] = "passed 8 exact necessary tests; proof required"
            candidate_output.write_text(json.dumps(payload, indent=2) + "\n")
            print(
                "TWO_PARAMETER CANDIDATE passed all exact samples; "
                f"wrote {candidate_output.relative_to(ROOT)}"
            )


if __name__ == "__main__":
    main()
