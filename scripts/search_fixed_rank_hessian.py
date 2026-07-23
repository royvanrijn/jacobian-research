#!/usr/bin/env python3
"""Search shared-factor BCW traces with generic Hessian rank as objective.

This deliberately differs from ``search_rank_aware_bcw.py``: that script
minimizes the essential dimension of the cubic map, whereas this experiment
minimizes the generic rank of the Hessian of the associated quartic
cotangent/symmetric lift.  Every rank is a deterministic good-prime lower
bound until a winning trace is frozen and given an exact syzygy certificate.
"""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path

from audit_fixed_rank_hessian_witness import (
    cotangent_hessian,
    specialization_profile,
)
from rank_compressed_bcw_homogenization import (
    constant_kernel_quotient,
    extract_quadratic_cubic,
    factor_cubic_output,
    rank_compressed_homogeneous_map,
)
from search_rank_aware_bcw import State, initial_state, score, state_key, transitions
from verify_shared_bcw_33_route import high_terms


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "artifacts" / "generated-results" / "fixed_rank_hessian_candidate.json"


def poly_components(components) -> list[dict[tuple[int, ...], object]]:
    return [
        {tuple(exponent): coefficient for exponent, coefficient in poly.terms() if coefficient}
        for poly in components
    ]


def hessian_profile(components, seed: int = 20_260_723) -> tuple[int, int, int]:
    jacobian, upper_left, hessian = cotangent_hessian(poly_components(components))
    return specialization_profile(jacobian, upper_left, hessian, seed)


def terminal_profile(state: State) -> dict[str, int]:
    quadratic, cubic = extract_quadratic_cubic(state.expressions, state.variables)
    factorization = factor_cubic_output(cubic)
    ambient_variables, ambient_h = rank_compressed_homogeneous_map(
        state.variables, quadratic, factorization
    )
    quotient = constant_kernel_quotient(ambient_variables, ambient_h)
    ambient = hessian_profile(ambient_h)
    descended = hessian_profile(quotient.quotient_h)
    return {
        "introduced_variables": state.introduced,
        "cubic_output_rank": len(factorization.c),
        "ambient_dimension": len(ambient_variables),
        "constant_kernel_dimension": quotient.kernel.cols,
        "quotient_dimension": len(quotient.quotient_variables),
        "ambient_jacobian_rank_mod_1000003": ambient[0],
        "ambient_hessian_rank_mod_1000003": ambient[2],
        "quotient_jacobian_rank_mod_1000003": descended[0],
        "quotient_hessian_rank_mod_1000003": descended[2],
    }


def encoded_plan(state: State, profile: dict[str, int]) -> dict[str, object]:
    return {
        "format": "fixed-rank-hessian-search-candidate-v1",
        "certification_status": (
            "good-prime search candidate; requires frozen collision replay and exact syzygy certificate"
        ),
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
    parser.add_argument("--width", type=int, default=64)
    parser.add_argument("--max-steps", type=int, default=24)
    parser.add_argument("--incumbent-rank", type=int, default=38)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if min(args.width, args.max_steps, args.incumbent_rank) < 1:
        parser.error("width, max-steps, and incumbent-rank must be positive")

    frontier = [initial_state()]
    best: tuple[tuple[object, ...], State, dict[str, int]] | None = None
    histogram: Counter[tuple[int, int]] = Counter()
    completed_count = 0
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
            profile = terminal_profile(candidate)
            completed_count += 1
            histogram[
                (
                    profile["ambient_hessian_rank_mod_1000003"],
                    profile["quotient_hessian_rank_mod_1000003"],
                )
            ] += 1
            key = (
                profile["quotient_hessian_rank_mod_1000003"],
                profile["ambient_hessian_rank_mod_1000003"],
                profile["quotient_dimension"],
                candidate.plan,
            )
            if best is None or key < best[0]:
                best = (key, candidate, profile)

        frontier = sorted(
            deduplicated.values(), key=lambda candidate: score(candidate, "rank-first")
        )[: args.width]
        best_rank = best[2]["quotient_hessian_rank_mod_1000003"] if best else None
        print(
            f"depth={depth} generated={generated} unique={len(deduplicated)} "
            f"kept={len(frontier)} completed_profiled={completed_count} best_rank={best_rank}",
            flush=True,
        )
        if not frontier:
            break

    print(
        "rank_histogram="
        + json.dumps(
            {
                f"ambient_{ambient}_quotient_{quotient}": count
                for (ambient, quotient), count in sorted(histogram.items())
            },
            sort_keys=True,
        )
    )
    if best is not None and best[2]["quotient_hessian_rank_mod_1000003"] < args.incumbent_rank:
        args.output.write_text(json.dumps(encoded_plan(best[1], best[2]), indent=2) + "\n")
        print(
            f"IMPROVEMENT rank={best[2]['quotient_hessian_rank_mod_1000003']}; "
            f"wrote {args.output.relative_to(ROOT)}"
        )
    else:
        print(
            f"NO IMPROVEMENT within width={args.width}, max_steps={args.max_steps}; "
            f"incumbent rank remains {args.incumbent_rank}"
        )


if __name__ == "__main__":
    main()
