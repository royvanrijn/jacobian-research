#!/usr/bin/env python3
"""Exact direction-collision boundary for degree-ten four-channel profiles.

The pairwise-distinct charts are treated by
``research_gvc3_degree10_distinct_four_channels.py``.  Here we check every
nonterminal set partition of the four channel directions for the four
genuinely new degree-ten profiles.  A cutoff is discovered by quotient
saturation modulo 101, replayed modulo 103 and 107, and promoted only when
msolve returns the literal reduced basis ``[1]`` over the rationals.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import sympy as sp

from research_gvc3_many_coherent_channels import (
    exact_unit,
    moment,
    modular_saturation_cutoff,
    primitive_polynomial,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "gvc3_degree10_four_channel_collisions.json"
)
ALL_PROFILES = tuple(itertools.combinations((2, 4, 6, 8, 10), 4))
PROFILES = tuple(profile for profile in ALL_PROFILES if profile != (2, 4, 6, 8))
PRIMES = (101, 103, 107)


def set_partitions(size: int) -> tuple[tuple[tuple[int, ...], ...], ...]:
    answer: list[tuple[tuple[int, ...], ...]] = []

    def visit(index: int, blocks: list[list[int]]) -> None:
        if index == size:
            answer.append(tuple(tuple(block) for block in blocks))
            return
        for block in blocks:
            block.append(index)
            visit(index + 1, blocks)
            block.pop()
        blocks.append([index])
        visit(index + 1, blocks)
        blocks.pop()

    visit(0, [])
    return tuple(answer)


PARTITIONS = tuple(
    groups for groups in set_partitions(4) if len(groups) not in (1, 4)
)


def compile_chart(
    degrees: tuple[int, ...],
    groups: tuple[tuple[int, ...], ...],
    max_order: int,
) -> dict[str, object]:
    coefficients = sp.symbols("a0:4")
    a0, a1, a2, a3 = coefficients
    variables = (a1, a2, a3)
    equations: dict[int, sp.Expr] = {}
    for order in range(2, max_order + 1):
        expression = moment(degrees, groups, order, coefficients, ())
        polynomial = primitive_polynomial(expression.subs(a0, 1), variables)
        if polynomial != 0:
            equations[order] = polynomial
    return {
        "degrees": degrees,
        "groups": groups,
        "variables": variables,
        "equations": equations,
        "saturation": a1 * a2 * a3,
        "moment_term_counts": {
            str(order): len(sp.Poly(equation, *variables).terms())
            for order, equation in equations.items()
        },
        "moment_sha256": {
            str(order): hashlib.sha256(str(equation).encode()).hexdigest()
            for order, equation in equations.items()
        },
    }


def certify_chart(
    data: dict[str, object],
    discovery_prime: int,
    replay_primes: tuple[int, ...],
    singular: str,
    modular_timeout: int,
    exact_timeout: int,
) -> dict[str, object]:
    equations = data["equations"]
    variables = data["variables"]
    saturation = data["saturation"]
    assert isinstance(equations, dict)
    assert isinstance(variables, tuple)
    candidates = sorted(equations)
    discovery_attempts: list[dict[str, object]] = []
    cutoff = 0
    selected: dict[int, sp.Expr] = {}
    for candidate in candidates:
        selected = {
            order: equation
            for order, equation in equations.items()
            if order <= candidate
        }
        result = modular_saturation_cutoff(
            selected,
            saturation,
            variables,
            discovery_prime,
            singular,
            modular_timeout,
        )
        discovery_attempts.append(result)
        if result.get("status") == "completed" and result.get("unit") == 1:
            cutoff = candidate
            break
    if cutoff == 0:
        raise RuntimeError(
            f"no modular unit through the compiled cutoff for "
            f"{data['degrees']} {data['groups']}"
        )

    modular_results = [discovery_attempts[-1]]
    replay_attempts: dict[str, list[dict[str, object]]] = {}
    modular_cutoffs = {str(discovery_prime): cutoff}
    for prime in replay_primes:
        attempts: list[dict[str, object]] = []
        replay_cutoff = 0
        replay: dict[str, object] | None = None
        for candidate in (order for order in candidates if order >= cutoff):
            replay_selected = {
                order: equation
                for order, equation in equations.items()
                if order <= candidate
            }
            replay = modular_saturation_cutoff(
                replay_selected,
                saturation,
                variables,
                prime,
                singular,
                modular_timeout,
            )
            attempts.append(replay)
            if replay.get("status") == "completed" and replay.get("unit") == 1:
                replay_cutoff = candidate
                break
        replay_attempts[str(prime)] = attempts
        if replay is None or replay_cutoff == 0:
            raise RuntimeError(
                f"modular replay found no unit at p={prime} for "
                f"{data['degrees']} {data['groups']}"
            )
        modular_cutoffs[str(prime)] = replay_cutoff
        modular_results.append(replay)

    exact = exact_unit(
        equations,
        saturation,
        variables,
        cutoff,
        exact_timeout,
        1,
    )
    if exact.get("unit") != 1:
        raise RuntimeError(
            f"exact promotion failed for {data['degrees']} {data['groups']}"
        )
    return {
        "harmonic_degrees": list(data["degrees"]),
        "groups": [list(block) for block in data["groups"]],
        "directions": len(data["groups"]),
        "cutoff": cutoff,
        "modular_cutoffs": modular_cutoffs,
        "nonzero_moment_orders": sorted(selected),
        "moment_term_counts": data["moment_term_counts"],
        "moment_sha256": data["moment_sha256"],
        "discovery_attempts": discovery_attempts,
        "replay_attempts": replay_attempts,
        "modular_unit_results": modular_results,
        "exact_Q_result": exact,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-order", type=int, default=12)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--modular-timeout", type=int, default=180)
    parser.add_argument("--exact-timeout", type=int, default=300)
    parser.add_argument("--singular", default="Singular")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()

    compiled = [
        compile_chart(profile, groups, arguments.max_order)
        for profile in PROFILES
        for groups in PARTITIONS
    ]
    results: list[dict[str, object]] = []
    with ThreadPoolExecutor(max_workers=arguments.workers) as executor:
        tasks = {
            executor.submit(
                certify_chart,
                data,
                PRIMES[0],
                PRIMES[1:],
                arguments.singular,
                arguments.modular_timeout,
                arguments.exact_timeout,
            ): data
            for data in compiled
        }
        for future in as_completed(tasks):
            result = future.result()
            results.append(result)
            print(
                "PASS",
                ",".join(map(str, result["harmonic_degrees"])),
                result["groups"],
                f"cutoff={result['cutoff']}",
                flush=True,
            )
    results.sort(key=lambda result: (result["harmonic_degrees"], result["groups"]))
    assert len(results) == len(PROFILES) * len(PARTITIONS) == 52
    assert all(result["exact_Q_result"]["unit"] == 1 for result in results)

    cutoff_distribution: dict[str, int] = {}
    for result in results:
        key = str(result["cutoff"])
        cutoff_distribution[key] = cutoff_distribution.get(key, 0) + 1
    artifact = {
        "format": "gvc3-degree10-four-channel-collisions-v1",
        "status": "exact characteristic-zero direction-collision obstruction",
        "balanced_degree": 10,
        "laplacian_power": 5,
        "profiles": [list(profile) for profile in PROFILES],
        "partitions_per_profile": len(PARTITIONS),
        "total_charts": len(results),
        "coefficient_saturation": "a1*a2*a3",
        "max_compiled_moment_order": arguments.max_order,
        "modular_discovery_primes": list(PRIMES),
        "cutoff_distribution": cutoff_distribution,
        "maximum_exact_cutoff": max(result["cutoff"] for result in results),
        "charts": results,
        "coefficient_boundaries": (
            "a zero coefficient reduces to the exact two/three-channel theorem GVC3IHC"
        ),
        "one_direction_terminal_reason": (
            "phase weight is one-sided, so every fixed multiplier is eventually killed"
        ),
        "not_in_scope": [
            "pairwise-distinct characteristic-zero promotion",
            "the degree-ten five-channel profile",
            "noncoherent or repeated harmonic profiles",
        ],
        "conclusion": (
            "all 52 nonterminal direction-collision coefficient-torus charts "
            "are empty over Q at their declared finite moment cutoffs"
        ),
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n"
    )
    print(
        "PASS 52 exact Q collision charts; cutoff distribution",
        cutoff_distribution,
    )


if __name__ == "__main__":
    main()
