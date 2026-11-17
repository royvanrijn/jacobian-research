#!/usr/bin/env python3
"""Beam-search BCW terminal traces by post-quotient Jacobian profile.

This is a bounded search diagnostic, not an exact rank certificate.  Every
terminal map is rank-compressed, homogenized, constant-kernel quotiented, and
evaluated at one deterministic point over a good prime.  Shortlisted maps
must be frozen and audited over QQ, as in generate_low_complexity_bcw_21.py.
"""

from __future__ import annotations

import argparse
from collections import Counter

import sympy as sp

from rank_compressed_bcw_homogenization import (
    constant_kernel_quotient,
    extract_quadratic_cubic,
    factor_cubic_output,
    rank_compressed_homogeneous_map,
)
from search_essential_bcw import replay_collision_points
from search_rank_aware_bcw import initial_state, score, state_key, transitions
from verify_shared_bcw_33_route import high_terms


def modular_rank(matrix: list[list[int]], prime: int) -> int:
    rows = [row.copy() for row in matrix]
    if not rows:
        return 0
    m, n = len(rows), len(rows[0])
    rank = 0
    for column in range(n):
        pivot = next((row for row in range(rank, m) if rows[row][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = pow(rows[rank][column], -1, prime)
        rows[rank] = [value * inverse % prime for value in rows[rank]]
        for row in range(rank + 1, m):
            if rows[row][column]:
                scale = rows[row][column]
                rows[row] = [
                    (left - scale * right) % prime
                    for left, right in zip(rows[row], rows[rank])
                ]
        rank += 1
    return rank


def modular_product(left: list[list[int]], right: list[list[int]], prime: int) -> list[list[int]]:
    n = len(left)
    answer = [[0] * n for _ in range(n)]
    nonzero_right = [[(j, value) for j, value in enumerate(row) if value] for row in right]
    for i, row in enumerate(left):
        for k, value in enumerate(row):
            if value:
                for j, right_value in nonzero_right[k]:
                    answer[i][j] = (answer[i][j] + value * right_value) % prime
    return answer


def residue(value, prime: int) -> int:
    return int(value.p) * pow(int(value.q), -1, prime) % prime


def terminal_profile(state, prime: int) -> tuple[int, int, int, int, int] | None:
    quadratic, cubic = extract_quadratic_cubic(state.expressions, state.variables)
    factorization = factor_cubic_output(cubic)
    variables, homogeneous = rank_compressed_homogeneous_map(
        state.variables, quadratic, factorization
    )
    quotient = constant_kernel_quotient(variables, homogeneous)
    projected_points = []
    for point in replay_collision_points(state):
        substitution = dict(zip(state.variables, point))
        ambient = sp.Matrix(
            point
            + [poly.eval(substitution) for poly in factorization.c]
            + [sp.Integer(1)]
        )
        projected_points.append(quotient.B * ambient)
    if len({tuple(point) for point in projected_points}) != 3:
        return None
    qvars, qh = quotient.quotient_variables, quotient.quotient_h
    point = {variable: (i * i + 3 * i + 5) % prime for i, variable in enumerate(qvars)}
    jacobian = [
        [residue(poly.diff(variable).eval(point), prime) for variable in qvars]
        for poly in qh
    ]
    power = jacobian
    power_ranks: list[int] = []
    for _ in range(len(qvars)):
        power_ranks.append(modular_rank(power, prime))
        if power_ranks[-1] == 0:
            break
        power = modular_product(power, jacobian, prime)
    terms = sum(sum(bool(coefficient) for _, coefficient in poly.terms()) for poly in qh)
    return (
        len(qvars),
        power_ranks[0],
        len(power_ranks),
        terms,
        quotient.kernel.cols,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--width", type=int, default=24)
    parser.add_argument("--max-steps", type=int, default=24)
    parser.add_argument("--prime", type=int, default=1_000_003)
    args = parser.parse_args()

    frontier = [initial_state()]
    profiles: Counter[tuple[int, int, int, int, int]] = Counter()
    for depth in range(1, args.max_steps + 1):
        deduplicated = {}
        completed = []
        for state in frontier:
            for candidate in transitions(state):
                signature, _ = high_terms(candidate.expressions, candidate.variables)
                if signature[0] <= 3:
                    completed.append(candidate)
                    continue
                key = state_key(candidate)
                previous = deduplicated.get(key)
                if previous is None or candidate.plan < previous.plan:
                    deduplicated[key] = candidate
        for candidate in completed:
            profile = terminal_profile(candidate, args.prime)
            if profile is not None:
                profiles[profile] += 1
        frontier = sorted(
            deduplicated.values(), key=lambda candidate: score(candidate, "rank-first")
        )[: args.width]
        print(
            f"depth={depth} unique={len(deduplicated)} kept={len(frontier)} "
            f"terminal_profiled={sum(profiles.values())}",
            flush=True,
        )
        if not frontier:
            break

    print("profile=(quotient_dimension,rank,index,terms,ambient_constant_kernel)")
    for profile, count in sorted(profiles.items()):
        print(f"count={count} profile={profile}")
    print(
        f"DIAGNOSTIC ONLY: ranks and indices are at x_i=i^2+3i+5 modulo {args.prime}; "
        "freeze and audit any improvement over QQ"
    )


if __name__ == "__main__":
    main()
