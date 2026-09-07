#!/usr/bin/env python3
"""Expose the rank term that must be changed to improve cotangent HN lifts.

For a cubic map ``H`` put ``N=JH`` and
``A=sum_i y_i Hess(H_i)``.  The cotangent Hessian is

    M = [ A  N^t ]
        [ N   0  ].

Over the fraction field, row/column reduction of ``N`` gives

    rank(M) = 2 rank(N) + rank(K^t A K),

where the columns of ``K`` span ``ker(N)``.  The last summand is called the
cotangent kernel excess here.  This script independently checks the identity
at several deterministic good-prime specializations for every certified
cubic source on the restricted-minima frontier.

The specializations are diagnostics for the excess.  Exact generic ranks
remain those certified in the source artifacts.
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

from audit_fixed_rank_hessian_witness import (
    PRIME,
    cotangent_hessian,
    decode_h,
    deterministic_point,
    evaluate_matrix,
    modular_rank,
)


ROOT = Path(__file__).resolve().parents[1]
ARTIFACTS = ROOT / "artifacts" / "generated-results"
OUTPUT = ARTIFACTS / "cotangent_kernel_excess_frontier.json"
SOURCES = (
    ("ambient-21", "low_complexity_bcw_21_counterexample.json"),
    ("index-18", "index_reduced_bcw_22_counterexample.json"),
    ("rank-17", "rank_reduced_bcw_24_counterexample.json"),
    ("hessian-rank-37", "hessian_rank_reduced_bcw_22_counterexample.json"),
)
SEEDS = (20_260_723, 20_260_730, 20_260_737)


def modular_nullspace(
    matrix: list[list[int]], prime: int
) -> list[list[int]]:
    """Return a matrix whose columns form a right-kernel basis."""

    rows = len(matrix)
    columns = len(matrix[0]) if matrix else 0
    work = [[value % prime for value in row] for row in matrix]
    pivots: list[int] = []
    pivot_row = 0
    for column in range(columns):
        pivot = next(
            (row for row in range(pivot_row, rows) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        inverse = pow(work[pivot_row][column], -1, prime)
        work[pivot_row] = [
            value * inverse % prime for value in work[pivot_row]
        ]
        for row in range(rows):
            if row != pivot_row and work[row][column]:
                factor = work[row][column]
                work[row] = [
                    (left - factor * right) % prime
                    for left, right in zip(work[row], work[pivot_row])
                ]
        pivots.append(column)
        pivot_row += 1
        if pivot_row == rows:
            break

    free = [column for column in range(columns) if column not in pivots]
    basis_columns: list[list[int]] = []
    for free_column in free:
        vector = [0] * columns
        vector[free_column] = 1
        for row, pivot_column in enumerate(pivots):
            vector[pivot_column] = -work[row][free_column] % prime
        basis_columns.append(vector)
    return [
        [basis_columns[column][row] for column in range(len(basis_columns))]
        for row in range(columns)
    ]


def transpose(matrix: list[list[int]]) -> list[list[int]]:
    return [list(column) for column in zip(*matrix)] if matrix else []


def multiply(
    left: list[list[int]], right: list[list[int]], prime: int
) -> list[list[int]]:
    if not left:
        return []
    columns = transpose(right)
    return [
        [
            sum(a * b for a, b in zip(row, column)) % prime
            for column in columns
        ]
        for row in left
    ]


def source_profile(label: str, filename: str) -> dict[str, object]:
    stored = json.loads((ARTIFACTS / filename).read_text())
    h = decode_h(stored)
    dimension = len(h)
    jacobian, upper_left, hessian = cotangent_hessian(h)
    samples = []
    for seed in SEEDS:
        point = deterministic_point(2 * dimension, PRIME, seed)
        n_matrix = evaluate_matrix(jacobian, point, PRIME)
        a_matrix = evaluate_matrix(upper_left, point, PRIME)
        hessian_matrix = evaluate_matrix(hessian, point, PRIME)
        kernel = modular_nullspace(n_matrix, PRIME)
        assert modular_rank(kernel, PRIME) == dimension - modular_rank(
            n_matrix, PRIME
        )
        restriction = multiply(
            multiply(transpose(kernel), a_matrix, PRIME),
            kernel,
            PRIME,
        )
        rank_n = modular_rank(n_matrix, PRIME)
        excess = modular_rank(restriction, PRIME)
        rank_hessian = modular_rank(hessian_matrix, PRIME)
        assert rank_hessian == 2 * rank_n + excess
        samples.append(
            {
                "seed": seed,
                "rank_JH": rank_n,
                "kernel_dimension": dimension - rank_n,
                "kernel_excess": excess,
                "cotangent_hessian_rank": rank_hessian,
                "identity_rhs": 2 * rank_n + excess,
            }
        )

    signatures = {
        (
            sample["rank_JH"],
            sample["kernel_dimension"],
            sample["kernel_excess"],
            sample["cotangent_hessian_rank"],
        )
        for sample in samples
    }
    assert len(signatures) == 1
    rank_n, kernel_dimension, excess, rank_hessian = signatures.pop()
    statistics = stored["statistics"]
    exact_rank = int(statistics["generic_rank_JH_over_QQ_x"])
    exact_hessian = statistics.get("generic_cotangent_hessian_rank_over_QQ_xy")
    assert rank_n == exact_rank
    if exact_hessian is not None:
        assert rank_hessian == int(exact_hessian)
    return {
        "label": label,
        "artifact": filename,
        "dimension": dimension,
        "certified_generic_rank_JH": exact_rank,
        "sampled_kernel_dimension": kernel_dimension,
        "sampled_kernel_excess": excess,
        "sampled_cotangent_hessian_rank": rank_hessian,
        "samples": samples,
    }


def main() -> None:
    profiles = [
        source_profile(label, filename) for label, filename in SOURCES
    ]
    by_label = {profile["label"]: profile for profile in profiles}
    assert by_label["ambient-21"]["sampled_kernel_excess"] == 2
    assert by_label["index-18"]["sampled_kernel_excess"] == 2
    assert by_label["rank-17"]["sampled_kernel_excess"] == 4
    assert by_label["hessian-rank-37"]["sampled_kernel_excess"] == 1
    rank_37 = by_label["hessian-rank-37"]
    assert rank_37["certified_generic_rank_JH"] == 18
    rank_37_stored = json.loads(
        (
            ARTIFACTS
            / "hessian_rank_reduced_bcw_22_counterexample.json"
        ).read_text()
    )
    rank_37_exact_hessian_rank = int(
        rank_37_stored["statistics"][
            "generic_cotangent_hessian_rank_over_QQ_xy"
        ]
    )
    exact_rank_37_excess = (
        rank_37_exact_hessian_rank
        - 2 * int(rank_37["certified_generic_rank_JH"])
    )
    assert exact_rank_37_excess == 1

    combined_search = json.loads(
        (
            ARTIFACTS
            / "restricted_bcw_circuit_search_all_w64.json"
        ).read_text()
    )
    excess_histogram: Counter[int] = Counter()
    terminal_count = 0
    for row in combined_search["terminal_histogram"]:
        excess = (
            int(row["quartic_hessian_rank"])
            - 2 * int(row["cubic_rank"])
        )
        excess_histogram[excess] += int(row["count"])
        terminal_count += int(row["count"])
    assert terminal_count == 140
    assert dict(sorted(excess_histogram.items())) == {
        1: 1,
        2: 8,
        3: 17,
        4: 31,
        5: 41,
        6: 23,
        7: 15,
        8: 4,
    }

    payload = {
        "format": "cotangent-kernel-excess-frontier-v1",
        "field_for_identity": "arbitrary field of characteristic not two",
        "identity": (
            "rank([[A,N^t],[N,0]]) = 2*rank(N) + "
            "rank(K^t*A*K), columns(K)=ker(N)"
        ),
        "proof": (
            "choose invertible row and column bases reducing N to "
            "diag(I_r,0); congruent block elimination splits off two "
            "rank-r hyperbolic blocks and leaves precisely the restriction "
            "of A to ker(N)"
        ),
        "specialization_prime": PRIME,
        "specialization_status": (
            "the block identity is exact; reported kernel-excess values are "
            "stable diagnostics at three good-prime points"
        ),
        "profiles": profiles,
        "exact_rank_37_excess": {
            "generic_rank_JH": 18,
            "generic_cotangent_hessian_rank": 37,
            "kernel_excess": exact_rank_37_excess,
            "exactness_reason": (
                "the two characteristic-zero generic ranks are independently "
                "certified, so their difference is exact"
            ),
        },
        "combined_nine_atom_search": {
            "artifact": "restricted_bcw_circuit_search_all_w64.json",
            "terminal_count": terminal_count,
            "sampled_excess_histogram": {
                str(excess): count
                for excess, count in sorted(excess_histogram.items())
            },
            "unique_excess_one_terminal": "qb+x2s",
            "excess_zero_terminals": 0,
            "status": "bounded finite-field diagnostic, not a lower bound",
        },
        "strict_improvement_targets": {
            "current_cotangent_hessian_rank": 37,
            "fixed_rank_18": (
                "kernel excess 0 gives Hessian rank 36; no larger "
                "improvement is possible without lowering rank(JH)"
            ),
            "rank_17": (
                "kernel excess at most 2 gives Hessian rank at most 36; "
                "the current rank-17 source has excess 4"
            ),
            "cotangent_dimension": (
                "beating quartic dimension 42 requires a cubic source in "
                "dimension at most 20"
            ),
        },
        "rank_index_coupling": {
            "identity": (
                "a nilpotent matrix of rank r has nilpotency index at most r+1"
            ),
            "proof": (
                "a Jordan block of length nu already contributes nu-1 to "
                "the rank"
            ),
            "frontier_consequence": (
                "a positive theorem for all cubic-homogeneous maps of index "
                "at most d implies nu_cub>=d+1 and r_cub>=d"
            ),
            "current_saturation": (
                "the rank-17 witness has index 18, so it saturates "
                "index=rank+1"
            ),
        },
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print("PASS cotangent excess: verified the block-rank identity at 12 points")
    print("PASS cotangent excess: rank-37 witness has exact generic excess 1")
    print("PASS cotangent excess: 140-terminal search has no sampled excess-zero map")
    print("TARGET cotangent excess: rank 18/excess 0 or rank 17/excess <=2")
    print(f"PASS cotangent excess: wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
