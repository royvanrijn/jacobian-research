#!/usr/bin/env python3
"""Exact low-degree census for rank-stratified balanced moments.

For

    V_d = Sym^d(k^2)^* tensor Sym^d(k^2)

write X_(d,r) for the determinantal locus of coefficient matrices of rank
at most r.  This script performs two finite-prefix tests:

* it computes the generic rank of the Jacobian of the first
  q_(d,r) = r(2(d+1)-r)-3 moments on a global factor chart, modulo a
  good prime;
* it computes the diagonal-SL_2 invariant Hilbert series of k[X_(d,r)]
  from the Cauchy decomposition and tests the candidate numerator for
  moment degrees 1,...,q_(d,r).

A negative numerator coefficient proves that the consecutive moments are
not a homogeneous system of parameters, hence that their common zero
fiber contains a semistable point.  It does not prove an all-order
moment-zero point or an SIC/GVC/GMC counterexample.
"""

from __future__ import annotations

import json
from functools import lru_cache
from math import factorial
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "rank_stratified_moment_census.json"
)
PRIME = 1_000_003
HILBERT_CUTOFF = 85
Polynomial = dict[tuple[int, int], int]
Character = dict[int, int]


def multiply_polynomials(
    left: Polynomial,
    right: Polynomial,
    prime: int,
) -> Polynomial:
    result: Polynomial = {}
    for (i, j), left_value in left.items():
        for (a, b), right_value in right.items():
            exponent = (i + a, j + b)
            result[exponent] = (
                result.get(exponent, 0) + left_value * right_value
            ) % prime
    return {exponent: value for exponent, value in result.items() if value}


def matrix_rank_mod(matrix: list[list[int]], prime: int) -> int:
    if not matrix:
        return 0
    work = [[value % prime for value in row] for row in matrix]
    row = 0
    for column in range(len(work[0])):
        pivot = next(
            (
                index
                for index in range(row, len(work))
                if work[index][column]
            ),
            None,
        )
        if pivot is None:
            continue
        work[row], work[pivot] = work[pivot], work[row]
        inverse = pow(work[row][column], prime - 2, prime)
        work[row] = [(value * inverse) % prime for value in work[row]]
        for index in range(len(work)):
            if index == row or not work[index][column]:
                continue
            scale = work[index][column]
            work[index] = [
                (left - scale * right) % prime
                for left, right in zip(
                    work[index],
                    work[row],
                    strict=True,
                )
            ]
        row += 1
        if row == len(work):
            break
    return row


def factor_chart(degree: int, rank: int, seed: int) -> tuple[
    list[list[int]],
    list[list[int]],
]:
    size = degree + 1
    state = seed

    def next_value() -> int:
        nonlocal state
        state = (48_271 * state) % 2_147_483_647
        return 1 + state % (PRIME - 1)

    u = [
        [
            (
                1
                if i == inner
                else 0
                if i < rank
                else next_value()
            )
            for inner in range(rank)
        ]
        for i in range(size)
    ]
    w = [
        [
            next_value()
            for j in range(size)
        ]
        for inner in range(rank)
    ]
    return u, w


def moment_jacobian_rank(
    degree: int,
    rank: int,
    orders: tuple[int, ...],
) -> tuple[int, int]:
    size = degree + 1
    chart_dimension = rank * (2 * size - rank)
    best_rank = 0
    best_seed = 0
    for seed in range(1, 9):
        u, w = factor_chart(degree, rank, seed)
        c = [
            [
                sum(
                    u[i][inner] * w[inner][j]
                    for inner in range(rank)
                )
                % PRIME
                for j in range(size)
            ]
            for i in range(size)
        ]
        coefficient_polynomial = {
            (i, j): c[i][j]
            for i in range(size)
            for j in range(size)
            if c[i][j]
        }
        previous_power: Polynomial = {(0, 0): 1}
        rows: list[list[int]] = []
        for order in range(1, max(orders) + 1):
            total_degree = degree * order
            gradient = [[0 for _ in range(size)] for _ in range(size)]
            for a in range(size):
                for b in range(size):
                    gradient[a][b] = (
                        order
                        * sum(
                            factorial(total_degree - index)
                            * factorial(index)
                            * previous_power.get(
                                (index - a, index - b),
                                0,
                            )
                            for index in range(total_degree + 1)
                        )
                    ) % PRIME
            if order in orders:
                row: list[int] = []
                for i in range(rank, size):
                    for inner in range(rank):
                        row.append(
                            sum(
                                gradient[i][j] * w[inner][j]
                                for j in range(size)
                            )
                            % PRIME
                        )
                for inner in range(rank):
                    for j in range(size):
                        row.append(
                            sum(
                                u[i][inner] * gradient[i][j]
                                for i in range(size)
                            )
                            % PRIME
                        )
                assert len(row) == chart_dimension
                rows.append(row)
            previous_power = multiply_polynomials(
                previous_power,
                coefficient_polynomial,
                PRIME,
            )
        current_rank = matrix_rank_mod(rows, PRIME)
        if current_rank > best_rank:
            best_rank = current_rank
            best_seed = seed
        if current_rank == len(orders):
            break
    return best_rank, best_seed


def multiply_by_binomial(
    polynomial: list[int],
    exponent: int,
) -> list[int]:
    result = polynomial + [0 for _ in range(exponent)]
    for index, coefficient in enumerate(polynomial):
        result[index + exponent] -= coefficient
    return result


def exact_polynomial_quotient(
    numerator: list[int],
    denominator: list[int],
) -> list[int]:
    quotient_degree = len(numerator) - len(denominator)
    assert quotient_degree >= 0
    assert denominator[0] == 1
    quotient = [0 for _ in range(quotient_degree + 1)]
    for degree in range(quotient_degree + 1):
        quotient[degree] = numerator[degree] - sum(
            denominator[index] * quotient[degree - index]
            for index in range(1, min(degree, len(denominator) - 1) + 1)
        )
    product = [0 for _ in range(len(numerator))]
    for i, left in enumerate(quotient):
        for j, right in enumerate(denominator):
            product[i + j] += left * right
    assert product == numerator
    return quotient


def principal_schur_character(
    degree: int,
    partition: tuple[int, ...],
) -> Character:
    """Principal-specialization form of S_partition(Sym^degree)."""
    size = degree + 1
    padded = partition + (0,) * (size - len(partition))
    numerator = [1]
    denominator = [1]
    for i in range(size):
        for j in range(i + 1, size):
            numerator = multiply_by_binomial(
                numerator,
                padded[i] - padded[j] + j - i,
            )
            denominator = multiply_by_binomial(
                denominator,
                j - i,
            )
    quotient = exact_polynomial_quotient(numerator, denominator)
    t_shift = sum(i * padded[i] for i in range(size))
    q_shift = -degree * sum(partition)
    return {
        q_shift + 2 * (t_shift + exponent): multiplicity
        for exponent, multiplicity in enumerate(quotient)
        if multiplicity
    }


def partitions_with_length_at_most(
    total: int,
    maximum_length: int,
) -> list[tuple[int, ...]]:
    result: list[tuple[int, ...]] = []

    def visit(remaining: int, ceiling: int, parts: list[int]) -> None:
        if not remaining:
            result.append(tuple(parts))
            return
        if len(parts) == maximum_length:
            return
        for part in range(min(ceiling, remaining), 0, -1):
            visit(remaining - part, part, parts + [part])

    visit(total, total, [])
    return result


def invariant_multiplicity(character: Character) -> int:
    weight_zero = sum(
        multiplicity * character.get(-weight, 0)
        for weight, multiplicity in character.items()
    )
    weight_two = sum(
        multiplicity * character.get(2 - weight, 0)
        for weight, multiplicity in character.items()
    )
    return weight_zero - weight_two


def hilbert_numerator(
    hilbert: list[int],
    parameter_degrees: tuple[int, ...],
) -> list[int]:
    result = list(hilbert)
    for degree in parameter_degrees:
        for index in range(len(result) - 1, degree - 1, -1):
            result[index] -= result[index - degree]
    return result


def hilbert_census(
    degree: int,
    maximum_rank: int,
) -> dict[int, list[int]]:
    coefficients = {
        rank: [0 for _ in range(HILBERT_CUTOFF + 1)]
        for rank in range(1, maximum_rank + 1)
    }

    @lru_cache(maxsize=None)
    def cached_schur(partition: tuple[int, ...]) -> tuple[
        tuple[int, int],
        ...,
    ]:
        return tuple(
            sorted(principal_schur_character(degree, partition).items())
        )

    for order in range(HILBERT_CUTOFF + 1):
        for partition in partitions_with_length_at_most(
            order,
            maximum_rank,
        ):
            character = dict(cached_schur(partition))
            contribution = invariant_multiplicity(character)
            for rank in range(
                max(1, len(partition)),
                maximum_rank + 1,
            ):
                coefficients[rank][order] += contribution
    return coefficients


def main() -> None:
    cases: dict[str, object] = {}
    for degree in (2, 3, 4):
        size = degree + 1
        maximum_rank = min(size, 4)
        hilbert_by_rank = hilbert_census(degree, maximum_rank)
        for rank in range(1, maximum_rank + 1):
            quotient_dimension = rank * (2 * size - rank) - 3
            if quotient_dimension <= 0:
                continue
            orders = tuple(range(1, quotient_dimension + 1))
            jacobian_rank, seed = moment_jacobian_rank(
                degree,
                rank,
                orders,
            )
            numerator = hilbert_numerator(
                hilbert_by_rank[rank],
                orders,
            )
            negative = next(
                (
                    (index, value)
                    for index, value in enumerate(numerator)
                    if value < 0
                ),
                None,
            )
            least_single_replacement = None
            if negative is not None:
                for replacement in range(
                    quotient_dimension + 1,
                    quotient_dimension + 6,
                ):
                    repaired_orders = (
                        tuple(range(1, quotient_dimension))
                        + (replacement,)
                    )
                    repaired_numerator = hilbert_numerator(
                        hilbert_by_rank[rank],
                        repaired_orders,
                    )
                    if any(value < 0 for value in repaired_numerator):
                        continue
                    repaired_rank, repaired_seed = moment_jacobian_rank(
                        degree,
                        rank,
                        repaired_orders,
                    )
                    if repaired_rank != quotient_dimension:
                        continue
                    predicted_top_degree = (
                        sum(repaired_orders) - size * rank
                    )
                    last_nonzero_degree = max(
                        index
                        for index, value in enumerate(repaired_numerator)
                        if value
                    )
                    assert last_nonzero_degree == predicted_top_degree
                    assert all(
                        repaired_numerator[index]
                        == repaired_numerator[
                            predicted_top_degree - index
                        ]
                        for index in range(predicted_top_degree + 1)
                    )
                    assert all(
                        value == 0
                        for value in repaired_numerator[
                            predicted_top_degree + 1:
                        ]
                    )
                    least_single_replacement = {
                        "orders": list(repaired_orders),
                        "jacobian_rank_mod_prime": repaired_rank,
                        "jacobian_point_seed": repaired_seed,
                        "numerator": {
                            "nonnegative_through": HILBERT_CUTOFF,
                            "last_nonzero_degree": last_nonzero_degree,
                            "predicted_gorenstein_top_degree": (
                                predicted_top_degree
                            ),
                            "palindromic": True,
                        },
                    }
                    break
            key = f"d{degree}_r{rank}"
            cases[key] = {
                "balanced_degree": degree,
                "rank_bound": rank,
                "matrix_size": size,
                "determinantal_dimension": (
                    rank * (2 * size - rank)
                ),
                "quotient_dimension": quotient_dimension,
                "consecutive_moment_orders": list(orders),
                "jacobian_rank_mod_prime": jacobian_rank,
                "jacobian_point_seed": seed,
                "first_negative_numerator_coefficient": (
                    list(negative) if negative is not None else None
                ),
                "least_hilbert_compatible_single_replacement": (
                    least_single_replacement
                ),
                "hilbert_cutoff": HILBERT_CUTOFF,
            }

    # Regressions against the two previously certified ambient calculations.
    assert cases["d4_r2"]["jacobian_rank_mod_prime"] == 13
    assert (
        cases["d4_r2"]["first_negative_numerator_coefficient"]
        == [69, -5266]
    )
    assert (
        cases["d3_r4"]["first_negative_numerator_coefficient"]
        == [63, -2186]
    )

    artifact = {
        "format": "rank-stratified-balanced-moment-census-v1",
        "field": "characteristic zero",
        "good_prime_for_jacobians": PRIME,
        "scope": {
            "balanced_degrees": [2, 3, 4],
            "rank_bounds": "all determinantal ranks in each degree",
            "hilbert_cutoff": HILBERT_CUTOFF,
        },
        "interpretation": {
            "full_jacobian_rank": (
                "proves algebraic independence of the displayed moments "
                "in characteristic zero"
            ),
            "negative_numerator": (
                "proves the consecutive moments are not a homogeneous "
                "system of parameters and their zero fiber has a "
                "semistable point"
            ),
            "warning": (
                "no finite-prefix semistable point is asserted to satisfy "
                "all moments; this census does not classify counterexamples"
            ),
        },
        "cases": cases,
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")
    print("PASS rank-stratified moment Jacobians in degrees two through four")
    print(
        "PASS invariant Hilbert numerators through degree "
        f"{HILBERT_CUTOFF}"
    )
    print("PASS recovered the certified d=4 rank-two and d=3 ambient obstructions")
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
