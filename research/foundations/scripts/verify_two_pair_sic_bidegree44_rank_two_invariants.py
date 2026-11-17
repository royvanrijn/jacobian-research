#!/usr/bin/env python3
"""Exact invariant-quotient tests on the rank-at-most-two V_4 locus.

For X_2={C in Mat_5: rank(C)<=2}, this checker proves:

* moments mu_1,...,mu_13 have Jacobian rank 13 on an exact rank-two chart;
* mu_1,...,mu_12,mu_14 also have Jacobian rank 13 there;
* the SL_2-invariant Hilbert series of Q[X_2] makes degrees 1,...,13
  impossible for a homogeneous system of parameters, with coefficient
  -5266 in degree 69;
* degrees 1,...,12,14 pass the same necessary test through degree 100,
  and the candidate numerator through its predicted top degree 82 is
  nonnegative and palindromic.

The Hilbert coefficients use the Cauchy decomposition

  Q[X_2]_n = direct_sum_(lambda partition n, length(lambda)<=2)
             S_lambda(Sym^4)^* tensor S_lambda(Sym^4)

and the SL_2 rule dim(M^SL2)=mult_M(0)-mult_M(2).
"""

from __future__ import annotations

import json
from math import factorial
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree44_rank_two_invariants.json"
)
PRIME = 1_000_003
HILBERT_CUTOFF = 100
WeightCharacter = dict[int, int]
Polynomial = dict[tuple[int, int], int]


def multiply_mod(
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


def moment_jacobian_rows(
    orders: range,
    prime: int,
) -> tuple[dict[int, list[int]], list[list[int]]]:
    # Global rank-two chart: pivot rows of U are the 2-by-2 identity.
    u = [
        [1, 0],
        [0, 1],
        [14, 17],
        [18, 4],
        [6, 13],
    ]
    w = [
        [8, 10, 1, 8, 4],
        [19, 1, 4, 6, 17],
    ]
    c = [
        [
            sum(u[i][inner] * w[inner][j] for inner in range(2)) % prime
            for j in range(5)
        ]
        for i in range(5)
    ]
    assert c == [
        [8, 10, 1, 8, 4],
        [19, 1, 4, 6, 17],
        [435, 157, 82, 214, 345],
        [220, 184, 34, 168, 140],
        [295, 73, 58, 126, 245],
    ]
    coefficient_polynomial = {
        (i, j): c[i][j]
        for i in range(5)
        for j in range(5)
        if c[i][j]
    }
    previous_power: Polynomial = {(0, 0): 1}
    rows: dict[int, list[int]] = {}
    nonpivot_rows = (2, 3, 4)
    for order in orders:
        degree = 4 * order
        gradient = [[0 for _ in range(5)] for _ in range(5)]
        for a in range(5):
            for b in range(5):
                gradient[a][b] = (
                    order
                    * sum(
                        factorial(degree - index)
                        * factorial(index)
                        * previous_power.get((index - a, index - b), 0)
                        for index in range(degree + 1)
                    )
                ) % prime

        # Chart variables are the six nonpivot entries of U followed by
        # all ten entries of W.  Apply the chain rule to C=UW.
        row: list[int] = []
        for i in nonpivot_rows:
            for inner in range(2):
                row.append(
                    sum(
                        gradient[i][j] * w[inner][j]
                        for j in range(5)
                    )
                    % prime
                )
        for inner in range(2):
            for j in range(5):
                row.append(
                    sum(
                        u[i][inner] * gradient[i][j]
                        for i in range(5)
                    )
                    % prime
                )
        assert len(row) == 16
        rows[order] = row
        previous_power = multiply_mod(
            previous_power,
            coefficient_polynomial,
            prime,
        )
    return rows, c


def rank_one_jacobian_rows(prime: int) -> list[list[int]]:
    # Gauge u_0=1 on C=u*w^T.  The chart has four remaining u variables
    # and five w variables.
    u = [1, 15, 12, 9, 5]
    w = [6, 1, 11, 17, 15]
    c = [[u[i] * w[j] % prime for j in range(5)] for i in range(5)]
    coefficient_polynomial = {
        (i, j): c[i][j] for i in range(5) for j in range(5)
    }
    previous_power: Polynomial = {(0, 0): 1}
    rows: list[list[int]] = []
    for order in range(1, 7):
        degree = 4 * order
        gradient = [[0 for _ in range(5)] for _ in range(5)]
        for a in range(5):
            for b in range(5):
                gradient[a][b] = (
                    order
                    * sum(
                        factorial(degree - index)
                        * factorial(index)
                        * previous_power.get((index - a, index - b), 0)
                        for index in range(degree + 1)
                    )
                ) % prime
        row = [
            sum(gradient[i][j] * w[j] for j in range(5)) % prime
            for i in range(1, 5)
        ]
        row.extend(
            sum(u[i] * gradient[i][j] for i in range(5)) % prime
            for j in range(5)
        )
        rows.append(row)
        previous_power = multiply_mod(
            previous_power,
            coefficient_polynomial,
            prime,
        )
    return rows


def add_characters(
    left: WeightCharacter,
    right: WeightCharacter,
    right_scale: int = 1,
) -> WeightCharacter:
    result = dict(left)
    for weight, multiplicity in right.items():
        result[weight] = (
            result.get(weight, 0) + right_scale * multiplicity
        )
        if not result[weight]:
            del result[weight]
    return result


def multiply_characters(
    left: WeightCharacter,
    right: WeightCharacter,
) -> WeightCharacter:
    result: WeightCharacter = {}
    for left_weight, left_multiplicity in left.items():
        for right_weight, right_multiplicity in right.items():
            weight = left_weight + right_weight
            result[weight] = (
                result.get(weight, 0)
                + left_multiplicity * right_multiplicity
            )
    return result


def complete_symmetric_characters(cutoff: int) -> list[WeightCharacter]:
    # Coefficients of product_(w=-4,-2,0,2,4) (1-t*q^w)^(-1).
    complete: list[WeightCharacter] = [{} for _ in range(cutoff + 2)]
    complete[0] = {0: 1}
    for weight in (-4, -2, 0, 2, 4):
        updated: list[WeightCharacter] = [{} for _ in range(cutoff + 2)]
        updated[0] = dict(complete[0])
        for degree in range(1, cutoff + 2):
            shifted = {
                exponent + weight: multiplicity
                for exponent, multiplicity in updated[degree - 1].items()
            }
            updated[degree] = add_characters(complete[degree], shifted)
        complete = updated
    return complete


def schur_two_row_character(
    first_row: int,
    second_row: int,
    complete: list[WeightCharacter],
) -> WeightCharacter:
    # Jacobi--Trudi: s_(a,b)=h_a*h_b-h_(a+1)*h_(b-1).
    result = multiply_characters(
        complete[first_row],
        complete[second_row],
    )
    if second_row:
        result = add_characters(
            result,
            multiply_characters(
                complete[first_row + 1],
                complete[second_row - 1],
            ),
            -1,
        )
    return result


def invariant_multiplicity(tensor_factor: WeightCharacter) -> int:
    # The character of S_lambda(V)^* tensor S_lambda(V) is the square.
    weight_zero = 0
    weight_two = 0
    for weight, multiplicity in tensor_factor.items():
        weight_zero += multiplicity * tensor_factor.get(-weight, 0)
        weight_two += multiplicity * tensor_factor.get(2 - weight, 0)
    return weight_zero - weight_two


def invariant_hilbert_coefficients(cutoff: int) -> list[int]:
    complete = complete_symmetric_characters(cutoff)
    coefficients: list[int] = []
    for degree in range(cutoff + 1):
        value = 0
        for second_row in range(degree // 2 + 1):
            first_row = degree - second_row
            value += invariant_multiplicity(
                schur_two_row_character(
                    first_row,
                    second_row,
                    complete,
                )
            )
        coefficients.append(value)
    return coefficients


def hilbert_numerator(
    hilbert: list[int],
    degrees: tuple[int, ...],
) -> list[int]:
    result = list(hilbert)
    for degree in degrees:
        for index in range(len(result) - 1, degree - 1, -1):
            result[index] -= result[index - degree]
    return result


def main() -> None:
    rows, point = moment_jacobian_rows(range(1, 15), PRIME)
    consecutive = tuple(range(1, 14))
    corrected = tuple(range(1, 13)) + (14,)
    consecutive_rank = matrix_rank_mod(
        [rows[order] for order in consecutive],
        PRIME,
    )
    corrected_rank = matrix_rank_mod(
        [rows[order] for order in corrected],
        PRIME,
    )
    assert consecutive_rank == corrected_rank == 13
    rank_one_moment_rank = matrix_rank_mod(
        rank_one_jacobian_rows(PRIME),
        PRIME,
    )
    assert rank_one_moment_rank == 6

    hilbert = invariant_hilbert_coefficients(HILBERT_CUTOFF)
    assert hilbert[:14] == [
        1,
        1,
        5,
        13,
        53,
        149,
        483,
        1274,
        3370,
        7994,
        18398,
        39472,
        81962,
        161896,
    ]
    consecutive_numerator = hilbert_numerator(hilbert, consecutive)
    assert consecutive_numerator[69] == -5266
    assert all(
        value >= 0
        for value in consecutive_numerator[:69]
    )

    corrected_numerator = hilbert_numerator(hilbert, corrected)
    assert all(value >= 0 for value in corrected_numerator)
    assert max(
        index
        for index, value in enumerate(corrected_numerator)
        if value
    ) == 82
    # The square 5-by-5 rank-at-most-two determinantal ring is Gorenstein
    # with a-invariant -5*2=-10.  Since SL_2 has no nontrivial character,
    # the invariant ring has the same a-invariant.  Thus an hsop of total
    # degree 1+...+12+14=92 would have numerator degree 92-10=82.
    corrected_predicted_top_degree = sum(corrected) - 10
    assert corrected_predicted_top_degree == 82
    assert all(
        corrected_numerator[index]
        == corrected_numerator[corrected_predicted_top_degree - index]
        for index in range(corrected_predicted_top_degree + 1)
    )

    complete = complete_symmetric_characters(HILBERT_CUTOFF)
    rank_one_hilbert = [
        invariant_multiplicity(complete[degree])
        for degree in range(HILBERT_CUTOFF + 1)
    ]
    assert rank_one_hilbert[:14] == [
        1,
        1,
        3,
        5,
        12,
        18,
        36,
        54,
        91,
        133,
        205,
        285,
        414,
        558,
    ]
    rank_one_numerator = hilbert_numerator(
        rank_one_hilbert,
        tuple(range(1, 7)),
    )
    assert all(value >= 0 for value in rank_one_numerator)
    assert [
        (index, value)
        for index, value in enumerate(rank_one_numerator)
        if value
    ] == [
        (0, 1),
        (2, 1),
        (3, 1),
        (4, 4),
        (5, 2),
        (6, 7),
        (7, 5),
        (8, 8),
        (9, 5),
        (10, 7),
        (11, 2),
        (12, 4),
        (13, 1),
        (14, 1),
        (16, 1),
    ]

    artifact = {
        "format": "two-pair-sic-bidegree44-rank-two-invariants-v1",
        "field": "characteristic zero",
        "determinantal_locus": "rank(C)<=2 in Mat_5",
        "affine_dimension": 16,
        "generic_invariant_quotient_dimension": 13,
        "jacobian_certificate": {
            "prime": PRIME,
            "rank_two_point": point,
            "degrees_1_through_13_rank": consecutive_rank,
            "degrees_1_through_12_and_14_rank": corrected_rank,
        },
        "hilbert_series": {
            "method": (
                "Cauchy decomposition for the rank-two determinantal "
                "ring and SL2 weight-zero-minus-weight-two multiplicity"
            ),
            "checked_through_degree": HILBERT_CUTOFF,
            "initial_coefficients_0_through_13": hilbert[:14],
            "degrees_1_through_13": {
                "first_negative_numerator_coefficient": {
                    "degree": 69,
                    "value": consecutive_numerator[69],
                },
                "hsop_status": "excluded exactly",
            },
            "degrees_1_through_12_and_14": {
                "negative_coefficients_through_cutoff": [],
                "last_nonzero_degree_through_cutoff": 82,
                "determinantal_a_invariant": -10,
                "predicted_hsop_numerator_degree": (
                    corrected_predicted_top_degree
                ),
                "palindromic_through_predicted_top_degree": True,
                "candidate_numerator_coefficient_sum": sum(
                    corrected_numerator[
                        : corrected_predicted_top_degree + 1
                    ]
                ),
                "hsop_status": "necessary Hilbert test only",
            },
        },
        "rank_one_boundary": {
            "affine_dimension": 9,
            "generic_invariant_quotient_dimension": 6,
            "moments_1_through_6_jacobian_rank": rank_one_moment_rank,
            "degrees_1_through_6_hilbert_numerator": [
                [index, value]
                for index, value in enumerate(rank_one_numerator)
                if value
            ],
            "numerator_coefficient_sum": sum(rank_one_numerator),
            "hsop_status": "necessary Hilbert test only",
            "remaining_geometry": (
                "finitely many exceptional squarefree quartic "
                "cross-ratios and the uniform low-root finite cutoff"
            ),
            "hilbert_mumford_observation": (
                "rank-one annihilator tensors A tensor ell^4 with "
                "A(ell)=0 are unstable: the ell-adapted one-parameter "
                "subgroup gives every surviving tensor coefficient "
                "strictly positive weight"
            ),
        },
        "consequence": (
            "the first thirteen moments have a semistable common zero "
            "on rank(C)<=2; its exact rank may be one or two"
        ),
        "written_source": (
            "extended-geometry/TWO_PAIR_SIC_BIDEGREE44_RANK_FRONTIER.md"
        ),
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")
    print("PASS rank<=2: moments 1,...,13 have exact Jacobian rank 13")
    print("PASS rank<=2: degrees 1,...,13 fail the Hilbert test at degree 69")
    print("PASS rank<=2: moments 1,...,12,14 have exact Jacobian rank 13")
    print(
        "PASS rank<=2: corrected numerator is nonnegative through 100 "
        "and palindromic through predicted degree 82"
    )
    print("PASS rank=1 boundary: moments 1,...,6 have Jacobian rank 6")
    print("PASS rank=1 boundary: degrees 1,...,6 pass the Hilbert test through 100")
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
