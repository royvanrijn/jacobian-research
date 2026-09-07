#!/usr/bin/env python3
"""Dependency-free exact checks for the uniform two-pair linear pivot.

The proof in TWO_PAIR_OPPOSITE_MONOMIAL_OBSTRUCTION.md reduces the pivot
matrix to a filtered endpoint calculation over Z_(2).  This script
independently reconstructs its entries from the terminating binomial
formulas and checks the predicted local Smith exponents in a finite range.
The finite range is a regression check, not the all-degree proof.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from math import comb
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_linear_pivot.json"
)


def valuation_two(value: Fraction | int) -> int:
    value = Fraction(value)
    if not value:
        raise ValueError("the 2-adic valuation of zero is infinite")
    numerator = abs(value.numerator)
    denominator = value.denominator
    result = 0
    while numerator % 2 == 0:
        numerator //= 2
        result += 1
    while denominator % 2 == 0:
        denominator //= 2
        result -= 1
    return result


def integral_power(height: int, exponent: int) -> Fraction:
    """Integral of t^height (1-t^2)^exponent on [0,1]."""

    result = Fraction(1, height + 1)
    for index in range(1, exponent + 1):
        result *= Fraction(2 * index, height + 2 * index + 1)
    return result


def positive_entry(degree: int, phase: int, order: int) -> Fraction:
    """2^(order-1) times the positive-side angular linear coefficient."""

    height = degree - phase
    power = order - 1
    if power < phase:
        return Fraction(0)
    return (
        integral_power(height, power)
        * generalized_binomial(power - height - 1, power - phase)
    )


def generalized_binomial(upper: int, lower: int) -> int:
    if lower < 0:
        return 0
    result = 1
    for index in range(1, lower + 1):
        result = result * (upper - index + 1) // index
    return result


def positive_entry_sum(
    degree: int,
    phase: int,
    order: int,
) -> Fraction:
    """The unsimplified terminating sum for positive_entry."""

    height = degree - phase
    power = order - 1
    result = Fraction(0)
    for index in range(power + 1):
        target = power - phase
        coefficient = (
            comb(power + 2 * index, target)
            if 0 <= target <= power + 2 * index
            else 0
        )
        result += Fraction(
            (-1) ** index
            * comb(power, index)
            * coefficient,
            height + 2 * index + 1,
        )
    return result


def negative_entry(degree: int, phase: int, order: int) -> Fraction:
    """2^(order-1+phase) times the negative-side angular coefficient."""

    height = degree - phase
    power = order - 1
    target = power + phase
    result = Fraction(0)
    for index in range(power + 1):
        coefficient = (
            comb(power + 2 * index, target)
            if target <= power + 2 * index
            else 0
        )
        if not coefficient:
            continue
        for endpoint_index in range(phase + 1):
            result += Fraction(
                (-1) ** (index + endpoint_index)
                * comb(power, index)
                * comb(phase, endpoint_index)
                * coefficient,
                height + 2 * index + 2 * endpoint_index + 1,
            )
    return result


def pivot_matrix(degree: int) -> list[list[Fraction]]:
    phases = [
        phase
        for phase in range(1, degree + 1)
        if (degree - phase) % 2 == 0
    ]
    columns = [
        (side, phase)
        for phase in phases
        for side in ("negative", "positive")
    ]
    return [
        [
            (
                negative_entry(degree, phase, order)
                if side == "negative"
                else positive_entry(degree, phase, order)
            )
            for side, phase in columns
        ]
        for order in range(2, 2 + len(columns))
    ]


def local_smith_exponents(
    matrix: list[list[Fraction]],
) -> list[int]:
    """Diagonal valuations from elimination over the DVR Z_(2)."""

    work = [row[:] for row in matrix]
    size = len(work)
    result: list[int] = []
    for pivot_index in range(size):
        candidates = [
            (valuation_two(work[row][column]), row, column)
            for row in range(pivot_index, size)
            for column in range(pivot_index, size)
            if work[row][column]
        ]
        assert candidates
        valuation, pivot_row, pivot_column = min(candidates)
        work[pivot_index], work[pivot_row] = (
            work[pivot_row],
            work[pivot_index],
        )
        for row in work:
            row[pivot_index], row[pivot_column] = (
                row[pivot_column],
                row[pivot_index],
            )
        pivot = work[pivot_index][pivot_index]
        result.append(valuation)
        for row in range(pivot_index + 1, size):
            multiplier = work[row][pivot_index] / pivot
            assert not multiplier or valuation_two(multiplier) >= 0
            for column in range(pivot_index, size):
                work[row][column] -= (
                    multiplier * work[pivot_index][column]
                )
    return result


def phase_ordered_pivot_exponents(
    matrix: list[list[Fraction]],
) -> list[int]:
    """Valuations of the written phase-by-phase elimination."""

    work = [row[:] for row in matrix]
    result: list[int] = []
    for pivot_index in range(len(work)):
        pivot = work[pivot_index][pivot_index]
        assert pivot
        result.append(valuation_two(pivot))
        for row in range(pivot_index + 1, len(work)):
            multiplier = work[row][pivot_index] / pivot
            for column in range(pivot_index, len(work)):
                work[row][column] -= (
                    multiplier * work[pivot_index][column]
                )
    return result


def factorial_valuation(index: int) -> int:
    result = 0
    for value in range(1, index + 1):
        result += valuation_two(value)
    return result


def contact_valuation(index: int) -> int:
    return index + factorial_valuation(index)


def predicted_exponents(degree: int) -> list[int]:
    half = (degree + 1) // 2
    if degree % 2 == 0:
        return [
            contact_valuation(2 * index)
            for index in range(1, half + 1)
            for _ in range(2)
        ]
    return (
        [contact_valuation(1)]
        + [
            contact_valuation(2 * index - 1)
            for index in range(2, half + 1)
            for _ in range(2)
        ]
        + [contact_valuation(2 * half + 1)]
    )


def original_determinant_valuation(
    degree: int,
    scaled_exponents: list[int],
) -> int:
    size = 2 * ((degree + 1) // 2)
    phases = [
        phase
        for phase in range(1, degree + 1)
        if (degree - phase) % 2 == 0
    ]
    row_order_factor = factorial_valuation(size + 1)
    row_power_factor = size * (size + 1) // 2
    negative_column_factor = sum(phases)
    return (
        sum(scaled_exponents)
        + row_order_factor
        - row_power_factor
        - negative_column_factor
    )


def certificate(degree: int) -> dict[str, object]:
    phases = [
        phase
        for phase in range(1, degree + 1)
        if (degree - phase) % 2 == 0
    ]
    for phase in phases:
        for order in range(2, 2 + 2 * len(phases)):
            assert positive_entry(
                degree,
                phase,
                order,
            ) == positive_entry_sum(degree, phase, order)
    matrix = pivot_matrix(degree)
    phase_block_valuations = {}
    for phase in phases:
        row_powers = (
            (phase - 1, phase)
            if phase % 2 == 0
            else (phase, phase + 1)
        )
        block = [
            [
                negative_entry(degree, phase, power + 1),
                positive_entry(degree, phase, power + 1),
            ]
            for power in row_powers
        ]
        determinant = (
            block[0][0] * block[1][1]
            - block[0][1] * block[1][0]
        )
        expected_block_valuation = (
            2 * contact_valuation(phase)
            if phase % 2 == 0
            else contact_valuation(phase)
            + contact_valuation(phase + 2)
        )
        assert valuation_two(determinant) == expected_block_valuation
        phase_block_valuations[str(phase)] = expected_block_valuation
    exponents = local_smith_exponents(matrix)
    ordered_exponents = phase_ordered_pivot_exponents(matrix)
    expected = predicted_exponents(degree)
    assert exponents == expected
    assert ordered_exponents == expected
    return {
        "degree": degree,
        "phases": phases,
        "phase_block_determinant_valuations": phase_block_valuations,
        "scaled_local_smith_exponents": exponents,
        "phase_ordered_pivot_exponents": ordered_exponents,
        "scaled_determinant_valuation": sum(exponents),
        "original_pivot_determinant_valuation": (
            original_determinant_valuation(degree, exponents)
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-degree", type=int, default=5)
    parser.add_argument("--max-degree", type=int, default=25)
    args = parser.parse_args()
    assert 5 <= args.min_degree <= args.max_degree
    certificates = [
        certificate(degree)
        for degree in range(args.min_degree, args.max_degree + 1)
    ]
    artifact = {
        "format": "two-pair-linear-pivot-v1",
        "status": (
            "finite exact regression for the written all-degree "
            "2-adic proof"
        ),
        "degree_range": [args.min_degree, args.max_degree],
        "endpoint_factorization": (
            "2p=-(1+x)(2+x)+(1-u)x^(-1)(1+x)^3"
        ),
        "contact_valuation": "delta_s=s+v_2(s!)",
        "degrees": certificates,
        "written_source": (
            "extended-geometry/"
            "TWO_PAIR_OPPOSITE_MONOMIAL_OBSTRUCTION.md"
        ),
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")
    print(
        "PASS exact terminating entry formulas and predicted local "
        f"Smith exponents in degrees {args.min_degree}..{args.max_degree}"
    )
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
