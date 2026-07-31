#!/usr/bin/env python3
"""Simultaneously reconstruct a normalized Ore operator over QQ.

Coefficientwise rational reconstruction needs a modulus roughly twice the
coefficient height because it reconstructs each numerator/denominator pair
in isolation.  A fitted Ore operator instead has one common projective
denominator.  This script uses an LLL lattice

    < M e_1, ..., M e_n, (a_1,...,a_n,1) >

to recover that denominator and selected numerators simultaneously, then
recovers every remaining numerator by centered reduction and tests the
whole operator at independent holdout primes.
"""

from __future__ import annotations

import argparse
from math import gcd
import json
from pathlib import Path

from flint import fmpz_mat

from research_two_pair_sic_bidegree33_rank_two_ore_reconstruct import (
    ROOT,
    crt_merge,
)


COMMON_CACHE = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_rank_two_ore_reconstruct_images.json"
)
RELATIVE_CACHE = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_rank_two_relative_ore_images.json"
)
COMPACT_RELATIVE_CACHE = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_rank_two_compact_relative_pf_images.json"
)
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_rank_two_relative_ore_lift.json"
)


def centered(value: int, modulus: int) -> int:
    value %= modulus
    return value - modulus if 2 * value > modulus else value


def nested(
    numerators: list[int],
    denominator: int,
    order: int,
    degree: int,
) -> list[list[list[int]]]:
    width = degree + 1
    return [
        [
            [numerators[shift * width + exponent], denominator]
            for exponent in range(width)
        ]
        for shift in range(order + 1)
    ]


def nested_integers(
    numerators: list[int],
    order: int,
    degree: int,
) -> list[list[int]]:
    width = degree + 1
    return [
        numerators[shift * width : (shift + 1) * width]
        for shift in range(order + 1)
    ]


def choose_indices(
    residues: list[int],
    modulus: int,
    count: int,
    offset: int,
) -> list[int]:
    eligible = [
        index
        for index, residue in enumerate(residues[:-1])
        if residue % modulus not in (0, 1, modulus - 1)
    ]
    if len(eligible) < count:
        raise ValueError("not enough nontrivial coefficients for lattice")
    return [
        eligible[(offset + position * len(eligible) // count) % len(eligible)]
        for position in range(count)
    ]


def holdout_matches(
    numerators: list[int],
    denominator: int,
    holdouts: list[tuple[int, list[int]]],
) -> dict[str, int]:
    matches = {}
    for prime, image in holdouts:
        if denominator % prime == 0:
            matches[str(prime)] = 0
            continue
        inverse = pow(denominator, -1, prime)
        matches[str(prime)] = sum(
            numerator * inverse % prime == residue
            for numerator, residue in zip(
                numerators,
                image,
                strict=True,
            )
        )
    return matches


def candidate_from_row(
    row: list[int],
    residues: list[int],
    modulus: int,
) -> tuple[list[int], int] | None:
    denominator = row[-1]
    if denominator == 0:
        return None
    numerators = [
        centered(residue * denominator, modulus)
        for residue in residues
    ]
    content = abs(denominator)
    for numerator in numerators:
        content = gcd(content, abs(numerator))
    if content > 1:
        denominator //= content
        numerators = [
            numerator // content for numerator in numerators
        ]
    if denominator < 0:
        denominator = -denominator
        numerators = [-numerator for numerator in numerators]
    if any(
        (numerator - residue * denominator) % modulus
        for numerator, residue in zip(
            numerators,
            residues,
            strict=True,
        )
    ):
        return None
    return numerators, denominator


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--kind",
        choices=("relative", "common", "compact"),
        default="relative",
    )
    parser.add_argument("--cache", type=Path)
    parser.add_argument("--prime-count", type=int)
    parser.add_argument("--holdout-count", type=int, default=2)
    parser.add_argument(
        "--dimensions",
        type=int,
        nargs="+",
        default=(4, 8, 12, 16, 24),
    )
    parser.add_argument("--offsets", type=int, default=4)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()

    if arguments.kind == "relative":
        cache_path = arguments.cache or RELATIVE_CACHE
        order, degree = 18, 18
    elif arguments.kind == "common":
        cache_path = arguments.cache or COMMON_CACHE
        order, degree = 14, 58
    else:
        cache_path = arguments.cache or COMPACT_RELATIVE_CACHE
        order, degree = 8, 72
    payload = json.loads(cache_path.read_text())
    images = sorted(
        (
            int(prime),
            [int(value) for value in image],
        )
        for prime, image in payload["images"].items()
    )
    if arguments.prime_count is not None:
        images = images[: arguments.prime_count]
    if not 1 <= arguments.holdout_count < len(images):
        raise ValueError("invalid holdout count")
    build = images[: -arguments.holdout_count]
    holdouts = images[-arguments.holdout_count :]
    residues: list[int] = []
    modulus = 1
    for prime, image in build:
        residues, modulus = crt_merge(residues, modulus, image, prime)

    width = len(residues)
    expected_width = (order + 1) * (degree + 1)
    if width != expected_width:
        raise ValueError(
            f"image width {width} does not match {expected_width}"
        )
    attempts = []
    solution = None
    for dimension in arguments.dimensions:
        for offset in range(arguments.offsets):
            indices = choose_indices(
                residues,
                modulus,
                dimension,
                offset,
            )
            rows = []
            for position in range(dimension):
                row = [0] * (dimension + 1)
                row[position] = modulus
                rows.append(row)
            rows.append(
                [residues[index] for index in indices] + [1]
            )
            reduced = fmpz_mat(rows).lll()
            for row_index in range(dimension + 1):
                row = [
                    int(reduced[row_index, column])
                    for column in range(dimension + 1)
                ]
                candidate = candidate_from_row(row, residues, modulus)
                if candidate is None:
                    continue
                numerators, denominator = candidate
                matches = holdout_matches(
                    numerators,
                    denominator,
                    holdouts,
                )
                maximum_bits = max(
                    [abs(denominator).bit_length()]
                    + [
                        abs(numerator).bit_length()
                        for numerator in numerators
                    ]
                )
                record = {
                    "dimension": dimension,
                    "offset": offset,
                    "row": row_index,
                    "denominator_bits": abs(denominator).bit_length(),
                    "maximum_coefficient_bits": maximum_bits,
                    "holdout_matches": matches,
                }
                attempts.append(record)
                if all(value == width for value in matches.values()):
                    solution = (
                        numerators,
                        denominator,
                        record,
                        indices,
                    )
                    break
            if solution is not None:
                break
        if solution is not None:
            break

    result = {
        "format": (
            "two-pair-sic-bidegree33-rank-two-"
            "simultaneous-ore-reconstruct-v1"
        ),
        "status": (
            "stable simultaneous rational reconstruction"
            if solution is not None
            else "no stable simultaneous rational reconstruction"
        ),
        "kind": arguments.kind,
        "point": int(payload["point"]),
        "operator": {
            "order": order,
            "m_degree": degree,
            "coefficient_count": width,
        },
        "build_primes": [prime for prime, _ in build],
        "holdout_primes": [prime for prime, _ in holdouts],
        "crt_modulus_bits": modulus.bit_length(),
        "attempts": attempts,
    }
    if solution is not None:
        numerators, denominator, record, indices = solution
        result.update(
            {
                "successful_lattice": record,
                "lattice_coefficient_indices": indices,
                "common_denominator": str(denominator),
                "primitive_integer_coefficients": nested_integers(
                    numerators,
                    order,
                    degree,
                ),
                "primitive_coefficients": nested(
                    numerators,
                    denominator,
                    order,
                    degree,
                ),
                "fully_stable": True,
            }
        )
    else:
        result["fully_stable"] = False
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(result, indent=2) + "\n")
    if solution is None:
        print("PASS no stable simultaneous reconstruction at this modulus")
    else:
        print(
            "PASS simultaneous reconstruction "
            f"at dimension {solution[2]['dimension']}"
        )
        print(
            "PASS maximum primitive coefficient height "
            f"{solution[2]['maximum_coefficient_bits']} bits"
        )
        print("PASS every coefficient matches every holdout")
    print(f"PASS wrote {arguments.output}")


if __name__ == "__main__":
    main()
