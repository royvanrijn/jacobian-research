#!/usr/bin/env python3
"""Bounded local screen for transverse two-section seeds in the max-200 census.

The input is the frozen exact six-root panel in
``elliptic_mestre_root_tuple_scale_max200.json``.  It retains only root tuples
whose affine-normalized leading invariant is a rational square.  For each
such seed, a small-prime exhaustive scan finds nonsingular projected affine
section residues, Hensel-lifts them in the two abscissa coordinates, and
checks the reconstructed rational sections in the recursive equations.  All
pairs of reconstructed sections then receive the exact seven-by-eight
Jacobian-rank calculation.

This is a bounded seed census, not a proof that a tuple has no other affine
sections when it yields none here.  It never expands the residual equations.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from math import isqrt
from pathlib import Path
from typing import Iterable

from probe_mestre_two_section_local_continuation import (
    Field,
    rational_mod,
    rational_reconstruction,
    residuals,
    row_reduce,
    solve_square,
)


Q = Fraction
DEFAULT_INPUT = Path(
    "archive/elliptic-curves/artifacts/generated-results/elliptic_mestre_root_tuple_scale_max200.json"
)
DEFAULT_PRIMES = (7, 11, 13, 17)


def rational_text(value: Fraction) -> str:
    value = Q(value)
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def normalized_moduli(roots: Iterable[int]) -> tuple[Fraction, ...]:
    roots = tuple(map(Q, roots))
    if len(roots) != 6 or roots[0] != 0 or roots[1] == 0:
        raise ValueError("the frozen census roots must begin with distinct 0,r1")
    scale = roots[1]
    normalized = tuple(root / scale for root in roots)
    h = [Q(1)]
    for root in normalized[2:]:
        updated = [Q(0)] * (len(h) + 1)
        for degree, coefficient in enumerate(h):
            updated[degree] -= root * coefficient
            updated[degree + 1] += coefficient
        h = updated
    return tuple(reversed(h[:-1]))


def leading_invariant(moduli: tuple[Fraction, ...]) -> Fraction:
    c1, c2, c3, c4 = moduli
    a1, a2, a3, a4 = c1 - 1, c2 - c1, c3 - c2, c4 - c3
    return 5 * a1**4 - 24 * a1**2 * a2 + 32 * a1 * a3 + 16 * a2**2 - 64 * a4


def rational_square(value: Fraction) -> bool:
    return (
        value >= 0
        and isqrt(value.numerator) ** 2 == value.numerator
        and isqrt(value.denominator) ** 2 == value.denominator
    )


def invertible_row_pair(jacobian: list[list[int]], prime: int) -> tuple[int, int] | None:
    for first in range(3):
        for second in range(first + 1, 3):
            determinant = (
                jacobian[first][0] * jacobian[second][1]
                - jacobian[first][1] * jacobian[second][0]
            ) % prime
            if determinant:
                return first, second
    return None


def lift_section(
    moduli: tuple[Fraction, ...],
    residue: tuple[int, int],
    prime: int,
    precision: int,
) -> tuple[Fraction, Fraction] | None:
    """Lift a nonsingular projected section and recognize a rational point."""

    reduced_moduli = [rational_mod(value, prime) for value in moduli]
    initial = residuals((*reduced_moduli, *residue, 0, 1), Field(prime))[1:4]
    jacobian = [list(value.gradient[4:6]) for value in initial]
    rows = invertible_row_pair(jacobian, prime)
    if rows is None:
        return None
    coordinates = list(residue)
    for exponent in range(1, precision):
        modulus = prime ** (exponent + 1)
        lifted_moduli = [rational_mod(value, modulus) for value in moduli]
        values = residuals((*lifted_moduli, *coordinates, 0, 1), Field(modulus))[1:4]
        right = [(-values[index].value // (prime**exponent)) % prime for index in rows]
        correction = solve_square([jacobian[index] for index in rows], right, prime)
        coordinates = [
            (value + prime**exponent * digit) % modulus
            for value, digit in zip(coordinates, correction)
        ]
        checked = residuals((*lifted_moduli, *coordinates, 0, 1), Field(modulus))[1:4]
        if any(value.value for value in checked):
            return None
    modulus = prime**precision
    section = tuple(rational_reconstruction(value, modulus) for value in coordinates)
    if any(value is None for value in section):
        return None
    exact = tuple(Q(value) for value in section)
    values = residuals((*moduli, *exact, 0, 1), Field())[1:4]
    return exact if not any(value.value for value in values) else None


def sections_at_prime(
    moduli: tuple[Fraction, ...], prime: int, precision: int
) -> tuple[tuple[Fraction, Fraction], ...]:
    """Enumerate the nonsingular mod-p solutions and lift every candidate."""

    reduced_moduli = [rational_mod(value, prime) for value in moduli]
    field = Field(prime)
    sections = set()
    for intercept in range(prime):
        for slope in range(prime):
            values = residuals((*reduced_moduli, intercept, slope, 0, 1), field)[1:4]
            if any(value.value for value in values):
                continue
            jacobian = [list(value.gradient[4:6]) for value in values]
            if invertible_row_pair(jacobian, prime) is None:
                continue
            lifted = lift_section(moduli, (intercept, slope), prime, precision)
            if lifted is not None:
                sections.add(lifted)
    return tuple(sorted(sections))


def sections_at_seed(
    moduli: tuple[Fraction, ...], primes: Iterable[int], precision: int
) -> tuple[tuple[Fraction, Fraction], ...]:
    """Use the first good listed prime; larger primes only cost more here."""

    for prime in primes:
        try:
            return sections_at_prime(moduli, prime, precision)
        except ValueError:
            continue
    raise ValueError("every supplied prime divided a normalized moduli denominator")


def sections_across_primes(
    moduli: tuple[Fraction, ...], primes: Iterable[int], precision: int
) -> tuple[tuple[Fraction, Fraction], ...]:
    """Union exact reconstructions over every supplied prime.

    A section can be singular in the projected abscissa variables at the first
    usable prime.  This remains a bounded reconstruction screen, but avoids
    treating that one-prime singularity as a chart-level miss.
    """

    sections = set()
    for prime in primes:
        try:
            sections.update(sections_at_prime(moduli, prime, precision))
        except ValueError:
            continue
    return tuple(sorted(sections))


def screen(
    input_path: Path,
    primes: Iterable[int],
    precision: int,
    *,
    start: int = 0,
    count: int | None = None,
    all_primes: bool = False,
) -> dict[str, object]:
    source = json.loads(input_path.read_text())
    records = source["complete_panel_screen"]["family_records"]
    candidates = []
    for record in records:
        roots = tuple(record["roots"])
        moduli = normalized_moduli(roots)
        invariant = leading_invariant(moduli)
        if not rational_square(invariant):
            continue
        candidates.append((roots, moduli, invariant))
    selected = candidates[start:] if count is None else candidates[start : start + count]
    output = []
    for roots, moduli, invariant in selected:
        sections = (
            sections_across_primes(moduli, primes, precision)
            if all_primes
            else sections_at_seed(moduli, primes, precision)
        )
        transverse_pairs = []
        for left in range(len(sections)):
            for right in range(left + 1, len(sections)):
                rank, pivots = row_reduce(
                    [
                        value.gradient
                        for value in residuals((*moduli, *sections[left], *sections[right]), Field())
                    ],
                    Field(),
                )
                if rank == 7:
                    transverse_pairs.append(
                        {
                            "sections": [
                                [rational_text(value) for value in section]
                                for section in (sections[left], sections[right])
                            ],
                            "rank": rank,
                            "pivot_columns": pivots,
                        }
                    )
        output.append(
            {
                "roots": list(roots),
                "leading_invariant": rational_text(invariant),
                "reconstructed_sections": [
                    [rational_text(value) for value in section] for section in sections
                ],
                "transverse_rank_seven_pairs": transverse_pairs,
            }
        )
    return {
        "status": "bounded square-leading two-section seed screen completed",
        "input": str(input_path),
        "candidate_count": len(candidates),
        "candidate_offset": start,
        "candidate_batch_count": len(output),
        "prime_order": list(primes),
        "hensel_precision": precision,
        "exact_transverse_pair_count": sum(
            len(row["transverse_rank_seven_pairs"]) for row in output
        ),
        "records_with_reconstructed_sections": sum(
            bool(row["reconstructed_sections"]) for row in output
        ),
        "records": output,
        "scope": (
            "only nonsingular projected sections recovered from the union over "
            "the listed small primes are covered; no absence statement is made"
            if all_primes
            else "only nonsingular projected sections recovered from the first usable "
            "small prime are covered; no absence statement is made"
        ),
    }


def merge_batches(paths: Iterable[Path]) -> dict[str, object]:
    """Merge disjoint bounded outputs into the canonical whole-screen JSON."""

    batches = [json.loads(path.read_text()) for path in paths]
    if not batches:
        raise ValueError("at least one batch is required")
    reference = batches[0]
    keys = ("input", "candidate_count", "prime_order", "hensel_precision", "scope")
    if any(any(batch[key] != reference[key] for key in keys) for batch in batches[1:]):
        raise ValueError("the supplied batches do not have the same screen parameters")
    records: dict[int, dict[str, object]] = {}
    for batch in batches:
        offset = int(batch["candidate_offset"])
        listed = batch["records"]
        if len(listed) != int(batch["candidate_batch_count"]):
            raise ValueError("a batch record count changed")
        for index, record in enumerate(listed, start=offset):
            if index in records:
                raise ValueError("the bounded batches overlap")
            records[index] = record
    count = int(reference["candidate_count"])
    if sorted(records) != list(range(count)):
        raise ValueError("the bounded batches do not cover the whole candidate interval")
    output = [records[index] for index in range(count)]
    return {
        "status": "bounded square-leading two-section seed screen completed",
        "input": reference["input"],
        "candidate_count": count,
        "prime_order": reference["prime_order"],
        "hensel_precision": reference["hensel_precision"],
        "exact_transverse_pair_count": sum(
            len(row["transverse_rank_seven_pairs"]) for row in output
        ),
        "records_with_reconstructed_sections": sum(
            bool(row["reconstructed_sections"]) for row in output
        ),
        "records": output,
        "scope": reference["scope"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--primes", nargs="+", type=int, default=DEFAULT_PRIMES)
    parser.add_argument("--precision", type=int, default=8)
    parser.add_argument(
        "--all-primes",
        action="store_true",
        help="union exact reconstructions over every listed small prime",
    )
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--count", type=int)
    parser.add_argument("--merge-batches", nargs="+", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = (
        merge_batches(args.merge_batches)
        if args.merge_batches is not None
        else screen(
            args.input, args.primes, args.precision,
            start=args.start, count=args.count, all_primes=args.all_primes,
        )
    )
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        args.output.write_text(rendered)


if __name__ == "__main__":
    main()
