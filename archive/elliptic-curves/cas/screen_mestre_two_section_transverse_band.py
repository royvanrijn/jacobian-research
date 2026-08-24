#!/usr/bin/env python3
"""Exact post-max-200 band screen for transverse Mestre two-section seeds."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

from probe_mestre_two_section_local_continuation import Field, residuals, row_reduce
from screen_mestre_two_section_transverse_seeds import (
    DEFAULT_PRIMES,
    leading_invariant,
    normalized_moduli,
    rational_square,
    rational_text,
    sections_at_seed,
)


Q = Fraction


def enumerate_band(first: int, last: int) -> tuple[tuple[tuple[int, ...], bool], ...]:
    compiler = shutil.which("c++")
    if compiler is None:
        raise FileNotFoundError("a C++17 compiler is required")
    source = Path(__file__).with_name("enumerate_mestre_root_tuples_scale_band.cpp")
    with tempfile.TemporaryDirectory(prefix="mestre-band-") as directory:
        binary = Path(directory) / "enumerator"
        subprocess.run([compiler, "-std=c++17", "-O3", "-DNDEBUG", str(source), "-o", str(binary)], check=True)
        output = subprocess.run([str(binary), str(first), str(last)], check=True, text=True, capture_output=True).stdout
    lines = output.splitlines()
    if not lines or lines[0] != "MESTRE_ROOT_TUPLES_BAND_V1":
        raise AssertionError("band enumerator header changed")
    records = []
    for line in lines[1:-1]:
        fields = line.split()
        if len(fields) != 8 or fields[0] != "R":
            raise AssertionError("malformed band root record")
        records.append((tuple(map(int, fields[1:7])), bool(int(fields[7]))))
    summary = lines[-1].split()
    if len(summary) != 6 or summary[:3] != ["S", str(first), str(last)]:
        raise AssertionError("band enumerator summary changed")
    if int(summary[4]) != len(records):
        raise AssertionError("band enumerator obstruction count changed")
    return tuple(records)


def sections_across_primes(
    moduli: tuple[Fraction, ...], precision: int
) -> tuple[tuple[Fraction, Fraction], ...]:
    """Union exact reconstructions over every declared small prime.

    A rational affine section need not have coordinates integral at the first
    prime where the normalized moduli are integral, so this avoids treating a
    single-prime miss as a chart-level negative result.
    """

    sections = set()
    for prime in DEFAULT_PRIMES:
        try:
            sections.update(sections_at_seed(moduli, (prime,), precision))
        except ValueError:
            continue
    return tuple(sorted(sections))


def screen(first: int, last: int, precision: int, *, all_primes: bool = False) -> dict[str, object]:
    records = enumerate_band(first, last)
    candidates = []
    for roots, symmetric in records:
        moduli = normalized_moduli(roots)
        if not symmetric and rational_square(leading_invariant(moduli)):
            candidates.append((roots, moduli))
    hits = []
    section_count_histogram: dict[int, int] = {}
    multi_companion_records = []
    multiple_section_records = []
    for roots, moduli in candidates:
        if all_primes:
            sections = sections_across_primes(moduli, precision)
        else:
            # Keep the original single-prime replay semantics available for
            # the pinned 201--205 artifact.  New exhaustive-in-list runs pass
            # ``all_primes=True`` below.
            sections = sections_at_seed(moduli, DEFAULT_PRIMES, precision)
        section_count_histogram[len(sections)] = section_count_histogram.get(len(sections), 0) + 1
        pairs = []
        rank_histogram: dict[int, int] = {}
        for left in range(len(sections)):
            for right in range(left + 1, len(sections)):
                rank, pivots = row_reduce(
                    [value.gradient for value in residuals((*moduli, *sections[left], *sections[right]), Field())],
                    Field(),
                )
                rank_histogram[rank] = rank_histogram.get(rank, 0) + 1
                if rank == 7:
                    pairs.append({"sections": [[rational_text(v) for v in item] for item in (sections[left], sections[right])], "pivot_columns": pivots})
        if pairs:
            hits.append({"roots": list(roots), "transverse_pairs": pairs})
        if len(sections) >= 2:
            multiple_section_records.append(
                {
                    "roots": list(roots),
                    "sections": [
                        [rational_text(value) for value in section]
                        for section in sections
                    ],
                    "pair_jacobian_rank_histogram": rank_histogram,
                }
            )
        if len(sections) >= 3:
            multi_companion_records.append(
                {
                    "roots": list(roots),
                    "section_count": len(sections),
                    "pair_jacobian_rank_histogram": rank_histogram,
                }
            )
    root_stream = "\n".join(",".join(map(str, roots)) for roots, _ in records)
    return {
        "status": "exact post-max-200 transverse seed band screen completed",
        "diameter_band": [first, last],
        "obstruction_zero_record_count": len(records),
        "nonreflection_square_leading_candidate_count": len(candidates),
        "candidate_root_sha256": hashlib.sha256(root_stream.encode()).hexdigest(),
        "hensel_primes": list(DEFAULT_PRIMES),
        "prime_strategy": (
            "union of exact reconstructed sections over every declared prime"
            if all_primes
            else "first declared prime with integral normalized moduli"
        ),
        "hensel_precision": precision,
        "rank_seven_pair_count": sum(len(item["transverse_pairs"]) for item in hits),
        "rank_seven_hits": hits,
        "reconstructed_section_count_histogram": {
            str(count): frequency for count, frequency in sorted(section_count_histogram.items())
        },
        "multi_companion_records": sorted(
            multi_companion_records,
            key=lambda record: (-int(record["section_count"]), record["roots"]),
        ),
        "multiple_section_records": sorted(
            multiple_section_records,
            key=lambda record: (-len(record["sections"]), record["roots"]),
        ),
        "scope": "a bounded new-diameter screen; no absence claim outside this band or section chart",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--first", type=int, default=201)
    parser.add_argument("--last", type=int, default=205)
    parser.add_argument("--precision", type=int, default=8)
    parser.add_argument(
        "--all-primes",
        action="store_true",
        help="union exact reconstructions over every declared small prime",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rendered = json.dumps(
        screen(args.first, args.last, args.precision, all_primes=args.all_primes),
        indent=2,
        sort_keys=True,
    ) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        args.output.write_text(rendered)


if __name__ == "__main__":
    main()
