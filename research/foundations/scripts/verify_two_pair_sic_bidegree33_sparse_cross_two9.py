#!/usr/bin/env python3
"""Classify every size-nine cross-plus-two support in bidegree (3,3).

A cross-plus-two support consists of one complete row, one complete column,
and any two further entries.  There are 576 supports in 156 orbits under
coefficient-matrix transpose and simultaneous reversal.

For every representative this checker constructs mu_1,...,mu_10 over QQ
on the normalized dense coefficient torus.  Exact msolve elimination gives
the unit ideal in every case.  All coordinate boundaries have support at
most eight and are covered by the existing sparse classification.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
from itertools import combinations
import json
from pathlib import Path
import shutil

from research_two_pair_sic_bidegree33_sparse_six_counterexample import (
    RESIDUAL_SYMBOL_POOL,
    msolve_expression,
    normalize_support,
    restricted_moment,
    screen_support,
    verify_restricted_formula,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_sparse_cross_two9.json"
)
MOMENT_ORDERS = tuple(range(1, 11))
POSITIONS = tuple((row, column) for row in range(4) for column in range(4))
Support = frozenset[tuple[int, int]]


def transpose(support: Support) -> Support:
    return frozenset((column, row) for row, column in support)


def reverse(support: Support) -> Support:
    return frozenset((3 - row, 3 - column) for row, column in support)


def symmetry_orbit(support: Support) -> frozenset[Support]:
    return frozenset(
        {
            support,
            transpose(support),
            reverse(support),
            transpose(reverse(support)),
        }
    )


def support_class() -> tuple[set[Support], set[Support], set[Support]]:
    all_crosses: set[Support] = set()
    row_arms: set[Support] = set()
    for full_row in range(4):
        for full_column in range(4):
            cross = {
                *((full_row, column) for column in range(4)),
                *((row, full_column) for row in range(4)),
            }
            for extras in combinations(set(POSITIONS) - cross, 2):
                all_crosses.add(frozenset(cross | set(extras)))
            for arm_row in range(4):
                if arm_row == full_row:
                    continue
                available = tuple(
                    column for column in range(4) if column != full_column
                )
                for arm_columns in combinations(available, 2):
                    row_arms.add(
                        frozenset(
                            cross
                            | {(arm_row, column) for column in arm_columns}
                        )
                    )
    aligned = row_arms | {transpose(support) for support in row_arms}
    matching = all_crosses - aligned
    if len(aligned) != 288 or len(matching) != 288 or len(all_crosses) != 576:
        raise AssertionError("unexpected cross-plus-two support count")
    return aligned, matching, all_crosses


def representatives(
    supports: set[Support],
) -> tuple[tuple[Support, frozenset[Support]], ...]:
    unseen = set(supports)
    answer: list[tuple[Support, frozenset[Support]]] = []
    while unseen:
        seed = min(unseen, key=lambda value: tuple(sorted(value)))
        orbit = symmetry_orbit(seed)
        if len(orbit) not in (2, 4) or not orbit <= supports:
            raise AssertionError("unexpected cross-plus-two symmetry orbit")
        unseen.difference_update(orbit)
        representative = min(orbit, key=lambda value: tuple(sorted(value)))
        answer.append((representative, orbit))
    answer.sort(key=lambda record: tuple(sorted(record[0])))
    orbit_sizes = {2: 0, 4: 0}
    for _, orbit in answer:
        orbit_sizes[len(orbit)] += 1
    if len(answer) != 156 or orbit_sizes != {2: 24, 4: 132}:
        raise AssertionError("unexpected cross-plus-two orbit count")
    return tuple(answer)


def equation_profile(
    support: tuple[tuple[int, int], ...],
) -> tuple[tuple[tuple[int, int], ...], list[int], str]:
    normalized = normalize_support(support)
    residuals = RESIDUAL_SYMBOL_POOL[:7]
    moments = [
        restricted_moment(order, normalized, residuals)
        for order in MOMENT_ORDERS
    ]
    serialization = "\n".join(
        f"mu_{order}={msolve_expression(moment)}"
        for order, moment in zip(MOMENT_ORDERS, moments, strict=True)
    )
    return (
        normalized,
        [len(moment.terms()) for moment in moments],
        hashlib.sha256(serialization.encode()).hexdigest(),
    )


def classify(
    representative: Support,
    orbit: frozenset[Support],
    msolve: str,
    timeout: int,
) -> dict[str, object]:
    support = tuple(sorted(representative))
    normalized, term_counts, equations_sha256 = equation_profile(support)
    result = screen_support(
        support,
        through=10,
        timeout=timeout,
        msolve=msolve,
        threads=1,
        rational_parametrization=False,
    )
    if (
        result["status"] != "excluded_on_coefficient_torus"
        or result["msolve"]["returncode"] != 0
        or result["msolve"]["result"] not in ("[-1]", "[-1]:")
    ):
        raise AssertionError(f"nonunit cross-plus-two system: {result}")
    return {
        "representative": [list(position) for position in support],
        "normalized_anchors": [list(normalized[0]), list(normalized[1])],
        "moment_term_counts": term_counts,
        "equations_sha256": equations_sha256,
        "dense_torus_moment_ideal_through_mu_10": "unit over QQ",
        "symmetry_orbit": [
            [list(position) for position in sorted(member)]
            for member in sorted(orbit, key=lambda value: tuple(sorted(value)))
        ],
    }


def updated_census(cross_arms: set[Support]) -> dict[str, object]:
    mixed = {
        frozenset(support)
        for support in combinations(POSITIONS, 9)
        if any(row > column for row, column in support)
        and any(row < column for row, column in support)
    }
    rectangles = {
        frozenset(
            (row, column)
            for row, column in POSITIONS
            if row != missing_row and column != missing_column
        )
        for missing_row in range(4)
        for missing_column in range(4)
    }
    fringes: set[Support] = set()
    for full_rows in combinations(range(4), 2):
        for fringe_row in range(4):
            if fringe_row in full_rows:
                continue
            for fringe_column in range(4):
                fringes.add(
                    frozenset(
                        {
                            *((row, column) for row in full_rows for column in range(4)),
                            (fringe_row, fringe_column),
                        }
                    )
                )
    fringes |= {transpose(support) for support in fringes}
    all_closed = rectangles | fringes | cross_arms
    closed = all_closed & mixed
    nonmixed_closed = all_closed - mixed
    if (
        len(mixed) != 11420
        or len(all_closed) != 688
        or len(closed) != 682
        or len(nonmixed_closed) != 6
    ):
        raise AssertionError("unexpected updated size-nine census")
    unseen = set(closed)
    orbit_sizes: dict[int, int] = {}
    while unseen:
        seed = unseen.pop()
        orbit = symmetry_orbit(seed)
        unseen.difference_update(orbit)
        orbit_sizes[len(orbit)] = orbit_sizes.get(len(orbit), 0) + 1
    if orbit_sizes != {2: 27, 4: 157}:
        raise AssertionError(f"unexpected closed orbit sizes: {orbit_sizes}")
    return {
        "mixed_support_count": 11420,
        "closed_coordinate_subspace_count": len(all_closed),
        "closed_nonmixed_support_count": len(nonmixed_closed),
        "closed_mixed_support_count": len(closed),
        "remaining_mixed_support_count": 10738,
        "closed_mixed_symmetry_orbit_count": 184,
        "closed_orbits_by_size": orbit_sizes,
        "remaining_mixed_symmetry_orbit_count": 2740,
        "remaining_orbits_by_size": {2: 111, 4: 2629},
        "scope": (
            "discrete standard-basis support orbits under transpose and "
            "simultaneous reversal, not continuous diagonal-SL2 orbits"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()
    if arguments.workers < 1 or arguments.timeout < 1:
        raise ValueError("--workers and --timeout must be positive")
    msolve = shutil.which("msolve")
    if msolve is None:
        raise RuntimeError("msolve is required")

    aligned, matching, supports = support_class()
    reps = representatives(supports)
    records: list[dict[str, object]] = []
    with ThreadPoolExecutor(max_workers=arguments.workers) as executor:
        futures = [
            executor.submit(
                classify,
                representative,
                orbit,
                msolve,
                arguments.timeout,
            )
            for representative, orbit in reps
        ]
        for future in as_completed(futures):
            records.append(future.result())
    records.sort(key=lambda record: record["representative"])

    artifact = {
        "format": "two-pair-sic-bidegree33-sparse-cross-plus-two9-v1",
        "field": "characteristic zero",
        "support_class": (
            "one complete row, one complete column, and any two extra "
            "entries outside their cross"
        ),
        "support_size": 9,
        "aligned_arm_support_count": len(aligned),
        "matching_extra_support_count": len(matching),
        "total_support_count": len(supports),
        "symmetry_orbit_count": len(reps),
        "symmetry_orbit_sizes": [len(orbit) for _, orbit in reps],
        "moment_orders": list(MOMENT_ORDERS),
        "exact_unit_system_count": len(records),
        "component_method": (
            "exact QQ moments with dense coefficient-torus Rabinowitsch "
            "localization; msolve returns the unit ideal through mu_10"
        ),
        "boundary_conclusion": (
            "every coordinate boundary has support at most eight and is "
            "SIC-safe by the complete support-eight theorem"
        ),
        "global_conclusion": (
            "all 576 cross-plus-two coordinate subspaces are SIC-safe"
        ),
        "updated_size_nine_census": updated_census(supports),
        "independent_formula_check": verify_restricted_formula(),
        "scope": (
            "complete exact classification of this 576-support class, not "
            "the full size-nine moment classification"
        ),
        "orbits": records,
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(artifact, indent=2) + "\n")

    print("PASS 288 aligned-arm and 288 matching-extra supports are distinct")
    print("PASS 156 symmetry orbits cover all 576 supports")
    print("PASS all 156 dense-torus ideals are units over QQ through mu_10")
    print("PASS 682 mixed size-nine supports closed; 10738 remain")


if __name__ == "__main__":
    main()
