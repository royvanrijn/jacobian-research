#!/usr/bin/env python3
"""Classify the 4+3+2 three-line size-nine supports in bidegree (3,3).

These are the mixed supports whose row counts, or column counts, are
0,2,3,4.  There are 1,148 supports in 287 four-element orbits under
transpose and simultaneous reversal.  Exact QQ moment equations through
mu_10 give the unit ideal on every normalized dense coefficient torus.
Together with the rectangle, fringe, and regular-three-line theorems this
closes every size-nine support having an empty row or column.
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
from verify_two_pair_sic_bidegree33_sparse_cross_two9 import (
    support_class as cross_support_class,
    symmetry_orbit,
    transpose,
)
from verify_two_pair_sic_bidegree33_sparse_three_line9 import (
    support_class as regular_three_line_support_class,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_sparse_three_line4329.json"
)
POSITIONS = tuple((row, column) for row in range(4) for column in range(4))
MOMENT_ORDERS = tuple(range(1, 11))
Support = frozenset[tuple[int, int]]


def is_mixed(support: Support) -> bool:
    return (
        any(row > column for row, column in support)
        and any(row < column for row, column in support)
    )


def line_counts(support: Support, axis: int) -> list[int]:
    return sorted(
        sum(position[axis] == value for position in support)
        for value in range(4)
    )


def support_class() -> set[Support]:
    answer = {
        frozenset(support)
        for support in combinations(POSITIONS, 9)
        if is_mixed(frozenset(support))
        and (
            line_counts(frozenset(support), 0) == [0, 2, 3, 4]
            or line_counts(frozenset(support), 1) == [0, 2, 3, 4]
        )
    }
    if len(answer) != 1148:
        raise AssertionError("unexpected 4+3+2 three-line support count")
    return answer


def representatives(
    supports: set[Support],
) -> tuple[tuple[Support, frozenset[Support]], ...]:
    unseen = set(supports)
    answer: list[tuple[Support, frozenset[Support]]] = []
    while unseen:
        seed = min(unseen, key=lambda value: tuple(sorted(value)))
        orbit = symmetry_orbit(seed)
        if len(orbit) != 4 or not orbit <= supports:
            raise AssertionError("unexpected 4+3+2 three-line symmetry orbit")
        unseen.difference_update(orbit)
        representative = min(orbit, key=lambda value: tuple(sorted(value)))
        answer.append((representative, orbit))
    answer.sort(key=lambda record: tuple(sorted(record[0])))
    if len(answer) != 287:
        raise AssertionError("unexpected 4+3+2 three-line orbit count")
    return tuple(answer)


def equation_profile(
    support: tuple[tuple[int, int], ...],
) -> tuple[list[int], str]:
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
        raise AssertionError(f"nonunit 4+3+2 three-line system: {result}")
    term_counts, equations_sha256 = equation_profile(support)
    return {
        "representative": [list(position) for position in support],
        "moment_term_counts": term_counts,
        "equations_sha256": equations_sha256,
        "dense_torus_moment_ideal_through_mu_10": "unit over QQ",
        "seconds": result["seconds"],
        "symmetry_orbit": [
            [list(position) for position in sorted(member)]
            for member in sorted(orbit, key=lambda value: tuple(sorted(value)))
        ],
    }


def updated_census(new_supports: set[Support]) -> dict[str, object]:
    mixed = {
        frozenset(support)
        for support in combinations(POSITIONS, 9)
        if is_mixed(frozenset(support))
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
    crosses = cross_support_class()[2]
    regular = regular_three_line_support_class()
    missing_line = rectangles | fringes | regular | new_supports
    if len(missing_line) != 1740:
        raise AssertionError("unexpected complete missing-line support count")
    closed = (missing_line | crosses) & mixed
    remaining = mixed - closed
    if len(closed) != 2310 or len(remaining) != 9110:
        raise AssertionError("unexpected updated mixed size-nine census")

    def orbit_sizes(supports: set[Support]) -> dict[int, int]:
        unseen = set(supports)
        answer: dict[int, int] = {}
        while unseen:
            seed = unseen.pop()
            orbit = symmetry_orbit(seed) & supports
            unseen.difference_update(orbit)
            answer[len(orbit)] = answer.get(len(orbit), 0) + 1
        return answer

    closed_orbits = orbit_sizes(closed)
    remaining_orbits = orbit_sizes(remaining)
    if closed_orbits != {2: 27, 4: 564}:
        raise AssertionError(f"unexpected closed orbit sizes: {closed_orbits}")
    if remaining_orbits != {2: 111, 4: 2222}:
        raise AssertionError(f"unexpected remaining orbit sizes: {remaining_orbits}")
    return {
        "mixed_support_count": len(mixed),
        "missing_row_or_column_support_count": len(missing_line),
        "closed_mixed_support_count": len(closed),
        "remaining_mixed_support_count": len(remaining),
        "closed_mixed_symmetry_orbit_count": sum(closed_orbits.values()),
        "closed_orbits_by_size": closed_orbits,
        "remaining_mixed_symmetry_orbit_count": sum(remaining_orbits.values()),
        "remaining_orbits_by_size": remaining_orbits,
        "scope": (
            "discrete standard-basis support orbits under transpose and "
            "simultaneous reversal, not continuous diagonal-SL2 orbits"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()
    if arguments.workers < 1 or arguments.timeout < 1:
        raise ValueError("--workers and --timeout must be positive")
    msolve = shutil.which("msolve")
    if msolve is None:
        raise RuntimeError("msolve is required")

    supports = support_class()
    reps = representatives(supports)
    records: list[dict[str, object]] = []
    with ThreadPoolExecutor(max_workers=arguments.workers) as executor:
        futures = [
            executor.submit(classify, representative, orbit, msolve, arguments.timeout)
            for representative, orbit in reps
        ]
        for future in as_completed(futures):
            records.append(future.result())
    records.sort(key=lambda record: record["representative"])
    if len(records) != 287:
        raise AssertionError("not every 4+3+2 orbit was classified")

    artifact = {
        "format": "two-pair-sic-bidegree33-sparse-three-line4329-v1",
        "field": "characteristic zero",
        "support_class": (
            "mixed nine-entry supports with row or column counts 0,2,3,4"
        ),
        "support_size": 9,
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
            "all 1148 4+3+2 three-line coordinate subspaces are SIC-safe; "
            "all 1740 size-nine supports with an empty row or column are now closed"
        ),
        "updated_size_nine_census": updated_census(supports),
        "independent_formula_check": verify_restricted_formula(),
        "scope": (
            "complete exact classification of this 1148-support class, not "
            "the full size-nine moment classification"
        ),
        "orbits": records,
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(artifact, indent=2) + "\n")

    print("PASS 287 four-element symmetry orbits cover all 1148 supports")
    print("PASS all 287 dense coefficient tori are units over QQ through mu_10")
    print("PASS every size-nine support with an empty row or column is closed")
    print("PASS 2310 mixed size-nine supports closed; 9110 remain")


if __name__ == "__main__":
    main()
