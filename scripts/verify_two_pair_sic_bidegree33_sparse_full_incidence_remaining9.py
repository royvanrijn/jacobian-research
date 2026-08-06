#!/usr/bin/env python3
"""Close the remaining full-incidence size-nine supports in bidegree (3,3).

After the (3,2,2,2)^2 class, exactly four unordered row/column partition
types remain.  They contain 7,050 mixed supports in 1,792 orbits under
transpose and simultaneous reversal.  Exact QQ moment
equations through mu_10 give the unit ideal on every normalized dense
coefficient torus.  Together with the earlier structured support
theorems, this closes the complete 11,420-support size-nine census.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import combinations
import json
from pathlib import Path
import shutil

from research_two_pair_sic_bidegree33_sparse_six_counterexample import (
    screen_support,
    verify_restricted_formula,
)
from verify_two_pair_sic_bidegree33_sparse_cross_two9 import symmetry_orbit


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_sparse_full_incidence_remaining9.json"
)
POSITIONS = tuple((row, column) for row in range(4) for column in range(4))
Support = frozenset[tuple[int, int]]
Partition = tuple[int, int, int, int]

P4221: Partition = (4, 2, 2, 1)
P3321: Partition = (3, 3, 2, 1)
P3222: Partition = (3, 2, 2, 2)
TARGET_TYPES = frozenset(
    {
        frozenset((P4221, P3222)),
        frozenset((P4221, P3321)),
        frozenset((P3321,)),
        frozenset((P3321, P3222)),
    }
)
TYPE_COUNTS = {
    "(4,2,2,1)x(3,2,2,2)": (1152, 288),
    "(4,2,2,1)x(3,3,2,1)": (1436, 359),
    "(3,3,2,1)x(3,3,2,1)": (1870, 497),
    "(3,3,2,1)x(3,2,2,2)": (2592, 648),
}


def is_mixed(support: Support) -> bool:
    return (
        any(row > column for row, column in support)
        and any(row < column for row, column in support)
    )


def line_counts(support: Support, axis: int) -> Partition:
    return tuple(
        sorted(
            (
                sum(position[axis] == value for position in support)
                for value in range(4)
            ),
            reverse=True,
        )
    )


def incidence_type(support: Support) -> frozenset[Partition]:
    return frozenset((line_counts(support, 0), line_counts(support, 1)))


def type_label(support: Support) -> str:
    rows = line_counts(support, 0)
    columns = line_counts(support, 1)
    ordered = (rows, columns)
    if P4221 in ordered:
        other = columns if rows == P4221 else rows
        return (
            "(4,2,2,1)x(3,2,2,2)"
            if other == P3222
            else "(4,2,2,1)x(3,3,2,1)"
        )
    if rows == columns == P3321:
        return "(3,3,2,1)x(3,3,2,1)"
    return "(3,3,2,1)x(3,2,2,2)"


def support_class() -> set[Support]:
    answer = {
        support
        for positions in combinations(POSITIONS, 9)
        if is_mixed(support := frozenset(positions))
        and incidence_type(support) in TARGET_TYPES
    }
    if len(answer) != 7050:
        raise AssertionError("unexpected remaining full-incidence support count")
    counts: dict[str, int] = {}
    for support in answer:
        label = type_label(support)
        counts[label] = counts.get(label, 0) + 1
    expected = {label: count for label, (count, _) in TYPE_COUNTS.items()}
    if counts != expected:
        raise AssertionError(f"unexpected incidence-type counts: {counts}")
    return answer


def representatives(
    supports: set[Support],
) -> tuple[tuple[Support, frozenset[Support]], ...]:
    unseen = set(supports)
    answer: list[tuple[Support, frozenset[Support]]] = []
    while unseen:
        seed = min(unseen, key=lambda value: tuple(sorted(value)))
        orbit = symmetry_orbit(seed)
        if len(orbit) not in (2, 4) or not orbit <= supports:
            raise AssertionError("unexpected remaining full-incidence symmetry orbit")
        unseen.difference_update(orbit)
        representative = min(orbit, key=lambda value: tuple(sorted(value)))
        answer.append((representative, orbit))
    answer.sort(key=lambda record: tuple(sorted(record[0])))
    if len(answer) != 1792:
        raise AssertionError("unexpected remaining full-incidence orbit count")
    counts: dict[str, int] = {}
    for representative, _ in answer:
        label = type_label(representative)
        counts[label] = counts.get(label, 0) + 1
    expected = {label: count for label, (_, count) in TYPE_COUNTS.items()}
    if counts != expected:
        raise AssertionError(f"unexpected incidence-type orbit counts: {counts}")
    return tuple(answer)


def full_census(new_supports: set[Support]) -> dict[str, object]:
    mixed = {
        support
        for positions in combinations(POSITIONS, 9)
        if is_mixed(support := frozenset(positions))
    }
    earlier = mixed - new_supports
    if len(mixed) != 11420 or len(earlier) != 4370:
        raise AssertionError("unexpected complete mixed size-nine census")

    def orbit_sizes(supports: set[Support]) -> dict[int, int]:
        unseen = set(supports)
        answer: dict[int, int] = {}
        while unseen:
            seed = unseen.pop()
            orbit = symmetry_orbit(seed) & supports
            unseen.difference_update(orbit)
            answer[len(orbit)] = answer.get(len(orbit), 0) + 1
        return answer

    earlier_orbits = orbit_sizes(earlier)
    new_orbits = orbit_sizes(new_supports)
    all_orbits = orbit_sizes(mixed)
    if earlier_orbits != {2: 79, 4: 1053}:
        raise AssertionError(f"unexpected earlier orbit sizes: {earlier_orbits}")
    if new_orbits != {2: 59, 4: 1733}:
        raise AssertionError(f"unexpected new orbit sizes: {new_orbits}")
    if all_orbits != {2: 138, 4: 2786}:
        raise AssertionError(f"unexpected complete orbit sizes: {all_orbits}")
    return {
        "mixed_support_count": len(mixed),
        "closed_mixed_support_count": len(mixed),
        "remaining_mixed_support_count": 0,
        "mixed_symmetry_orbit_count": sum(all_orbits.values()),
        "orbits_by_size": all_orbits,
        "new_support_count": len(new_supports),
        "new_symmetry_orbit_count": sum(new_orbits.values()),
        "scope": (
            "discrete standard-basis support orbits under transpose and "
            "simultaneous reversal, not continuous diagonal-SL2 orbits"
        ),
    }


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
        raise AssertionError(f"nonunit remaining full-incidence system: {result}")
    return {
        "incidence_type": type_label(representative),
        "representative": [list(position) for position in support],
        "dense_torus_moment_ideal_through_mu_10": "unit over QQ",
        "seconds": result["seconds"],
        "symmetry_orbit": [
            [list(position) for position in sorted(member)]
            for member in sorted(orbit, key=lambda value: tuple(sorted(value)))
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--timeout", type=int, default=600)
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
    with ProcessPoolExecutor(max_workers=arguments.workers) as executor:
        futures = [
            executor.submit(classify, representative, orbit, msolve, arguments.timeout)
            for representative, orbit in reps
        ]
        for completed, future in enumerate(as_completed(futures), start=1):
            records.append(future.result())
            if completed % 64 == 0 or completed == len(futures):
                print(f"classified {completed}/{len(futures)}", flush=True)
    records.sort(key=lambda record: record["representative"])
    if len(records) != 1792:
        raise AssertionError("not every remaining full-incidence orbit was classified")

    artifact = {
        "format": "two-pair-sic-bidegree33-sparse-full-incidence-remaining9-v1",
        "field": "characteristic zero",
        "support_class": (
            "the four mixed full-incidence row/column partition types "
            "remaining after the (3,2,2,2)^2 class"
        ),
        "support_size": 9,
        "incidence_type_counts": {
            label: {"supports": counts[0], "symmetry_orbits": counts[1]}
            for label, counts in TYPE_COUNTS.items()
        },
        "total_support_count": len(supports),
        "symmetry_orbit_count": len(reps),
        "symmetry_orbit_sizes": [len(orbit) for _, orbit in reps],
        "moment_orders": [1, 10],
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
            "all 7050 remaining full-incidence coordinate subspaces are "
            "SIC-safe, completing the mixed size-nine census"
        ),
        "complete_size_nine_census": full_census(supports),
        "independent_formula_check": verify_restricted_formula(),
        "scope": (
            "complete exact standard-basis size-nine support classification, "
            "not a continuous diagonal-SL2 orbit classification"
        ),
        "orbits": records,
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(artifact, indent=2) + "\n")

    print("PASS 1792 symmetry orbits cover all 7050 supports")
    print("PASS all 1792 dense coefficient tori are units over QQ through mu_10")
    print("PASS all 11420 mixed size-nine supports are SIC-safe")


if __name__ == "__main__":
    main()
