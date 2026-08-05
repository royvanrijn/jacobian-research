#!/usr/bin/env python3
"""Classify all complete-two-row/column fringe supports in bidegree (3,3).

A row-fringe support consists of any two complete rows, together with one
entry in a third row.  Column-fringe supports are their transposes.  There
are 48 supports of each type and 24 four-element orbits under coefficient-
matrix transpose and simultaneous row/column reversal.

For one representative of every orbit, this checker normalizes the dense
coefficient torus, constructs the exact moments mu_1,...,mu_14 over QQ,
and asks msolve for the saturated moment scheme.  Every representative is
the unit ideal.  Coordinate boundaries have support at most eight and are
covered by the existing complete support-eight classification.

This is one structured size-nine class, not the full size-nine census.
"""

from __future__ import annotations

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
    / "two_pair_sic_bidegree33_sparse_two_row_fringe9.json"
)
MOMENT_ORDERS = tuple(range(1, 15))
Support = frozenset[tuple[int, int]]
POSITIONS = tuple((row, column) for row in range(4) for column in range(4))


def transpose(support: Support) -> Support:
    return frozenset((column, row) for row, column in support)


def reverse(support: Support) -> Support:
    return frozenset((3 - row, 3 - column) for row, column in support)


def symmetry_orbit(support: Support) -> frozenset[Support]:
    return frozenset({
        support,
        transpose(support),
        reverse(support),
        transpose(reverse(support)),
    })


def row_fringe_supports() -> set[Support]:
    answer: set[Support] = set()
    for full_rows in combinations(range(4), 2):
        fringe_rows = tuple(
            row for row in range(4) if row not in full_rows
        )
        for fringe_row in fringe_rows:
            for fringe_column in range(4):
                answer.add(frozenset({
                    *((row, column) for row in full_rows for column in range(4)),
                    (fringe_row, fringe_column),
                }))
    return answer


def support_class() -> tuple[set[Support], set[Support], set[Support]]:
    rows = row_fringe_supports()
    columns = {transpose(support) for support in rows}
    if len(rows) != 48 or len(columns) != 48:
        raise AssertionError("unexpected row/column fringe count")
    if rows & columns:
        raise AssertionError("row and column fringe classes must be disjoint")
    return rows, columns, rows | columns


def orbit_size_distribution(supports: set[Support]) -> dict[int, int]:
    unseen = set(supports)
    distribution: dict[int, int] = {}
    while unseen:
        seed = unseen.pop()
        orbit = symmetry_orbit(seed)
        unseen.difference_update(orbit)
        distribution[len(orbit)] = distribution.get(len(orbit), 0) + 1
    return distribution


def mixed_size_nine_census(fringe_supports: set[Support]) -> dict[str, object]:
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
    closed = rectangles | fringe_supports
    remaining = mixed - closed
    if len(mixed) != 11420 or len(rectangles) != 16 or len(closed) != 112:
        raise AssertionError("unexpected mixed size-nine support census")
    distributions = {
        "all": orbit_size_distribution(mixed),
        "rectangles": orbit_size_distribution(rectangles),
        "fringes": orbit_size_distribution(fringe_supports),
        "closed_union": orbit_size_distribution(closed),
        "remaining": orbit_size_distribution(remaining),
    }
    expected = {
        "all": {2: 138, 4: 2786},
        "rectangles": {2: 4, 4: 2},
        "fringes": {4: 24},
        "closed_union": {2: 4, 4: 26},
        "remaining": {2: 134, 4: 2760},
    }
    if distributions != expected:
        raise AssertionError(f"unexpected symmetry census: {distributions}")
    return {
        "mixed_support_count": len(mixed),
        "mixed_symmetry_orbit_count": sum(distributions["all"].values()),
        "mixed_orbits_by_size": distributions["all"],
        "closed_support_count": len(closed),
        "closed_symmetry_orbit_count": sum(
            distributions["closed_union"].values()
        ),
        "closed_orbits_by_size": distributions["closed_union"],
        "remaining_support_count": len(remaining),
        "remaining_symmetry_orbit_count": sum(
            distributions["remaining"].values()
        ),
        "remaining_orbits_by_size": distributions["remaining"],
        "scope": (
            "discrete standard-basis support orbits under transpose and "
            "simultaneous reversal, not continuous diagonal-SL2 orbits"
        ),
    }


def orbit_representatives(
    supports: set[Support],
) -> tuple[tuple[Support, frozenset[Support]], ...]:
    unseen = set(supports)
    records: list[tuple[Support, frozenset[Support]]] = []
    while unseen:
        seed = min(unseen, key=lambda value: tuple(sorted(value)))
        orbit = symmetry_orbit(seed)
        if len(orbit) != 4 or not orbit <= supports:
            raise AssertionError("the fringe class must have four-element orbits")
        unseen.difference_update(orbit)
        representative = min(
            orbit,
            key=lambda value: tuple(sorted(value)),
        )
        records.append((representative, orbit))
    records.sort(key=lambda record: tuple(sorted(record[0])))
    if len(records) != 24:
        raise AssertionError("unexpected fringe symmetry orbit count")
    return tuple(records)


def support_shape(support: Support) -> str:
    row_degrees = sorted(
        (sum((row, column) in support for column in range(4)) for row in range(4)),
        reverse=True,
    )
    column_degrees = sorted(
        (sum((row, column) in support for row in range(4)) for column in range(4)),
        reverse=True,
    )
    if row_degrees == [4, 4, 1, 0] and column_degrees == [3, 2, 2, 2]:
        return "row fringe"
    if row_degrees == [3, 2, 2, 2] and column_degrees == [4, 4, 1, 0]:
        return "column fringe"
    raise AssertionError("support is not in the declared fringe class")


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
        hashlib.sha256(serialization.encode("utf-8")).hexdigest(),
    )


def classify_representative(
    representative: Support,
    orbit: frozenset[Support],
    msolve: str,
) -> dict[str, object]:
    support = tuple(sorted(representative))
    normalized, term_counts, equations_sha256 = equation_profile(support)
    record = screen_support(
        support,
        through=14,
        timeout=120,
        msolve=msolve,
        threads=4,
        rational_parametrization=False,
    )
    if (
        record["status"] != "excluded_on_coefficient_torus"
        or record["msolve"]["returncode"] != 0
        or record["msolve"]["result"] not in ("[-1]", "[-1]:")
    ):
        raise AssertionError(
            f"nonunit or failed exact fringe system on {support}: {record}"
        )
    return {
        "representative": [list(position) for position in support],
        "representative_shape": support_shape(representative),
        "normalized_anchors": [
            list(normalized[0]),
            list(normalized[1]),
        ],
        "moment_term_counts": term_counts,
        "equations_sha256": equations_sha256,
        "dense_torus_moment_ideal": "unit over QQ",
        "symmetry_orbit": [
            [list(position) for position in sorted(member)]
            for member in sorted(orbit, key=lambda value: tuple(sorted(value)))
        ],
    }


def main() -> None:
    msolve = shutil.which("msolve")
    if msolve is None:
        raise RuntimeError("msolve is required")

    rows, columns, supports = support_class()
    representatives = orbit_representatives(supports)
    classifications = [
        classify_representative(representative, orbit, msolve)
        for representative, orbit in representatives
    ]
    size_nine_census = mixed_size_nine_census(supports)

    artifact = {
        "format": "two-pair-sic-bidegree33-sparse-two-row-fringe9-v2",
        "field": "characteristic zero",
        "support_class": (
            "any two complete rows plus one entry in a third row, "
            "together with the transposed any-two-complete-columns class"
        ),
        "support_size": 9,
        "row_fringe_support_count": len(rows),
        "column_fringe_support_count": len(columns),
        "total_support_count": len(supports),
        "symmetry_group": [
            "coefficient-matrix transpose",
            "simultaneous row/column reversal",
        ],
        "symmetry_orbit_count": len(representatives),
        "symmetry_orbit_sizes": [len(orbit) for _, orbit in representatives],
        "full_mixed_size_nine_support_census": size_nine_census,
        "moment_orders": list(MOMENT_ORDERS),
        "component_method": (
            "exact QQ moment equations with dense coefficient-torus "
            "Rabinowitsch localization; msolve returns the unit ideal "
            "for every orbit representative"
        ),
        "exact_unit_system_count": len(classifications),
        "dense_torus_conclusion": (
            "no support in this class has a pure-moment point with all "
            "nine displayed coefficients nonzero"
        ),
        "boundary_conclusion": (
            "every coordinate boundary has support at most eight and is "
            "SIC-safe by the existing complete support-eight theorem"
        ),
        "global_conclusion_for_this_support_class": (
            "all 96 complete-two-row/column fringe coordinate "
            "subspaces are SIC-safe"
        ),
        "scope": (
            "complete exact classification of one 96-support structured "
            "size-nine class; not the full 11420 mixed size-nine census or "
            "the dense bidegree-(3,3) orbit classification"
        ),
        "independent_formula_check": verify_restricted_formula(),
        "orbits": classifications,
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")

    print("PASS 48 row and 48 column fringe supports are distinct")
    print("PASS 24 four-element symmetry orbits cover all 96 supports")
    print("PASS all 24 dense-torus moment ideals are units over QQ")
    print("PASS all complete-two-row/column fringe supports are SIC-safe")
    print("PASS 2924 mixed size-nine support orbits; 2894 remain open")


if __name__ == "__main__":
    main()
