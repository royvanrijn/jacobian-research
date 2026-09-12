#!/usr/bin/env python3
"""Classify the full-incidence (3,2,2,2)^2 size-nine supports.

There are 816 mixed supports whose row and column count partitions are
both 3+2+2+2.  They form 230 orbits under transpose and simultaneous
reversal.  Exact QQ elimination through mu_10 gives 228 unit dense-torus
systems and two systems with one rational point each.  This checker proves
that both points have coefficient rank three and enter the fixed chamber
i>j under an explicit contraction-preserving flag change.  They are
therefore one-sided and have mixed cutoff m>e.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
from itertools import combinations
import json
from pathlib import Path
import shutil

import sympy as sp

from research_two_pair_sic_bidegree33_sparse_six_counterexample import (
    RESIDUAL_SYMBOL_POOL,
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
from verify_two_pair_sic_bidegree33_sparse_full_line43119 import (
    support_class as full_line_4311_support_class,
)
from verify_two_pair_sic_bidegree33_sparse_three_line9 import (
    certify_rur_box,
    coefficient_data,
    flag_transform,
    mixed_contraction,
    support_class as regular_three_line_support_class,
)
from verify_two_pair_sic_bidegree33_sparse_three_line4329 import (
    support_class as three_line_432_support_class,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_sparse_full_incidence32229.json"
)
POSITIONS = tuple((row, column) for row in range(4) for column in range(4))
Support = frozenset[tuple[int, int]]
X, Y = sp.symbols("x y")


SURVIVORS: dict[
    Support,
    tuple[tuple[sp.Rational, ...], sp.Rational],
] = {
    frozenset(
        ((0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (2, 0), (2, 3), (3, 2), (3, 3))
    ): (
        (
            sp.Rational(2, 9),
            sp.Rational(-9, 4),
            sp.Rational(-3, 2),
            sp.Rational(9, 8),
            sp.Rational(1, 3),
            sp.Rational(-3, 4),
            sp.Rational(-1, 2),
        ),
        sp.Rational(2, 3),
    ),
    frozenset(
        ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 0), (2, 3), (3, 0), (3, 1))
    ): (
        (
            sp.Rational(2, 9),
            sp.Rational(-3),
            sp.Rational(-1),
            sp.Rational(27),
            sp.Rational(1),
            sp.Rational(-81),
            sp.Rational(-27),
        ),
        sp.Rational(1, 3),
    ),
}


def is_mixed(support: Support) -> bool:
    return (
        any(row > column for row, column in support)
        and any(row < column for row, column in support)
    )


def line_counts(support: Support, axis: int) -> tuple[int, ...]:
    return tuple(
        sorted(
            (
                sum(position[axis] == value for position in support)
                for value in range(4)
            ),
            reverse=True,
        )
    )


def support_class() -> set[Support]:
    answer = {
        frozenset(support)
        for support in combinations(POSITIONS, 9)
        if is_mixed(frozenset(support))
        and line_counts(frozenset(support), 0) == (3, 2, 2, 2)
        and line_counts(frozenset(support), 1) == (3, 2, 2, 2)
    }
    if len(answer) != 816:
        raise AssertionError("unexpected (3,2,2,2)^2 support count")
    return answer


def representatives(
    supports: set[Support],
) -> tuple[tuple[Support, frozenset[Support]], ...]:
    unseen = set(supports)
    answer: list[tuple[Support, frozenset[Support]]] = []
    orbit_sizes: dict[int, int] = {}
    while unseen:
        seed = min(unseen, key=lambda value: tuple(sorted(value)))
        orbit = symmetry_orbit(seed)
        if len(orbit) not in (2, 4) or not orbit <= supports:
            raise AssertionError("unexpected (3,2,2,2)^2 symmetry orbit")
        unseen.difference_update(orbit)
        representative = min(orbit, key=lambda value: tuple(sorted(value)))
        answer.append((representative, orbit))
        orbit_sizes[len(orbit)] = orbit_sizes.get(len(orbit), 0) + 1
    answer.sort(key=lambda record: tuple(sorted(record[0])))
    if len(answer) != 230 or orbit_sizes != {2: 52, 4: 178}:
        raise AssertionError("unexpected (3,2,2,2)^2 orbit census")
    if set(SURVIVORS) - {representative for representative, _ in answer}:
        raise AssertionError("a recorded survivor is not an orbit representative")
    return tuple(answer)


def certify_survivor(
    support: tuple[tuple[int, int], ...],
    result: str,
) -> dict[str, object]:
    residuals, flag = SURVIVORS[frozenset(support)]
    if normalize_support(support) != support:
        raise AssertionError("survivor anchors differ from the recorded chart")
    h = sp.prod(residuals) ** -1
    certify_rur_box(result, (h, *residuals))

    symbols = RESIDUAL_SYMBOL_POOL[:7]
    substitution = dict(zip(symbols, residuals, strict=True))
    moments = [
        restricted_moment(order, support, symbols)
        for order in range(1, 15)
    ]
    if any(moment.as_expr().subs(substitution) != 0 for moment in moments):
        raise AssertionError("recorded survivor fails a pure moment through mu_14")

    matrix, polynomial = coefficient_data(support, residuals)
    if matrix.rank() != 3:
        raise AssertionError("a full-incidence survivor is not rank three")
    transformed = flag_transform(polynomial, flag)
    transformed_poly = sp.Poly(transformed, X, Y, domain=sp.QQ)
    transformed_support = [
        exponent
        for exponent, coefficient in transformed_poly.terms()
        if coefficient
    ]
    if not transformed_support or not all(
        row > column for row, column in transformed_support
    ):
        raise AssertionError("flag transform did not enter the positive chamber")

    coefficients = (sp.Rational(1), sp.Rational(1), *residuals)
    mixed_checks: list[dict[str, object]] = []
    low_nonzero = False
    for degree in (1, 2):
        for order in range(1, degree + 3):
            for dual_index in range(degree + 1):
                for coordinate_index in range(degree + 1):
                    value = mixed_contraction(
                        support,
                        coefficients,
                        order,
                        degree,
                        dual_index,
                        coordinate_index,
                    )
                    if order > degree and value != 0:
                        raise AssertionError("fixed-flag mixed cutoff failed")
                    if order <= degree and value != 0:
                        low_nonzero = True
                    mixed_checks.append(
                        {
                            "degree": degree,
                            "order": order,
                            "dual_index": dual_index,
                            "coordinate_index": coordinate_index,
                            "value": str(value),
                        }
                    )
    if not low_nonzero:
        raise AssertionError("survivor has no nonzero low-degree mixed value")

    return {
        "normalized_residual_coordinates": [str(value) for value in residuals],
        "coefficient_matrix": [
            [str(matrix[row, column]) for column in range(4)]
            for row in range(4)
        ],
        "coefficient_rank": 3,
        "dehomogenized_factorization": str(sp.factor(polynomial)),
        "flag_parameter": str(flag),
        "transformed_positive_weight_form": str(sp.factor(transformed)),
        "transformed_support": [list(exponent) for exponent in transformed_support],
        "pure_recurrence": "nu_(m+1)=0 for m>=0, with nu_1=0",
        "mixed_cutoff": "all bidegree-(e,e) multiplier contractions vanish for m>e",
        "mixed_checks": mixed_checks,
        "rur_sha256": hashlib.sha256(" ".join(result.split()).encode()).hexdigest(),
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
    closed = (
        rectangles
        | fringes
        | cross_support_class()[2]
        | regular_three_line_support_class()
        | three_line_432_support_class()
        | full_line_4311_support_class()
        | new_supports
    ) & mixed
    remaining = mixed - closed
    if len(closed) != 4370 or len(remaining) != 7050:
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
    if closed_orbits != {2: 79, 4: 1053}:
        raise AssertionError(f"unexpected closed orbit sizes: {closed_orbits}")
    if remaining_orbits != {2: 59, 4: 1733}:
        raise AssertionError(f"unexpected remaining orbit sizes: {remaining_orbits}")
    return {
        "mixed_support_count": len(mixed),
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
        rational_parametrization=True,
    )
    record: dict[str, object] = {
        "representative": [list(position) for position in support],
        "symmetry_orbit": [
            [list(position) for position in sorted(member)]
            for member in sorted(orbit, key=lambda value: tuple(sorted(value)))
        ],
        "seconds": result["seconds"],
    }
    if result["status"] == "excluded_on_coefficient_torus":
        record["classification"] = "unit ideal over QQ through mu_10"
        return record
    if result["status"] != "msolve_nonempty" or representative not in SURVIVORS:
        raise AssertionError(f"unexpected full-incidence result: {result}")
    rur = str(result["msolve"]["result_head"])
    record["classification"] = "unique rational fixed-flag rank-three point"
    record["survivor_certificate"] = certify_survivor(support, rur)
    return record


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
    with ThreadPoolExecutor(max_workers=arguments.workers) as executor:
        futures = [
            executor.submit(classify, representative, orbit, msolve, arguments.timeout)
            for representative, orbit in reps
        ]
        for future in as_completed(futures):
            records.append(future.result())
    records.sort(key=lambda record: record["representative"])
    unit_count = sum(
        record["classification"] == "unit ideal over QQ through mu_10"
        for record in records
    )
    survivor_count = len(records) - unit_count
    if unit_count != 228 or survivor_count != 2:
        raise AssertionError("unexpected full-incidence classification counts")

    artifact = {
        "format": "two-pair-sic-bidegree33-sparse-full-incidence32229-v1",
        "field": "characteristic zero",
        "support_class": (
            "mixed size-nine supports with row and column partitions 3,2,2,2"
        ),
        "support_size": 9,
        "total_support_count": len(supports),
        "symmetry_orbit_count": len(reps),
        "symmetry_orbit_sizes": [len(orbit) for _, orbit in reps],
        "moment_orders": [1, 10],
        "exact_unit_system_count": unit_count,
        "exact_fixed_flag_system_count": survivor_count,
        "component_method": (
            "exact QQ coefficient-torus localization and msolve RURs; the "
            "two nonunit systems each have one rational rank-three point"
        ),
        "relative_period": (
            "after the certified flag change every monomial has weight "
            "i-j>=1, so CT_u P(u,t)^m=0"
        ),
        "recurrence_certificate": (
            "over each survivor residue field the normalized pure sequence "
            "satisfies nu_(m+1)=0 with initial value nu_1=0"
        ),
        "mixed_multiplier_conclusion": (
            "nonzero low-degree mixed values occur, but weight gives the "
            "uniform cutoff m>e for every bidegree-(e,e) multiplier"
        ),
        "boundary_conclusion": (
            "every coordinate boundary has support at most eight and is "
            "SIC-safe by the complete support-eight theorem"
        ),
        "global_conclusion": (
            "all 816 full-incidence (3,2,2,2)^2 coordinate subspaces are SIC-safe"
        ),
        "updated_size_nine_census": updated_census(supports),
        "independent_formula_check": verify_restricted_formula(),
        "scope": (
            "complete exact classification of this 816-support class, not "
            "the full size-nine moment classification"
        ),
        "orbits": records,
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(artifact, indent=2) + "\n")

    print("PASS 230 symmetry orbits cover all 816 supports")
    print("PASS 228 dense coefficient tori are units over QQ through mu_10")
    print("PASS two unique residual points are fixed-flag rank-three nullcone points")
    print("PASS 4370 mixed size-nine supports closed; 7050 remain")


if __name__ == "__main__":
    main()
