#!/usr/bin/env python3
"""Classify the regular three-line size-nine supports in bidegree (3,3).

The class consists of the mixed supports with three occupied rows carrying
three entries each, together with the transposed class, after removing the
sixteen 3-by-3 rectangles.  It has 480 supports in 120 four-element orbits
under transpose and simultaneous reversal.

Exact QQ elimination through mu_10 makes 114 representative coefficient
tori empty.  The other six have one rational point each.  This checker
certifies their exact coordinates, rank two, and a contraction-preserving
flag change that puts every term in i>j.  They are therefore fixed-flag
one-sided-nullcone points, with the explicit mixed cutoff m>e for a
balanced multiplier of bidegree (e,e).
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from fractions import Fraction
import hashlib
from itertools import combinations
import json
from math import factorial
from pathlib import Path
import re
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


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_sparse_three_line9.json"
)
POSITIONS = tuple((row, column) for row in range(4) for column in range(4))
Support = frozenset[tuple[int, int]]
X, Y = sp.symbols("x y")


SURVIVORS: dict[
    Support,
    tuple[tuple[sp.Rational, ...], sp.Rational],
] = {
    frozenset(
        ((0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 3), (2, 0), (2, 2), (2, 3))
    ): (
        (
            sp.Rational(1, 4),
            sp.Rational(-8, 3),
            sp.Rational(-2),
            sp.Rational(1, 6),
            sp.Rational(4, 3),
            sp.Rational(-1),
            sp.Rational(-1, 3),
        ),
        sp.Rational(1, 2),
    ),
    frozenset(
        ((0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 3), (3, 1), (3, 2), (3, 3))
    ): (
        (
            sp.Rational(1, 4),
            sp.Rational(-2),
            sp.Rational(-3, 2),
            sp.Rational(1, 8),
            sp.Rational(-2),
            sp.Rational(-2),
            sp.Rational(-1, 2),
        ),
        sp.Rational(1, 2),
    ),
    frozenset(
        ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 0), (2, 2), (3, 0), (3, 1))
    ): (
        (
            sp.Rational(3, 16),
            sp.Rational(-6),
            sp.Rational(-3, 2),
            sp.Rational(-48),
            sp.Rational(3),
            sp.Rational(128),
            sp.Rational(32),
        ),
        sp.Rational(1, 4),
    ),
    frozenset(
        ((0, 0), (0, 1), (0, 3), (1, 0), (1, 2), (1, 3), (2, 1), (2, 2), (2, 3))
    ): (
        (
            sp.Rational(-4, 27),
            sp.Rational(-3, 2),
            sp.Rational(2),
            sp.Rational(8, 9),
            sp.Rational(-9, 4),
            sp.Rational(-3),
            sp.Rational(-1),
        ),
        sp.Rational(2, 3),
    ),
    frozenset(
        ((0, 1), (0, 2), (0, 3), (1, 0), (1, 2), (1, 3), (2, 0), (2, 1), (2, 3))
    ): (
        (
            sp.Rational(1, 4),
            sp.Rational(8, 3),
            sp.Rational(-2),
            sp.Rational(-2, 3),
            sp.Rational(-16, 3),
            sp.Rational(-4),
            sp.Rational(1, 3),
        ),
        sp.Rational(1, 2),
    ),
    frozenset(
        ((0, 1), (0, 2), (0, 3), (1, 0), (1, 2), (1, 3), (3, 0), (3, 1), (3, 2))
    ): (
        (
            sp.Rational(1, 4),
            sp.Rational(2),
            sp.Rational(-3, 2),
            sp.Rational(-1, 2),
            sp.Rational(-8),
            sp.Rational(-8),
            sp.Rational(-2),
        ),
        sp.Rational(1, 2),
    ),
}


def is_mixed(support: Support) -> bool:
    return (
        any(row > column for row, column in support)
        and any(row < column for row, column in support)
    )


def regular_on_axis(support: Support, axis: int) -> bool:
    counts = [
        sum(position[axis] == value for position in support)
        for value in range(4)
    ]
    return sorted(counts) == [0, 3, 3, 3]


def support_class() -> set[Support]:
    mixed = {
        frozenset(support)
        for support in combinations(POSITIONS, 9)
        if is_mixed(frozenset(support))
    }
    rectangles = {
        support
        for support in mixed
        if regular_on_axis(support, 0) and regular_on_axis(support, 1)
    }
    answer = {
        support
        for support in mixed
        if regular_on_axis(support, 0) or regular_on_axis(support, 1)
    } - rectangles
    if len(rectangles) != 16 or len(answer) != 480:
        raise AssertionError("unexpected regular three-line support census")
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
            raise AssertionError("unexpected regular three-line symmetry orbit")
        unseen.difference_update(orbit)
        representative = min(orbit, key=lambda value: tuple(sorted(value)))
        answer.append((representative, orbit))
    answer.sort(key=lambda record: tuple(sorted(record[0])))
    if len(answer) != 120:
        raise AssertionError("unexpected regular three-line orbit count")
    if set(SURVIVORS) - {representative for representative, _ in answer}:
        raise AssertionError("a recorded survivor is not an orbit representative")
    return tuple(answer)


def fraction_token(token: str) -> Fraction:
    token = token.strip()
    power_of_two = re.fullmatch(r"(-?\d+)\s*/\s*2\^(\d+)", token)
    if power_of_two:
        return Fraction(
            int(power_of_two.group(1)),
            2 ** int(power_of_two.group(2)),
        )
    rational = re.fullmatch(r"(-?\d+)\s*/\s*(\d+)", token)
    if rational:
        return Fraction(int(rational.group(1)), int(rational.group(2)))
    if re.fullmatch(r"-?\d+", token):
        return Fraction(int(token))
    raise ValueError(f"unrecognized msolve interval endpoint: {token!r}")


def certify_rur_box(result: str, expected: tuple[sp.Rational, ...]) -> None:
    header = re.match(r"\[0,\s*\[0,\s*(\d+),\s*(\d+),", result)
    if header is None:
        raise AssertionError("survivor result is not zero-dimensional")
    if re.search(
        r"\[1,\s*\[\s*\[1,\s*\[-?\d+,\s*-?\d+\]\]",
        result,
    ) is None:
        raise AssertionError("survivor RUR does not have a linear eliminant")
    pairs = re.findall(
        r"\[\s*([^,\[\]]+)\s*,\s*([^,\[\]]+)\s*\]",
        result,
    )
    if len(pairs) < len(expected):
        raise AssertionError("survivor RUR omitted isolating intervals")
    intervals = [
        (fraction_token(left), fraction_token(right))
        for left, right in pairs[-len(expected):]
    ]
    for value, (left, right) in zip(expected, intervals, strict=True):
        rational = Fraction(int(value.p), int(value.q))
        if not left <= rational <= right:
            raise AssertionError(
                f"expected coordinate {value} is outside [{left},{right}]"
            )


def coefficient_data(
    support: tuple[tuple[int, int], ...],
    residuals: tuple[sp.Rational, ...],
) -> tuple[sp.Matrix, sp.Expr]:
    coefficients = (sp.Integer(1), sp.Integer(1), *residuals)
    matrix = sp.zeros(4)
    polynomial = sp.Integer(0)
    for coefficient, (row, column) in zip(
        coefficients,
        support,
        strict=True,
    ):
        matrix[row, column] = coefficient
        polynomial += coefficient * X**row * Y**column
    return matrix, sp.expand(polynomial)


def flag_transform(polynomial: sp.Expr, flag: sp.Rational) -> sp.Expr:
    transformed = sp.cancel(
        (1 - flag * Y) ** 3
        * polynomial.subs({X: X + flag, Y: Y / (1 - flag * Y)})
    )
    return sp.expand(transformed)


def mixed_contraction(
    support: tuple[tuple[int, int], ...],
    coefficients: tuple[sp.Rational, ...],
    order: int,
    degree: int,
    dual_index: int,
    coordinate_index: int,
) -> sp.Rational:
    power: dict[tuple[int, int], sp.Rational] = {(0, 0): sp.Rational(1)}
    for _ in range(order):
        updated: dict[tuple[int, int], sp.Rational] = {}
        for (left, right), scalar in power.items():
            for (row, column), coefficient in zip(
                support,
                coefficients,
                strict=True,
            ):
                exponent = (left + row, right + column)
                updated[exponent] = (
                    updated.get(exponent, sp.Rational(0))
                    + scalar * coefficient
                )
        power = updated
    total_degree = 3 * order + degree
    answer = sp.Rational(0)
    for (row, column), coefficient in power.items():
        if row + dual_index != column + coordinate_index:
            continue
        diagonal = row + dual_index
        answer += (
            coefficient
            * factorial(diagonal)
            * factorial(total_degree - diagonal)
        )
    return answer


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
    moments = [restricted_moment(order, support, symbols) for order in range(1, 11)]
    if any(moment.as_expr().subs(substitution) != 0 for moment in moments):
        raise AssertionError("recorded survivor does not satisfy mu_1,...,mu_10")

    matrix, polynomial = coefficient_data(support, residuals)
    if matrix.rank() != 2:
        raise AssertionError("a regular three-line survivor is not rank two")
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
        "coefficient_rank": 2,
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
    crosses = cross_support_class()[2]
    closed = (rectangles | fringes | crosses | new_supports) & mixed
    remaining = mixed - closed
    if len(closed) != 1162 or len(remaining) != 10258:
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
    if closed_orbits != {2: 27, 4: 277}:
        raise AssertionError(f"unexpected closed orbit sizes: {closed_orbits}")
    if remaining_orbits != {2: 111, 4: 2509}:
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
        raise AssertionError(f"unexpected regular three-line result: {result}")
    rur = str(result["msolve"]["result_head"])
    record["classification"] = "unique rational fixed-flag rank-two point"
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
    if unit_count != 114 or survivor_count != 6:
        raise AssertionError("unexpected regular three-line classification counts")

    artifact = {
        "format": "two-pair-sic-bidegree33-sparse-three-line9-v1",
        "field": "characteristic zero",
        "support_class": (
            "three occupied rows with three entries each, or transpose, "
            "excluding the 3-by-3 rectangles"
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
            "six nonunit systems each have one rational rank-two point"
        ),
        "relative_period": (
            "after the certified flag change every monomial has u-weight "
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
            "all 480 regular three-line coordinate subspaces are SIC-safe"
        ),
        "updated_size_nine_census": updated_census(supports),
        "independent_formula_check": verify_restricted_formula(),
        "scope": (
            "complete exact classification of this 480-support class, not "
            "the full size-nine moment classification"
        ),
        "orbits": records,
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(artifact, indent=2) + "\n")

    print("PASS 120 four-element symmetry orbits cover all 480 supports")
    print("PASS 114 dense coefficient tori are units over QQ through mu_10")
    print("PASS six unique residual points are fixed-flag rank-two nullcone points")
    print("PASS 1162 mixed size-nine supports closed; 10258 remain")


if __name__ == "__main__":
    main()
