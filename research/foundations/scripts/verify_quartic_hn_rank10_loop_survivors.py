#!/usr/bin/env python3
"""Verify the frozen rank-ten Gale-loop census and trace sieve.

The external-catalogue extraction is frozen in
``artifacts/generated-results/quartic_hn_rank10_loop_survivors.json``.
This checker independently replays every returned rankline and colour,
constructs its normalized characteristic-zero realization ideal, verifies
the three unit ideals and rational realizations for the other 63
simplifications, and evaluates the one-loop first-two-trace self-square
scheme at one exact rational point of every realizable coloured type.

The one-loop calculation proves an excluded Zariski-open neighbourhood of
each displayed point.  It does not exclude special closed realization
strata.  Types with at least two loops are closed by the written
complementary-splitting and low-dimensional Hesse argument.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import random
from collections import Counter
from itertools import combinations
from pathlib import Path

import sympy as sp

from verify_quartic_hn_rank10_matroid_survivors import (
    bases_from_rankline,
    hyperplanes,
    normalized_realization,
    rank_table,
)


CENSUS_PATH = Path(
    "artifacts/generated-results/quartic_hn_rank10_loop_survivors.json"
)
CENSUS_SHA256 = "a7f5c04a05cdd6d32c723ae4ecfdfbb8c2018cf2072f00254f7c3a50c15e71b9"
EMPTY_Q_KEYS = {(7, 13), (8, 269), (9, 187398)}
SELF_SQUARE_SATURATION_CLOSED = {
    (6, 1, "111123"),
    (6, 10, "111123"),
    (6, 11, "111123"),
    (7, 4, "1111311"),
    (7, 6, "1111113"),
    (7, 8, "1111113"),
    (7, 10, "1111113"),
    (7, 15, "1111113"),
    (7, 15, "1131111"),
    (7, 15, "1221111"),
    (7, 19, "1111113"),
    (7, 19, "1111131"),
    (7, 26, "1111113"),
    (7, 26, "1111122"),
    (7, 26, "2111112"),
    (7, 54, "1111122"),
    (7, 72, "1111113"),
    (8, 586, "11121111"),
    (8, 586, "21111111"),
    (8, 589, "11111112"),
    (9, 188849, "111111111"),
}


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def actual_basis_set(matrix: sp.Matrix) -> set[int]:
    return {
        sum(1 << element for element in subset)
        for subset in combinations(range(matrix.cols), 4)
        if matrix[:, subset].det() != 0
    }


def assert_weighted_constraints(
    loops: int,
    n: int,
    bases: list[int],
    weights: tuple[int, ...],
) -> None:
    if loops + sum(weights) != 10:
        raise AssertionError((loops, weights))
    if max(weights) > 3 or sum(weight == 3 for weight in weights) > 1:
        raise AssertionError(weights)
    ranks = rank_table(n, bases)
    for plane in hyperplanes(n, ranks):
        cocircuit_size = sum(
            weights[element]
            for element in range(n)
            if not plane & (1 << element)
        )
        if cocircuit_size < 3:
            raise AssertionError((loops, weights, plane, cocircuit_size))

    full = (1 << n) - 1
    for basis in bases:
        complement_support = full
        for element in range(n):
            if basis & (1 << element) and weights[element] == 1:
                complement_support ^= 1 << element
        complement_rank = ranks[complement_support]
        has_coloop = any(
            weights[element] - bool(basis & (1 << element)) == 1
            and ranks[complement_support ^ (1 << element)] < complement_rank
            for element in range(n)
        )
        if not has_coloop:
            raise AssertionError((loops, weights, basis))


def linear_parametrization(
    variables: tuple[sp.Symbol, ...], equations: list[sp.Expr]
) -> tuple[dict[sp.Symbol, sp.Expr], tuple[sp.Symbol, ...], list[sp.Expr]]:
    substitutions: dict[sp.Symbol, sp.Expr] = {}
    remaining = [sp.factor(equation) for equation in equations if equation != 0]
    free = list(variables)
    while remaining:
        choices: list[tuple[int, sp.Symbol, sp.Expr, sp.Expr]] = []
        for equation in remaining:
            numerator = sp.together(equation).as_numer_denom()[0]
            for variable in free:
                if sp.Poly(numerator, variable).degree() != 1:
                    continue
                solutions = sp.solve(
                    numerator, variable, dict=False, simplify=False
                )
                if solutions:
                    choices.append(
                        (
                            sp.count_ops(solutions[0]),
                            variable,
                            solutions[0],
                            equation,
                        )
                    )
        if not choices:
            break
        _, variable, value, chosen = min(choices, key=lambda item: item[0])
        substitutions = {
            key: sp.cancel(expression.subs(variable, value))
            for key, expression in substitutions.items()
        }
        substitutions[variable] = sp.cancel(value)
        free.remove(variable)
        updated = []
        for equation in remaining:
            if equation == chosen:
                continue
            numerator = sp.together(
                equation.subs(variable, value)
            ).as_numer_denom()[0]
            numerator = sp.factor(numerator)
            if numerator != 0:
                updated.append(numerator)
        remaining = updated
    return substitutions, tuple(free), remaining


def special_rational_point(
    key: tuple[int, int],
    matrix: sp.Matrix,
    variables: tuple[sp.Symbol, ...],
) -> sp.Matrix | None:
    if key == (8, 111):
        values = (1, -2, -14, -7, -35, -7, -35)
        return matrix.subs(dict(zip(variables, values, strict=True)))
    if key == (8, 214):
        values = (1, -1, 7, -1, -7)
        return matrix.subs(dict(zip(variables, values, strict=True)))
    return None


def rational_realization(
    key: tuple[int, int],
    matrix: sp.Matrix,
    variables: tuple[sp.Symbol, ...],
    equations: list[sp.Expr],
    bases: list[int],
) -> sp.Matrix:
    special = special_rational_point(key, matrix, variables)
    if special is not None:
        if actual_basis_set(special) != set(bases):
            raise AssertionError((key, "bad special point"))
        return special

    substitutions, free, remaining = linear_parametrization(
        variables, equations
    )
    if remaining:
        raise AssertionError((key, "unresolved realization equations", remaining))
    generator = random.Random(1000 * matrix.cols + key[1])
    for _ in range(5000):
        values: dict[sp.Symbol, sp.Expr] = {
            variable: generator.randint(-30, 30) or 1 for variable in free
        }
        pending = dict(substitutions)
        for _ in range(len(pending) + 1):
            progress = False
            for variable, expression in list(pending.items()):
                value = sp.cancel(expression.subs(values))
                if not value.free_symbols:
                    values[variable] = value
                    del pending[variable]
                    progress = True
            if not progress:
                break
        if pending:
            continue
        if any(sp.together(equation.subs(values)) != 0 for equation in equations):
            continue
        candidate = matrix.subs(values)
        if actual_basis_set(candidate) == set(bases):
            return candidate
    raise AssertionError((key, "rational point search failed"))


def coordinate_square_data(
    gale: sp.Matrix,
) -> tuple[sp.Matrix, sp.Matrix]:
    kernel = sp.Matrix.hstack(*gale.nullspace())
    variables = sp.symbols(f"u0:{kernel.cols}")
    forms = kernel * sp.Matrix(variables)
    monomials: list[tuple[int, ...]] = []
    columns: list[dict[tuple[int, ...], sp.Expr]] = []
    for form in forms:
        terms = dict(sp.Poly(sp.expand(form**2), *variables).terms())
        columns.append(terms)
        for monomial in terms:
            if monomial not in monomials:
                monomials.append(monomial)
    matrix = sp.zeros(len(monomials), gale.cols)
    for column, terms in enumerate(columns):
        for monomial, coefficient in terms.items():
            matrix[monomials.index(monomial), column] = coefficient
    return kernel, matrix


def self_square_chart_is_empty(
    kernel: sp.Matrix, square_matrix: sp.Matrix
) -> bool:
    parameters = sp.symbols(f"t0:{kernel.cols}")
    active_row = kernel * sp.Matrix(parameters)
    _, pivot_rows = square_matrix.T.rref()
    independent = square_matrix[list(pivot_rows), :]
    equations = [
        sp.factor(equation)
        for equation in independent
        * sp.Matrix([entry**2 for entry in active_row])
        if equation != 0
    ]
    for parameter in parameters:
        basis = sp.groebner(
            equations + [parameter - 1], *parameters, order="grevlex"
        )
        if not (
            len(basis.polys) == 1
            and basis.polys[0].as_expr() == 1
        ):
            return False
    return True


def universal_self_square_is_empty(
    gale: sp.Matrix,
    realization_variables: tuple[sp.Symbol, ...],
    realization_equations: list[sp.Expr],
    open_minors: list[sp.Expr],
) -> bool:
    """Check the projective self-square incidence on the realization open set."""

    kernel = sp.Matrix.hstack(*gale.nullspace())
    if kernel.shape != (9, 5) or gale * kernel != sp.zeros(4, 5):
        raise AssertionError("bad symbolic Gale kernel")
    if any(
        sp.together(entry).as_numer_denom()[1] != 1
        for entry in kernel
    ):
        raise AssertionError("nonpolynomial symbolic Gale kernel")

    parameters = sp.symbols("qhnw10_t0:5")
    test_variables = sp.symbols("qhnw10_u0:5")
    active_row = kernel * sp.Matrix(parameters)
    kernel_forms = kernel * sp.Matrix(test_variables)
    quadratic = sp.Poly(
        sp.expand(
            sum(
                active_row[index] ** 2 * kernel_forms[index] ** 2
                for index in range(9)
            )
        ),
        *test_variables,
    )
    self_square_equations: list[sp.Expr] = []
    for coefficient in quadratic.coeffs():
        coefficient = sp.factor(coefficient)
        if coefficient != 0 and coefficient not in self_square_equations:
            self_square_equations.append(coefficient)

    open_factors: list[sp.Expr] = []
    for minor in open_minors:
        _, factorization = sp.factor_list(minor, *realization_variables)
        for factor, _ in factorization:
            factor = sp.Poly(
                factor, *realization_variables
            ).monic().as_expr()
            if factor not in open_factors:
                open_factors.append(factor)
    open_inverse = sp.Symbol("qhnw10_open_inverse")
    open_product = sp.prod(open_factors)
    if open_product == 0:
        raise AssertionError("zero realization-open product")
    for parameter in parameters:
        basis = sp.groebner(
            realization_equations
            + self_square_equations
            + [parameter - 1, open_inverse * open_product - 1],
            *parameters,
            open_inverse,
            *realization_variables,
            order="grevlex",
        )
        if not (
            len(basis.polys) == 1
            and basis.polys[0].as_expr() == 1
        ):
            return False
    return True


def has_matroidal_self_square_support(
    n: int, bases: list[int], weights: tuple[int, ...]
) -> bool:
    ranks = rank_table(n, bases)
    copies = [
        element
        for element, multiplicity in enumerate(weights)
        for _ in range(multiplicity)
    ]
    all_copies = (1 << len(copies)) - 1

    def direction_mask(copy_mask: int) -> int:
        result = 0
        for copy, element in enumerate(copies):
            if copy_mask & (1 << copy):
                result |= 1 << element
        return result

    for support in range(1, 1 << len(copies)):
        support_rank = ranks[direction_mask(support)]
        cyclic = all(
            ranks[direction_mask(support ^ (1 << copy))] == support_rank
            for copy in range(len(copies))
            if support & (1 << copy)
        )
        if not cyclic:
            continue
        projected_waring_rank = (
            support.bit_count()
            - 4
            + ranks[direction_mask(all_copies ^ support)]
        )
        if projected_waring_rank <= support.bit_count() // 2:
            return True
    return False


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--details", action="store_true")
    arguments = parser.parse_args()

    census_path = CENSUS_PATH.resolve()
    if file_sha256(census_path) != CENSUS_SHA256:
        raise AssertionError("stale loop census artifact")
    census = json.loads(census_path.read_text())
    if census["counts"] != {
        "catalogue_keys": 102,
        "coloured_isomorphism_types_by_loops": {
            "1": 47,
            "2": 49,
            "3": 18,
            "4": 1,
        },
        "labelled_survivors_by_loops": {
            "1": 126,
            "2": 84,
            "3": 24,
            "4": 1,
        },
        "underlying_simplifications": 66,
    }:
        raise AssertionError(census["counts"])

    schemes: dict[
        tuple[int, int],
        tuple[
            sp.Matrix,
            tuple[sp.Symbol, ...],
            list[sp.Expr],
            list[int],
            list[sp.Expr],
        ],
    ] = {}
    for record in census["records"]:
        n = record["simple_points"]
        key = (n, record["catalogue_index"])
        rankline = record["rankline"]
        bases = bases_from_rankline(n, rankline)
        for encoded in record["weight_representatives"]:
            assert_weighted_constraints(
                record["loops"],
                n,
                bases,
                tuple(int(value) for value in encoded),
            )
        if key not in schemes:
            matrix, variables, equations, open_minors = normalized_realization(
                n, bases
            )
            schemes[key] = (
                matrix,
                variables,
                equations,
                bases,
                open_minors,
            )

    rational_points: dict[tuple[int, int], sp.Matrix] = {}
    empty_keys: set[tuple[int, int]] = set()
    for key, (matrix, variables, equations, bases, _) in schemes.items():
        if key in EMPTY_Q_KEYS:
            basis = sp.groebner(equations, *variables, order="grevlex")
            if not (
                len(basis.polys) == 1
                and basis.polys[0].as_expr() == 1
            ):
                raise AssertionError((key, "expected unit ideal"))
            empty_keys.add(key)
            continue
        rational_points[key] = rational_realization(
            key, matrix, variables, equations, bases
        )

    if empty_keys != EMPTY_Q_KEYS or len(rational_points) != 63:
        raise AssertionError((empty_keys, len(rational_points)))

    qbar_counts: Counter[int] = Counter()
    unit_colour_counts: Counter[int] = Counter()
    one_loop_full_square = 0
    one_loop_empty_self_square = 0
    one_loop_matroid_support_closed = 0
    one_loop_special_candidates: set[tuple[int, int, str]] = set()
    for record in census["records"]:
        key = (record["simple_points"], record["catalogue_index"])
        colour_count = len(record["weight_representatives"])
        if key in empty_keys:
            unit_colour_counts[record["loops"]] += colour_count
            continue
        simplification = rational_points[key]
        for encoded in record["weight_representatives"]:
            qbar_counts[record["loops"]] += 1
            if record["loops"] != 1:
                continue
            weights = tuple(int(value) for value in encoded)
            if not has_matroidal_self_square_support(
                record["simple_points"],
                schemes[key][3],
                weights,
            ):
                one_loop_matroid_support_closed += 1
            else:
                one_loop_special_candidates.add((*key, encoded))
            gale = sp.Matrix.hstack(
                *[
                    simplification[:, element]
                    for element, multiplicity in enumerate(weights)
                    for _ in range(multiplicity)
                ]
            )
            kernel, square_matrix = coordinate_square_data(gale)
            if square_matrix.rank() == gale.cols:
                one_loop_full_square += 1
                status = "coordinate-square-full-rank"
            elif self_square_chart_is_empty(kernel, square_matrix):
                one_loop_empty_self_square += 1
                status = "projective-self-square-empty"
            else:
                raise AssertionError((key, encoded, "trace survivor"))
            if arguments.details:
                print(
                    "QHNW10_ONE_LOOP_TRACE_RECORD",
                    key,
                    encoded,
                    status,
                )

    if qbar_counts != Counter({1: 45, 2: 48, 3: 17, 4: 1}):
        raise AssertionError(qbar_counts)
    if unit_colour_counts != Counter({1: 2, 2: 1, 3: 1}):
        raise AssertionError(unit_colour_counts)
    if (one_loop_full_square, one_loop_empty_self_square) != (10, 35):
        raise AssertionError(
            (one_loop_full_square, one_loop_empty_self_square)
        )
    if one_loop_matroid_support_closed != 16:
        raise AssertionError(one_loop_matroid_support_closed)
    if len(one_loop_special_candidates) != 29:
        raise AssertionError(one_loop_special_candidates)
    if not SELF_SQUARE_SATURATION_CLOSED <= one_loop_special_candidates:
        raise AssertionError(
            SELF_SQUARE_SATURATION_CLOSED - one_loop_special_candidates
        )

    saturation_closed = 0
    records_by_key = {
        (record["simple_points"], record["catalogue_index"], record["loops"]): record
        for record in census["records"]
    }
    for n, catalogue_index, encoded in sorted(SELF_SQUARE_SATURATION_CLOSED):
        key = (n, catalogue_index)
        record = records_by_key[(*key, 1)]
        matrix, variables, equations, _, open_minors = schemes[key]
        weights = tuple(int(value) for value in encoded)
        gale = sp.Matrix.hstack(
            *[
                matrix[:, element]
                for element, multiplicity in enumerate(weights)
                for _ in range(multiplicity)
            ]
        )
        if not universal_self_square_is_empty(
            gale, variables, equations, open_minors
        ):
            raise AssertionError((key, encoded, "nonempty saturated incidence"))
        saturation_closed += 1
        if arguments.details:
            print(
                "QHNW10_ONE_LOOP_SATURATION_RECORD",
                key,
                encoded,
                "unit-in-all-five-projective-charts",
            )
    if saturation_closed != 21:
        raise AssertionError(saturation_closed)

    print("QHNW10_LOOP_ABSTRACT_COLOURED_TYPES=115")
    print("QHNW10_LOOP_QBAR_COLOURED_TYPES=111")
    print("QHNW10_LOOP_EMPTY_Q_IDEAL_TYPES=4")
    print("QHNW10_MULTI_LOOP_QBAR_TYPES_CLOSED=66")
    print("QHNW10_ONE_LOOP_QBAR_TYPES=45")
    print("QHNW10_ONE_LOOP_SQUARE_FULL_RANK=10")
    print("QHNW10_ONE_LOOP_PROJECTIVE_SELF_SQUARE_EMPTY=35")
    print("QHNW10_ONE_LOOP_MATROID_SUPPORT_CLOSED=16")
    print("QHNW10_ONE_LOOP_SATURATION_CLOSED=21")
    print("QHNW10_LOOP_QBAR_TYPES_UNIVERSALLY_CLOSED=103")
    print("QHNW10_ONE_LOOP_SPECIAL_LOCUS_TYPES=8")
    print("QHNW10_ONE_LOOP_GENERIC_TRACE_SURVIVORS=0")
    print("QUARTIC_HN_RANK10_LOOP_SURVIVOR_AUDIT_PASS")


if __name__ == "__main__":
    main()
