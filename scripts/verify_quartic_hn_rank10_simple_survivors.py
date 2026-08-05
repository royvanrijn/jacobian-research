#!/usr/bin/env python3
"""Verify the five simple OP-QHNW10 survivors and close their Gram branch.

The companion enumerator constructs the complete catalogue-relative abstract
census.  This script freezes that census, independently rechecks every
ten-point matroid, constructs its normalized realization scheme over QQ,
exhibits an exact rational point, and verifies the uniform five-point-line
trace obstruction.

For the distinguished rank-two Gale flat F of size six, ``ker(K_F)`` has
dimension four.  Its six coordinate squares and fifteen off-diagonal products
of coordinate squares are independent on the realization open set.  The
first two Hessian traces therefore make the Gram block G_FF zero.  Gale
duality says that the six corresponding Waring vectors span at least a
four-plane.  This contradicts Witt index three in a nondegenerate six-space
and excludes every rank-six Gram matrix.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import sys
from pathlib import Path

import sympy as sp

from enumerate_quartic_hn_rank10_simple_survivors import assert_simple_survivor
from verify_quartic_hn_rank10_matroid_survivors import (
    normalized_realization,
    rank_table,
)


CENSUS_PATH = Path(
    "artifacts/generated-results/quartic_hn_rank10_simple_survivors.json"
)
CENSUS_SHA256 = "9304dff4987e85e9e26694a62ea38fd5549441689cb5a297b3659f5616aac725"
EXPECTED_BASIS_COUNTS = (114, 110, 108, 105, 90)
EXPECTED_SCHEME_SIZES = (
    (7, 1, 73),
    (5, 0, 41),
    (5, 0, 50),
    (4, 0, 38),
    (4, 0, 40),
)
RATIONAL_VALUES = (
    (-2, -1, 2, 2, 3, -2, -1),
    (-2, -2, -1, 2, -1),
    (-2, -1, 2, 3, -2),
    (-2, -2, -1, 2),
    (-2, -1, 2, -2),
)
SIX_POINT_FLAT = (0, 1, 2, 3, 4, 9)
SQUARE_MINOR_ROWS = (0, 1, 2, 3, 5, 9)
SQUARE_PAIR_MINOR_ROWS = (0, 1, 2, 3, 4, 5, 9, 10, 11, 12, 25, 26, 27, 28, 30)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coefficient_matrix(
    forms: list[sp.Expr], variables: tuple[sp.Symbol, ...], degree: int
) -> sp.Matrix:
    monomials = [
        exponents
        for exponents in itertools.product(range(degree + 1), repeat=len(variables))
        if sum(exponents) == degree
    ]
    result = sp.zeros(len(monomials), len(forms))
    for column, form in enumerate(forms):
        polynomial = sp.Poly(sp.expand(form), *variables)
        for row, monomial in enumerate(monomials):
            result[row, column] = polynomial.coeff_monomial(monomial)
    return result


def normal_form(
    expression: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    equations: list[sp.Expr],
) -> sp.Expr:
    expression = sp.cancel(expression)
    if not equations:
        return expression
    basis = sp.groebner(equations, *variables, order="grevlex", domain=sp.QQ)
    numerator, denominator = sp.fraction(expression)
    numerator = basis.reduce(sp.expand(numerator))[1]
    denominator = basis.reduce(sp.expand(denominator))[1]
    return sp.cancel(numerator / denominator)


def factor_keys(
    expression: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    equations: list[sp.Expr],
) -> set[str]:
    expression = normal_form(expression, variables, equations)
    numerator, denominator = sp.fraction(sp.together(expression))
    result: set[str] = set()
    for value in (numerator, denominator):
        _, factors = sp.factor_list(value, *variables)
        for factor, _ in factors:
            polynomial = sp.Poly(factor, *variables, domain=sp.QQ)
            if polynomial.total_degree() == 0:
                continue
            monic = sp.expand(polynomial.as_expr() / polynomial.LC())
            result.add(sp.sstr(monic))
    return result


def open_factor_keys(
    open_minors: list[sp.Expr],
    variables: tuple[sp.Symbol, ...],
    equations: list[sp.Expr],
) -> set[str]:
    result: set[str] = set()
    for minor in open_minors:
        result.update(factor_keys(minor, variables, equations))
    return result


def actual_bases(matrix: sp.Matrix) -> set[int]:
    n = matrix.cols
    return {
        sum(1 << element for element in subset)
        for subset in itertools.combinations(range(n), 4)
        if matrix[:, subset].det() != 0
    }


def matrix_as_lists(matrix: sp.Matrix) -> list[list[str]]:
    return [
        [sp.sstr(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def verify_record(
    record: dict[str, object], rational_values: tuple[int, ...]
) -> dict[str, object]:
    survivor_index = int(record["survivor_index"])
    bases = [int(value) for value in record["bases"]]
    assert_simple_survivor(10, bases)
    ranks = rank_table(10, bases)
    six_flat = sum(1 << element for element in SIX_POINT_FLAT)
    complement = ((1 << 10) - 1) ^ six_flat
    waring_flat_rank = 6 - 4 + ranks[complement]
    if ranks[six_flat] != 2 or waring_flat_rank <= 3:
        raise AssertionError((survivor_index, ranks[six_flat], ranks[complement]))
    if any(
        ranks[six_flat | (1 << element)] == 2
        for element in range(10)
        if element not in SIX_POINT_FLAT
    ):
        raise AssertionError((survivor_index, "six-set is not a flat"))

    gale, variables, equations, open_minors = normalized_realization(10, bases)
    expected_sizes = EXPECTED_SCHEME_SIZES[survivor_index]
    if (len(variables), len(equations), len(open_minors)) != expected_sizes:
        raise AssertionError(
            (
                survivor_index,
                (len(variables), len(equations), len(open_minors)),
                expected_sizes,
            )
        )
    substitution = dict(zip(variables, rational_values, strict=True))
    if any(equation.subs(substitution) != 0 for equation in equations):
        raise AssertionError((survivor_index, "rational point misses ideal"))
    if any(minor.subs(substitution) == 0 for minor in open_minors):
        raise AssertionError((survivor_index, "rational point misses open set"))
    rational_gale = gale.subs(substitution)
    if actual_bases(rational_gale) != set(bases):
        raise AssertionError((survivor_index, "wrong rational matroid"))

    restricted_gale = gale[:, SIX_POINT_FLAT]
    kernel_columns = restricted_gale.nullspace()
    if len(kernel_columns) != 4:
        raise AssertionError((survivor_index, len(kernel_columns)))
    value_matrix = sp.Matrix.hstack(*kernel_columns).T
    if restricted_gale * value_matrix.T != sp.zeros(4, 4):
        raise AssertionError((survivor_index, "bad restricted kernel"))
    test_variables = sp.symbols("qhnw10_u0:4")
    values = value_matrix.T * sp.Matrix(test_variables)

    square_matrix = coefficient_matrix(
        [values[index] ** 2 for index in range(6)], test_variables, 2
    )
    square_pair_matrix = coefficient_matrix(
        [
            values[left] ** 2 * values[right] ** 2
            for left, right in itertools.combinations(range(6), 2)
        ],
        test_variables,
        4,
    )
    square_determinant = sp.factor(
        square_matrix.extract(SQUARE_MINOR_ROWS, range(6)).det()
    )
    square_pair_determinant = sp.factor(
        square_pair_matrix.extract(SQUARE_PAIR_MINOR_ROWS, range(15)).det()
    )
    if square_determinant.subs(substitution) == 0:
        raise AssertionError((survivor_index, "bad square minor"))
    if square_pair_determinant.subs(substitution) == 0:
        raise AssertionError((survivor_index, "bad square-pair minor"))

    allowed_factors = open_factor_keys(open_minors, variables, equations)
    square_factors = factor_keys(square_determinant, variables, equations)
    pair_factors = factor_keys(square_pair_determinant, variables, equations)
    if not square_factors <= allowed_factors:
        raise AssertionError((survivor_index, square_factors - allowed_factors))
    if not pair_factors <= allowed_factors:
        raise AssertionError((survivor_index, pair_factors - allowed_factors))

    return {
        "survivor_index": survivor_index,
        "basis_count": len(bases),
        "normalized_gale_matrix": matrix_as_lists(gale),
        "realization_variables": [str(variable) for variable in variables],
        "realization_ideal_generators": [sp.sstr(value) for value in equations],
        "open_minor_count": len(open_minors),
        "rational_values": list(rational_values),
        "rational_gale_matrix": matrix_as_lists(rational_gale),
        "six_point_flat": list(SIX_POINT_FLAT),
        "six_point_complement_rank": ranks[complement],
        "six_point_waring_span_rank": waring_flat_rank,
        "square_minor": sp.sstr(square_determinant),
        "square_pair_minor": sp.sstr(square_pair_determinant),
        "trace_status": "first-two-traces-force-zero-six-by-six-Gram-block",
        "rank_six_gram_exists": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--details",
        action="store_true",
        help="print exact realization and trace records as JSON lines",
    )
    arguments = parser.parse_args()
    census_path = CENSUS_PATH.resolve()
    if file_sha256(census_path) != CENSUS_SHA256:
        raise AssertionError("stale simple-survivor census artifact")
    census = json.loads(census_path.read_text())
    if census["counts"]["abstract_isomorphism_types"] != 5:
        raise AssertionError(census["counts"])
    if census["counts"]["labelled_extension_models"] != 23:
        raise AssertionError(census["counts"])
    records = census["records"]
    if tuple(record["basis_count"] for record in records) != EXPECTED_BASIS_COUNTS:
        raise AssertionError("unexpected simple survivor records")

    details = [
        verify_record(record, rational_values)
        for record, rational_values in zip(records, RATIONAL_VALUES, strict=True)
    ]
    if arguments.details:
        for detail in details:
            print("QHNW10_SIMPLE_REALIZATION_RECORD=" + json.dumps(detail, sort_keys=True))
    print(f"QHNW10_SIMPLE_CENSUS_SHA256={CENSUS_SHA256}")
    print("QHNW10_SIMPLE_RESIDUAL_DELETIONS=13")
    print("QHNW10_SIMPLE_LABELLED_EXTENSIONS=23")
    print("QHNW10_SIMPLE_ABSTRACT_SURVIVORS=5")
    print("QHNW10_SIMPLE_QBAR_SURVIVORS=5")
    print("QHNW10_SIMPLE_SIX_POINT_TRACE_CLOSED=5")
    print("QHNW10_SIMPLE_RANK6_GRAM_SURVIVORS=0")
    print("QUARTIC_HN_RANK10_SIMPLE_SURVIVOR_AUDIT_PASS")


if __name__ == "__main__":
    main()
