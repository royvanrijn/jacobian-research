#!/usr/bin/env python3
"""Classify every nine-entry 3-by-3 rectangle in bidegree (3,3).

A rectangle support is obtained by deleting one row and one column from
the 4-by-4 coefficient matrix.  There are sixteen such supports and six
orbits under transpose and simultaneous row/column reversal, both of which
preserve the contraction moments.

For one representative of each orbit, this checker works over QQ on the
dense coefficient torus.  It computes the exact finite moment scheme
through order fourteen and then localizes it at each of the nine 2-by-2
minors.  Every one of the resulting 54 systems is the unit ideal.  Hence
the reduced dense moment scheme has coefficient rank one.  The balanced
cubic rank-one theorem makes those points one-sided and SIC-safe.  Every
coordinate boundary has support at most eight and is covered by the
existing complete sparse classification.

This is an exact characteristic-zero component classification, not a
modular fiber screen and not a classification of all nine-entry supports.
"""

from __future__ import annotations

import hashlib
from itertools import combinations
import json
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

import sympy as sp

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
    / "two_pair_sic_bidegree33_sparse_rectangle9.json"
)
POSITIONS = tuple((i, j) for i in range(4) for j in range(4))
MOMENT_ORDERS = tuple(range(1, 15))


def rectangle_support(
    missing_row: int,
    missing_column: int,
) -> tuple[tuple[int, int], ...]:
    return tuple(
        (row, column)
        for row in range(4)
        for column in range(4)
        if row != missing_row and column != missing_column
    )


def symmetry_orbit(pair: tuple[int, int]) -> tuple[tuple[int, int], ...]:
    row, column = pair
    return tuple(sorted({
        (row, column),
        (column, row),
        (3 - row, 3 - column),
        (3 - column, 3 - row),
    }))


def missing_pair_orbits() -> tuple[tuple[tuple[int, int], ...], ...]:
    unseen = {(row, column) for row in range(4) for column in range(4)}
    answer = []
    while unseen:
        orbit = symmetry_orbit(min(unseen))
        unseen.difference_update(orbit)
        answer.append(orbit)
    return tuple(answer)


def normalized_coefficient_matrix(
    support: tuple[tuple[int, int], ...],
) -> tuple[
    tuple[tuple[int, int], ...],
    tuple[sp.Symbol, ...],
    sp.Matrix,
]:
    normalized = normalize_support(support)
    residuals = RESIDUAL_SYMBOL_POOL[:7]
    coefficients: dict[tuple[int, int], sp.Expr] = {
        normalized[0]: sp.Integer(1),
        normalized[1]: sp.Integer(1),
    }
    coefficients.update(
        {
            position: residual
            for position, residual in zip(
                normalized[2:],
                residuals,
                strict=True,
            )
        }
    )
    matrix = sp.Matrix(
        [
            [coefficients.get((row, column), 0) for column in range(4)]
            for row in range(4)
        ]
    )
    return normalized, residuals, matrix


def finite_scheme_profile(
    support: tuple[tuple[int, int], ...],
    msolve: str,
) -> tuple[int, int, str]:
    record = screen_support(
        support,
        through=14,
        timeout=120,
        msolve=msolve,
        threads=4,
        rational_parametrization=True,
    )
    if record["status"] != "msolve_nonempty":
        raise AssertionError(f"unexpected rectangle status: {record}")
    text = str(record["msolve"]["result_head"])
    header = re.match(r"\[0,\s*\[0,\s*(\d+),\s*(\d+),", text)
    eliminant = re.search(r"\[1,\s*\[\s*\[\s*(\d+),", text)
    if header is None or eliminant is None:
        raise AssertionError("could not parse the exact rectangle RUR")
    variable_count, scheme_degree = map(int, header.groups())
    if variable_count != 9:
        raise AssertionError("unexpected normalized rectangle variable count")
    compact = " ".join(text.split())
    return (
        scheme_degree,
        int(eliminant.group(1)),
        hashlib.sha256(compact.encode("utf-8")).hexdigest(),
    )


def minor_localization_is_unit(
    normalized: tuple[tuple[int, int], ...],
    residuals: tuple[sp.Symbol, ...],
    minor: sp.Expr,
    moments: list[sp.Poly],
    msolve: str,
) -> bool:
    h = sp.Symbol("h")
    inverse = sp.Symbol("rinv")
    saturation_monomial = sp.prod(residuals)
    inverse_equation = sp.Poly(
        sp.expand(inverse * minor - 1),
        inverse,
        *residuals,
        domain=sp.ZZ,
    )
    equations = [
        msolve_expression(
            sp.Poly(
                h * saturation_monomial - 1,
                h,
                inverse,
                *residuals,
                domain=sp.ZZ,
            )
        ),
        msolve_expression(inverse_equation),
        *(msolve_expression(moment) for moment in moments if not moment.is_zero),
    ]
    with tempfile.TemporaryDirectory(
        prefix="sic33-rectangle9-minor-"
    ) as directory:
        input_path = Path(directory) / "system.ms"
        output_path = Path(directory) / "result.ms"
        input_path.write_text(
            ",".join((str(h), str(inverse), *(str(value) for value in residuals)))
            + "\n0\n"
            + ",\n".join(equations)
            + "\n",
            encoding="utf-8",
        )
        completed = subprocess.run(
            [
                msolve,
                "-f",
                str(input_path),
                "-o",
                str(output_path),
                "-t",
                "4",
                "-l",
                "2",
                "-v",
                "0",
            ],
            text=True,
            capture_output=True,
            check=False,
            timeout=120,
        )
        result = (
            output_path.read_text(encoding="utf-8").strip()
            if output_path.exists()
            else ""
        )
    if completed.returncode != 0:
        raise AssertionError(
            f"msolve failed on {normalized}, minor {minor}: "
            f"{completed.stderr[-1000:]}"
        )
    return result in ("[-1]", "[-1]:")


def classify_representative(
    missing_pair: tuple[int, int],
    msolve: str,
) -> dict[str, object]:
    support = rectangle_support(*missing_pair)
    normalized, residuals, matrix = normalized_coefficient_matrix(support)
    moments = [
        restricted_moment(order, normalized, residuals)
        for order in MOMENT_ORDERS
    ]
    scheme_degree, reduced_degree, rur_hash = finite_scheme_profile(
        support,
        msolve,
    )

    active_rows = tuple(row for row in range(4) if row != missing_pair[0])
    active_columns = tuple(
        column for column in range(4) if column != missing_pair[1]
    )
    minor_records = []
    for rows in combinations(active_rows, 2):
        for columns in combinations(active_columns, 2):
            minor = sp.expand(matrix.extract(rows, columns).det())
            if not minor:
                raise AssertionError("an active rectangle minor vanished identically")
            is_unit = minor_localization_is_unit(
                normalized,
                residuals,
                minor,
                moments,
                msolve,
            )
            if not is_unit:
                raise AssertionError(
                    f"non-rank-one rectangle component at {missing_pair}, "
                    f"minor rows={rows}, columns={columns}"
                )
            minor_records.append({
                "rows": list(rows),
                "columns": list(columns),
                "minor": str(minor),
                "localized_moment_ideal": "unit over QQ",
            })

    if len(minor_records) != 9:
        raise AssertionError("a 3-by-3 rectangle must have nine 2-by-2 minors")
    return {
        "missing_row_column": list(missing_pair),
        "support": [list(position) for position in support],
        "normalized_anchors": [
            list(normalized[0]),
            list(normalized[1]),
        ],
        "moment_scheme_degree": scheme_degree,
        "squarefree_eliminant_degree": reduced_degree,
        "rur_head_sha256": rur_hash,
        "minor_localizations": minor_records,
        "all_minor_localizations_unit": True,
        "reduced_dense_scheme_rank": 1,
    }


def main() -> None:
    msolve = shutil.which("msolve")
    if msolve is None:
        raise RuntimeError("msolve is required")
    formula_check = verify_restricted_formula()

    orbits = missing_pair_orbits()
    if len(orbits) != 6 or sum(len(orbit) for orbit in orbits) != 16:
        raise AssertionError("unexpected rectangle symmetry orbit count")

    classifications = []
    for orbit in orbits:
        representative = min(orbit)
        classification = classify_representative(representative, msolve)
        classification["symmetry_orbit"] = [
            list(pair) for pair in orbit
        ]
        classifications.append(classification)

    artifact = {
        "format": "two-pair-sic-bidegree33-sparse-rectangle9-v1",
        "field": "characteristic zero",
        "support_class": (
            "nine-entry 3-by-3 rectangles in the standard 4-by-4 "
            "coefficient matrix"
        ),
        "rectangle_support_count": 16,
        "symmetry_group": [
            "coefficient-matrix transpose",
            "simultaneous row/column reversal",
        ],
        "symmetry_orbit_count": 6,
        "moment_orders": list(MOMENT_ORDERS),
        "component_method": (
            "exact QQ rational-univariate finite schemes, followed by "
            "Rabinowitsch localization at every 2-by-2 minor"
        ),
        "minor_localized_system_count": 54,
        "all_minor_localizations_unit": True,
        "dense_torus_conclusion": (
            "the reduced moment scheme on every rectangle has coefficient "
            "rank one and is SIC-safe by the balanced cubic rank-one theorem"
        ),
        "boundary_conclusion": (
            "every coordinate boundary has support at most eight and is "
            "SIC-safe by SIC2B33SP8"
        ),
        "global_conclusion_for_this_support_class": (
            "all sixteen nine-entry rectangle coordinate subspaces are "
            "SIC-safe"
        ),
        "scope": (
            "complete exact classification of one structured size-nine "
            "support class; not the full 11420-support size-nine census or "
            "the dense bidegree-(3,3) orbit classification"
        ),
        "independent_formula_check": formula_check,
        "orbits": classifications,
        "written_source": (
            "extended-geometry/TWO_PAIR_SIC_BIDEGREE33_FRONTIER.md"
        ),
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")

    print("PASS six exact symmetry orbits cover all 16 rectangle supports")
    print("PASS 54 characteristic-zero minor localizations are unit ideals")
    print("PASS every dense rectangle moment component has rank one")
    print("PASS all nine-entry 3-by-3 rectangle supports are SIC-safe")


if __name__ == "__main__":
    main()
