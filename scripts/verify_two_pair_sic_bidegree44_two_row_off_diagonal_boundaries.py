#!/usr/bin/env python3
"""Classify all coordinate boundaries of the off-diagonal two-row charts.

For every row pair r<s, start from

    {r} x ({0,...,4}-{r}) union {s} x ({0,...,4}-{s})

and take every proper coefficient support.  Simultaneous reversal acts on
this family.  Zero, fixed-flag one-sided, and rank-one coefficient tori
are certified without moment elimination.  A two-row torus whose row
column supports differ has exact rank two automatically.  When the row
column supports agree, the rank-one closed locus is SIC-safe and every
2-by-2-minor open is checked separately.  All remaining systems use
exact QQ moments and coefficient-torus Rabinowitsch localization through
mu_10.  A nonunit result is rerun with an exact rational parametrization
and aborts the theorem artifact for component analysis.

No ambient rank-at-most-two determinantal ideal is used.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from hashlib import sha256
from itertools import combinations
import json
from math import factorial, gcd
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

import sympy as sp

from explore_two_pair_sic_bidegree44_two_row_off_diagonal import (
    compositions,
    support_for,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree44_two_row_off_diagonal_boundaries.json"
)
Support = frozenset[tuple[int, int]]


@dataclass(frozen=True)
class NormalizedChart:
    support: tuple[tuple[int, int], ...]
    anchors: tuple[tuple[int, int], tuple[int, int]]
    residuals: tuple[sp.Symbol, ...]
    coefficients: tuple[sp.Expr, ...]


def reversal(support: Support) -> Support:
    return frozenset((4 - row, 4 - column) for row, column in support)


def support_key(support: Support) -> tuple[tuple[int, int], ...]:
    return tuple(sorted(support))


def canonical_support(support: Support) -> tuple[tuple[int, int], ...]:
    return min(support_key(support), support_key(reversal(support)))


def boundary_orbits() -> tuple[Support, ...]:
    representatives: set[tuple[tuple[int, int], ...]] = set()
    all_boundaries: set[Support] = set()
    for rows in combinations(range(5), 2):
        full = support_for(rows)
        for size in range(len(full)):
            for subset in combinations(full, size):
                support = frozenset(subset)
                all_boundaries.add(support)
                representatives.add(canonical_support(support))
    answer = tuple(
        frozenset(representative)
        for representative in sorted(
            representatives,
            key=lambda value: (len(value), value),
        )
    )
    covered: set[Support] = set()
    for representative in answer:
        covered.add(representative)
        covered.add(reversal(representative))
    if covered != all_boundaries:
        raise AssertionError("boundary reversal orbits do not cover the family")
    if len(answer) != 1174:
        raise AssertionError("unexpected boundary orbit count")
    return answer


def weights(support: Support) -> tuple[int, ...]:
    return tuple(row - column for row, column in support)


def row_column_sets(support: Support) -> tuple[frozenset[int], ...]:
    rows = sorted({row for row, _ in support})
    return tuple(
        frozenset(column for row_value, column in support if row_value == row)
        for row in rows
    )


def stratum(support: Support) -> str:
    if not support:
        return "rank_zero"
    support_weights = weights(support)
    if all(weight > 0 for weight in support_weights) or all(
        weight < 0 for weight in support_weights
    ):
        return "fixed_flag_one_sided"
    rows = {row for row, _ in support}
    columns = {column for _, column in support}
    if len(rows) == 1 or len(columns) == 1:
        return "rank_one"
    column_sets = row_column_sets(support)
    if len(column_sets) != 2:
        raise AssertionError("a boundary uses more than two rows")
    if column_sets[0] == column_sets[1]:
        return "rank_one_closed_and_exact_rank_two_open"
    return "exact_rank_two"


def normalize(support: Support) -> NormalizedChart:
    negative = sorted(
        (position for position in support if position[0] - position[1] < 0),
        key=lambda position: (position[0] - position[1], position),
    )
    positive = sorted(
        (position for position in support if position[0] - position[1] > 0),
        key=lambda position: (position[0] - position[1], position),
    )
    if not negative or not positive:
        raise ValueError("only mixed supports have quotient charts")
    anchors = (negative[0], positive[-1])
    tail = tuple(sorted(support - frozenset(anchors)))
    ordered = (*anchors, *tail)
    residuals = sp.symbols(f"z0:{len(tail)}")
    coefficients = (sp.Integer(1), sp.Integer(1), *residuals)
    return NormalizedChart(ordered, anchors, residuals, coefficients)


def primitive_moment(chart: NormalizedChart, order: int) -> sp.Poly | int:
    terms: dict[tuple[int, ...], int] = {}
    order_factorial = factorial(order)
    for counts in compositions(order, len(chart.support)):
        row_degree = sum(
            count * position[0]
            for count, position in zip(counts, chart.support, strict=True)
        )
        column_degree = sum(
            count * position[1]
            for count, position in zip(counts, chart.support, strict=True)
        )
        if row_degree != column_degree:
            continue
        denominator = 1
        for count in counts:
            denominator *= factorial(count)
        scalar = (
            order_factorial
            // denominator
            * factorial(row_degree)
            * factorial(4 * order - row_degree)
        )
        exponent = counts[2:]
        terms[exponent] = terms.get(exponent, 0) + scalar

    content = 0
    for scalar in terms.values():
        content = gcd(content, abs(scalar))
    if content:
        terms = {exponent: scalar // content for exponent, scalar in terms.items()}
    if not chart.residuals:
        return terms.get((), 0)
    return sp.Poly.from_dict(terms, chart.residuals, domain=sp.ZZ)


def coefficient_map(chart: NormalizedChart) -> dict[tuple[int, int], sp.Expr]:
    return dict(zip(chart.support, chart.coefficients, strict=True))


def row_minors(chart: NormalizedChart) -> tuple[tuple[tuple[int, int], sp.Expr], ...]:
    support = frozenset(chart.support)
    rows = sorted({row for row, _ in support})
    column_sets = row_column_sets(support)
    if len(rows) != 2 or column_sets[0] != column_sets[1]:
        return ()
    coefficient = coefficient_map(chart)
    answer = []
    for left, right in combinations(sorted(column_sets[0]), 2):
        minor = sp.expand(
            coefficient[(rows[0], left)] * coefficient[(rows[1], right)]
            - coefficient[(rows[0], right)] * coefficient[(rows[1], left)]
        )
        answer.append(((left, right), minor))
    return tuple(answer)


def msolve_result(
    chart: NormalizedChart,
    moments: tuple[sp.Poly | int, ...],
    timeout: int,
    *,
    invert: sp.Expr | None,
    rational_parametrization: bool,
) -> dict[str, object]:
    nonzero_constants = [
        int(moment)
        for moment in moments
        if isinstance(moment, int) and moment
    ]
    if not chart.residuals:
        if not nonzero_constants:
            raise AssertionError("a zero-variable mixed chart has no nonzero moment")
        return {
            "returncode": 0,
            "status": "unit_ideal",
            "method": "nonzero normalized constant moment",
            "first_nonzero_order": next(
                order
                for order, moment in enumerate(moments, start=1)
                if isinstance(moment, int) and moment
            ),
            "result_sha256": sha256(b"constant-unit").hexdigest(),
            "result_bytes": 0,
        }
    if nonzero_constants:
        raise AssertionError("a positive-variable chart has a constant unit moment")

    msolve = shutil.which("msolve")
    if msolve is None:
        raise RuntimeError("msolve is required")
    h = sp.Symbol("h")
    localizer = sp.prod(chart.residuals)
    if invert is not None:
        localizer *= invert
    equations = [
        sp.expand(h * localizer - 1),
        *(
            moment.as_expr()
            for moment in moments
            if isinstance(moment, sp.Poly) and not moment.is_zero
        ),
    ]
    with tempfile.TemporaryDirectory(prefix="sic44-two-row-offdiag-boundary-") as directory:
        input_path = Path(directory) / "system.ms"
        output_path = Path(directory) / "result.ms"
        input_path.write_text(
            ",".join(str(value) for value in (h, *chart.residuals))
            + "\n0\n"
            + ",\n".join(str(value).replace("**", "^") for value in equations)
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
                "1",
                "-l",
                "2",
                "-v",
                "0",
                *(["-P", "1"] if rational_parametrization else []),
            ],
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
        result = (
            output_path.read_text(encoding="utf-8").strip()
            if output_path.exists()
            else ""
        )
    compact = " ".join(result.split())
    header = re.match(r"\[0,\s*\[0,\s*(\d+),\s*(\d+),", result)
    if result in ("[-1]", "[-1]:"):
        status = "unit_ideal"
    elif re.match(r"\[0,\s*\[1,", result):
        status = "nonempty_approximate"
    elif result.startswith("[1,"):
        status = "positive_dimensional"
    elif header is not None:
        status = "zero_dimensional"
    else:
        status = "unparsed"
    return {
        "returncode": completed.returncode,
        "status": status,
        "variables": int(header.group(1)) if header is not None else None,
        "scheme_degree": int(header.group(2)) if header is not None else None,
        "result_sha256": sha256(compact.encode()).hexdigest(),
        "result_bytes": len(result.encode()),
        "result_head": result[:3000],
        "stderr_tail": completed.stderr[-1000:],
    }


def compact_result(result: dict[str, object]) -> dict[str, object]:
    return {
        key: value
        for key, value in result.items()
        if key not in {"result_head", "stderr_tail"}
    }


def polynomial_metadata(moment: sp.Poly | int) -> dict[str, object]:
    expression = moment.as_expr() if isinstance(moment, sp.Poly) else sp.Integer(moment)
    text = str(expression)
    return {
        "is_zero": expression == 0,
        "total_degree": 0 if expression == 0 else sp.Poly(expression).total_degree(),
        "term_count": len(sp.Poly(expression).terms()) if expression != 0 else 0,
        "sha256": sha256(text.encode()).hexdigest(),
    }


def solve_open(
    chart: NormalizedChart,
    moments: tuple[sp.Poly | int, ...],
    timeout: int,
    invert: sp.Expr | None,
) -> dict[str, object]:
    result = msolve_result(
        chart,
        moments,
        timeout,
        invert=invert,
        rational_parametrization=False,
    )
    if result["returncode"] != 0 or result["status"] == "unparsed":
        raise RuntimeError(f"msolve failed: {result}")
    if result["status"] != "unit_ideal":
        rur = msolve_result(
            chart,
            moments,
            timeout,
            invert=invert,
            rational_parametrization=True,
        )
        raise AssertionError(
            "nonunit boundary component requires analysis: "
            + json.dumps(rur, sort_keys=True)
        )
    return compact_result(result)


def classify_support(support: Support, timeout: int) -> dict[str, object]:
    support_stratum = stratum(support)
    base = {
        "representative": [list(position) for position in sorted(support)],
        "reversal": [list(position) for position in sorted(reversal(support))],
        "support_size": len(support),
        "stratum": support_stratum,
    }
    if support_stratum in {"rank_zero", "fixed_flag_one_sided", "rank_one"}:
        base["classification"] = {
            "rank_zero": "zero form is SIC-safe",
            "fixed_flag_one_sided": "fixed-standard-flag one-sided and SIC-safe",
            "rank_one": "coefficient rank one and SIC-safe by split-symbol theorem",
        }[support_stratum]
        return base

    chart = normalize(support)
    moments = tuple(primitive_moment(chart, order) for order in range(1, 11))
    base["normalization_anchors"] = [list(position) for position in chart.anchors]
    base["residual_coordinate_count"] = len(chart.residuals)
    base["moment_term_counts"] = [
        0 if isinstance(moment, int) and not moment
        else 1 if isinstance(moment, int)
        else 0 if moment.is_zero
        else len(moment.terms())
        for moment in moments
    ]

    if support_stratum == "exact_rank_two":
        base["classification"] = "exact-rank-two coefficient torus is unit through mu_10"
        try:
            base["unit_certificate"] = solve_open(chart, moments, timeout, None)
        except AssertionError as error:
            raise AssertionError(
                f"support {sorted(support)}: {error}"
            ) from error
        return base

    minors = row_minors(chart)
    if not minors:
        raise AssertionError("a mixed-rank support has no minor opens")
    base["classification"] = (
        "rank-one closed locus is SIC-safe; every exact-rank-two minor open "
        "is unit through mu_10"
    )
    minor_certificates = []
    for columns, minor in minors:
        try:
            certificate = solve_open(chart, moments, timeout, minor)
        except AssertionError as error:
            raise AssertionError(
                f"support {sorted(support)}, minor columns {columns}: {error}"
            ) from error
        minor_certificates.append(
            {
                "columns": list(columns),
                "minor": str(minor),
                "certificate": certificate,
            }
        )
    base["minor_open_certificates"] = minor_certificates
    return base


def delayed_fibre_audit(timeout: int) -> list[dict[str, object]]:
    specifications = (
        {
            "support": frozenset({(0, 1), (0, 4), (2, 0), (2, 3)}),
            "invert_columns": None,
            "survives_through": 8,
            "killed_by": 9,
            "scheme_degree": 8,
        },
        {
            "support": frozenset(
                {(1, 0), (1, 2), (1, 4), (3, 0), (3, 2), (3, 4)}
            ),
            "invert_columns": (0, 2),
            "survives_through": 9,
            "killed_by": 10,
            "scheme_degree": 252,
        },
    )
    answer = []
    for specification in specifications:
        support = specification["support"]
        chart = normalize(support)
        invert_columns = specification["invert_columns"]
        invert = None
        if invert_columns is not None:
            invert = dict(row_minors(chart))[invert_columns]
        moments = tuple(
            primitive_moment(chart, order)
            for order in range(1, int(specification["killed_by"]) + 1)
        )
        surviving = msolve_result(
            chart,
            moments[: int(specification["survives_through"])],
            timeout,
            invert=invert,
            rational_parametrization=True,
        )
        killed = msolve_result(
            chart,
            moments,
            timeout,
            invert=invert,
            rational_parametrization=False,
        )
        if (
            surviving["returncode"] != 0
            or surviving["status"] != "zero_dimensional"
            or surviving["scheme_degree"] != specification["scheme_degree"]
        ):
            raise AssertionError(f"unexpected delayed fibre: {surviving}")
        if killed["returncode"] != 0 or killed["status"] != "unit_ideal":
            raise AssertionError(f"delayed fibre was not killed: {killed}")
        answer.append(
            {
                "support": [list(position) for position in sorted(support)],
                "normalization_anchors": [
                    list(position) for position in chart.anchors
                ],
                "localized_minor_columns": (
                    list(invert_columns) if invert_columns is not None else None
                ),
                "survives_through": specification["survives_through"],
                "exact_scheme_degree": specification["scheme_degree"],
                "rational_parametrization_certificate": compact_result(surviving),
                "killed_by": specification["killed_by"],
                "killing_moment": polynomial_metadata(moments[-1]),
                "unit_certificate": compact_result(killed),
            }
        )
    return answer


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()
    if arguments.workers < 1 or arguments.timeout < 1:
        raise ValueError("--workers and --timeout must be positive")

    supports = boundary_orbits()
    expected_strata = {
        "rank_zero": 1,
        "fixed_flag_one_sided": 222,
        "rank_one": 19,
        "rank_one_closed_and_exact_rank_two_open": 16,
        "exact_rank_two": 916,
    }
    census: dict[str, int] = {}
    for support in supports:
        name = stratum(support)
        census[name] = census.get(name, 0) + 1
    if census != expected_strata:
        raise AssertionError(f"unexpected boundary stratum census: {census}")

    records: list[dict[str, object]] = []
    with ProcessPoolExecutor(max_workers=arguments.workers) as executor:
        futures = [
            executor.submit(classify_support, support, arguments.timeout)
            for support in supports
        ]
        for completed, future in enumerate(as_completed(futures), start=1):
            records.append(future.result())
            if completed % 64 == 0 or completed == len(futures):
                print(f"classified {completed}/{len(futures)}", flush=True)
    records.sort(key=lambda record: (record["support_size"], record["representative"]))

    unit_system_count = sum(
        record["stratum"] == "exact_rank_two"
        for record in records
    ) + sum(
        len(record.get("minor_open_certificates", []))
        for record in records
    )
    delayed_fibres = delayed_fibre_audit(arguments.timeout)
    artifact = {
        "format": "two-pair-sic-bidegree44-two-row-off-diagonal-boundaries-v1",
        "field": "characteristic zero",
        "full_chart_family": "ten off-diagonal two-row supports",
        "proper_boundary_reversal_orbit_count": len(supports),
        "stratum_census": census,
        "moment_orders": [1, 10],
        "exact_QQ_unit_system_count": unit_system_count,
        "delayed_fibres": delayed_fibres,
        "rank_method": (
            "support combinatorics, with every 2-by-2 minor open localized "
            "separately when rank one and rank two share a coefficient torus"
        ),
        "global_conclusion": (
            "every coordinate boundary of all ten off-diagonal two-row charts "
            "is SIC-safe"
        ),
        "scope": (
            "complete coordinate-boundary classification for this chart family; "
            "non-coordinate U charts remain open"
        ),
        "orbits": records,
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")

    print("PASS 1174 reversal orbits cover every proper coordinate boundary")
    print(f"PASS {unit_system_count} exact QQ rank-two opens are units through mu_10")
    print("PASS all ten complete off-diagonal two-row coordinate subspaces are SIC-safe")


if __name__ == "__main__":
    main()
