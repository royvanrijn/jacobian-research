#!/usr/bin/env python3
"""Probe every off-diagonal two-row rank-two quartic chart orbit.

For rows r<s, use the eight-position support

    {r} x ({0,...,4} - {r})  union  {s} x ({0,...,4} - {s}).

Fixing U=(e_r,e_s) removes the internal GL_2 gauge.  On the dense
coefficient torus the minor in columns r,s is nonzero, so every point has
coefficient rank exactly two.  Overall scaling and the diagonal torus
normalize two coefficients of distinct weights, leaving six residual
coordinates.  The script asks msolve whether the exact QQ moment ideal
through a requested order is the unit ideal.  Any nonunit system is rerun
with a rational parametrization so that candidates are components over
QQ rather than isolated modular fibres.

This is an exploratory chart probe.  Its output is not a theorem until
every nonunit result and every coordinate boundary is separately audited.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
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


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree44_two_row_off_diagonal_probe.json"
)
ROW_PAIR_REPRESENTATIVES = ((0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3))


def compositions(total: int, parts: int):
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in compositions(total - first, parts - 1):
            yield (first, *tail)


def support_for(rows: tuple[int, int]) -> tuple[tuple[int, int], ...]:
    support = tuple(
        sorted(
            (row, column)
            for row in rows
            for column in range(5)
            if column != row
        )
    )
    if len(support) != 8:
        raise AssertionError("unexpected off-diagonal two-row support")
    return support


def normalized_chart(
    rows: tuple[int, int],
) -> tuple[
    tuple[tuple[int, int], ...],
    tuple[tuple[int, int], tuple[int, int]],
    tuple[sp.Symbol, ...],
    tuple[sp.Expr, ...],
]:
    support = support_for(rows)
    anchors = support[:2]
    weights = tuple(row - column for row, column in anchors)
    if weights[0] == weights[1]:
        raise AssertionError("normalization anchors have dependent torus weights")
    residuals = sp.symbols("z0:6")
    coefficients = (sp.Integer(1), sp.Integer(1), *residuals)
    return support, anchors, residuals, coefficients


def primitive_moment(
    support: tuple[tuple[int, int], ...],
    residuals: tuple[sp.Symbol, ...],
    order: int,
) -> sp.Poly:
    terms: dict[tuple[int, ...], int] = {}
    order_factorial = factorial(order)
    for counts in compositions(order, len(support)):
        row_degree = sum(
            count * position[0]
            for count, position in zip(counts, support, strict=True)
        )
        column_degree = sum(
            count * position[1]
            for count, position in zip(counts, support, strict=True)
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
    return sp.Poly.from_dict(terms, residuals, domain=sp.ZZ)


def run_msolve(
    residuals: tuple[sp.Symbol, ...],
    moments: tuple[sp.Poly, ...],
    timeout: int,
    *,
    rational_parametrization: bool,
) -> dict[str, object]:
    msolve = shutil.which("msolve")
    if msolve is None:
        raise RuntimeError("msolve is required")
    h = sp.Symbol("h")
    equations = [
        sp.expand(h * sp.prod(residuals) - 1),
        *(moment.as_expr() for moment in moments if not moment.is_zero),
    ]
    with tempfile.TemporaryDirectory(prefix="sic44-two-row-offdiag-") as directory:
        input_path = Path(directory) / "system.ms"
        output_path = Path(directory) / "result.ms"
        input_path.write_text(
            ",".join(str(value) for value in (h, *residuals))
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
                *( ["-P", "1"] if rational_parametrization else [] ),
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
        "result_head": result[:4000],
        "stderr_tail": completed.stderr[-2000:],
    }


def classify(rows: tuple[int, int], through: int, timeout: int) -> dict[str, object]:
    support, anchors, residuals, coefficients = normalized_chart(rows)
    del coefficients
    moments = tuple(
        primitive_moment(support, residuals, order)
        for order in range(1, through + 1)
    )
    first = run_msolve(
        residuals,
        moments,
        timeout,
        rational_parametrization=False,
    )
    if first["returncode"] != 0 or first["status"] == "unparsed":
        raise RuntimeError(f"msolve failed on rows {rows}: {first}")
    component = None
    if first["status"] != "unit_ideal":
        component = run_msolve(
            residuals,
            moments,
            timeout,
            rational_parametrization=True,
        )
        if component["returncode"] != 0 or component["status"] == "unparsed":
            raise RuntimeError(f"RUR failed on rows {rows}: {component}")

    left, right = rows
    exact_minor = ((left, left), (left, right), (right, left), (right, right))
    return {
        "row_pair": list(rows),
        "reversal_partner": [4 - right, 4 - left],
        "support": [list(position) for position in support],
        "normalization_anchors": [list(position) for position in anchors],
        "normalization_weights": [row - column for row, column in anchors],
        "residual_coordinate_count": len(residuals),
        "internal_GL2_gauge": "removed by U=(e_r,e_s)",
        "exact_rank_two_minor_positions": [list(position) for position in exact_minor],
        "exact_rank_two_minor": "-c_rs*c_sr, nonzero on the coefficient torus",
        "moment_profiles": [
            {
                "order": order,
                "total_degree": None if moment.is_zero else int(moment.total_degree()),
                "terms": len(moment.terms()),
            }
            for order, moment in enumerate(moments, start=1)
        ],
        "classification_through_requested_order": first,
        "component_rational_parametrization": component,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--through", type=int, default=8)
    parser.add_argument("--timeout", type=int, default=600)
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()
    if arguments.through < 1 or arguments.timeout < 1 or arguments.workers < 1:
        raise ValueError("--through, --timeout, and --workers must be positive")

    records: list[dict[str, object]] = []
    with ProcessPoolExecutor(max_workers=arguments.workers) as executor:
        futures = {
            executor.submit(classify, rows, arguments.through, arguments.timeout): rows
            for rows in ROW_PAIR_REPRESENTATIVES
        }
        for future in as_completed(futures):
            record = future.result()
            records.append(record)
            print(
                f"rows={tuple(record['row_pair'])} "
                f"status={record['classification_through_requested_order']['status']}",
                flush=True,
            )
    records.sort(key=lambda record: record["row_pair"])

    artifact = {
        "format": "two-pair-sic-bidegree44-two-row-off-diagonal-probe-v1",
        "field": "characteristic zero",
        "degree": 4,
        "moment_orders": [1, arguments.through],
        "row_pair_orbit_count": len(records),
        "row_pair_orbits": records,
        "scope": (
            "exact dense coefficient-torus probes on six direct rank-two "
            "two-row quotient-chart orbits; boundaries and generic U charts remain open"
        ),
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
