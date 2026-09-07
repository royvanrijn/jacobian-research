#!/usr/bin/env python3
"""Exact component probe on a two-separated-row rank-two SIC chart.

The coefficient support is

    {0} x {1,2,3,4}  union  {4} x {0,1,2,3}.

Overall scaling and the contraction-preserving diagonal torus normalize
c_01=c_43=1 over the algebraic closure.  The remaining six coefficients
are saturated, so c_04*c_40 is nonzero and every point on the chart has
coefficient rank exactly two.  The script forms the bidegree-(4,4) pure
moments directly in these quotient coordinates and asks msolve for the
exact characteristic-zero coefficient-torus scheme.

This is a chart probe.  Its output is not promoted to a theorem unless a
separate checker records and verifies every resulting component or proves
the localized system is the unit ideal.
"""

from __future__ import annotations

import argparse
import hashlib
from math import factorial, gcd
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
from typing import Any

import sympy as sp


SUPPORT = ((0, 1), (4, 3), (0, 2), (0, 3), (0, 4), (4, 0), (4, 1), (4, 2))
COEFFICIENTS = (sp.Integer(1), sp.Integer(1), *sp.symbols("z0:6"))
RESIDUALS = tuple(COEFFICIENTS[2:])


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--through", type=int, default=14)
    parser.add_argument("--timeout", type=int, default=600)
    parser.add_argument("--threads", type=int, default=4)
    return parser.parse_args()


def compositions(total: int, parts: int):
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in compositions(total - first, parts - 1):
            yield (first, *tail)


def moment(order: int) -> sp.Poly:
    terms: dict[tuple[int, ...], int] = {}
    order_factorial = factorial(order)
    for counts in compositions(order, len(SUPPORT)):
        row_degree = sum(
            count * position[0]
            for count, position in zip(counts, SUPPORT, strict=True)
        )
        column_degree = sum(
            count * position[1]
            for count, position in zip(counts, SUPPORT, strict=True)
        )
        if row_degree != column_degree:
            continue
        denominator = 1
        for count in counts:
            denominator *= factorial(count)
        scalar = (
            (order_factorial // denominator)
            * factorial(row_degree)
            * factorial(4 * order - row_degree)
        )
        exponent = counts[2:]
        terms[exponent] = terms.get(exponent, 0) + scalar

    content = 0
    for scalar in terms.values():
        content = gcd(content, abs(scalar))
    if content:
        terms = {
            exponent: scalar // content
            for exponent, scalar in terms.items()
        }
    return sp.Poly.from_dict(terms, RESIDUALS, domain=sp.ZZ)


def solve(
    through: int,
    timeout: int,
    threads: int,
    *,
    rational_parametrization: bool = True,
) -> dict[str, Any]:
    msolve = shutil.which("msolve")
    if msolve is None:
        raise RuntimeError("msolve is required")
    moments = [moment(order) for order in range(1, through + 1)]
    h = sp.Symbol("h")
    equations = [
        str(sp.expand(h * sp.prod(RESIDUALS) - 1)).replace("**", "^")
    ] + [
        str(value.as_expr()).replace("**", "^")
        for value in moments
        if not value.is_zero
    ]

    with tempfile.TemporaryDirectory(prefix="sic44-two-row-") as directory:
        input_path = Path(directory) / "system.ms"
        output_path = Path(directory) / "result.ms"
        input_path.write_text(
            ",".join((str(h), *(str(value) for value in RESIDUALS)))
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
                str(threads),
                "-l",
                "2",
                "-v",
                "0",
            ]
            + (["-P", "1"] if rational_parametrization else []),
            text=True,
            capture_output=True,
            check=False,
            timeout=timeout if timeout else None,
        )
        result = (
            output_path.read_text(encoding="utf-8").strip()
            if output_path.exists()
            else ""
        )

    compact = " ".join(result.split())
    header = re.match(r"\[0,\s*\[0,\s*(\d+),\s*(\d+),", result)
    status = "unparsed"
    if result in ("[-1]", "[-1]:"):
        status = "unit_ideal"
    elif result.startswith("[1,"):
        status = "positive_dimensional"
    elif header is not None:
        status = "zero_dimensional"

    return {
        "returncode": completed.returncode,
        "through": through,
        "profiles": [
            (order, value.total_degree(), len(value.terms()))
            for order, value in enumerate(moments, start=1)
        ],
        "status": status,
        "variables": int(header.group(1)) if header is not None else None,
        "scheme_degree": int(header.group(2)) if header is not None else None,
        "result_sha256": hashlib.sha256(compact.encode()).hexdigest(),
        "result_bytes": len(result.encode()),
        "result": result,
        "stdout_tail": completed.stdout[-2000:],
        "stderr_tail": completed.stderr[-2000:],
    }


def main() -> None:
    arguments = parse_arguments()
    record = solve(arguments.through, arguments.timeout, arguments.threads)
    result = str(record["result"])
    print(f"returncode={record['returncode']}")
    print(f"orders=1..{arguments.through}")
    print(f"profiles={record['profiles']}")
    if record["status"] == "unit_ideal":
        print("status=unit_ideal")
        return
    if record["status"] == "positive_dimensional":
        print("status=positive_dimensional")
        return
    if record["status"] == "zero_dimensional":
        print("status=zero_dimensional")
        print(f"variables={record['variables']}")
        print(f"scheme_degree={record['scheme_degree']}")
        print(f"result_sha256={record['result_sha256']}")
        print(f"result_bytes={record['result_bytes']}")
        return
    print("status=unparsed")
    print(result[:2000])


if __name__ == "__main__":
    main()
