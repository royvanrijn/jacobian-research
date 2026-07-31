#!/usr/bin/env python3
"""Sharded exact coefficient-torus screen for five-entry supports.

After normalizing two distinct-weight coefficients, a mixed-sign
five-entry support has three residual coordinates.  The equation
``h*z0*z1*z2-1`` restricts to the dense coefficient torus.  Incremental
characteristic-zero Groebner bases test the contraction moments there.
The support boundary is covered by the completed support-at-most-four
screen.
"""

from __future__ import annotations

import argparse
from itertools import combinations
import json
from math import factorial, gcd
from pathlib import Path
import shutil
import subprocess
import tempfile
import time

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_sparse_five_support_screen.json"
)
POSITIONS = tuple((i, j) for i in range(4) for j in range(4))
Z0, Z1, Z2, H = sp.symbols("z0 z1 z2 h")
RESIDUAL_SYMBOLS = (Z0, Z1, Z2)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--through", type=int, default=14)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--per-support-seconds", type=float, default=30.0)
    parser.add_argument(
        "--backend",
        choices=("sympy", "msolve"),
        default="sympy",
    )
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def normalize_support(
    support: tuple[tuple[int, int], ...],
) -> tuple[tuple[int, int], ...]:
    for first, second in combinations(range(5), 2):
        if (
            support[first][0] - support[first][1]
            != support[second][0] - support[second][1]
        ):
            return (
                support[first],
                support[second],
                *(
                    support[index]
                    for index in range(5)
                    if index not in (first, second)
                ),
            )
    raise ValueError("support has only one weight")


def compositions(total: int, parts: int):
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in compositions(total - first, parts - 1):
            yield (first, *tail)


def restricted_moment(
    order: int,
    normalized: tuple[tuple[int, int], ...],
) -> sp.Poly:
    answer: dict[tuple[int, int, int], int] = {}
    order_factorial = factorial(order)
    for counts in compositions(order, 5):
        x_degree = sum(
            count * position[0]
            for count, position in zip(counts, normalized)
        )
        y_degree = sum(
            count * position[1]
            for count, position in zip(counts, normalized)
        )
        if x_degree != y_degree:
            continue
        denominator = 1
        for count in counts:
            denominator *= factorial(count)
        coefficient = (
            order_factorial
            // denominator
            * factorial(3 * order - x_degree)
            * factorial(x_degree)
        )
        exponent = counts[2:]
        answer[exponent] = answer.get(exponent, 0) + coefficient
    content = 0
    for coefficient in answer.values():
        content = gcd(content, abs(coefficient))
    if content:
        answer = {
            exponent: coefficient // content
            for exponent, coefficient in answer.items()
        }
    return sp.Poly.from_dict(answer, RESIDUAL_SYMBOLS, domain=sp.ZZ)


def verify_restricted_formula() -> dict[str, object]:
    """Compare the count recursion with direct bivariate expansion."""

    support = ((0, 1), (1, 2), (2, 0), (2, 3), (3, 1))
    normalized = normalize_support(support)
    x_symbol, y_symbol = sp.symbols("x y")
    source = sum(
        coefficient
        * x_symbol ** position[0]
        * y_symbol ** position[1]
        for coefficient, position in zip(
            (1, 1, Z0, Z1, Z2),
            normalized,
        )
    )
    for order in range(1, 7):
        expanded = sp.expand(source**order)
        direct = 0
        for diagonal_degree in range(3 * order + 1):
            direct += (
                factorial(3 * order - diagonal_degree)
                * factorial(diagonal_degree)
                * expanded.coeff(x_symbol, diagonal_degree).coeff(
                    y_symbol,
                    diagonal_degree,
                )
            )
        direct_primitive = sp.Poly(
            direct,
            *RESIDUAL_SYMBOLS,
            domain=sp.ZZ,
        ).primitive()[1]
        if direct_primitive != restricted_moment(order, normalized):
            raise AssertionError(
                f"restricted moment mismatch on {support}, mu{order}"
            )
    return {
        "support": [list(position) for position in support],
        "orders": [1, 6],
        "passed": True,
    }


def is_unit_basis(basis: sp.GroebnerBasis) -> bool:
    return (
        len(basis.polys) == 1
        and basis.polys[0].total_degree() == 0
        and basis.polys[0].LC() != 0
    )


def screen_support(
    support: tuple[tuple[int, int], ...],
    through: int,
    soft_seconds: float,
) -> dict[str, object]:
    normalized = normalize_support(support)
    equations = [H * Z0 * Z1 * Z2 - 1]
    started = time.monotonic()
    decisive_order = None
    basis_size = None
    profiles = []
    for order in range(1, through + 1):
        moment = restricted_moment(order, normalized)
        profiles.append(
            {
                "order": order,
                "degree": (
                    -1 if moment.is_zero else int(moment.total_degree())
                ),
                "terms": len(moment.terms()),
            }
        )
        if moment.is_zero:
            continue
        equations.append(moment.as_expr())
        basis = sp.groebner(
            equations,
            H,
            Z0,
            Z1,
            Z2,
            order="grevlex",
            domain=sp.QQ,
        )
        basis_size = len(basis.polys)
        if is_unit_basis(basis):
            decisive_order = order
            break
        if time.monotonic() - started > soft_seconds:
            break
    seconds = time.monotonic() - started
    return {
        "support": [list(position) for position in support],
        "weights": [i - j for i, j in support],
        "anchors": [list(normalized[0]), list(normalized[1])],
        "residual": [
            list(position)
            for position in normalized[2:]
        ],
        "profiles": profiles,
        "decisive_order": decisive_order,
        "last_basis_size": basis_size,
        "seconds": round(seconds, 6),
        "status": (
            "excluded_on_coefficient_torus"
            if decisive_order is not None
            else (
                "soft_time_limit"
                if seconds > soft_seconds
                else "survives_tested_orders"
            )
        ),
    }


def msolve_expression(polynomial: sp.Poly) -> str:
    return str(polynomial.as_expr()).replace("**", "^")


def screen_support_msolve(
    support: tuple[tuple[int, int], ...],
    through: int,
    timeout: int,
) -> dict[str, object]:
    """Submit one saturated dense-support ideal exactly over QQ."""

    normalized = normalize_support(support)
    moments = []
    profiles = []
    started = time.monotonic()
    for order in range(1, through + 1):
        moment = restricted_moment(order, normalized)
        profiles.append(
            {
                "order": order,
                "degree": (
                    -1 if moment.is_zero else int(moment.total_degree())
                ),
                "terms": len(moment.terms()),
            }
        )
        if not moment.is_zero:
            moments.append(moment)
    msolve = shutil.which("msolve")
    if msolve is None:
        return {
            "support": [list(position) for position in support],
            "status": "missing_msolve",
        }
    equations = ["h*z0*z1*z2-1", *(
        msolve_expression(moment)
        for moment in moments
    )]
    with tempfile.TemporaryDirectory(
        prefix="sic33-sparse-five-msolve-",
    ) as directory:
        input_path = Path(directory) / "system.ms"
        output_path = Path(directory) / "result.ms"
        input_path.write_text(
            "h,z0,z1,z2\n0\n"
            + ",\n".join(equations)
            + "\n",
            encoding="utf-8",
        )
        try:
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
                    "1",
                ],
                text=True,
                capture_output=True,
                timeout=timeout if timeout else None,
                check=False,
            )
        except subprocess.TimeoutExpired as error:
            return {
                "support": [list(position) for position in support],
                "weights": [i - j for i, j in support],
                "anchors": [list(normalized[0]), list(normalized[1])],
                "residual": [
                    list(position)
                    for position in normalized[2:]
                ],
                "profiles": profiles,
                "seconds": round(time.monotonic() - started, 6),
                "status": "timeout",
                "stdout_tail": (error.stdout or "")[-4000:],
                "stderr_tail": (error.stderr or "")[-4000:],
            }
        result = (
            output_path.read_text(encoding="utf-8").strip()
            if output_path.exists()
            else ""
        )
    unit = result in ("[-1]:", "[-1]")
    record: dict[str, object] = {
        "support": [list(position) for position in support],
        "weights": [i - j for i, j in support],
        "anchors": [list(normalized[0]), list(normalized[1])],
        "residual": [
            list(position)
            for position in normalized[2:]
        ],
        "decisive_order": through if unit else None,
        "seconds": round(time.monotonic() - started, 6),
        "status": (
            "excluded_on_coefficient_torus"
            if unit
            else (
                "msolve_nonempty"
                if completed.returncode == 0 and result
                else "msolve_failed"
            )
        ),
        "msolve": (
            {
                "returncode": completed.returncode,
                "result": result,
            }
            if unit
            else {
            "returncode": completed.returncode,
            "result_head": result[:2000],
            "result_tail": result[-2000:],
            "stdout_tail": completed.stdout[-4000:],
            "stderr_tail": completed.stderr[-4000:],
            }
        ),
    }
    if not unit:
        record["profiles"] = profiles
    return record


def main() -> None:
    arguments = parse_arguments()
    formula_check = verify_restricted_formula()
    mixed = [
        support
        for support in combinations(POSITIONS, 5)
        if min(i - j for i, j in support) < 0
        and max(i - j for i, j in support) > 0
    ]
    stop = (
        len(mixed)
        if not arguments.limit
        else min(len(mixed), arguments.start + arguments.limit)
    )
    records = [
        (
            screen_support(
                support,
                arguments.through,
                arguments.per_support_seconds,
            )
            if arguments.backend == "sympy"
            else screen_support_msolve(
                support,
                arguments.through,
                arguments.timeout,
            )
        )
        for support in mixed[arguments.start:stop]
    ]
    counts: dict[str, int] = {}
    for record in records:
        status = str(record["status"])
        counts[status] = counts.get(status, 0) + 1
    payload = {
        "calculation": "two_pair_sic_bidegree33_sparse_five_support_screen",
        "scope": (
            "exact characteristic-zero coefficient-torus Groebner screen; "
            "sharded enumeration with boundaries covered by the exact "
            "support-at-most-four calculation"
        ),
        "through": arguments.through,
        "backend": arguments.backend,
        "mixed_support_count": len(mixed),
        "start": arguments.start,
        "stop": stop,
        "counts": counts,
        "complete_mixed_enumeration": (
            arguments.start == 0
            and stop == len(mixed)
            and all(
                record["status"] == "excluded_on_coefficient_torus"
                for record in records
            )
        ),
        "independent_formula_check": formula_check,
        "records": records,
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({key: payload[key] for key in (
        "through",
        "mixed_support_count",
        "start",
        "stop",
        "counts",
        "complete_mixed_enumeration",
    )}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
