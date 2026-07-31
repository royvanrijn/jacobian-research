#!/usr/bin/env python3
"""Exact/sharded screen of four-entry bidegree-(3,3) supports.

This is the two-parameter continuation of
``research_two_pair_sic_bidegree33_sparse_counterexample.py``.  On every
mixed-sign four-entry support, overall scaling and the diagonal torus
normalize two distinct-weight coefficients to one.  The remaining
coefficient-torus chart has coordinates ``z,w``.  Adding
``h*z*w-1`` saturates by their product, and an exact Groebner basis equal
to ``[1]`` excludes that dense support.  Boundary points have support at
most three and are handled by the companion exhaustive screen.

The script is deliberately sharded.  Some supports can be much harder
than others, and a bounded shard is an exact partial result rather than a
claim that an unfinished enumeration was exhaustive.
"""

from __future__ import annotations

import argparse
from itertools import combinations
import json
from math import factorial, gcd
from pathlib import Path
import time

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_sparse_four_support_screen.json"
)
POSITIONS = tuple((i, j) for i in range(4) for j in range(4))
Z, W, H = sp.symbols("z w h")


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--through", type=int, default=14)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument(
        "--per-support-seconds",
        type=float,
        default=30.0,
        help="soft elapsed-time ledger; SymPy calls cannot be interrupted",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def normalize_support(
    support: tuple[tuple[int, int], ...],
) -> tuple[
    tuple[int, int],
    tuple[int, int],
    tuple[int, int],
    tuple[int, int],
]:
    for first, second in combinations(range(4), 2):
        if (
            support[first][0] - support[first][1]
            != support[second][0] - support[second][1]
        ):
            residual = [
                support[index]
                for index in range(4)
                if index not in (first, second)
            ]
            return support[first], support[second], residual[0], residual[1]
    raise ValueError("support has only one weight")


def restricted_moment(
    order: int,
    normalized: tuple[
        tuple[int, int],
        tuple[int, int],
        tuple[int, int],
        tuple[int, int],
    ],
) -> sp.Poly:
    answer: dict[tuple[int, int], int] = {}
    order_factorial = factorial(order)
    first, second, third, fourth = normalized
    for first_count in range(order + 1):
        for second_count in range(order - first_count + 1):
            for third_count in range(
                order - first_count - second_count + 1
            ):
                fourth_count = (
                    order - first_count - second_count - third_count
                )
                counts = (
                    first_count,
                    second_count,
                    third_count,
                    fourth_count,
                )
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
                multinomial = order_factorial // (
                    factorial(first_count)
                    * factorial(second_count)
                    * factorial(third_count)
                    * factorial(fourth_count)
                )
                coefficient = (
                    multinomial
                    * factorial(3 * order - x_degree)
                    * factorial(x_degree)
                )
                exponent = (third_count, fourth_count)
                answer[exponent] = answer.get(exponent, 0) + coefficient
    content = 0
    for coefficient in answer.values():
        content = gcd(content, abs(coefficient))
    if content:
        answer = {
            exponent: coefficient // content
            for exponent, coefficient in answer.items()
        }
    return sp.Poly.from_dict(answer, (Z, W), domain=sp.ZZ)


def verify_restricted_formula() -> list[dict[str, object]]:
    """Compare the count recursion with direct bivariate expansion."""

    x_symbol, y_symbol = sp.symbols("x y")
    checks = []
    for support in (
        ((0, 1), (1, 2), (2, 3), (3, 0)),
        ((0, 3), (1, 0), (2, 1), (3, 2)),
    ):
        normalized = normalize_support(support)
        coefficients = (1, 1, Z, W)
        source = sum(
            coefficient
            * x_symbol ** position[0]
            * y_symbol ** position[1]
            for coefficient, position in zip(coefficients, normalized)
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
                Z,
                W,
                domain=sp.ZZ,
            ).primitive()[1]
            recursive = restricted_moment(order, normalized)
            if direct_primitive != recursive:
                raise AssertionError(
                    f"restricted moment mismatch on {support}, mu{order}"
                )
        checks.append(
            {
                "support": [list(position) for position in support],
                "orders": [1, 6],
                "passed": True,
            }
        )
    return checks


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
    equations = [H * Z * W - 1]
    profiles = []
    started = time.monotonic()
    decisive_order = None
    basis_size = None
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
            Z,
            W,
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
        "residual": [list(normalized[2]), list(normalized[3])],
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


def main() -> None:
    arguments = parse_arguments()
    formula_checks = verify_restricted_formula()
    mixed = [
        support
        for support in combinations(POSITIONS, 4)
        if min(i - j for i, j in support) < 0
        and max(i - j for i, j in support) > 0
    ]
    stop = (
        len(mixed)
        if not arguments.limit
        else min(len(mixed), arguments.start + arguments.limit)
    )
    selected = mixed[arguments.start:stop]
    records = [
        screen_support(
            support,
            arguments.through,
            arguments.per_support_seconds,
        )
        for support in selected
    ]
    counts: dict[str, int] = {}
    for record in records:
        status = str(record["status"])
        counts[status] = counts.get(status, 0) + 1
    payload = {
        "calculation": "two_pair_sic_bidegree33_sparse_four_support_screen",
        "scope": (
            "exact characteristic-zero coefficient-torus Groebner screen; "
            "sharded enumeration, with support-at-most-three boundaries "
            "handled by the companion exact screen"
        ),
        "through": arguments.through,
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
        "independent_formula_checks": formula_checks,
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
