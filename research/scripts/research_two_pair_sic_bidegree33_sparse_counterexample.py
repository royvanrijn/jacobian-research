#!/usr/bin/env python3
"""Exhaustive exact screen of three-entry bidegree-(3,3) supports.

For

    F_C(x,y) = sum_{0 <= i,j <= 3} c_ij x^i y^j,

the raw contraction moment is

    mu_m = sum_I (3m-I)! I! [x^I y^I] F_C(x,y)^m.

On a three-entry support containing both positive and negative diagonal
weights ``i-j``, overall scaling and the residual diagonal torus normalize
two entries of distinct weights to one.  The third entry is a single
parameter ``z``.  This script computes the resulting univariate moments
exactly and their saturated gcd.  A constant gcd excludes a common
moment-zero point with all three displayed coefficients nonzero.

Supports contained in a strict weight half-space are Hilbert--Mumford
one-sided.  Supports involving weight zero but only one strict sign reduce
to the already elementary one- or two-variable diagonal moment problem.
The all-weight-zero case is covered by the exact diagonal-slice theorem
in ``TWO_PAIR_SIC_BIDEGREE33_FRONTIER.md`` and is reported separately.
"""

from __future__ import annotations

import argparse
from itertools import combinations
import json
from math import factorial
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_sparse_three_support_screen.json"
)
POSITIONS = tuple((i, j) for i in range(4) for j in range(4))
Z = sp.Symbol("z")


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--through", type=int, default=14)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def normalized_support(
    support: tuple[tuple[int, int], ...],
) -> tuple[tuple[int, int], tuple[int, int], tuple[int, int]]:
    """Choose two distinct-weight anchors and one residual coordinate."""

    for first, second in combinations(range(3), 2):
        if (
            support[first][0] - support[first][1]
            != support[second][0] - support[second][1]
        ):
            residual = 3 - first - second
            return support[first], support[second], support[residual]
    raise ValueError("support has only one weight")


def restricted_moment(
    order: int,
    anchors_and_residual: tuple[
        tuple[int, int],
        tuple[int, int],
        tuple[int, int],
    ],
) -> sp.Poly:
    """Return the exact one-parameter raw moment."""

    first, second, residual = anchors_and_residual
    answer = 0
    factorial_order = factorial(order)
    for first_count in range(order + 1):
        for second_count in range(order - first_count + 1):
            residual_count = order - first_count - second_count
            diagonal_degree_x = (
                first_count * first[0]
                + second_count * second[0]
                + residual_count * residual[0]
            )
            diagonal_degree_y = (
                first_count * first[1]
                + second_count * second[1]
                + residual_count * residual[1]
            )
            if diagonal_degree_x != diagonal_degree_y:
                continue
            multinomial = factorial_order // (
                factorial(first_count)
                * factorial(second_count)
                * factorial(residual_count)
            )
            answer += (
                multinomial
                * factorial(3 * order - diagonal_degree_x)
                * factorial(diagonal_degree_x)
                * Z**residual_count
            )
    return sp.Poly(answer, Z, domain=sp.ZZ)


def remove_zero_root(polynomial: sp.Poly) -> sp.Poly:
    """Saturate a univariate polynomial by the residual coordinate."""

    while polynomial.degree() > 0 and polynomial.eval(0) == 0:
        polynomial = sp.exquo(polynomial, sp.Poly(Z, Z, domain=sp.ZZ))
    return polynomial.primitive()[1]


def screen_mixed_support(
    support: tuple[tuple[int, int], ...],
    through: int,
) -> dict[str, object]:
    normalized = normalized_support(support)
    common: sp.Poly | None = None
    profiles = []
    decisive_order = None
    for order in range(1, through + 1):
        moment = restricted_moment(order, normalized)
        if moment.is_zero:
            profiles.append({"order": order, "degree": -1, "terms": 0})
            continue
        primitive = moment.primitive()[1]
        profiles.append(
            {
                "order": order,
                "degree": int(primitive.degree()),
                "terms": len(primitive.terms()),
            }
        )
        common = primitive if common is None else sp.gcd(common, primitive)
        common = remove_zero_root(common)
        if common.degree() == 0:
            decisive_order = order
            break
    if common is None:
        common = sp.Poly(0, Z, domain=sp.ZZ)
    return {
        "support": [list(position) for position in support],
        "weights": [i - j for i, j in support],
        "anchors": [list(normalized[0]), list(normalized[1])],
        "residual": list(normalized[2]),
        "profiles": profiles,
        "decisive_order": decisive_order,
        "saturated_gcd_degree": (
            -1 if common.is_zero else int(common.degree())
        ),
        "saturated_gcd": str(common.as_expr()),
        "status": (
            "excluded_on_coefficient_torus"
            if not common.is_zero and common.degree() == 0
            else "survives_tested_orders"
        ),
    }


def diagonal_reduction_status(weights: tuple[int, int, int]) -> str:
    zeros = weights.count(0)
    if zeros == 3:
        return "covered_by_exact_diagonal_slice"
    if zeros <= 2 and (all(weight >= 0 for weight in weights) or all(
        weight <= 0 for weight in weights
    )):
        return "reduces_to_at_most_two_diagonal_entries"
    raise ValueError(weights)


def main() -> None:
    arguments = parse_arguments()
    if arguments.through < 1:
        raise ValueError("--through must be positive")

    records = []
    counts: dict[str, int] = {}
    for support in combinations(POSITIONS, 3):
        weights = tuple(i - j for i, j in support)
        if all(weight > 0 for weight in weights) or all(
            weight < 0 for weight in weights
        ):
            status = "strict_half_space_nullcone"
            record = {
                "support": [list(position) for position in support],
                "weights": list(weights),
                "status": status,
            }
        elif min(weights) < 0 < max(weights):
            record = screen_mixed_support(support, arguments.through)
            status = str(record["status"])
        elif len(set(weights)) == 1 and weights[0] != 0:
            status = "strict_half_space_nullcone"
            record = {
                "support": [list(position) for position in support],
                "weights": list(weights),
                "status": status,
            }
        else:
            status = diagonal_reduction_status(weights)
            record = {
                "support": [list(position) for position in support],
                "weights": list(weights),
                "status": status,
            }
        counts[status] = counts.get(status, 0) + 1
        records.append(record)

    survivors = [
        record
        for record in records
        if record["status"] == "survives_tested_orders"
    ]
    payload = {
        "calculation": "two_pair_sic_bidegree33_sparse_three_support_screen",
        "through": arguments.through,
        "support_size": 3,
        "support_count": len(records),
        "counts": counts,
        "survivors": survivors,
        "all_mixed_supports_excluded": not survivors,
        "scope": (
            "exact characteristic-zero coefficient-torus screen for every "
            "three-entry support; strict half-space and diagonal reductions "
            "are classified separately"
        ),
        "records": records,
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({key: payload[key] for key in (
        "through",
        "support_count",
        "counts",
        "all_mixed_supports_excluded",
        "survivors",
    )}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
