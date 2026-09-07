#!/usr/bin/env python3
"""Rank-seven local continuation at a new five-companion Mestre seed.

The normalized roots ``(0,1,7,8,9,11)`` satisfy the Mestre condition and
have D=48^2.  The recursive ordinate-eliminated equations from
``probe_mestre_two_section_local_continuation`` vanish for

    x_1=(61+7*T)/5,  x_2=33/5.

Unlike the two previously known six-companion fibers, this labelled pair has
Jacobian rank seven in the eight coordinates ``(c1,...,c4,x01,x11,x02,x12)``.
The script gives a bounded regular Hensel continuation along a free coordinate
at several good primes.  The separate
``verify_mestre_transverse_two_section_component.py`` script proves the
recognized rational component.  This local probe itself does not address
section independence, a height pairing, or rank at least fourteen.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path
from typing import Iterable

from probe_mestre_two_section_local_continuation import (
    Field,
    FormalSeries,
    VARIABLES,
    rational_mod,
    rational_reconstruction,
    residuals,
    residuals_from_jets,
    row_reduce,
    solve_square,
    solve_square_over_q,
)


Q = Fraction
MODULI = (Q(-35), Q(455), Q(-2605), Q(5544))
SECTIONS = ((Q(61, 5), Q(7, 5)), (Q(33, 5), Q(0)))
DEFAULT_PRIMES = (17, 19, 23)


def transverse_coordinate(jacobian: list[list[Fraction]]) -> int:
    """Prefer a moduli coordinate as the one-dimensional local parameter."""

    for index in (3, 0, 1, 2, 4, 5, 6, 7):
        minor = [[value for column, value in enumerate(row) if column != index] for row in jacobian]
        if row_reduce(minor, Field())[0] == 7:
            return index
    raise AssertionError("the rank-seven point has no transverse coordinate")


def hensel_branch(prime: int, precision: int, free: int) -> dict[str, object]:
    """Lift after moving the free coordinate by p in Z_p."""

    seed = MODULI + SECTIONS[0] + SECTIONS[1]
    coordinates = [rational_mod(value, prime) for value in seed]
    jacobian = [list(map(int, value.gradient)) for value in residuals(coordinates, Field(prime))]
    for candidate in (free, *range(len(VARIABLES))):
        columns = [index for index in range(len(VARIABLES)) if index != candidate]
        minor = [[row[index] for index in columns] for row in jacobian]
        if row_reduce(minor, Field(prime))[0] == 7:
            free = candidate
            break
    else:
        raise AssertionError(f"no transverse coordinate exists modulo {prime}")
    target = seed[free] + prime
    for exponent in range(1, precision):
        modulus = prime ** (exponent + 1)
        trial = list(coordinates)
        trial[free] = rational_mod(target, modulus)
        values = residuals(trial, Field(modulus, tangent=False))
        right = [(-int(value.value) // (prime**exponent)) % prime for value in values]
        correction = solve_square(minor, right, prime)
        for column, digit in zip(columns, correction):
            trial[column] = (trial[column] + prime**exponent * digit) % modulus
        coordinates = trial
    modulus = prime**precision
    if any(value.value for value in residuals(coordinates, Field(modulus, tangent=False))):
        raise AssertionError("regular Hensel lift failed")
    reconstructions = [rational_reconstruction(value, modulus) for value in coordinates]
    exact_reconstruction = False
    if all(value is not None for value in reconstructions):
        exact_reconstruction = all(value.value == 0 for value in residuals(reconstructions, Field()))
    return {
        "prime": prime,
        "precision": precision,
        "modulus": modulus,
        "free_coordinate": VARIABLES[free],
        "free_coordinate_target": str(target),
        "coordinates_mod_prime_power": dict(zip(VARIABLES, coordinates)),
        "all_seven_residuals_zero_mod_prime_power": True,
        "low_height_rational_reconstruction": {
            "all_coordinates_reconstructed": all(value is not None for value in reconstructions),
            "exactly_satisfies_the_seven_Q_equations": exact_reconstruction,
        },
    }


def formal_transverse_series(order: int = 12) -> dict[str, object]:
    """Implicitly solve all seven rows in Q[[t]] with c4=c4(0)+t."""

    seed = MODULI + SECTIONS[0] + SECTIONS[1]
    exact = residuals(seed, Field())
    jacobian = [list(value.gradient) for value in exact]
    free = transverse_coordinate(jacobian)
    columns = [index for index in range(len(VARIABLES)) if index != free]
    minor = [[row[index] for index in columns] for row in jacobian]
    if row_reduce(minor, Field())[0] != 7:
        raise AssertionError("the characteristic-zero transverse minor is singular")
    coefficients = [[Q(0)] * (order + 1) for _ in VARIABLES]
    for index, value in enumerate(seed):
        coefficients[index][0] = value
    coefficients[free][1] = Q(1)
    field = Field(series_order=order)
    for degree in range(1, order + 1):
        values = residuals_from_jets(
            [FormalSeries.seed(field, row) for row in coefficients]
        )
        correction = solve_square_over_q(
            minor, [-value.coefficients[degree] for value in values]
        )
        for column, value in zip(columns, correction):
            coefficients[column][degree] = value
        checked = residuals_from_jets(
            [FormalSeries.seed(field, row) for row in coefficients]
        )
        if any(value.coefficients[degree] != 0 for value in checked):
            raise AssertionError("formal rank-seven implicit solve failed")
    return {
        "parameter": "c4=5544+t",
        "order": order,
        "all_seven_coefficients_through_order_vanish": True,
        "coordinate_series": {
            name: [str(value) for value in row]
            for name, row in zip(VARIABLES, coefficients)
        },
    }


def run(primes: Iterable[int], precision: int) -> dict[str, object]:
    seed = MODULI + SECTIONS[0] + SECTIONS[1]
    exact = residuals(seed, Field())
    if any(value.value for value in exact):
        raise AssertionError("the new transverse seed is not an exact residual zero")
    jacobian = [list(value.gradient) for value in exact]
    rank, pivots = row_reduce(jacobian, Field())
    if rank != 7:
        raise AssertionError(f"expected exact rank 7, found {rank}")
    free = transverse_coordinate(jacobian)
    checks = []
    branches = []
    for prime in primes:
        reduced = residuals([rational_mod(value, prime) for value in seed], Field(prime))
        reduced_rank, reduced_pivots = row_reduce(
            [list(value.gradient) for value in reduced], Field(prime)
        )
        if reduced_rank != 7:
            raise AssertionError(f"rank dropped modulo {prime}")
        checks.append({
            "prime": prime,
            "jacobian_rank": reduced_rank,
            "pivot_columns": [VARIABLES[index] for index in reduced_pivots],
        })
        branches.append(hensel_branch(prime, precision, free))
    return {
        "status": "rank-seven bounded local two-section continuation recorded",
        "normalized_roots": [0, 1, 7, 8, 9, 11],
        "leading_invariant_D": "2304=48^2",
        "sections": [[str(value) for value in section] for section in SECTIONS],
        "equations": "M plus (E2,E1,E0) for each section, evaluated recursively",
        "expanded_residuals_materialized": False,
        "exact_jacobian_rank": rank,
        "exact_pivot_columns": [VARIABLES[index] for index in pivots],
        "transverse_free_coordinate": VARIABLES[free],
        "finite_field_tangent_checks": checks,
        "hensel_continuations": branches,
        "exact_component_verifier": "verify_mestre_transverse_two_section_component.py",
        "not_established": [
            "pair intersections at infinity or a Shioda Gram matrix",
            "saturation or independence from the generic rank-13 subgroup",
            "generic rank at least 14",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primes", nargs="+", type=int, default=DEFAULT_PRIMES)
    parser.add_argument("--precision", type=int, default=5)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rendered = json.dumps(run(args.primes, args.precision), indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        args.output.write_text(rendered)


if __name__ == "__main__":
    main()
