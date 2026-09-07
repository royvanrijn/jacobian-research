#!/usr/bin/env python3
"""Collision-first finite-field search for sparse HC(4) potentials.

This is a bounded experiment, not a proof of HC(4) and not a
characteristic-zero lifting theorem.

Normalize a prospective gradient collision to the points 0 and e_0 and use

    q = x_0*x_1 + x_2*x_3,

whose Hessian determinant is one.  For a degree bound d, start with

    psi_base = q - x_1*x_0**(d-1).

The displayed higher term cancels the gradient difference of q at 0 and
e_0.  The script constructs a basis of *every* polynomial of degrees 3
through d whose gradient has the same value at those two points.  It then
exhausts affine perturbations of psi_base supported on at most two basis
directions, with arbitrary nonzero coefficients over GF(p).

Candidates are rejected by determinant evaluations at deterministic points.
This rejection is exact: one failed evaluation proves that the determinant
is not the constant one.  Any evaluation survivor is checked by a full sparse
polynomial expansion of det(Hess(psi))-1 over GF(p).

Use primes greater than the degree bound.  This avoids vanishing derivative
coefficients and makes the collision-kernel basis have the same support in
characteristic zero and in every searched characteristic.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from itertools import combinations, permutations
from pathlib import Path
from time import perf_counter
from typing import Iterable


Exponent = tuple[int, int, int, int]
RationalTerm = tuple[Exponent, Fraction]
Direction = tuple[RationalTerm, ...]
Matrix = tuple[int, ...]
Polynomial = dict[Exponent, int]

VARIABLE_COUNT = 4
ZERO_EXPONENT: Exponent = (0, 0, 0, 0)


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            return False
        divisor += 1
    return True


def weak_compositions(total: int, slots: int) -> Iterable[tuple[int, ...]]:
    if slots == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in weak_compositions(total - first, slots - 1):
            yield (first,) + tail


def carrier(degree: int, coordinate: int) -> Exponent:
    exponents = [0] * VARIABLE_COUNT
    if coordinate == 0:
        exponents[0] = degree
    else:
        exponents[0] = degree - 1
        exponents[coordinate] = 1
    return tuple(exponents)  # type: ignore[return-value]


def collision_kernel_directions(degree_bound: int) -> list[Direction]:
    """Basis for the degree-[3,d] perturbations with zero gradient jump."""

    top_carriers = {
        coordinate: carrier(degree_bound, coordinate)
        for coordinate in range(VARIABLE_COUNT)
    }
    directions: list[Direction] = []
    for degree in range(3, degree_bound + 1):
        for exponent in weak_compositions(degree, VARIABLE_COUNT):
            visible_coordinate = next(
                (
                    coordinate
                    for coordinate in range(VARIABLE_COUNT)
                    if exponent == carrier(degree, coordinate)
                ),
                None,
            )
            if visible_coordinate is None:
                directions.append(((exponent, Fraction(1)),))
                continue
            if degree == degree_bound:
                continue
            if visible_coordinate == 0:
                top_scale = Fraction(-degree, degree_bound)
            else:
                top_scale = Fraction(-1)
            directions.append(
                (
                    (exponent, Fraction(1)),
                    (top_carriers[visible_coordinate], top_scale),
                )
            )
    expected = (
        sum(
            len(list(weak_compositions(degree, VARIABLE_COUNT)))
            for degree in range(3, degree_bound + 1)
        )
        - VARIABLE_COUNT
    )
    assert len(directions) == expected
    return directions


def fraction_mod(value: Fraction, prime: int) -> int:
    return value.numerator * pow(value.denominator, -1, prime) % prime


def hessian_matrix_at(
    terms: Iterable[RationalTerm],
    point: tuple[int, int, int, int],
    prime: int,
) -> Matrix:
    result = [0] * (VARIABLE_COUNT * VARIABLE_COUNT)
    for exponent, rational_coefficient in terms:
        coefficient = fraction_mod(rational_coefficient, prime)
        for row in range(VARIABLE_COUNT):
            for column in range(VARIABLE_COUNT):
                multiplier = exponent[row] * (
                    exponent[column] - (1 if row == column else 0)
                )
                if multiplier == 0:
                    continue
                derived_exponent = list(exponent)
                derived_exponent[row] -= 1
                derived_exponent[column] -= 1
                value = coefficient * multiplier
                for coordinate, power in zip(
                    point, derived_exponent, strict=True
                ):
                    value = value * pow(coordinate, power, prime) % prime
                index = VARIABLE_COUNT * row + column
                result[index] = (result[index] + value) % prime
    return tuple(result)


def determinant_four(matrix: Matrix, prime: int) -> int:
    """Unrolled 4-by-4 determinant."""

    (
        a,
        b,
        c,
        d,
        e,
        f,
        g,
        h,
        i,
        j,
        k,
        ell,
        m,
        n,
        o,
        p,
    ) = matrix
    minor_00 = f * (k * p - ell * o) - g * (j * p - ell * n) + h * (
        j * o - k * n
    )
    minor_01 = e * (k * p - ell * o) - g * (i * p - ell * m) + h * (
        i * o - k * m
    )
    minor_02 = e * (j * p - ell * n) - f * (i * p - ell * m) + h * (
        i * n - j * m
    )
    minor_03 = e * (j * o - k * n) - f * (i * o - k * m) + g * (
        i * n - j * m
    )
    return (
        a * minor_00 - b * minor_01 + c * minor_02 - d * minor_03
    ) % prime


def add_scaled_matrices(
    base: Matrix,
    left_scale: int,
    left: Matrix,
    right_scale: int,
    right: Matrix,
    prime: int,
) -> Matrix:
    return tuple(
        (
            base[index]
            + left_scale * left[index]
            + right_scale * right[index]
        )
        % prime
        for index in range(VARIABLE_COUNT * VARIABLE_COUNT)
    )


def add_polynomial_term(
    polynomial: Polynomial,
    exponent: Exponent,
    coefficient: int,
    prime: int,
) -> None:
    value = (polynomial.get(exponent, 0) + coefficient) % prime
    if value:
        polynomial[exponent] = value
    else:
        polynomial.pop(exponent, None)


def multiply_polynomials(
    left: Polynomial, right: Polynomial, prime: int
) -> Polynomial:
    result: Polynomial = {}
    for left_exponent, left_coefficient in left.items():
        for right_exponent, right_coefficient in right.items():
            exponent = tuple(
                left_exponent[index] + right_exponent[index]
                for index in range(VARIABLE_COUNT)
            )
            add_polynomial_term(
                result,
                exponent,  # type: ignore[arg-type]
                left_coefficient * right_coefficient,
                prime,
            )
    return result


def permutation_sign(permutation: tuple[int, ...]) -> int:
    inversions = sum(
        permutation[left] > permutation[right]
        for left in range(VARIABLE_COUNT)
        for right in range(left + 1, VARIABLE_COUNT)
    )
    return -1 if inversions % 2 else 1


PERMUTATIONS = tuple(
    (permutation, permutation_sign(permutation))
    for permutation in permutations(range(VARIABLE_COUNT))
)


def exact_hessian_determinant(
    terms: Iterable[tuple[Exponent, int]], prime: int
) -> Polynomial:
    hessian: list[list[Polynomial]] = [
        [{} for _ in range(VARIABLE_COUNT)]
        for _ in range(VARIABLE_COUNT)
    ]
    for exponent, coefficient in terms:
        for row in range(VARIABLE_COUNT):
            for column in range(VARIABLE_COUNT):
                multiplier = exponent[row] * (
                    exponent[column] - (1 if row == column else 0)
                )
                if multiplier == 0:
                    continue
                derived_exponent = list(exponent)
                derived_exponent[row] -= 1
                derived_exponent[column] -= 1
                add_polynomial_term(
                    hessian[row][column],
                    tuple(derived_exponent),  # type: ignore[arg-type]
                    coefficient * multiplier,
                    prime,
                )
    determinant: Polynomial = {}
    for permutation, sign in PERMUTATIONS:
        product: Polynomial = {ZERO_EXPONENT: sign % prime}
        for row, column in enumerate(permutation):
            product = multiply_polynomials(
                product, hessian[row][column], prime
            )
            if not product:
                break
        for exponent, coefficient in product.items():
            add_polynomial_term(determinant, exponent, coefficient, prime)
    return determinant


def deterministic_points(
    prime: int, count: int
) -> list[tuple[int, int, int, int]]:
    points: list[tuple[int, int, int, int]] = [
        (1 % prime, 2 % prime, 3 % prime, 4 % prime),
        (2 % prime, 1 % prime, 4 % prime, 3 % prime),
        (1 % prime, 1 % prime, 2 % prime, 3 % prime),
    ]
    state = 1
    while len(points) < count:
        coordinates = []
        for _ in range(VARIABLE_COUNT):
            state = (1103515245 * state + 12345) & 0x7FFFFFFF
            coordinates.append(1 + state % (prime - 1))
        point = tuple(coordinates)
        if point not in points:
            points.append(point)  # type: ignore[arg-type]
    return points[:count]


def base_terms(degree_bound: int) -> list[RationalTerm]:
    return [
        ((1, 1, 0, 0), Fraction(1)),
        ((0, 0, 1, 1), Fraction(1)),
        (carrier(degree_bound, 1), Fraction(-1)),
    ]


def candidate_terms_mod(
    degree_bound: int,
    directions: list[Direction],
    indices: tuple[int, ...],
    coefficients: tuple[int, ...],
    prime: int,
) -> list[tuple[Exponent, int]]:
    combined: Polynomial = {}
    for exponent, coefficient in base_terms(degree_bound):
        add_polynomial_term(
            combined, exponent, fraction_mod(coefficient, prime), prime
        )
    for index, scalar in zip(indices, coefficients, strict=True):
        for exponent, coefficient in directions[index]:
            add_polynomial_term(
                combined,
                exponent,
                scalar * fraction_mod(coefficient, prime),
                prime,
            )
    return sorted(combined.items())


def gradient_jump(
    terms: Iterable[tuple[Exponent, int]], prime: int
) -> tuple[int, int, int, int]:
    jump = [0] * VARIABLE_COUNT
    for exponent, coefficient in terms:
        for coordinate in range(VARIABLE_COUNT):
            derived = exponent[coordinate]
            if derived == 0:
                continue
            remaining = list(exponent)
            remaining[coordinate] -= 1
            value_at_e0 = derived * coefficient
            if any(remaining[index] for index in range(1, VARIABLE_COUNT)):
                value_at_e0 = 0
            value_at_zero = value_at_e0 if sum(remaining) == 0 else 0
            jump[coordinate] += value_at_e0 - value_at_zero
    return tuple(value % prime for value in jump)  # type: ignore[return-value]


def direction_record(direction: Direction) -> list[dict[str, object]]:
    return [
        {
            "exponents": list(exponent),
            "coefficient": (
                str(coefficient.numerator)
                if coefficient.denominator == 1
                else f"{coefficient.numerator}/{coefficient.denominator}"
            ),
        }
        for exponent, coefficient in direction
    ]


def search_prime_degree(
    prime: int,
    degree_bound: int,
    support_bound: int,
    point_count: int,
) -> dict[str, object]:
    assert prime > degree_bound
    assert support_bound in (1, 2)
    directions = collision_kernel_directions(degree_bound)
    points = deterministic_points(prime, point_count)
    base_hessians = [
        hessian_matrix_at(base_terms(degree_bound), point, prime)
        for point in points
    ]
    direction_hessians = [
        [
            hessian_matrix_at(direction, point, prime)
            for direction in directions
        ]
        for point in points
    ]
    coefficient_values = range(1, prime)
    exact_candidates: list[dict[str, object]] = []
    evaluated = 0
    evaluation_survivors = 0
    started = perf_counter()

    def test(
        indices: tuple[int, ...], coefficients: tuple[int, ...]
    ) -> None:
        nonlocal evaluated, evaluation_survivors
        evaluated += 1
        for point_index in range(point_count):
            left_index = indices[0]
            left_coefficient = coefficients[0]
            if len(indices) == 1:
                right_index = left_index
                right_coefficient = 0
            else:
                right_index = indices[1]
                right_coefficient = coefficients[1]
            matrix = add_scaled_matrices(
                base_hessians[point_index],
                left_coefficient,
                direction_hessians[point_index][left_index],
                right_coefficient,
                direction_hessians[point_index][right_index],
                prime,
            )
            if determinant_four(matrix, prime) != 1:
                return
        evaluation_survivors += 1
        terms = candidate_terms_mod(
            degree_bound, directions, indices, coefficients, prime
        )
        determinant = exact_hessian_determinant(terms, prime)
        if determinant != {ZERO_EXPONENT: 1}:
            return
        assert gradient_jump(terms, prime) == (0, 0, 0, 0)
        exact_candidates.append(
            {
                "direction_indices": list(indices),
                "direction_coefficients_mod_p": list(coefficients),
                "directions": [
                    direction_record(directions[index]) for index in indices
                ],
                "potential_terms_mod_p": [
                    {
                        "exponents": list(exponent),
                        "coefficient": coefficient,
                    }
                    for exponent, coefficient in terms
                ],
            }
        )

    for index in range(len(directions)):
        for coefficient in coefficient_values:
            test((index,), (coefficient,))
    if support_bound >= 2:
        for left_index, right_index in combinations(
            range(len(directions)), 2
        ):
            for left_coefficient in coefficient_values:
                for right_coefficient in coefficient_values:
                    test(
                        (left_index, right_index),
                        (left_coefficient, right_coefficient),
                    )

    return {
        "prime": prime,
        "degree_bound": degree_bound,
        "collision_points": [[0, 0, 0, 0], [1, 0, 0, 0]],
        "direction_count": len(directions),
        "support_bound": support_bound,
        "deterministic_points": [list(point) for point in points],
        "potentials_evaluated": evaluated,
        "evaluation_survivors": evaluation_survivors,
        "exact_candidate_count": len(exact_candidates),
        "exact_candidates": exact_candidates,
        "elapsed_seconds": round(perf_counter() - started, 6),
    }


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--degrees",
        type=int,
        nargs="+",
        default=[5, 6, 7, 8],
        help="degree bounds to search (default: 5 6 7 8)",
    )
    parser.add_argument(
        "--primes",
        type=int,
        nargs="+",
        default=[11],
        help="prime characteristics, each larger than every degree",
    )
    parser.add_argument(
        "--support-bound",
        type=int,
        choices=(1, 2),
        default=2,
        help="number of collision-kernel directions (default: 2)",
    )
    parser.add_argument(
        "--points",
        type=int,
        default=6,
        help="deterministic rejection points before exact expansion",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="optional JSON output path",
    )
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    if any(degree < 3 for degree in arguments.degrees):
        raise SystemExit("every degree bound must be at least 3")
    if any(
        prime <= max(arguments.degrees) for prime in arguments.primes
    ):
        raise SystemExit("every prime must exceed every searched degree")
    if any(not is_prime(prime) for prime in arguments.primes):
        raise SystemExit("--primes accepts prime numbers only")

    records = []
    for prime in arguments.primes:
        for degree in arguments.degrees:
            record = search_prime_degree(
                prime,
                degree,
                arguments.support_bound,
                arguments.points,
            )
            elapsed_seconds = record.pop("elapsed_seconds")
            records.append(record)
            print(
                "HC4_FF_SEARCH"
                f" p={prime}"
                f" degree={degree}"
                f" directions={record['direction_count']}"
                f" evaluated={record['potentials_evaluated']}"
                f" point_survivors={record['evaluation_survivors']}"
                f" exact_candidates={record['exact_candidate_count']}"
                f" seconds={elapsed_seconds}"
            )

    common_supports: dict[str, list[int]] = {}
    for record in records:
        prime = int(record["prime"])
        for candidate in record["exact_candidates"]:
            signature = json.dumps(
                [
                    record["degree_bound"],
                    candidate["direction_indices"],
                ],
                separators=(",", ":"),
            )
            common_supports.setdefault(signature, []).append(prime)
    repeated_supports = [
        {
            "signature": signature,
            "primes": sorted(set(primes)),
        }
        for signature, primes in sorted(common_supports.items())
        if len(set(primes)) >= 2
    ]

    payload = {
        "status": "bounded finite-field experiment; not a proof",
        "normalization": {
            "quadratic_part": "x0*x1 + x2*x3",
            "collision": "grad(psi)(0) = grad(psi)(1,0,0,0)",
            "constant_hessian_determinant": 1,
        },
        "records": records,
        "supports_occurring_at_multiple_primes": repeated_supports,
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    payload["content_sha256_before_hash_field"] = hashlib.sha256(
        canonical.encode()
    ).hexdigest()
    if arguments.output is not None:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n"
        )
        print(f"WROTE {arguments.output}")


if __name__ == "__main__":
    main()
