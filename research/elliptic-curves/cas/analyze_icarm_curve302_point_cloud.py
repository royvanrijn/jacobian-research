#!/usr/bin/env python3
"""Exact, bounded reconstruction probes for the ICARM curve-302 point cloud.

The calculations deliberately do not assume an R17 parent.  Curve 273 is a
second high-rank target and curve 245 is a negative control with a known
Fermigier--Mestre generic-rank-12 parent.

The finite-reduction row spaces, rational interpolation tests, elementary
coordinate identities and fixed-x deformation tests are exact.  Every
negative conclusion is limited to the explicitly declared bounds and models.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
import hashlib
from itertools import combinations
import json
from math import isqrt
from pathlib import Path
import platform
from typing import Iterable, Sequence

from sympy import Matrix, __version__ as sympy_version

import icarm_curve245
import icarm_curve273
import icarm_curve302
from mod2_reduction_independence import (
    finite_add,
    finite_curve_points,
    finite_multiply,
    finite_subtract,
)


Q = Fraction
FinitePoint = tuple[int, int] | None
RationalPoint = tuple[Fraction, Fraction]
ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "icarm_curve302_point_cloud_v1.json"
)
FERMIGIER_CONTROL_PATH = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "record_rank17_fingerprint_calibration_v1.json"
)
PRIME_BOUND = 1000
SMALL_SQUARECLASS_PRIME_BOUND = 97
INTERPOLATION_TOTAL_DEGREE = 6


@dataclass(frozen=True)
class CurveData:
    label: str
    role: str
    coefficients: tuple[Fraction, ...]
    points: tuple[RationalPoint, ...]
    short_coefficients: tuple[Fraction, ...]
    short_points: tuple[RationalPoint, ...]
    first_block_size: int | None


CURVES = (
    CurveData(
        "curve302",
        "primary rank-at-least-31 reconstruction target",
        icarm_curve302.GENERAL_WEIERSTRASS_COEFFICIENTS,
        icarm_curve302.POINTS,
        icarm_curve302.short_coefficients(),
        icarm_curve302.SHORT_POINTS,
        17,
    ),
    CurveData(
        "curve273",
        "rank-at-least-30 comparison target",
        icarm_curve273.GENERAL_WEIERSTRASS_COEFFICIENTS,
        icarm_curve273.POINTS,
        icarm_curve273.short_coefficients(),
        icarm_curve273.SHORT_POINTS,
        17,
    ),
    CurveData(
        "curve245-negative-control",
        "known Fermigier--Mestre generic-rank-12 negative control",
        icarm_curve245.GENERAL_WEIERSTRASS_COEFFICIENTS,
        icarm_curve245.POINTS,
        icarm_curve245.short_coefficients(),
        icarm_curve245.SHORT_POINTS,
        None,
    ),
)


def primes_up_to(bound: int) -> tuple[int, ...]:
    primes: list[int] = []
    for candidate in range(2, bound + 1):
        if all(candidate % prime for prime in primes if prime * prime <= candidate):
            primes.append(candidate)
    return tuple(primes)


def reduce_rational(value: Fraction, prime: int) -> int:
    value = Q(value)
    if value.denominator % prime == 0:
        raise ValueError("point denominator is not invertible")
    return value.numerator * pow(value.denominator, -1, prime) % prime


def rref_mod(rows: Iterable[Sequence[int]], width: int, prime: int) -> list[list[int]]:
    """Canonical reduced row-echelon basis over a prime field."""

    matrix = [
        [int(value) % prime for value in row]
        for row in rows
        if any(int(value) % prime for value in row)
    ]
    if any(len(row) != width for row in matrix):
        raise ValueError("matrix row width changed")
    rank = 0
    for column in range(width):
        pivot = next(
            (index for index in range(rank, len(matrix)) if matrix[index][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        inverse = pow(matrix[rank][column], -1, prime)
        matrix[rank] = [(value * inverse) % prime for value in matrix[rank]]
        for index in range(len(matrix)):
            if index == rank or not matrix[index][column]:
                continue
            multiplier = matrix[index][column]
            matrix[index] = [
                (left - multiplier * right) % prime
                for left, right in zip(matrix[index], matrix[rank])
            ]
        rank += 1
        if rank == len(matrix):
            break
    return matrix[:rank]


def nullspace_mod(rows: Iterable[Sequence[int]], width: int, prime: int) -> list[list[int]]:
    """Canonical basis for coefficient relations annihilated by ``rows``."""

    rref = rref_mod(rows, width, prime)
    pivots = [next(index for index, value in enumerate(row) if value) for row in rref]
    free = [column for column in range(width) if column not in pivots]
    basis = []
    for column in free:
        vector = [0] * width
        vector[column] = 1
        for row, pivot in zip(rref, pivots):
            vector[pivot] = (-row[column]) % prime
        basis.append(vector)
    return basis


def row_rank_restricted(rows: Sequence[Sequence[int]], indices: Sequence[int], prime: int) -> int:
    return len(rref_mod(([row[index] for index in indices] for row in rows), len(indices), prime))


def quotient_signature(
    curve: CurveData, prime: int, relation_prime: int
) -> dict[str, object] | None:
    """Return the canonical point-image row space in E(F_p)/ell E(F_p)."""

    if prime == relation_prime or prime == 2:
        return None
    _, _, _, coefficient_a_q, coefficient_b_q = curve.short_coefficients
    coefficient_a = reduce_rational(coefficient_a_q, prime)
    coefficient_b = reduce_rational(coefficient_b_q, prime)
    discriminant = -16 * (4 * coefficient_a**3 + 27 * coefficient_b**2)
    if discriminant % prime == 0:
        return None
    try:
        reduced_points = [
            (reduce_rational(x, prime), reduce_rational(y, prime))
            for x, y in curve.short_points
        ]
    except ValueError:
        return None

    finite_points = finite_curve_points(coefficient_a, coefficient_b, prime)
    multiples = {
        finite_multiply(point, relation_prime, coefficient_a, prime)
        for point in finite_points
    }

    # Deterministic quotient basis.  The returned RREF below removes dependence
    # on this temporary basis and on the enumeration of quotient generators.
    span: list[tuple[FinitePoint, tuple[int, ...]]] = [(None, ())]
    quotient_dimension = 0
    for point in finite_points:
        if any(
            finite_subtract(point, representative, coefficient_a, prime) in multiples
            for representative, _coordinates in span
        ):
            continue
        old_span = tuple(span)
        span = []
        for scalar in range(relation_prime):
            scalar_point = finite_multiply(point, scalar, coefficient_a, prime)
            for representative, coordinates in old_span:
                span.append(
                    (
                        finite_add(representative, scalar_point, coefficient_a, prime),
                        coordinates + (scalar,),
                    )
                )
        quotient_dimension += 1
    if len(span) * len(multiples) != len(finite_points):
        raise ArithmeticError("finite quotient representatives do not cover the group")
    if len(span) != relation_prime**quotient_dimension:
        raise ArithmeticError("finite quotient is not an elementary prime group")

    rows = [[0] * len(curve.points) for _ in range(quotient_dimension)]
    for point_index, reduced in enumerate(reduced_points):
        coordinates = next(
            (
                coordinates
                for representative, coordinates in span
                if finite_subtract(reduced, representative, coefficient_a, prime)
                in multiples
            ),
            None,
        )
        if coordinates is None:
            raise ArithmeticError("a reduced point missed every quotient coset")
        for coordinate_index, value in enumerate(coordinates):
            rows[coordinate_index][point_index] = value
    canonical_rows = rref_mod(rows, len(curve.points), relation_prime)
    if len(canonical_rows) != quotient_dimension:
        raise ArithmeticError("the submitted points do not span the local quotient")
    return {
        "prime": prime,
        "group_order": len(finite_points),
        "multiple_subgroup_order": len(multiples),
        "ambient_quotient_dimension": quotient_dimension,
        "canonical_image_row_space_rref": canonical_rows,
        "canonical_relation_space_dimension": len(curve.points) - len(canonical_rows),
    }


def pair_separation(rows: Sequence[Sequence[int]], left: int, right: int) -> bool:
    return any(row[left] != row[right] for row in rows)


def kummer_analysis(
    curve: CurveData, relation_prime: int, prime_bound: int
) -> dict[str, object]:
    local_blocks = []
    cumulative_rows: list[list[int]] = []
    selected_primes = []
    rank_trajectory = []
    point_count = len(curve.points)
    pair_counts = {
        pair: {"separated": 0, "equal": 0} for pair in combinations(range(point_count), 2)
    }
    for prime in primes_up_to(prime_bound):
        block = quotient_signature(curve, prime, relation_prime)
        if block is None or not block["ambient_quotient_dimension"]:
            continue
        rows = block["canonical_image_row_space_rref"]
        for pair, counts in pair_counts.items():
            key = "separated" if pair_separation(rows, *pair) else "equal"
            counts[key] += 1
        new_rows = rref_mod(
            [*cumulative_rows, *rows], point_count, relation_prime
        )
        rank_before = len(cumulative_rows)
        rank_after = len(new_rows)
        block["submitted_image_dimension"] = len(rows)
        block["cumulative_image_dimension"] = rank_after
        block["increases_cumulative_rank"] = rank_after > rank_before
        local_blocks.append(block)
        if rank_after > rank_before:
            selected_primes.append(prime)
            rank_trajectory.append(
                {
                    "prime": prime,
                    "rank_before": rank_before,
                    "rank_after": rank_after,
                }
            )
        cumulative_rows = new_rows

    all_indices = list(range(point_count))
    if curve.first_block_size is None:
        blocks = {"all_points": all_indices}
    else:
        split = curve.first_block_size
        blocks = {
            "all_points": all_indices,
            "first_17": list(range(split)),
            "remaining_points": list(range(split, point_count)),
        }
    subset_dimensions = {
        label: row_rank_restricted(cumulative_rows, indices, relation_prime)
        for label, indices in blocks.items()
    }
    if curve.first_block_size is not None:
        subset_dimensions["remaining_dimension_mod_first_17"] = (
            subset_dimensions["all_points"] - subset_dimensions["first_17"]
        )

    nearest_pairs = sorted(
        (
            {
                "point_indices": [left + 1, right + 1],
                **counts,
            }
            for (left, right), counts in pair_counts.items()
        ),
        key=lambda item: (-item["equal"], item["point_indices"]),
    )[:12]

    def pair_profile(pairs: Iterable[tuple[int, int]]) -> dict[str, object]:
        pairs = list(pairs)
        separated = [pair_counts[pair]["separated"] for pair in pairs]
        local_count = len(local_blocks)
        return {
            "pair_count": len(pairs),
            "informative_local_prime_count": local_count,
            "mean_separation_fraction": str(
                Q(sum(separated), len(pairs) * local_count) if pairs and local_count else Q(0)
            ),
            "minimum_separating_prime_count": min(separated, default=0),
            "maximum_separating_prime_count": max(separated, default=0),
        }

    pair_profiles = {"all_pairs": pair_profile(pair_counts)}
    if curve.first_block_size is not None:
        split = curve.first_block_size
        pair_profiles.update(
            {
                "within_first_17": pair_profile(combinations(range(split), 2)),
                "within_remaining_points": pair_profile(
                    combinations(range(split, point_count), 2)
                ),
                "cross_block": pair_profile(
                    (left, right)
                    for left in range(split)
                    for right in range(split, point_count)
                ),
            }
        )
    return {
        "relation_prime": relation_prime,
        "finite_reduction_prime_bound": prime_bound,
        "canonicalization": (
            "RREF of the submitted-point image row space; invariant under a "
            "change of basis in each local quotient"
        ),
        "informative_local_prime_count": len(local_blocks),
        "greedy_rank_increase_primes": selected_primes,
        "rank_trajectory": rank_trajectory,
        "combined_image_dimension": len(cumulative_rows),
        "combined_relation_dimension": point_count - len(cumulative_rows),
        "combined_canonical_image_row_space_rref": cumulative_rows,
        "combined_canonical_relation_space_rref": nullspace_mod(
            cumulative_rows, point_count, relation_prime
        ),
        "subset_dimensions": subset_dimensions,
        "pair_separation_profiles": pair_profiles,
        "nearest_pairs_by_number_of_equal_local_classes": nearest_pairs,
        "local_blocks": local_blocks,
    }


def valuation(value: Fraction, prime: int) -> int:
    value = Q(value)
    numerator = abs(value.numerator)
    denominator = value.denominator
    result = 0
    while numerator % prime == 0:
        numerator //= prime
        result += 1
    while denominator % prime == 0:
        denominator //= prime
        result -= 1
    return result


def local_squareclass(value: Fraction, prime: int) -> str:
    """Canonical Q_p squareclass label for a nonzero rational number."""

    value = Q(value)
    exponent = valuation(value, prime)
    numerator = value.numerator
    denominator = value.denominator
    while numerator % prime == 0:
        numerator //= prime
    while denominator % prime == 0:
        denominator //= prime
    if prime == 2:
        unit = numerator * pow(denominator, -1, 8) % 8
        return f"valuation_parity={exponent % 2};odd_unit_mod_8={unit}"
    unit = numerator * pow(denominator, -1, prime) % prime
    legendre = pow(unit, (prime - 1) // 2, prime)
    sign = 1 if legendre == 1 else -1
    return f"valuation_parity={exponent % 2};unit_legendre={sign}"


def factor_small_integer(value: int) -> dict[str, int]:
    value = abs(int(value))
    factors: dict[str, int] = {}
    for prime in primes_up_to(isqrt(value)):
        if prime * prime > value:
            break
        while value % prime == 0:
            factors[str(prime)] = factors.get(str(prime), 0) + 1
            value //= prime
    if value > 1:
        factors[str(value)] = factors.get(str(value), 0) + 1
    return factors


def is_square_fraction(value: Fraction) -> bool:
    value = Q(value)
    if value < 0:
        return False
    numerator_root = isqrt(value.numerator)
    denominator_root = isqrt(value.denominator)
    return (
        numerator_root * numerator_root == value.numerator
        and denominator_root * denominator_root == value.denominator
    )


def pair_collision_records(
    values: Sequence[Fraction], operation: str
) -> list[dict[str, object]]:
    buckets: dict[Fraction, list[tuple[int, int]]] = {}
    for left, right in combinations(range(len(values)), 2):
        if operation == "sum":
            key = values[left] + values[right]
        elif operation == "difference":
            key = values[left] - values[right]
        elif operation == "product":
            key = values[left] * values[right]
        else:  # pragma: no cover - internal misuse guard
            raise ValueError(operation)
        buckets.setdefault(key, []).append((left + 1, right + 1))
    return [
        {
            "value": str(value),
            "point_pairs": [list(pair) for pair in pairs],
        }
        for value, pairs in sorted(buckets.items(), key=lambda item: str(item[0]))
        if len(pairs) > 1
    ]


def coordinate_patterns(curve: CurveData, squareclass_prime_bound: int) -> dict[str, object]:
    x_values = [Q(point[0]) for point in curve.points]
    short_x_values = [Q(point[0]) for point in curve.short_points]
    denominator_records = []
    denominator_clusters: dict[int, list[int]] = {}
    for index, x_value in enumerate(x_values, 1):
        root = isqrt(x_value.denominator)
        if root * root != x_value.denominator:
            raise ArithmeticError("an x-coordinate denominator ceased to be a square")
        denominator_clusters.setdefault(root, []).append(index)
        denominator_records.append(
            {
                "point_index": index,
                "denominator_root": root,
                "denominator_root_factorization": factor_small_integer(root),
            }
        )

    exact_relations = {
        "arithmetic_progressions": [],
        "geometric_progressions": [],
        "equal_squareclasses_of_x": [],
        "absolute_square_differences": [],
    }
    for left, middle, right in combinations(range(len(x_values)), 3):
        triple = (x_values[left], x_values[middle], x_values[right])
        for center in range(3):
            ends = [position for position in range(3) if position != center]
            if triple[ends[0]] + triple[ends[1]] == 2 * triple[center]:
                indices = [left + 1, middle + 1, right + 1]
                exact_relations["arithmetic_progressions"].append(
                    {
                        "point_indices": indices,
                        "middle_point_index": indices[center],
                    }
                )
            if triple[ends[0]] * triple[ends[1]] == triple[center] ** 2:
                indices = [left + 1, middle + 1, right + 1]
                exact_relations["geometric_progressions"].append(
                    {
                        "point_indices": indices,
                        "middle_point_index": indices[center],
                    }
                )
    for left, right in combinations(range(len(x_values)), 2):
        if is_square_fraction(x_values[left] / x_values[right]):
            exact_relations["equal_squareclasses_of_x"].append([left + 1, right + 1])
        if is_square_fraction(abs(x_values[left] - x_values[right])):
            exact_relations["absolute_square_differences"].append([left + 1, right + 1])

    _, _, _, coefficient_a_q, coefficient_b_q = curve.short_coefficients
    local_clusters = []
    numerator_divisibility_clusters = []
    for prime in primes_up_to(squareclass_prime_bound):
        classes: dict[str, list[int]] = {}
        divisible = []
        for index, value in enumerate(x_values, 1):
            classes.setdefault(local_squareclass(value, prime), []).append(index)
            if value.numerator % prime == 0:
                divisible.append(index)
        repeated_classes = {
            label: indices for label, indices in classes.items() if len(indices) >= 3
        }
        if repeated_classes:
            local_clusters.append(
                {
                    "prime": prime,
                    "good_reduction_on_short_model": (
                        -16
                        * (
                            4 * int(coefficient_a_q) ** 3
                            + 27 * int(coefficient_b_q) ** 2
                        )
                    )
                    % prime
                    != 0,
                    "classes_of_size_at_least_3": repeated_classes,
                }
            )
        if 3 <= len(divisible) < len(x_values):
            numerator_divisibility_clusters.append(
                {
                    "prime": prime,
                    "good_reduction_on_short_model": (
                        -16
                        * (
                            4 * int(coefficient_a_q) ** 3
                            + 27 * int(coefficient_b_q) ** 2
                        )
                    )
                    % prime
                    != 0,
                    "point_indices": divisible,
                }
            )

    return {
        "coordinate_choice": (
            "public minimal-model x for denominator and squareclass diagnostics; "
            "canonical integral short-model X for the deformation test"
        ),
        "denominator_records": denominator_records,
        "repeated_denominator_root_clusters": {
            str(root): indices
            for root, indices in denominator_clusters.items()
            if len(indices) > 1
        },
        "integral_x_point_indices": [
            index for index, value in enumerate(x_values, 1) if value.denominator == 1
        ],
        "small_prime_bound_for_local_squareclasses": squareclass_prime_bound,
        "small_prime_numerator_divisibility_clusters": numerator_divisibility_clusters,
        "local_squareclass_clusters": local_clusters,
        "exact_low_complexity_relations": exact_relations,
        "repeated_pair_sums": pair_collision_records(x_values, "sum"),
        "repeated_oriented_pair_differences": pair_collision_records(
            x_values, "difference"
        ),
        "repeated_pair_products": pair_collision_records(x_values, "product"),
        "short_x_values_sha256": hashlib.sha256(
            json.dumps([str(value) for value in short_x_values], separators=(",", ":")).encode()
        ).hexdigest(),
    }


def normalize_null_vector(vector: Sequence[object]) -> list[Fraction]:
    values = [Q(int(value.p), int(value.q)) for value in vector]
    first = next(value for value in values if value)
    return [value / first for value in values]


def fit_rational_function(
    predictors: Sequence[int],
    targets: Sequence[Fraction],
    numerator_degree: int,
    denominator_degree: int,
) -> tuple[list[Fraction], list[Fraction]] | None:
    rows = []
    for predictor, target in zip(predictors, targets):
        rows.append(
            [
                *[Q(predictor) ** degree for degree in range(numerator_degree + 1)],
                *[
                    -Q(target) * Q(predictor) ** degree
                    for degree in range(denominator_degree + 1)
                ],
            ]
        )
    matrix = Matrix(rows)
    nullspace = matrix.nullspace()
    if len(nullspace) != 1:
        return None
    vector = normalize_null_vector(list(nullspace[0]))
    numerator = vector[: numerator_degree + 1]
    denominator = vector[numerator_degree + 1 :]
    if all(not value for value in denominator):
        return None
    for predictor in predictors:
        if not sum(
            coefficient * Q(predictor) ** degree
            for degree, coefficient in enumerate(denominator)
        ):
            return None
    return numerator, denominator


def evaluate_rational_function(
    model: tuple[Sequence[Fraction], Sequence[Fraction]], predictor: int
) -> Fraction | None:
    numerator, denominator = model
    top = sum(
        coefficient * Q(predictor) ** degree
        for degree, coefficient in enumerate(numerator)
    )
    bottom = sum(
        coefficient * Q(predictor) ** degree
        for degree, coefficient in enumerate(denominator)
    )
    return None if not bottom else top / bottom


def interpolation_group(
    curve: CurveData, label: str, indices: Sequence[int], total_degree: int
) -> dict[str, object]:
    if len(indices) < 8:
        raise ValueError("held-out interpolation groups need at least eight points")
    # Seven alternating points are frozen training data.  Degree allocations
    # summing to six are then minimally determined and evaluated only on the
    # remaining, genuinely held-out points.
    training = list(indices[::2][: total_degree + 1])
    held_out = [index for index in indices if index not in training]
    if len(training) != total_degree + 1:
        raise ValueError("interpolation training set has the wrong size")
    predictors = [index + 1 for index in training]
    targets_by_coordinate = {
        "x": [curve.short_points[index][0] for index in training],
        "y": [curve.short_points[index][1] for index in training],
    }
    models: dict[str, list[dict[str, object]]] = {"x": [], "y": []}
    for coordinate, targets in targets_by_coordinate.items():
        for degree_sum in range(total_degree + 1):
            for numerator_degree in range(degree_sum + 1):
                denominator_degree = degree_sum - numerator_degree
                model = fit_rational_function(
                    predictors,
                    targets,
                    numerator_degree,
                    denominator_degree,
                )
                if model is None:
                    continue
                evaluations = [
                    evaluate_rational_function(model, index + 1) for index in held_out
                ]
                matches = [
                    index + 1
                    for index, prediction in zip(held_out, evaluations)
                    if prediction is not None
                    and prediction
                    == curve.short_points[index][0 if coordinate == "x" else 1]
                ]
                coefficient_payload = [
                    [str(value) for value in model[0]],
                    [str(value) for value in model[1]],
                ]
                models[coordinate].append(
                    {
                        "numerator_degree": numerator_degree,
                        "denominator_degree": denominator_degree,
                        "held_out_exact_match_count": len(matches),
                        "held_out_exact_match_point_indices": matches,
                        "held_out_pole_point_indices": [
                            index + 1
                            for index, prediction in zip(held_out, evaluations)
                            if prediction is None
                        ],
                        "normalized_coefficient_sha256": hashlib.sha256(
                            json.dumps(
                                coefficient_payload, separators=(",", ":")
                            ).encode()
                        ).hexdigest(),
                    }
                )
    best = {
        coordinate: max(
            (item["held_out_exact_match_count"] for item in coordinate_models),
            default=0,
        )
        for coordinate, coordinate_models in models.items()
    }
    joint_best = 0
    joint_indices: list[int] = []
    for x_model in models["x"]:
        x_matches = set(x_model["held_out_exact_match_point_indices"])
        for y_model in models["y"]:
            common = sorted(
                x_matches.intersection(y_model["held_out_exact_match_point_indices"])
            )
            if len(common) > joint_best:
                joint_best = len(common)
                joint_indices = common
    return {
        "label": label,
        "point_indices": [index + 1 for index in indices],
        "predictor": "one-based public submission index",
        "coordinate_model": "rational function P(i)/Q(i) on integral short X,Y",
        "training_point_indices": [index + 1 for index in training],
        "held_out_point_indices": [index + 1 for index in held_out],
        "degree_protocol": f"deg(P)+deg(Q)<={total_degree}",
        "best_held_out_exact_match_count": best,
        "best_joint_xy_held_out_exact_match_count": joint_best,
        "best_joint_xy_held_out_point_indices": joint_indices,
        "models": models,
    }


def interpolation_analysis(curve: CurveData, total_degree: int) -> dict[str, object]:
    groups = [("all_points", list(range(len(curve.points))))]
    if curve.first_block_size is not None:
        groups.extend(
            (
                ("first_17", list(range(curve.first_block_size))),
                (
                    "remaining_points",
                    list(range(curve.first_block_size, len(curve.points))),
                ),
            )
        )
    return {
        "status": "exact training/held-out evaluation; no approximate tolerance",
        "groups": [
            interpolation_group(curve, label, indices, total_degree)
            for label, indices in groups
        ],
    }


def quadratic_roots(coefficients: Sequence[Fraction]) -> set[Fraction]:
    constant, linear, quadratic = map(Q, coefficients)
    if quadratic:
        discriminant = linear * linear - 4 * quadratic * constant
        if not is_square_fraction(discriminant):
            return set()
        numerator_root = isqrt(discriminant.numerator)
        denominator_root = isqrt(discriminant.denominator)
        root = Q(numerator_root, denominator_root)
        return {
            (-linear + root) / (2 * quadratic),
            (-linear - root) / (2 * quadratic),
        }
    if linear:
        return {-constant / linear}
    return set()


def polynomial_subtract(left: Sequence[Fraction], right: Sequence[Fraction]) -> list[Fraction]:
    return [Q(a) - Q(b) for a, b in zip(left, right)]


def polynomial_scale(values: Sequence[Fraction], scalar: Fraction) -> list[Fraction]:
    return [Q(scalar) * Q(value) for value in values]


def deformation_analysis(curve: CurveData) -> dict[str, object]:
    """Search exact fixed-X quadratic pencils with linear point sections.

    For short ``Y^2=X^3+A*X+B``, put

      A(t)=A+a*t+c*t^2, B(t)=B+b*t+d*t^2,
      Y_i(t)=Y_i+h_i*t, X_i(t)=X_i.

    The first-order equation gives ``h_i=(a*X_i+b)/(2*Y_i)``.  The second
    order closes exactly iff the points ``(X_i,h_i^2)`` are collinear, with
    line ``c*X+d``.  A direction [a:b] inferred from a triple is tested
    against every other point, so points beyond the seed triple are held-out
    exact hits.
    """

    x_values = [Q(point[0]) for point in curve.short_points]
    y_values = [Q(point[1]) for point in curve.short_points]
    if len(set(x_values)) != len(x_values) or any(not value for value in y_values):
        raise ArithmeticError("the fixed-X deformation test needs distinct X and nonzero Y")

    q_polynomials = [
        [Q(1, 4 * y * y), Q(x, 2 * y * y), Q(x * x, 4 * y * y)]
        for x, y in zip(x_values, y_values)
    ]
    directions: set[tuple[int, int]] = set()
    seed_triples_tested = 0
    for left, middle, right in combinations(range(len(x_values)), 3):
        seed_triples_tested += 1
        first = polynomial_scale(
            polynomial_subtract(q_polynomials[middle], q_polynomials[left]),
            x_values[right] - x_values[left],
        )
        second = polynomial_scale(
            polynomial_subtract(q_polynomials[right], q_polynomials[left]),
            x_values[middle] - x_values[left],
        )
        condition = polynomial_subtract(first, second)
        for root in quadratic_roots(condition):
            directions.add((root.numerator, root.denominator))
        if condition[2] == 0:
            directions.add((1, 0))

    candidates = []
    maximum_size = 2
    for a, b in sorted(directions):
        squares = [
            Q((a * x + b) ** 2, 4 * y * y) for x, y in zip(x_values, y_values)
        ]
        line_sets: dict[tuple[Fraction, Fraction], set[int]] = {}
        for left, right in combinations(range(len(x_values)), 2):
            slope = (squares[right] - squares[left]) / (
                x_values[right] - x_values[left]
            )
            intercept = squares[left] - slope * x_values[left]
            line_sets.setdefault((slope, intercept), set()).update((left, right))
        for (slope, intercept), indices in line_sets.items():
            complete = {
                index
                for index, (x, square) in enumerate(zip(x_values, squares))
                if square == slope * x + intercept
            }
            if len(complete) < 3 or complete != indices:
                continue
            maximum_size = max(maximum_size, len(complete))
            if len(complete) >= 4:
                candidates.append(
                    {
                        "direction_a_b": [a, b],
                        "quadratic_coefficient_c": str(slope),
                        "quadratic_coefficient_d": str(intercept),
                        "preserved_point_indices": [index + 1 for index in sorted(complete)],
                        "seed_size": 3,
                        "held_out_exact_hit_count": len(complete) - 3,
                    }
                )
    unique = {
        (
            tuple(candidate["direction_a_b"]),
            tuple(candidate["preserved_point_indices"]),
        ): candidate
        for candidate in candidates
    }
    candidates = sorted(
        unique.values(),
        key=lambda item: (
            -len(item["preserved_point_indices"]),
            item["preserved_point_indices"],
            item["direction_a_b"],
        ),
    )
    return {
        "model": (
            "fixed short-model X; A(t)=A+a*t+c*t^2, B(t)=B+b*t+d*t^2; "
            "Y_i(t)=Y_i+h_i*t"
        ),
        "first_order_warning": (
            "unrestricted first-order formal lifting is automatic at every supplied "
            "point and therefore has no reconstruction value"
        ),
        "exact_closure_condition": (
            "h_i=(a*X_i+b)/(2*Y_i), and all (X_i,h_i^2) in a preserved "
            "subset lie on the line c*X+d"
        ),
        "seed_triples_tested": seed_triples_tested,
        "rational_projective_directions_from_seed_triples": len(directions),
        "maximum_preserved_subset_size": maximum_size,
        "exact_bounded_conclusion": (
            "no rational direction preserves three supplied points"
            if maximum_size == 2
            else "one or more rational directions preserve at least three supplied points"
        ),
        "candidates_with_at_least_one_held_out_hit": candidates,
    }


def curve_digest(curve: CurveData) -> str:
    payload = {
        "coefficients": [str(value) for value in curve.coefficients],
        "points": [[str(x), str(y)] for x, y in curve.points],
    }
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def mestre_two_cover_coset_calibration() -> dict[str, object]:
    """Verify the common mod-2 coset of the known quartic control points."""

    payload = json.loads(FERMIGIER_CONTROL_PATH.read_text(encoding="utf-8"))
    control = payload["negative_control"]["exact_true_generic_subgroup"]
    coordinates = [
        [int(value) for value in row]
        for row in control["coordinates_by_generic_point_in_public_basis"]
    ]
    if len(coordinates) != 13 or any(len(row) != 20 for row in coordinates):
        raise ArithmeticError("the Fermigier control coordinate matrix changed")
    parities = [[value % 2 for value in row] for row in coordinates]
    common_parity = parities[0]
    if any(row != common_parity for row in parities):
        raise ArithmeticError("the soluble-quartic common-coset control failed")
    visible_sum = [sum(row[column] for row in coordinates[:12]) for column in range(20)]
    if any(visible_sum):
        raise ArithmeticError("the twelve visible Mestre points ceased to sum to zero")
    return {
        "general_fact": (
            "after choosing a rational point on a soluble 2-covering quartic, "
            "its map to the Jacobian has image in one affine coset of 2E(Q)"
        ),
        "control": "curve245 reconstructed Fermigier--Mestre quartic",
        "quartic_point_count": 13,
        "twelve_visible_points_sum_to_zero_in_public_basis": True,
        "all_thirteen_transport_to_one_mod2_coset": True,
        "common_public_basis_parity_vector": common_parity,
        "source_artifact": str(FERMIGIER_CONTROL_PATH.relative_to(ROOT)),
        "source_artifact_sha256": hashlib.sha256(
            FERMIGIER_CONTROL_PATH.read_bytes()
        ).hexdigest(),
        "curve302_raw_point_implication": (
            "the full-rank mod-2 code separates every pair of submitted points, "
            "so no two raw submitted points are images from one soluble quartic "
            "2-cover; a hidden quartic construction would require nontrivial "
            "Mordell--Weil combinations or a different map"
        ),
    }


def analyze_curve(
    curve: CurveData, prime_bound: int, squareclass_prime_bound: int
) -> dict[str, object]:
    on_curve = {
        "curve302": icarm_curve302.on_curve,
        "curve273": icarm_curve273.on_curve,
        "curve245-negative-control": icarm_curve245.on_curve,
    }[curve.label]
    if any(not on_curve(point) for point in curve.points):
        raise ArithmeticError(f"{curve.label}: a pinned point is off the curve")
    return {
        "label": curve.label,
        "role": curve.role,
        "point_count": len(curve.points),
        "input_sha256": curve_digest(curve),
        "finite_reduction_kummer_codes": {
            "mod_2": kummer_analysis(curve, 2, prime_bound),
            "mod_3": kummer_analysis(curve, 3, prime_bound),
        },
        "coordinate_patterns": coordinate_patterns(curve, squareclass_prime_bound),
        "held_out_low_degree_interpolation": interpolation_analysis(
            curve, INTERPOLATION_TOTAL_DEGREE
        ),
        "fixed_x_quadratic_deformations": deformation_analysis(curve),
    }


def build_payload(prime_bound: int, squareclass_prime_bound: int) -> dict[str, object]:
    return {
        "schema": "elliptic-curves.icarm-curve302-point-cloud.v1",
        "status": (
            "exact bounded reconstruction probes; no construction provenance or "
            "Mordell--Weil upper bound is claimed"
        ),
        "primary_question": (
            "Does the 31-point configuration itself expose a low-complexity parent, "
            "without assuming a 17+14 or R17 decomposition?"
        ),
        "bounds": {
            "finite_good_reduction_prime_bound": prime_bound,
            "local_squareclass_prime_bound": squareclass_prime_bound,
            "interpolation_total_degree": INTERPOLATION_TOTAL_DEGREE,
            "deformation_model": "fixed-X quadratic coefficient pencil with linear Y-sections",
        },
        "software": {
            "python": platform.python_version(),
            "sympy": sympy_version,
        },
        "generation": {
            "command": (
                "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
                "elliptic-curves/cas/analyze_icarm_curve302_point_cloud.py"
            ),
            "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "control_input_sha256": hashlib.sha256(
                FERMIGIER_CONTROL_PATH.read_bytes()
            ).hexdigest(),
        },
        "mestre_two_cover_coset_calibration": mestre_two_cover_coset_calibration(),
        "curves": [
            analyze_curve(curve, prime_bound, squareclass_prime_bound)
            for curve in CURVES
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime-bound", type=int, default=PRIME_BOUND)
    parser.add_argument(
        "--squareclass-prime-bound",
        type=int,
        default=SMALL_SQUARECLASS_PRIME_BOUND,
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.prime_bound < 50 or args.squareclass_prime_bound < 11:
        raise SystemExit("declared bounds are too small for the pinned analysis")
    payload = build_payload(args.prime_bound, args.squareclass_prime_bound)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if args.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"FAIL: {args.output} differs from recomputation")
        print(
            f"PASS|{args.output}|sha256="
            f"{hashlib.sha256(rendered.encode()).hexdigest()}"
        )
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered, encoding="utf-8")
    print(
        f"WROTE|{args.output}|sha256="
        f"{hashlib.sha256(rendered.encode()).hexdigest()}"
    )


if __name__ == "__main__":
    main()
