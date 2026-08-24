#!/usr/bin/env python3
"""Exact local geometry search around two low-conductor six-root fibers.

This experiment has two deliberately small lanes.

* The four accidental mod-3 pivot directions at each pinned fiber define the
  genus-one slices ``x = +/-T+n``.  A bounded exact search on those slices
  supplies nearby rational parameters carrying at least one prescribed
  non-generic quartic point.
* The complete diameter-at-most-300 Mestre census supplies the closest
  admissible integer root configurations in affine-normalized L1 distance.
  The pinned value of ``T/diameter`` is transported to each neighbor.

Promoted fibers receive exact conductor/root-number data, a bounded quartic
point search, exact quartic-to-Jacobian mapping, and a finite-reduction mod-3
subgroup-rank certificate.  A bounded negative result is not a rank upper
bound.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
from decimal import Decimal
from fractions import Fraction
import hashlib
from itertools import combinations, permutations
import json
from math import gcd, isqrt, lcm
from pathlib import Path
import platform
import re
import subprocess
import sys
from typing import Any, Iterable, Sequence

from mestre_root_tuples import SixRootMestreConstruction
from nagao_linear_sections import LINEAR_COMPANION_SECTIONS
from search_mestre_root_tuple_scale import (
    bounded_quartic_points,
    canonical_signless_points,
    capped_minimal_curve_data,
    point_digest,
    primitive_visible_points,
    quartic_point_to_jacobian,
    sha256_file,
)
from search_mestre_root_tuple_scale_max200 import (
    gf_l_rank_and_pivots,
    mod3_independence_certificate,
)


Q = Fraction
TARGET_LOG_CONDUCTOR = Decimal("182.72")
STACK_BYTES = 512_000_000
DEFAULT_OUTPUT_DIRECTORY = Path(
    "artifacts/local/elliptic-curves/six-root-low-conductor-centers-v1"
)
ROOT = Path(__file__).resolve().parents[2]


@dataclass(frozen=True)
class Seed:
    label: str
    roots: tuple[int, ...]
    parameter: Fraction
    anchor_height: int
    certificate_path: str
    expected_rank: int
    forced_generic_dimension: int


SEEDS = (
    Seed(
        "nagao-rank13-roots-low-conductor-r16",
        (0, 25, 57, 104, 116, 148),
        Q(62, 35),
        1_000_000,
        "archive/elliptic-curves/artifacts/generated-results/"
        "elliptic_mestre_02557104116148_t62_35_rank16_certificate.json",
        16,
        12,
    ),
    Seed(
        "max300-low-conductor-r15",
        (0, 2, 136, 217, 261, 290),
        Q(2),
        5_000,
        "archive/elliptic-curves/artifacts/generated-results/"
        "elliptic_mestre_02136217261290_t2_rank15_certificate.json",
        15,
        11,
    ),
)


def rational_text(value: Fraction) -> str:
    value = Q(value)
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def fraction_record(value: Fraction) -> dict[str, Any]:
    value = Q(value)
    return {
        "value": rational_text(value),
        "numerator": value.numerator,
        "denominator": value.denominator,
    }


def point_record(point: tuple[Fraction, Fraction]) -> dict[str, str]:
    return {"x": rational_text(point[0]), "y": rational_text(point[1])}


def parse_fraction(text: str) -> Fraction:
    return Q(text)


def canonical_digest(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def atomic_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def poly_trim(values: Sequence[Fraction]) -> tuple[Fraction, ...]:
    answer = [Q(value) for value in values]
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return tuple(answer)


def poly_add(left: Sequence[Fraction], right: Sequence[Fraction]) -> tuple[Fraction, ...]:
    length = max(len(left), len(right))
    return poly_trim(
        tuple(
            (Q(left[index]) if index < len(left) else Q(0))
            + (Q(right[index]) if index < len(right) else Q(0))
            for index in range(length)
        )
    )


def poly_multiply(
    left: Sequence[Fraction], right: Sequence[Fraction]
) -> tuple[Fraction, ...]:
    answer = [Q(0)] * (len(left) + len(right) - 1)
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            answer[left_index + right_index] += Q(left_value) * Q(right_value)
    return poly_trim(answer)


def poly_scale(values: Sequence[Fraction], scalar: Fraction) -> tuple[Fraction, ...]:
    return poly_trim(tuple(Q(scalar) * Q(value) for value in values))


def poly_evaluate(values: Sequence[Fraction], argument: Fraction) -> Fraction:
    answer = Q(0)
    for value in reversed(values):
        answer = answer * Q(argument) + Q(value)
    return answer


def interpolate(points: Sequence[tuple[Fraction, Fraction]]) -> tuple[Fraction, ...]:
    """Exact Lagrange interpolation in ascending power order."""

    answer = (Q(0),)
    for index, (x_value, y_value) in enumerate(points):
        basis = (Q(1),)
        denominator = Q(1)
        for other_index, (other_x, _) in enumerate(points):
            if other_index == index:
                continue
            basis = poly_multiply(basis, (-Q(other_x), Q(1)))
            denominator *= Q(x_value) - Q(other_x)
        answer = poly_add(answer, poly_scale(basis, Q(y_value) / denominator))
    return poly_trim(answer)


def bivariate_quartic_coefficients(
    construction: SixRootMestreConstruction,
) -> tuple[tuple[Fraction, ...], ...]:
    """Return coefficients by x-power and then T-power, exactly."""

    answer = []
    for x_power in range(5):
        samples = tuple(
            (Q(parameter), construction.primitive_quartic_coefficients(Q(parameter))[x_power])
            for parameter in range(1, 8)
        )
        polynomial = interpolate(samples)
        check = construction.primitive_quartic_coefficients(Q(8))[x_power]
        if poly_evaluate(polynomial, Q(8)) != check:
            raise AssertionError("quartic coefficient interpolation failed")
        answer.append(polynomial)
    return tuple(answer)


def slice_polynomial(
    bivariate: Sequence[Sequence[Fraction]], slope: Fraction, intercept: Fraction
) -> tuple[Fraction, ...]:
    linear = (Q(intercept), Q(slope))
    power = (Q(1),)
    answer = (Q(0),)
    for coefficient_polynomial in bivariate:
        answer = poly_add(answer, poly_multiply(coefficient_polynomial, power))
        power = poly_multiply(power, linear)
    return poly_trim(answer)


def integer_factorization(value: int) -> Counter[int]:
    value = abs(int(value))
    if value <= 1:
        return Counter()
    if isqrt(value) ** 2 == value:
        root = isqrt(value)
        inner = integer_factorization(root)
        return Counter({prime: 2 * exponent for prime, exponent in inner.items()})
    process = subprocess.run(
        ("factor", "--", str(value)),
        check=True,
        text=True,
        capture_output=True,
        timeout=30,
    )
    factors = [int(item) for item in process.stdout.split(":", 1)[1].split()]
    answer = Counter(factors)
    product = 1
    for prime, exponent in answer.items():
        product *= prime**exponent
    if product != value:
        raise AssertionError("GNU factor output did not multiply back")
    return answer


def square_normalize_rational_polynomial(
    coefficients: Sequence[Fraction],
) -> dict[str, Any]:
    """Clear denominators by a square and remove constant square content."""

    raw = poly_trim(coefficients)
    denominator = 1
    for value in raw:
        denominator = lcm(denominator, Q(value).denominator)
    integer = tuple(int(Q(value) * denominator**2) for value in raw)
    content = gcd(*(abs(value) for value in integer))
    factors = integer_factorization(content)
    square_root = 1
    for prime, exponent in factors.items():
        square_root *= prime ** (exponent // 2)
    square_content = square_root**2
    normalized = tuple(value // square_content for value in integer)
    if tuple(Q(value) * Q(square_root, denominator) ** 2 for value in normalized) != raw:
        raise AssertionError("slice square normalization failed")
    return {
        "raw_coefficients": [rational_text(value) for value in raw],
        "normalized_integer_coefficients": list(normalized),
        "ordinate_scale_original_per_normalized": rational_text(Q(square_root, denominator)),
        "raw_degree": len(raw) - 1,
        "normalized_degree": len(poly_trim(tuple(Q(value) for value in normalized))) - 1,
        "integer_content": str(content),
        "removed_square_content": str(square_content),
        "residual_content": str(content // square_content),
    }


def evaluate_quartic(coefficients: Sequence[Fraction], x_value: Fraction) -> Fraction:
    return poly_evaluate(coefficients, x_value)


def pool_with_sources(
    construction: SixRootMestreConstruction,
    parameter: Fraction,
    searched: Sequence[tuple[Fraction, Fraction]],
    forced: Sequence[tuple[Fraction, Fraction]] = (),
) -> tuple[
    tuple[tuple[Fraction, Fraction], ...],
    tuple[dict[str, Any], ...],
    int,
    int,
]:
    visible = primitive_visible_points(construction, parameter)
    tagged = [
        (f"visible-{index:02d}", point) for index, point in enumerate(visible)
    ]
    tagged.extend((f"forced-{index:02d}", point) for index, point in enumerate(forced))
    tagged.extend((f"searched-{index:04d}", point) for index, point in enumerate(searched))
    by_jacobian_x: dict[Fraction, tuple[tuple[Fraction, Fraction], dict[str, Any]]] = {}
    for label, quartic_point in tagged:
        jacobian_point = quartic_point_to_jacobian(construction, parameter, quartic_point)
        by_jacobian_x.setdefault(
            jacobian_point[0],
            (
                jacobian_point,
                {
                    "source": label,
                    "quartic_point": point_record(quartic_point),
                    "jacobian_point": point_record(jacobian_point),
                },
            ),
        )
    items = tuple(by_jacobian_x.values())
    sources = tuple(item[1] for item in items)
    visible_columns = sum(item["source"].startswith("visible-") for item in sources)
    prescribed_columns = visible_columns + sum(
        item["source"].startswith("forced-") for item in sources
    )
    return (
        tuple(item[0] for item in items),
        sources,
        visible_columns,
        prescribed_columns,
    )


def rank_restricted_columns(certificate: dict[str, Any], column_count: int) -> int:
    rows = []
    for signature in certificate["signatures"]:
        rows.extend(tuple(row[:column_count]) for row in signature["rows"])
    rank, _ = gf_l_rank_and_pivots(rows, column_count, 3)
    return rank


def generic_x_values(seed: Seed, parameter: Fraction) -> set[Fraction]:
    construction = SixRootMestreConstruction(tuple(Q(value) for value in seed.roots))
    answer = {point[0] for point in primitive_visible_points(construction, parameter)}
    if seed.label.startswith("nagao-rank13"):
        answer.update(section.point(parameter)[0] for section in LINEAR_COMPANION_SECTIONS)
    return answer


def anchor_analysis(seed: Seed) -> dict[str, Any]:
    construction = SixRootMestreConstruction(tuple(Q(value) for value in seed.roots))
    parameter = seed.parameter
    raw = bounded_quartic_points(
        construction.primitive_quartic_coefficients(parameter),
        height_bound=seed.anchor_height,
        timeout=180,
        stack_bytes=STACK_BYTES,
    )
    searched = canonical_signless_points(raw)
    points, sources, visible_count, _ = pool_with_sources(construction, parameter, searched)
    certificate = mod3_independence_certificate(
        construction.primitive_jacobian_coefficients(parameter),
        points,
        prime_bound=1000 if seed.expected_rank == 16 else 499,
    )
    rank = certificate["combined_exact_rank_over_F3"]
    if rank != seed.expected_rank:
        raise AssertionError(f"{seed.label} anchor rank replay changed")
    pivot_indices = [index - 1 for index in certificate["independent_subset_indices_one_based"]]
    generic_x = generic_x_values(seed, parameter)
    pivots = []
    accidentals = []
    for index in pivot_indices:
        source = sources[index]
        quartic_point = (
            parse_fraction(source["quartic_point"]["x"]),
            parse_fraction(source["quartic_point"]["y"]),
        )
        generic = quartic_point[0] in generic_x
        record = {"pool_index_one_based": index + 1, **source, "generic_abscissa": generic}
        pivots.append(record)
        if not generic:
            accidentals.append(quartic_point)
    expected_accidental = seed.expected_rank - seed.forced_generic_dimension
    if len(accidentals) != expected_accidental:
        raise AssertionError(
            f"{seed.label} recovered {len(accidentals)} accidental pivots, expected {expected_accidental}"
        )
    certificate_path = ROOT / seed.certificate_path
    frozen = json.loads(certificate_path.read_text())
    return {
        "label": seed.label,
        "roots": list(seed.roots),
        "parameter": rational_text(parameter),
        "parameter_height": max(abs(parameter.numerator), parameter.denominator),
        "anchor_search_height": seed.anchor_height,
        "searched_signless_quartic_point_count": len(searched),
        "pool_point_count_modulo_inverse": len(points),
        "pool_point_sha256": point_digest(points),
        "forced_generic_dimension": seed.forced_generic_dimension,
        "exact_specialization_rank_lower_bound": rank,
        "accidental_pivot_count": len(accidentals),
        "accidental_pivot_points": [point_record(point) for point in accidentals],
        "pivot_records": pivots,
        "finite_reduction_certificate": certificate,
        "frozen_certificate": {
            "path": seed.certificate_path,
            "sha256": sha256_file(certificate_path),
            "result_sha256": frozen.get("result_sha256"),
            "log_conductor": frozen["curve"]["log_conductor"],
            "conductor": frozen["curve"]["conductor"],
            "root_number": frozen["curve"]["root_number"],
        },
    }


def differences(roots: Sequence[int]) -> Counter[int]:
    return Counter(right - left for left, right in combinations(roots, 2))


def perfect_matchings(values: tuple[int, ...]) -> tuple[tuple[tuple[int, int], ...], ...]:
    if not values:
        return ((),)
    first = values[0]
    answer = []
    for index in range(1, len(values)):
        second = values[index]
        remainder = values[1:index] + values[index + 1 :]
        for tail in perfect_matchings(remainder):
            answer.append(((first, second), *tail))
    return tuple(answer)


def determinant3(matrix: Sequence[Sequence[Fraction]]) -> Fraction:
    a, b, c = matrix
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def mobius_matrix(
    source: Sequence[Fraction], target: Sequence[Fraction]
) -> tuple[Fraction, Fraction, Fraction, Fraction]:
    rows = tuple(
        (Q(x_value), Q(1), -Q(y_value) * Q(x_value), -Q(y_value))
        for x_value, y_value in zip(source, target)
    )
    vector = []
    for omitted in range(4):
        minor = tuple(tuple(row[index] for index in range(4) if index != omitted) for row in rows)
        vector.append((Q(-1) if omitted % 2 else Q(1)) * determinant3(minor))
    if all(value == 0 for value in vector):
        raise AssertionError("three point pairs did not determine a Mobius map")
    return tuple(vector)  # type: ignore[return-value]


def normalized_matrix(matrix: Sequence[Fraction]) -> tuple[int, ...]:
    denominator = lcm(*(Q(value).denominator for value in matrix))
    integers = [int(Q(value) * denominator) for value in matrix]
    common = gcd(*(abs(value) for value in integers))
    integers = [value // common for value in integers]
    first = next(value for value in integers if value)
    if first < 0:
        integers = [-value for value in integers]
    return tuple(integers)


def projective_automorphisms(roots: Sequence[int]) -> tuple[tuple[int, ...], ...]:
    root_values = tuple(Q(value) for value in roots)
    root_set = set(root_values)
    maps = set()
    source = root_values[:3]
    for target in permutations(root_values, 3):
        matrix = mobius_matrix(source, target)
        a_value, b_value, c_value, d_value = matrix
        images = []
        valid = True
        for root in root_values:
            denominator = c_value * root + d_value
            if denominator == 0:
                valid = False
                break
            images.append((a_value * root + b_value) / denominator)
        if valid and set(images) == root_set:
            maps.add(normalized_matrix(matrix))
    return tuple(sorted(maps))


def cross_ratio_orbit(value: Fraction) -> tuple[Fraction, ...]:
    value = Q(value)
    return tuple(
        sorted(
            {
                value,
                1 / value,
                1 - value,
                1 / (1 - value),
                value / (value - 1),
                (value - 1) / value,
            }
        )
    )


def root_invariants(roots: tuple[int, ...]) -> dict[str, Any]:
    diameter = roots[-1] - roots[0]
    mean = Q(sum(roots), len(roots))
    difference_counts = differences(roots)
    matching_records = []
    for matching in perfect_matchings(tuple(roots)):
        centers = tuple(Q(left + right, 2) for left, right in matching)
        lengths = tuple(right - left for left, right in matching)
        center_mean = sum(centers, Q(0)) / 3
        dispersion = sum((center - center_mean) ** 2 for center in centers)
        matching_records.append(
            {
                "pairs": [list(pair) for pair in matching],
                "centers": [rational_text(value) for value in centers],
                "lengths": list(lengths),
                "center_spread": rational_text(max(centers) - min(centers)),
                "center_variance_sum": rational_text(dispersion),
            }
        )
    matching_records.sort(
        key=lambda item: (
            Q(item["center_spread"]),
            Q(item["center_variance_sum"]),
            item["pairs"],
        )
    )
    cross_ratios = []
    for quadruple in combinations(tuple(Q(value) for value in roots), 4):
        a_value, b_value, c_value, d_value = quadruple
        value = (a_value - c_value) * (b_value - d_value) / (
            (a_value - d_value) * (b_value - c_value)
        )
        orbit = cross_ratio_orbit(value)
        cross_ratios.append(rational_text(orbit[0]))
    automorphisms = projective_automorphisms(roots)
    construction = SixRootMestreConstruction(tuple(Q(value) for value in roots))
    return {
        "roots": list(roots),
        "diameter": diameter,
        "mean": rational_text(mean),
        "centered_roots": [rational_text(Q(root) - mean) for root in roots],
        "normalized_roots": [rational_text(Q(root - roots[0], diameter)) for root in roots],
        "successive_gaps": [right - left for left, right in zip(roots, roots[1:])],
        "difference_multiplicities": [
            {"difference": difference, "multiplicity": multiplicity}
            for difference, multiplicity in sorted(difference_counts.items())
        ],
        "reflection_symmetric": construction.is_reflection_symmetric,
        "mestre_quartic_obstruction": rational_text(construction.quartic_condition),
        "quartic_content": rational_text(construction.quartic_content),
        "collision_parameters": [rational_text(value) for value in construction.collision_parameters()],
        "five_most_balanced_pairings": matching_records[:5],
        "projective_cross_ratio_orbit_representatives": sorted(cross_ratios, key=Q),
        "projective_automorphism_group_order": len(automorphisms),
        "projective_automorphism_matrices": [list(matrix) for matrix in automorphisms],
    }


def build_slice_records(seed: Seed, analysis: dict[str, Any], *, height_bound: int) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    construction = SixRootMestreConstruction(tuple(Q(value) for value in seed.roots))
    bivariate = bivariate_quartic_coefficients(construction)
    slices = []
    by_parameter: dict[Fraction, dict[str, Any]] = {}
    for accidental_index, point_data in enumerate(analysis["accidental_pivot_points"]):
        source_point = (Q(point_data["x"]), Q(point_data["y"]))
        for slope in (-1, 1):
            intercept = source_point[0] - slope * seed.parameter
            polynomial = slice_polynomial(bivariate, Q(slope), intercept)
            normalized = square_normalize_rational_polynomial(polynomial)
            if normalized["normalized_degree"] != 4:
                raise AssertionError("an x=+/-T accidental slice ceased to be genus one")
            normalized_coefficients = tuple(
                Q(value) for value in normalized["normalized_integer_coefficients"]
            )
            raw = bounded_quartic_points(
                normalized_coefficients,
                height_bound=height_bound,
                timeout=90,
                stack_bytes=STACK_BYTES,
            )
            signless = canonical_signless_points(raw)
            scale = Q(normalized["ordinate_scale_original_per_normalized"])
            mapped = []
            for parameter, normalized_ordinate in signless:
                if parameter == 0:
                    continue
                x_value = Q(slope) * parameter + intercept
                ordinate = scale * normalized_ordinate
                quartic = construction.primitive_quartic_coefficients(parameter)
                if ordinate**2 != evaluate_quartic(quartic, x_value):
                    raise AssertionError("a mapped genus-one-slice point missed the Mestre quartic")
                canonical_parameter = abs(parameter)
                canonical_x = x_value
                if ordinate**2 != evaluate_quartic(
                    construction.primitive_quartic_coefficients(canonical_parameter),
                    canonical_x,
                ):
                    raise AssertionError("T -> -T canonicalization changed the quartic point")
                mapped.append(
                    {
                        "parameter": rational_text(canonical_parameter),
                        "quartic_point": point_record((canonical_x, ordinate)),
                    }
                )
                record = by_parameter.setdefault(
                    canonical_parameter,
                    {
                        "seed_label": seed.label,
                        "roots": list(seed.roots),
                        "parameter": rational_text(canonical_parameter),
                        "parameter_height": max(
                            abs(canonical_parameter.numerator), canonical_parameter.denominator
                        ),
                        "distance_from_anchor": rational_text(abs(canonical_parameter - seed.parameter)),
                        "sources": [],
                        "forced_quartic_points": {},
                    },
                )
                record["sources"].append(
                    {
                        "accidental_index_zero_based": accidental_index,
                        "slope": slope,
                        "intercept": rational_text(intercept),
                    }
                )
                record["forced_quartic_points"][rational_text(canonical_x)] = point_record(
                    (canonical_x, ordinate)
                )
            slices.append(
                {
                    "identifier": f"{seed.label}-a{accidental_index + 1:02d}-m{slope:+d}",
                    "accidental_index_zero_based": accidental_index,
                    "source_point": point_record(source_point),
                    "slope": slope,
                    "intercept": rational_text(intercept),
                    "normalization": normalized,
                    "height_bound": height_bound,
                    "signed_points_returned": len(raw),
                    "signless_points_returned": len(signless),
                    "mapped_points": mapped,
                }
            )
    candidates = []
    for parameter, record in by_parameter.items():
        if parameter == seed.parameter:
            continue
        record["forced_quartic_points"] = list(record["forced_quartic_points"].values())
        record["slice_source_count"] = len(record["sources"])
        candidates.append(record)
    candidates.sort(
        key=lambda item: (
            item["parameter_height"],
            Q(item["distance_from_anchor"]),
            Q(item["parameter"]),
        )
    )
    return slices, candidates


def load_census_roots() -> tuple[tuple[int, ...], ...]:
    max200_path = ROOT / "archive/elliptic-curves/artifacts/generated-results/elliptic_mestre_root_tuple_scale_max200_census.json"
    max300_path = ROOT / "archive/elliptic-curves/artifacts/generated-results/elliptic_mestre_root_tuple_scale_max300_census.json"
    max200 = json.loads(max200_path.read_text())
    max300 = json.loads(max300_path.read_text())
    roots = [
        tuple(values)
        for values in max200["tuple_populations"]["generically_nonsingular_nonreflection_roots"]
    ]
    roots.extend(
        tuple(values)
        for values in max300["tuple_populations"]["genuinely_new_nonsingular_roots"]
    )
    if len(roots) != 2329 or len(set(roots)) != 2329:
        raise AssertionError("the complete max-root-300 nonsingular census changed")
    return tuple(roots)


def affine_root_distance(left: Sequence[int], right: Sequence[int]) -> Fraction:
    return sum(
        (
            abs(Q(left[index] - left[0], left[-1] - left[0]) - Q(right[index] - right[0], right[-1] - right[0]))
            for index in range(1, 5)
        ),
        Q(0),
    )


def neighbor_candidates(seed: Seed, census: Sequence[tuple[int, ...]], keep: int) -> list[dict[str, Any]]:
    ordered = sorted(
        (roots for roots in census if roots != seed.roots),
        key=lambda roots: (affine_root_distance(seed.roots, roots), roots[-1], roots),
    )[:keep]
    answer = []
    for roots in ordered:
        parameter = seed.parameter * Q(roots[-1] - roots[0], seed.roots[-1] - seed.roots[0])
        construction = SixRootMestreConstruction(tuple(Q(value) for value in roots))
        if construction.quartic_condition != 0:
            raise AssertionError("a census neighbor failed Mestre's obstruction")
        degeneracy = construction.visible_point_degeneracy(parameter)
        answer.append(
            {
                "seed_label": seed.label,
                "roots": list(roots),
                "parameter": rational_text(parameter),
                "parameter_height": max(abs(parameter.numerator), parameter.denominator),
                "distance_from_anchor_root_shape": rational_text(
                    affine_root_distance(seed.roots, roots)
                ),
                "transported_normalized_parameter": rational_text(
                    Q(parameter, roots[-1] - roots[0])
                ),
                "source_kind": "affine-nearest admissible root tuple with T/diameter transported exactly",
                "forced_quartic_points": [],
                "admissible": (
                    construction.primitive_quartic_discriminant(parameter) != 0
                    and degeneracy.collision_loss == 0
                    and degeneracy.zero_ordinates == 0
                ),
            }
        )
    return answer


def candidate_key(record: dict[str, Any]) -> str:
    return "r" + "_".join(str(value) for value in record["roots"]) + "_t" + record["parameter"].replace("/", "_")


def candidate_worker(payload: dict[str, Any]) -> dict[str, Any]:
    roots = tuple(int(value) for value in payload["roots"])
    parameter = Q(payload["parameter"])
    construction = SixRootMestreConstruction(tuple(Q(value) for value in roots))
    if construction.primitive_quartic_discriminant(parameter) == 0:
        return {**payload, "status": "excluded singular specialization"}
    degeneracy = construction.visible_point_degeneracy(parameter)
    if degeneracy.collision_loss or degeneracy.zero_ordinates:
        return {
            **payload,
            "status": "excluded visible-point degeneracy",
            "collision_loss": degeneracy.collision_loss,
            "zero_ordinates": degeneracy.zero_ordinates,
        }
    coefficients = construction.primitive_jacobian_coefficients(parameter)
    try:
        global_data = capped_minimal_curve_data(
            coefficients, timeout=20, stack_bytes=STACK_BYTES
        )
        below_target: bool | None = (
            Decimal(global_data["log_conductor"]) < TARGET_LOG_CONDUCTOR
        )
    except Exception as error:
        # Global reduction can be substantially harder than the bounded point
        # and finite-reduction stages.  Preserve that timeout honestly without
        # turning it into a rank-coverage gap.
        global_data = {"status": "error", "error": repr(error)}
        below_target = None
    raw = bounded_quartic_points(
        construction.primitive_quartic_coefficients(parameter),
        height_bound=5_000,
        timeout=60,
        stack_bytes=STACK_BYTES,
    )
    searched = canonical_signless_points(raw)
    forced = tuple(
        (Q(point["x"]), Q(point["y"])) for point in payload.get("forced_quartic_points", [])
    )
    quartic = construction.primitive_quartic_coefficients(parameter)
    if any(y_value**2 != evaluate_quartic(quartic, x_value) for x_value, y_value in forced):
        raise AssertionError("a promoted forced point missed its specialized quartic")
    points, sources, visible_count, prescribed_count = pool_with_sources(
        construction, parameter, searched, forced
    )
    certificate = mod3_independence_certificate(
        coefficients, points, prime_bound=499
    )
    visible_rank = rank_restricted_columns(certificate, visible_count)
    prescribed_rank = rank_restricted_columns(certificate, prescribed_count)
    return {
        **payload,
        "status": "completed exact bounded specialization triage",
        "global_curve": global_data,
        "below_strict_log_conductor_182_72": below_target,
        "height_bound": 5_000,
        "signed_quartic_points_returned": len(raw),
        "signless_quartic_points_returned": len(searched),
        "forced_quartic_point_count": len(forced),
        "pool_point_count_modulo_inverse": len(points),
        "pool_point_sha256": point_digest(points),
        "pool_sources": list(sources),
        "exact_visible_subgroup_dimension": visible_rank,
        "exact_prescribed_subgroup_dimension": prescribed_rank,
        "exact_specialization_rank_lower_bound": certificate["combined_exact_rank_over_F3"],
        "finite_reduction_certificate": certificate,
        "bounded_search_is_not_a_rank_upper_bound": True,
    }


def pareto_frontier(records: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    usable = [
        record
        for record in records
        if record.get("exact_specialization_rank_lower_bound") is not None
        and record.get("log_conductor") is not None
    ]
    answer = []
    for candidate in usable:
        dominated = False
        for other in usable:
            if other is candidate:
                continue
            weak = (
                other["exact_specialization_rank_lower_bound"] >= candidate["exact_specialization_rank_lower_bound"]
                and Decimal(other["log_conductor"]) <= Decimal(candidate["log_conductor"])
                and other["parameter_height"] <= candidate["parameter_height"]
                and other["root_height"] <= candidate["root_height"]
            )
            strict = (
                other["exact_specialization_rank_lower_bound"] > candidate["exact_specialization_rank_lower_bound"]
                or Decimal(other["log_conductor"]) < Decimal(candidate["log_conductor"])
                or other["parameter_height"] < candidate["parameter_height"]
                or other["root_height"] < candidate["root_height"]
            )
            if weak and strict:
                dominated = True
                break
        if not dominated:
            answer.append(candidate)
    return sorted(
        answer,
        key=lambda item: (
            -item["exact_specialization_rank_lower_bound"],
            Decimal(item["log_conductor"]),
            item["parameter_height"],
            item["root_height"],
        ),
    )


def markdown_pareto(records: Sequence[dict[str, Any]]) -> str:
    lines = [
        "# Six-root low-conductor center search: exact Pareto frontier",
        "",
        "A finite-reduction lower bound is unconditional. The bounded H=5,000 search is not a rank upper bound.",
        "",
        "| source | roots | T | forced directions | exact LB | ln N | W | parameter height | root height |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for record in records:
        lines.append(
            "| {source} | {roots} | {parameter} | {forced} | {rank} | {logn} | {rootno} | {ph} | {rh} |".format(
                source=record["source"],
                roots=",".join(str(value) for value in record["roots"]),
                parameter=record["parameter"],
                forced=record["forced_directions"],
                rank=record["exact_specialization_rank_lower_bound"],
                logn=record["log_conductor"],
                rootno=record["root_number"],
                ph=record["parameter_height"],
                rh=record["root_height"],
            )
        )
    lines.append("")
    return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-directory", type=Path, default=DEFAULT_OUTPUT_DIRECTORY)
    parser.add_argument("--slice-height", type=int, default=200_000)
    parser.add_argument("--slice-height-keep", type=int, default=16)
    parser.add_argument("--slice-distance-keep", type=int, default=8)
    parser.add_argument("--neighbor-keep", type=int, default=12)
    parser.add_argument("--workers", type=int, default=8)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if not 1 <= args.workers <= 12:
        raise SystemExit("--workers must lie in [1,12]")
    output = ROOT / args.output_directory
    output.mkdir(parents=True, exist_ok=True)
    certificates_directory = output / "candidate-certificates"
    certificates_directory.mkdir(exist_ok=True)

    seed_analyses = []
    all_slices = []
    slice_candidates_by_seed = {}
    invariant_records = []
    for seed in SEEDS:
        print(f"anchor {seed.label}", flush=True)
        analysis = anchor_analysis(seed)
        seed_analyses.append(analysis)
        invariant_records.append(root_invariants(seed.roots))
        slices, candidates = build_slice_records(seed, analysis, height_bound=args.slice_height)
        all_slices.extend(slices)
        slice_candidates_by_seed[seed.label] = candidates
        print(
            f"slices {seed.label}: {len(slices)} slices, {len(candidates)} non-anchor parameters",
            flush=True,
        )

    census = load_census_roots()
    promoted: dict[tuple[tuple[int, ...], Fraction], dict[str, Any]] = {}
    selection_audit = []
    for seed in SEEDS:
        candidates = slice_candidates_by_seed[seed.label]
        height_selected = candidates[: args.slice_height_keep]
        distance_selected = sorted(
            candidates,
            key=lambda item: (
                Q(item["distance_from_anchor"]), item["parameter_height"], Q(item["parameter"])
            ),
        )[: args.slice_distance_keep]
        chosen = {Q(item["parameter"]): item for item in (*height_selected, *distance_selected)}
        for item in chosen.values():
            item = {**item, "source_kind": "accidental genus-one x=+/-T+n slice"}
            promoted[(tuple(item["roots"]), Q(item["parameter"]))] = item
        neighbors = neighbor_candidates(seed, census, args.neighbor_keep)
        for item in neighbors:
            if item["admissible"]:
                promoted[(tuple(item["roots"]), Q(item["parameter"]))] = item
        selection_audit.append(
            {
                "seed_label": seed.label,
                "all_nonanchor_slice_parameters": len(candidates),
                "height_selected": len(height_selected),
                "distance_selected": len(distance_selected),
                "slice_union_selected": len(chosen),
                "neighbor_selected": sum(item["admissible"] for item in neighbors),
                "neighbor_inadmissible": sum(not item["admissible"] for item in neighbors),
            }
        )
    inputs = sorted(
        promoted.values(),
        key=lambda item: (
            item["seed_label"],
            item["source_kind"],
            tuple(item["roots"]),
            Q(item["parameter"]),
        ),
    )
    input_payload = {
        "scope": "frozen targeted promotion population",
        "selection": selection_audit,
        "candidate_count": len(inputs),
        "candidates": inputs,
    }
    input_payload["result_sha256"] = canonical_digest(input_payload)
    atomic_json(output / "candidate-input.json", input_payload)
    slice_payload = {
        "slice_height_bound": args.slice_height,
        "slices": all_slices,
    }
    slice_payload["result_sha256"] = canonical_digest(slice_payload)
    atomic_json(output / "slice-search.json", slice_payload)

    results = []
    pending = []
    for item in inputs:
        path = certificates_directory / f"{candidate_key(item)}.json"
        if path.exists():
            cached = json.loads(path.read_text())
            if tuple(cached["roots"]) != tuple(item["roots"]) or Q(cached["parameter"]) != Q(item["parameter"]):
                raise AssertionError("a cached candidate certificate key collided")
            if (
                cached.get("status", "").startswith("completed")
                and "exact_prescribed_subgroup_dimension" not in cached
            ):
                source_labels = [item["source"] for item in cached["pool_sources"]]
                visible_count = sum(label.startswith("visible-") for label in source_labels)
                prescribed_count = visible_count + sum(
                    label.startswith("forced-") for label in source_labels
                )
                cached["exact_visible_subgroup_dimension"] = rank_restricted_columns(
                    cached["finite_reduction_certificate"], visible_count
                )
                cached["exact_prescribed_subgroup_dimension"] = rank_restricted_columns(
                    cached["finite_reduction_certificate"], prescribed_count
                )
                cached["result_sha256"] = canonical_digest(
                    {key: value for key, value in cached.items() if key != "result_sha256"}
                )
                atomic_json(path, cached)
            if cached.get("status") == "error":
                pending.append(item)
            else:
                results.append(cached)
        else:
            pending.append(item)
    print(f"candidate exact stage: cached={len(results)} pending={len(pending)}", flush=True)
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(candidate_worker, item): item for item in pending}
        for completed, future in enumerate(as_completed(futures), start=1):
            item = futures[future]
            try:
                record = future.result()
            except Exception as error:
                record = {**item, "status": "error", "error": repr(error)}
            record["result_sha256"] = canonical_digest(record)
            atomic_json(certificates_directory / f"{candidate_key(item)}.json", record)
            results.append(record)
            rank = record.get("exact_specialization_rank_lower_bound")
            print(
                f"exact {completed}/{len(pending)} {candidate_key(item)} "
                f"rank={rank} lnN={record.get('global_curve', {}).get('log_conductor')} "
                f"status={record['status']}",
                flush=True,
            )
            if rank is not None and rank > 16:
                print(f"ALERT exact LB>{16}: {candidate_key(item)} rank={rank}", flush=True)

    results.sort(key=lambda item: (tuple(item["roots"]), Q(item["parameter"])))
    result_summary = {
        "completed": sum(item["status"].startswith("completed") for item in results),
        "errors": sum(item["status"] == "error" for item in results),
        "rank_lower_bound_distribution": dict(
            sorted(
                Counter(
                    item.get("exact_specialization_rank_lower_bound")
                    for item in results
                    if item.get("exact_specialization_rank_lower_bound") is not None
                ).items()
            )
        ),
        "rank_above_16": [
            candidate_key(item)
            for item in results
            if item.get("exact_specialization_rank_lower_bound", -1) > 16
        ],
        "target_qualified_rank_above_16": [
            candidate_key(item)
            for item in results
            if item.get("exact_specialization_rank_lower_bound", -1) > 16
            and item.get("below_strict_log_conductor_182_72")
        ],
    }

    pareto_rows = []
    for seed, analysis in zip(SEEDS, seed_analyses):
        frozen = analysis["frozen_certificate"]
        pareto_rows.append(
            {
                "source": "pinned seed",
                "roots": list(seed.roots),
                "parameter": rational_text(seed.parameter),
                "forced_directions": seed.forced_generic_dimension,
                "exact_specialization_rank_lower_bound": seed.expected_rank,
                "log_conductor": frozen["log_conductor"],
                "root_number": frozen["root_number"],
                "parameter_height": max(abs(seed.parameter.numerator), seed.parameter.denominator),
                "root_height": seed.roots[-1] - seed.roots[0],
            }
        )
    for item in results:
        if (
            item.get("exact_specialization_rank_lower_bound") is None
            or item.get("global_curve", {}).get("log_conductor") is None
        ):
            continue
        pareto_rows.append(
            {
                "source": item["source_kind"],
                "roots": item["roots"],
                "parameter": item["parameter"],
                "forced_directions": item["exact_prescribed_subgroup_dimension"],
                "exact_specialization_rank_lower_bound": item["exact_specialization_rank_lower_bound"],
                "log_conductor": item["global_curve"]["log_conductor"],
                "root_number": item["global_curve"]["root_number"],
                "parameter_height": item["parameter_height"],
                "root_height": item["roots"][-1] - item["roots"][0],
            }
        )
    pareto = pareto_frontier(pareto_rows)
    (output / "pareto.md").write_text(markdown_pareto(pareto))

    artifact = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "completed bounded exact two-seed center search",
        "scope_warning": "bounded slice and H=5000 searches are not rank upper bounds",
        "seed_analyses": seed_analyses,
        "root_invariants": invariant_records,
        "slice_search": {
            "path": str((output / "slice-search.json").relative_to(ROOT)),
            "sha256": sha256_file(output / "slice-search.json"),
            "result_sha256": slice_payload["result_sha256"],
        },
        "candidate_input": {
            "path": str((output / "candidate-input.json").relative_to(ROOT)),
            "sha256": sha256_file(output / "candidate-input.json"),
            "result_sha256": input_payload["result_sha256"],
        },
        "candidate_results": results,
        "summary": result_summary,
        "pareto_frontier": pareto,
        "provenance": {
            "script_path": str(Path(__file__).resolve().relative_to(ROOT)),
            "script_sha256": sha256_file(Path(__file__).resolve()),
            "reproducing_command": "PYTHONPATH=elliptic-curves/cas python3 " + " ".join(sys.argv),
            "max200_census_sha256": sha256_file(
                ROOT / "archive/elliptic-curves/artifacts/generated-results/elliptic_mestre_root_tuple_scale_max200_census.json"
            ),
            "max300_census_sha256": sha256_file(
                ROOT / "archive/elliptic-curves/artifacts/generated-results/elliptic_mestre_root_tuple_scale_max300_census.json"
            ),
        },
        "software": {"python": platform.python_version(), "platform": platform.platform()},
    }
    artifact["result_sha256"] = canonical_digest(
        {key: value for key, value in artifact.items() if key != "generated_at_utc"}
    )
    atomic_json(output / "summary.json", artifact)
    print(
        f"complete candidates={len(results)} distribution={result_summary['rank_lower_bound_distribution']} "
        f"rank>16={result_summary['rank_above_16']} result_sha256={artifact['result_sha256']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
