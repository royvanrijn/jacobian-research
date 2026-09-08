#!/usr/bin/env sage -python
"""Audit the coordinate-height mechanism in completed half-lattice searches.

The pointed quartic attached to ``Q`` is a degree-two coordinate chart on the
same elliptic curve.  This script reconstructs the horizontal PARI
minimization/reduction map for every chart in the detailed blind ledgers and
records quantities that were available before point search: lattice depth,
quartic coefficient heights, exact binary-quartic invariants, map size, and a
known-basis calibration of reduced-coordinate versus centered canonical
height.

The rank-28 ledger has a complete post-freeze labeling of every non-generic
point returned by the selected charts.  For that ledger only, the script also
performs an explicitly posthoc oracle audit.  It verifies that reduced
coordinate height exactly reconstructs point visibility at the declared
bound and compares it with the Neron--Tate midpoint quantity

    (1/4) * hhat(2*R-Q) = hhat(R-Q/2).

Nothing here assigns Selmer or covering meaning to a chart, and a bounded miss
has no rank, saturation, or point-absence consequence.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import gzip
from hashlib import sha256
import json
import math
from pathlib import Path
import statistics
import sys
from typing import Any, Iterable, Sequence

from cypari2 import Pari


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
ARTIFACTS = ROOT / "artifacts/generated-results/elliptic-curves"
OUTPUT = ARTIFACTS / "half_lattice_height_compression_analysis_v1.json.gz"

RANK28_BLIND = ARTIFACTS / "half_lattice_fake_descent_rank28_blind_v1.json"
RANK28_VERIFY = ARTIFACTS / "half_lattice_fake_descent_rank28_verification_v1.json"
CURVE385_ITERATED = ARTIFACTS / "curve385_iterated_half_lattice_blind_v1.json"
CURVE385_PRIMARY = ARTIFACTS / "curve385_sparse_quotient_rank32_primary_ledger_v1.json.gz"
CURVE398_ADAPTIVE = ARTIFACTS / "curve398_mw16_adaptive_half_lattice_blind_v1.json"
R17_MATRIX_BLIND = ARTIFACTS / "half_lattice_fake_descent_r17_matrix_blind_v1.json"
R17_MATRIX_VERIFY = ARTIFACTS / "half_lattice_fake_descent_r17_matrix_verification_v1.json"
RANK29_CONTROLS_BLIND = ARTIFACTS / "half_lattice_rank29_controls_blind_v1.json"
RANK29_CONTROLS_VERIFY = ARTIFACTS / "half_lattice_rank29_controls_verification_v1.json"

sys.path.insert(0, str(CAS))

from alternate_quartic_covers import short_add  # noqa: E402


Q = Fraction
Point = tuple[Fraction, Fraction]


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    if path.suffix == ".gz":
        with gzip.open(path, "rt") as handle:
            return json.load(handle)
    return json.loads(path.read_text())


def point(record: dict[str, str]) -> Point:
    return Q(record["x"]), Q(record["y"])


def rational_string(value: Fraction | int) -> str:
    value = Q(value)
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def rational_height(value: Fraction | int) -> int:
    value = Q(value)
    return max(abs(value.numerator), value.denominator)


def rational_bits(value: Fraction | int) -> int:
    return rational_height(value).bit_length()


def integer_bits(value: int) -> int:
    return abs(int(value)).bit_length()


def lcm(left: int, right: int) -> int:
    return abs(left * right) // math.gcd(left, right)


def median(values: Sequence[float]) -> float:
    if not values:
        raise ValueError("the median of an empty sequence is undefined")
    return float(statistics.median(values))


def quantiles(values: Sequence[float]) -> dict[str, str]:
    if not values:
        raise ValueError("quantiles require a nonempty sequence")
    ordered = sorted(float(value) for value in values)

    def at(fraction: float) -> float:
        return ordered[round(fraction * (len(ordered) - 1))]

    return {
        "minimum": decimal(at(0.0)),
        "q25": decimal(at(0.25)),
        "median": decimal(at(0.5)),
        "q75": decimal(at(0.75)),
        "maximum": decimal(at(1.0)),
    }


def decimal(value: float) -> str:
    return format(float(value), ".15g")


def binary_rank(values: Iterable[int]) -> int:
    pivots: dict[int, int] = {}
    for item in values:
        value = int(item)
        while value:
            pivot = value.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = value
                break
            value ^= pivots[pivot]
    return len(pivots)


def auc(rows: Sequence[dict[str, Any]], field: str, label: str) -> dict[str, Any] | None:
    positive = [float(row[field]) for row in rows if row[label]]
    negative = [float(row[field]) for row in rows if not row[label]]
    if not positive or not negative:
        return None
    wins = sum(left > right for left in positive for right in negative)
    ties = sum(left == right for left in positive for right in negative)
    raw = (wins + ties / 2) / (len(positive) * len(negative))
    return {
        "auc_larger_predicts_event": decimal(raw),
        "best_orientation_auc": decimal(max(raw, 1 - raw)),
        "better_orientation": "larger" if raw >= 0.5 else "smaller",
        "event_mean": decimal(statistics.mean(positive)),
        "nonevent_mean": decimal(statistics.mean(negative)),
    }


def binary_quartic_invariants(coefficients_ascending: Sequence[Fraction | int]) -> tuple[Fraction, Fraction]:
    """Return classical I,J for a binary quartic in ascending affine order."""

    if len(coefficients_ascending) > 5:
        raise ValueError("a binary quartic has at most five coefficients")
    coefficients = list(map(Q, coefficients_ascending)) + [Q(0)] * (
        5 - len(coefficients_ascending)
    )
    e_value, d_value, c_value, b_value, a_value = coefficients
    invariant_i = 12 * a_value * e_value - 3 * b_value * d_value + c_value**2
    invariant_j = (
        72 * a_value * c_value * e_value
        + 9 * b_value * c_value * d_value
        - 27 * a_value * d_value**2
        - 27 * b_value**2 * e_value
        - 2 * c_value**3
    )
    return invariant_i, invariant_j


def matrix_product(left: Sequence[Sequence[int]], right: Sequence[Sequence[int]]) -> list[list[int]]:
    return [
        [sum(int(left[i][k]) * int(right[k][j]) for k in range(2)) for j in range(2)]
        for i in range(2)
    ]


def primitive_matrix(matrix: Sequence[Sequence[int]]) -> list[list[int]]:
    entries = [int(matrix[i][j]) for i in range(2) for j in range(2)]
    common = 0
    for entry in entries:
        common = math.gcd(common, abs(entry))
    if common == 0:
        raise ArithmeticError("the horizontal model map vanished")
    entries = [entry // common for entry in entries]
    first = next(entry for entry in entries if entry)
    if first < 0:
        entries = [-entry for entry in entries]
    return [entries[:2], entries[2:]]


def discriminant_quartic(
    p_coefficients: Sequence[int], q_coefficients: Sequence[int]
) -> list[int]:
    answer = []
    for index in range(5):
        p_value = p_coefficients[index] if index < len(p_coefficients) else 0
        q_square = sum(
            q_coefficients[j] * q_coefficients[index - j]
            for j in range(len(q_coefficients))
            if 0 <= index - j < len(q_coefficients)
        )
        answer.append(4 * p_value + q_square)
    return answer


def transform_binary_quartic(
    coefficients_ascending: Sequence[int], matrix: Sequence[Sequence[int]]
) -> list[int]:
    """Return ``F(aX+bZ,cX+dZ)`` in ascending affine powers of X."""

    a_value, b_value = map(int, matrix[0])
    c_value, d_value = map(int, matrix[1])
    answer = [0] * 5
    for degree_x, coefficient in enumerate(coefficients_ascending):
        degree_z = 4 - degree_x
        for left_x in range(degree_x + 1):
            left = (
                math.comb(degree_x, left_x)
                * a_value**left_x
                * b_value ** (degree_x - left_x)
            )
            for right_x in range(degree_z + 1):
                right = (
                    math.comb(degree_z, right_x)
                    * c_value**right_x
                    * d_value ** (degree_z - right_x)
                )
                answer[left_x + right_x] += int(coefficient) * left * right
    return answer


def small_unimodular_matrices(bound: int = 2) -> list[list[list[int]]]:
    candidates = []
    seen = set()
    for a_value in range(-bound, bound + 1):
        for b_value in range(-bound, bound + 1):
            for c_value in range(-bound, bound + 1):
                for d_value in range(-bound, bound + 1):
                    if abs(a_value * d_value - b_value * c_value) != 1:
                        continue
                    matrix = primitive_matrix(
                        [[a_value, b_value], [c_value, d_value]]
                    )
                    key = tuple(entry for row in matrix for entry in row)
                    if key not in seen:
                        seen.add(key)
                        candidates.append(matrix)
    candidates.sort(
        key=lambda matrix: (
            max(abs(entry) for row in matrix for entry in row),
            sum(abs(entry) for row in matrix for entry in row),
            tuple(entry for row in matrix for entry in row),
        )
    )
    return candidates


SMALL_UNIMODULAR_MATRICES = small_unimodular_matrices()


def pari_polynomial_coefficients(polynomial: Any, length: int = 5) -> list[int]:
    return [int(polynomial.polcoef(index)) for index in range(length)]


def pari_transform(
    pari: Pari,
    raw_coefficients: Sequence[Fraction],
    recorded_reduced: dict[str, Any],
) -> dict[str, Any]:
    denominator = 1
    for coefficient in raw_coefficients:
        denominator = lcm(denominator, coefficient.denominator)
    integral = [int(coefficient * denominator * denominator) for coefficient in raw_coefficients]
    polynomial = "+".join(
        f"({coefficient})*x^{index}" for index, coefficient in enumerate(integral)
    )
    result = pari(
        "my(m1,m2,C0=["
        + polynomial
        + ",0],C1,C2);"
        + "C1=hyperellminimalmodel(C0,&m1);"
        + "C2=hyperellred(C1,&m2);[C2,m1,m2]"
    )
    curve, first, second = result[0], result[1], result[2]
    reduced_p = pari_polynomial_coefficients(curve[0])
    reduced_q = pari_polynomial_coefficients(curve[1])
    while reduced_p and reduced_p[-1] == 0:
        reduced_p.pop()
    while reduced_q and reduced_q[-1] == 0:
        reduced_q.pop()
    recorded_p = list(map(int, recorded_reduced["P_coefficients_ascending"]))
    recorded_q = list(map(int, recorded_reduced["Q_coefficients_ascending"]))
    exact_recorded_model_match = reduced_p == recorded_p and reduced_q == recorded_q
    recomputed_bits = max(integer_bits(value) for value in reduced_p + reduced_q)

    def horizontal(transformation: Any) -> list[list[int]]:
        matrix = transformation[1]
        return [[int(matrix[i, j]) for j in range(2)] for i in range(2)]

    first_horizontal = horizontal(first)
    second_horizontal = horizontal(second)
    composed = primitive_matrix(matrix_product(first_horizontal, second_horizontal))
    recomputed_quartic = discriminant_quartic(reduced_p, reduced_q)
    recorded_quartic = discriminant_quartic(recorded_p, recorded_q)
    horizontal_tie_change = None
    for candidate in SMALL_UNIMODULAR_MATRICES:
        if transform_binary_quartic(recomputed_quartic, candidate) == recorded_quartic:
            horizontal_tie_change = candidate
            break
    if horizontal_tie_change is None:
        raise ArithmeticError(
            "the recorded reduced model is not a small unimodular horizontal tie of the recomputation"
        )
    # candidate maps the recorded reduced coordinate to the recomputed one.
    # Compose it on the right to recover the horizontal map actually used by
    # the sealed ledger, even when hyperellred chose the other LLL tie here.
    composed = primitive_matrix(matrix_product(composed, horizontal_tie_change))
    determinant = composed[0][0] * composed[1][1] - composed[0][1] * composed[1][0]
    if determinant == 0:
        raise ArithmeticError("the composed horizontal model map is singular")

    stage_sizes = []
    for transformation in (first, second):
        matrix = horizontal(transformation)
        correction = pari_polynomial_coefficients(transformation[2])
        entries = [int(transformation[0])] + [entry for row in matrix for entry in row] + correction
        stage_sizes.append(max(integer_bits(entry) for entry in entries))

    reduced_i, reduced_j = binary_quartic_invariants(recorded_quartic)
    if reduced_i.denominator != 1 or reduced_j.denominator != 1:
        raise ArithmeticError("the reduced integral model acquired fractional invariants")
    return {
        "denominator": denominator,
        "integral_coefficients": integral,
        "horizontal_matrix": composed,
        "horizontal_determinant": determinant,
        "horizontal_map_bits": max(integer_bits(entry) for row in composed for entry in row),
        "horizontal_determinant_bits": integer_bits(determinant),
        "stage_map_bits": stage_sizes,
        "reduced_invariants": (int(reduced_i), int(reduced_j)),
        "exact_recorded_model_match": exact_recorded_model_match,
        "recomputed_reduced_bits": recomputed_bits,
        "horizontal_tie_change": horizontal_tie_change,
    }


def reduced_coordinate(parameter: Fraction, matrix: Sequence[Sequence[int]]) -> Fraction | None:
    """Map raw ``t`` to reduced ``s`` when ``t=M(s)``."""

    a_value, b_value = map(int, matrix[0])
    c_value, d_value = map(int, matrix[1])
    numerator = d_value * parameter.numerator - b_value * parameter.denominator
    denominator = -c_value * parameter.numerator + a_value * parameter.denominator
    if denominator == 0:
        return None
    return Q(numerator, denominator)


def reduced_coordinate_data(
    parameter: Fraction, matrix: Sequence[Sequence[int]]
) -> tuple[Fraction | None, int]:
    """Return reduced coordinate and the exact homogeneous preimage gcd."""

    a_value, b_value = map(int, matrix[0])
    c_value, d_value = map(int, matrix[1])
    numerator = d_value * parameter.numerator - b_value * parameter.denominator
    denominator = -c_value * parameter.numerator + a_value * parameter.denominator
    cancellation = math.gcd(abs(numerator), abs(denominator))
    if denominator == 0:
        return None, cancellation
    return Q(numerator, denominator), cancellation


def raw_parameter(curve_a: Fraction, base_point: Point, target: Point) -> Fraction | None:
    if target == base_point:
        return None
    x_value, y_value = target
    x_base, y_base = base_point
    if x_value == x_base:
        if y_value != -y_base or y_base == 0:
            raise ArithmeticError("unexpected exceptional point in the degree-two chart")
        return -(3 * x_base**2 + curve_a) / (2 * y_base)
    return (y_value + y_base) / (x_value - x_base)


def raw_parameter_data(
    curve_a: Fraction, base_point: Point, target: Point
) -> tuple[Fraction | None, int]:
    """Return raw slope and its exact cross-multiplication cancellation gcd."""

    if target == base_point:
        return None, 0
    x_value, y_value = target
    x_base, y_base = base_point
    if x_value == x_base:
        return raw_parameter(curve_a, base_point, target), 1
    y_numerator = (
        y_value.numerator * y_base.denominator
        + y_base.numerator * y_value.denominator
    )
    y_denominator = y_value.denominator * y_base.denominator
    x_numerator = (
        x_value.numerator * x_base.denominator
        - x_base.numerator * x_value.denominator
    )
    x_denominator = x_value.denominator * x_base.denominator
    numerator = y_numerator * x_denominator
    denominator = y_denominator * x_numerator
    cancellation = math.gcd(abs(numerator), abs(denominator))
    return Q(numerator, denominator), cancellation


def quadratic(gram: Sequence[Sequence[float]], vector: Sequence[int]) -> float:
    return sum(
        int(vector[i]) * float(gram[i][j]) * int(vector[j])
        for i in range(len(vector))
        for j in range(len(vector))
    )


def canonical_height_gram(pari: Pari, model: Sequence[Fraction], basis: Sequence[Point]) -> list[list[float]]:
    a_value, b_value = model[3], model[4]
    curve = pari(
        f"ellinit([0,0,0,{rational_string(a_value)},{rational_string(b_value)}])"
    )
    points = pari(
        "["
        + ",".join(
            f"[{rational_string(item[0])},{rational_string(item[1])}]" for item in basis
        )
        + "]"
    )
    raw = curve.ellheightmatrix(points)
    return [[float(raw[i, j]) for j in range(len(basis))] for i in range(len(basis))]


def event_priority(source: str) -> int:
    fields = source.split(":")
    return int(fields[fields.index("priority") + 1])


def detailed_datasets() -> tuple[list[dict[str, Any]], dict[str, Any], dict[str, Any]]:
    rank28 = load_json(RANK28_BLIND)
    rank28_verify = load_json(RANK28_VERIFY)
    curve385 = load_json(CURVE385_ITERATED)
    curve385_primary = load_json(CURVE385_PRIMARY)
    curve398 = load_json(CURVE398_ADAPTIVE)

    rank28_events = {
        int(row["half_lattice_mask"])
        for row in rank28_verify["half_lattice_class_summary"]["productive_centers"]
        if int(row["incremental_quotient_gain"]) > 0
    }
    curve385_iteration = curve385["iterations"][0]
    curve385_events = {
        event_priority(event["sources"][0])
        for event in curve385_iteration["discovered_group_saturation"]["events"]
        if event["type"] == "NEW_Q_INDEPENDENT_DIRECTION"
    }
    curve398_initial_events = {
        event_priority(event["sources"][0])
        for event in curve398["initial_search"]["discovered_group_saturation"]["events"]
        if event["type"] == "NEW_Q_INDEPENDENT_DIRECTION"
    }
    curve398_iteration = curve398["iterations"][0]
    curve398_iteration_events = {
        event_priority(event["sources"][0])
        for event in curve398_iteration["discovered_group_saturation"]["events"]
        if event["type"] == "NEW_Q_INDEPENDENT_DIRECTION"
    }

    datasets = [
        {
            "id": "rank28-selected-union",
            "model": list(map(Q, rank28["fibre"]["short_model"])),
            "basis": [point(item) for item in rank28["fibre"]["generic_points"]],
            "records": rank28["search_records"],
            "event_ids": rank28_events,
            "record_id": lambda row: int(row["mask"]),
            "depth": lambda row: float(row["specialized_depth"]),
            "search": lambda row: row,
        },
        {
            "id": "curve385-first-quotient-iteration",
            "model": list(map(Q, curve385["curve"]["short_model"])),
            "basis": [point(item) for item in curve385_iteration["basis_before"]],
            "records": curve385_iteration["cover_records"],
            "event_ids": curve385_events,
            "record_id": lambda row: int(row["priority"]),
            "depth": lambda row: float(row["canonical_depth"]),
            "search": lambda row: row["search"],
        },
        {
            "id": "curve398-initial-deepest12",
            "model": list(map(Q, curve398["curve"]["short_model"])),
            "basis": [point(item) for item in curve398["curve"]["generic_points"]],
            "records": curve398["initial_search"]["cover_records"],
            "event_ids": curve398_initial_events,
            "record_id": lambda row: int(row["priority"]),
            "depth": lambda row: float(row["specialized_depth"]),
            "search": lambda row: row["search"],
        },
        {
            "id": "curve398-first-quotient-iteration",
            "model": list(map(Q, curve398["curve"]["short_model"])),
            "basis": [point(item) for item in curve398["current_basis"][:21]],
            "records": curve398_iteration["cover_records"],
            "event_ids": curve398_iteration_events,
            "record_id": lambda row: int(row["priority"]),
            "depth": lambda row: float(row["canonical_depth"]),
            "search": lambda row: row["search"],
        },
    ]
    primary_state = curve385_primary["lattice_states"][0]
    for stage in primary_state["stages"]:
        datasets.append(
            {
                "id": f"curve385-primary-{stage['id']}",
                "model": list(map(Q, curve385_primary["curve"]["short_model"])),
                "basis": [point(item) for item in primary_state["basis"]],
                "records": stage["cover_records"],
                "event_ids": set(),
                "record_id": lambda row: int(row["priority"]),
                "depth": lambda row: float(row["canonical_depth"]),
                "search": lambda row: row["search"],
            }
        )
    return datasets, rank28, rank28_verify


def chart_census(pari: Pari, datasets: Sequence[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[tuple[str, int], dict[str, Any]]]:
    all_rows: list[dict[str, Any]] = []
    summaries: list[dict[str, Any]] = []
    runtime: dict[tuple[str, int], dict[str, Any]] = {}
    for dataset in datasets:
        model = dataset["model"]
        basis = dataset["basis"]
        gram = canonical_height_gram(pari, model, basis)
        expected_i = -48 * model[3]
        expected_j = -1728 * model[4]
        reduced_invariants: set[tuple[int, int]] = set()
        rows = []
        maximum_depth_error = 0.0
        for source_record in dataset["records"]:
            search = dataset["search"](source_record)
            identifier = dataset["record_id"](source_record)
            raw_coefficients = list(map(Q, search["raw_quartic_coefficients_ascending"]))
            invariant_i, invariant_j = binary_quartic_invariants(raw_coefficients)
            if invariant_i != expected_i or invariant_j != expected_j:
                raise ArithmeticError(
                    f"{dataset['id']} chart {identifier}: raw quartic invariants changed"
                )
            transformation = pari_transform(pari, raw_coefficients, search["reduced_model"])
            reduced_invariants.add(transformation["reduced_invariants"])
            representative = list(map(int, search["representative"]))
            if len(representative) != len(basis):
                raise ArithmeticError("the representative and basis dimensions disagree")
            depth = float(dataset["depth"](source_record))
            recomputed_depth = quadratic(gram, representative) / 4
            maximum_depth_error = max(maximum_depth_error, abs(depth - recomputed_depth))
            base = point(search["base_point"])

            coordinate_logs = []
            distortions = []
            raw_slope_cancellation_bits = []
            horizontal_cancellation_bits = []
            horizontal_compression = []
            for index, basis_point in enumerate(basis):
                for sign in (1, -1):
                    target = basis_point if sign == 1 else (basis_point[0], -basis_point[1])
                    parameter, raw_cancellation = raw_parameter_data(
                        model[3], base, target
                    )
                    if parameter is None:
                        continue
                    reduced, horizontal_cancellation = reduced_coordinate_data(
                        parameter, transformation["horizontal_matrix"]
                    )
                    if reduced is None:
                        continue
                    log_height = math.log(rational_height(reduced))
                    vector = [-entry for entry in representative]
                    vector[index] += 2 * sign
                    centered_degree_two_height = quadratic(gram, vector) / 2
                    coordinate_logs.append(log_height)
                    distortions.append(log_height - centered_degree_two_height)
                    raw_slope_cancellation_bits.append(integer_bits(raw_cancellation))
                    horizontal_cancellation_bits.append(
                        integer_bits(horizontal_cancellation)
                    )
                    horizontal_compression.append(
                        math.log(rational_height(parameter)) - log_height
                    )
            if not coordinate_logs:
                raise ArithmeticError("a chart sent every signed known basis point to infinity")

            raw_bits = int(search["raw_rational_coefficient_maximum_bits"])
            integral_bits = int(search["integral_model_maximum_coefficient_bits"])
            reduced_bits = int(search["reduced_model"]["maximum_coefficient_bits"])
            denominator_bits = int(search["denominator_clearing_factor_bits"])
            density = float(search["local_stage"]["joint_independent_density_product"])
            base_bits = max(rational_bits(base[0]), rational_bits(base[1]))
            row = {
                "dataset": dataset["id"],
                "chart_id": identifier,
                "chart_hex": search.get("hex"),
                "depth": decimal(depth),
                "first_independent_event": identifier in dataset["event_ids"],
                "finite_point_occurrences": len(search.get("finite_curve_points", [])),
                "has_finite_point": bool(search.get("finite_curve_points", [])),
                "raw_coefficient_bits": raw_bits,
                "denominator_clearing_bits": denominator_bits,
                "integral_coefficient_bits": integral_bits,
                "reduced_coefficient_bits": reduced_bits,
                "coefficient_reduction_bits": integral_bits - reduced_bits,
                "base_point_coordinate_bits": base_bits,
                "modular_density_product": decimal(density),
                "horizontal_reduced_to_raw_matrix": [
                    [str(entry) for entry in line]
                    for line in transformation["horizontal_matrix"]
                ],
                "horizontal_map_bits": transformation["horizontal_map_bits"],
                "horizontal_determinant": str(transformation["horizontal_determinant"]),
                "horizontal_determinant_bits": transformation["horizontal_determinant_bits"],
                "stage_map_bits": transformation["stage_map_bits"],
                "exact_recorded_reduced_model_match": transformation[
                    "exact_recorded_model_match"
                ],
                "recomputed_reduced_coefficient_bits": transformation[
                    "recomputed_reduced_bits"
                ],
                "raw_slope_cancellation_median_bits": median(
                    raw_slope_cancellation_bits
                ),
                "horizontal_cancellation_median_bits": median(
                    horizontal_cancellation_bits
                ),
                "raw_to_reduced_median_log_compression": median(
                    horizontal_compression
                ),
                "known_basis_reduced_log_height": {
                    "minimum": decimal(min(coordinate_logs)),
                    "median": decimal(median(coordinate_logs)),
                    "maximum": decimal(max(coordinate_logs)),
                },
                "known_basis_height_distortion": {
                    "definition": "log H(s(P)) - hhat(2P-Q)/2 for signed displayed basis points",
                    "minimum": decimal(min(distortions)),
                    "median": decimal(median(distortions)),
                    "maximum": decimal(max(distortions)),
                },
                "known_basis_exact_cancellation": {
                    "raw_slope_cross_multiplication_gcd_bits": {
                        "minimum": min(raw_slope_cancellation_bits),
                        "median": median(raw_slope_cancellation_bits),
                        "maximum": max(raw_slope_cancellation_bits),
                    },
                    "horizontal_preimage_gcd_bits": {
                        "minimum": min(horizontal_cancellation_bits),
                        "median": median(horizontal_cancellation_bits),
                        "maximum": max(horizontal_cancellation_bits),
                    },
                    "raw_to_reduced_log_height_compression": {
                        "minimum": decimal(min(horizontal_compression)),
                        "median": decimal(median(horizontal_compression)),
                        "maximum": decimal(max(horizontal_compression)),
                    },
                },
            }
            rows.append(row)
            all_rows.append(row)
            runtime[(dataset["id"], identifier)] = {
                "base_point": base,
                "horizontal_matrix": transformation["horizontal_matrix"],
                "model": model,
            }

        if len(reduced_invariants) != 1:
            raise ArithmeticError(
                f"{dataset['id']}: reduced binary-quartic invariants vary by chart"
            )
        predictor_fields = (
            "depth",
            "raw_coefficient_bits",
            "integral_coefficient_bits",
            "reduced_coefficient_bits",
            "coefficient_reduction_bits",
            "base_point_coordinate_bits",
            "modular_density_product",
            "horizontal_map_bits",
            "horizontal_determinant_bits",
            "raw_slope_cancellation_median_bits",
            "horizontal_cancellation_median_bits",
            "raw_to_reduced_median_log_compression",
        )
        event_predictors = {
            field: auc(rows, field, "first_independent_event") for field in predictor_fields
        }
        finite_point_predictors = {
            field: auc(rows, field, "has_finite_point") for field in predictor_fields
        }
        summaries.append(
            {
                "id": dataset["id"],
                "chart_count": len(rows),
                "first_independent_event_chart_count": sum(
                    row["first_independent_event"] for row in rows
                ),
                "finite_point_chart_count": sum(row["has_finite_point"] for row in rows),
                "exact_recorded_reduced_model_match_count": sum(
                    row["exact_recorded_reduced_model_match"] for row in rows
                ),
                "depth_quantiles": quantiles([float(row["depth"]) for row in rows]),
                "maximum_recomputed_depth_error": decimal(maximum_depth_error),
                "raw_binary_quartic_invariants": {
                    "I": rational_string(expected_i),
                    "J": rational_string(expected_j),
                    "formula": "I=-48*A; J=-1728*B",
                    "constant_across_every_chart": True,
                },
                "reduced_discriminant_quartic_invariants": {
                    "I": str(next(iter(reduced_invariants))[0]),
                    "J": str(next(iter(reduced_invariants))[1]),
                    "constant_across_every_chart": True,
                },
                "presearch_predictors_of_first_independent_event": event_predictors,
                "presearch_predictors_of_any_finite_point": finite_point_predictors,
            }
        )
    return all_rows, summaries, runtime


def prefix_gain_labels(
    ordered_masks: Sequence[int], quotient_masks_by_chart: dict[int, list[int]]
) -> tuple[list[int], int]:
    """Label the charts that first enlarge the accumulated F_2 quotient span."""

    accumulated: list[int] = []
    previous_rank = 0
    gains = []
    for mask in ordered_masks:
        accumulated.extend(quotient_masks_by_chart.get(int(mask), []))
        current_rank = binary_rank(accumulated)
        gains.append(current_rank - previous_rank)
        previous_rank = current_rank
    return gains, previous_rank


def compact_control_predictor_panel() -> dict[str, Any]:
    """Add the earlier positive controls whose ledgers store compact chart metrics."""

    r17_blind = load_json(R17_MATRIX_BLIND)
    r17_verify = load_json(R17_MATRIX_VERIFY)
    rank29_blind = load_json(RANK29_CONTROLS_BLIND)
    rank29_verify = load_json(RANK29_CONTROLS_VERIFY)
    predictor_fields = (
        "depth",
        "integral_coefficient_bits",
        "reduced_coefficient_bits",
        "modular_density_product",
    )
    cases = []

    r17_fibres = {item["id"]: item for item in r17_blind["fibres"]}
    for verification in r17_verify["positive_controls"]:
        identifier = verification["id"]
        fibre = r17_fibres[identifier]
        quotient_masks_by_chart = {
            int(item["half_lattice_mask"]): list(
                map(int, item["recovered_quotient_masks"])
            )
            for item in verification["productive_half_lattice_centers"]
        }
        ordered_masks = [int(item["mask"]) for item in fibre["cover_records"]]
        gains, final_rank = prefix_gain_labels(ordered_masks, quotient_masks_by_chart)
        expected_rank = int(verification["blind_recovered_quotient_dimension"])
        if final_rank != expected_rank:
            raise ArithmeticError(f"{identifier}: compact prefix rank changed")
        rows = [
            {
                "chart_mask": int(record["mask"]),
                "prefix_quotient_gain": gain,
                "first_independent_event": gain > 0,
                "depth": float(record["specialized_depth"]),
                "integral_coefficient_bits": int(
                    record["integral_model_maximum_coefficient_bits"]
                ),
                "reduced_coefficient_bits": int(
                    record["reduced_model_maximum_coefficient_bits"]
                ),
                "modular_density_product": float(
                    record["independent_modular_density_product"]
                ),
            }
            for record, gain in zip(fibre["cover_records"], gains)
        ]
        cases.append(
            {
                "id": identifier,
                "lineage": "published-MW17-positive-control",
                "displayed_exceptional_quotient_dimension": int(
                    verification["public_exceptional_quotient_dimension"]
                ),
                "blind_recovered_quotient_dimension": final_rank,
                "chart_count": len(rows),
                "first_independent_event_chart_count": sum(
                    row["first_independent_event"] for row in rows
                ),
                "presearch_predictors_of_first_independent_event": {
                    field: auc(rows, field, "first_independent_event")
                    for field in predictor_fields
                },
                "rows": rows,
            }
        )

    rank29_results = {item["label"]: item for item in rank29_blind["results"]}
    for verification in rank29_verify["results"]:
        identifier = verification["label"]
        result = rank29_results[identifier]
        quotient_masks_by_chart: dict[int, list[int]] = {}
        for relation in verification["relations"]:
            for source_mask in relation["source_half_class_masks"]:
                quotient_masks_by_chart.setdefault(int(source_mask), []).append(
                    int(relation["quotient_mask"])
                )
        ordered_masks = [int(item["mask"]) for item in result["cover_records"]]
        gains, final_rank = prefix_gain_labels(ordered_masks, quotient_masks_by_chart)
        expected_rank = int(verification["blind_selected_exact_quotient_rank"])
        if final_rank != expected_rank:
            raise ArithmeticError(f"{identifier}: compact prefix rank changed")
        rows = [
            {
                "chart_mask": int(record["mask"]),
                "prefix_quotient_gain": gain,
                "first_independent_event": gain > 0,
                "depth": float(record["specialized_depth"]),
                "integral_coefficient_bits": int(record["integral_coefficient_bits"]),
                "reduced_coefficient_bits": int(record["reduced_coefficient_bits"]),
                "modular_density_product": float(record["modular_density_product"]),
            }
            for record, gain in zip(result["cover_records"], gains)
        ]
        cases.append(
            {
                "id": identifier,
                "lineage": result["lineage"],
                "displayed_exceptional_quotient_dimension": int(
                    verification["displayed_exceptional_quotient_rank"]
                ),
                "blind_recovered_quotient_dimension": final_rank,
                "chart_count": len(rows),
                "first_independent_event_chart_count": sum(
                    row["first_independent_event"] for row in rows
                ),
                "presearch_predictors_of_first_independent_event": {
                    field: auc(rows, field, "first_independent_event")
                    for field in predictor_fields
                },
                "rows": rows,
            }
        )

    return {
        "boundary": (
            "These sealed compact ledgers expose depth, integral/reduced coefficient "
            "height, and modular density, but not the full horizontal maps reconstructed "
            "in the detailed chart census. Outcome labels are added from their separate "
            "verification ledgers."
        ),
        "case_count": len(cases),
        "chart_count": sum(item["chart_count"] for item in cases),
        "cases": cases,
    }


def rank28_posthoc_oracle(
    pari: Pari,
    rank28: dict[str, Any],
    verification: dict[str, Any],
    rows: Sequence[dict[str, Any]],
    runtime: dict[tuple[str, int], dict[str, Any]],
) -> dict[str, Any]:
    dataset = "rank28-selected-union"
    chart_rows = [row for row in rows if row["dataset"] == dataset]
    row_by_mask = {int(row["chart_id"]): row for row in chart_rows}
    targets = [
        {
            "point": point(item["point"]),
            "quotient_mask": int(item["exceptional_quotient_mod2_mask"]),
            "source_masks": set(map(int, item["source_half_lattice_masks"])),
        }
        for item in verification["labeled_blind_points"]
    ]
    productive = {
        int(item["half_lattice_mask"])
        for item in verification["half_lattice_class_summary"]["productive_centers"]
    }
    incremental = {
        int(item["half_lattice_mask"]): int(item["incremental_quotient_gain"])
        for item in verification["half_lattice_class_summary"]["productive_centers"]
    }
    model = list(map(Q, rank28["fibre"]["short_model"]))
    curve = pari(
        f"ellinit([0,0,0,{rational_string(model[3])},{rational_string(model[4])}])"
    )
    height_cache: dict[Point, float] = {}

    def canonical_height(item: Point) -> float:
        if item not in height_cache:
            pari_point = pari(
                f"[{rational_string(item[0])},{rational_string(item[1])}]"
            )
            height_cache[item] = float(curve.ellheight(pari_point))
        return height_cache[item]

    confusion = {"true_positive": 0, "false_positive": 0, "false_negative": 0, "true_negative": 0}
    chart_oracle_rows = []
    prefix_masks: list[int] = []
    predicted_rank = 0
    prefix_matches = True
    height_bound = int(rank28["declared_search_budget"]["height_bound_each"])
    curve_coefficients = (Q(0), Q(0), Q(0), model[3], model[4])
    for source_record in rank28["search_records"]:
        mask = int(source_record["mask"])
        state = runtime[(dataset, mask)]
        base = state["base_point"]
        matrix = state["horizontal_matrix"]
        base_height = canonical_height(base)
        minimum_centered = math.inf
        minimum_energy = math.inf
        visible_masks = []
        visible_count = 0
        for target in targets:
            target_point = target["point"]
            parameter = raw_parameter(model[3], base, target_point)
            reduced = None if parameter is None else reduced_coordinate(parameter, matrix)
            visible = reduced is None or rational_height(reduced) <= height_bound
            recorded_source = mask in target["source_masks"]
            if visible and recorded_source:
                confusion["true_positive"] += 1
            elif visible:
                confusion["false_positive"] += 1
            elif recorded_source:
                confusion["false_negative"] += 1
            else:
                confusion["true_negative"] += 1
            if visible:
                visible_count += 1
                visible_masks.append(target["quotient_mask"])

            doubled = short_add(curve_coefficients, target_point, target_point)
            if doubled is None:
                raise ArithmeticError("a non-torsion target doubled to infinity")
            delta = short_add(curve_coefficients, doubled, (base[0], -base[1]))
            if delta is None:
                centered = 0.0
            else:
                centered = canonical_height(delta) / 4
            energy = base_height / 2 + 2 * centered
            minimum_centered = min(minimum_centered, centered)
            minimum_energy = min(minimum_energy, energy)

        old_rank = predicted_rank
        prefix_masks.extend(visible_masks)
        predicted_rank = binary_rank(prefix_masks)
        predicted_gain = predicted_rank - old_rank
        recorded_gain = incremental.get(mask, 0)
        prefix_matches &= predicted_gain == recorded_gain
        row_by_mask[mask]["posthoc_rank28_target_oracle"] = {
            "minimum_centered_half_height": decimal(minimum_centered),
            "minimum_fiber_energy": decimal(minimum_energy),
            "visible_candidate_count_at_declared_bound": visible_count,
            "predicted_incremental_quotient_gain": predicted_gain,
        }
        chart_oracle_rows.append(
            {
                "productive": mask in productive,
                "minimum_centered_half_height": minimum_centered,
                "minimum_fiber_energy": minimum_energy,
            }
        )

    productive_minima = [
        row["minimum_centered_half_height"] for row in chart_oracle_rows if row["productive"]
    ]
    nonproductive_minima = [
        row["minimum_centered_half_height"] for row in chart_oracle_rows if not row["productive"]
    ]
    centered_auc = auc(chart_oracle_rows, "minimum_centered_half_height", "productive")
    energy_auc = auc(chart_oracle_rows, "minimum_fiber_energy", "productive")
    exact_visibility = confusion["false_positive"] == 0 and confusion["false_negative"] == 0
    if not exact_visibility:
        raise ArithmeticError("reduced-coordinate visibility no longer matches the sealed source ledger")
    if not prefix_matches or predicted_rank != 11:
        raise ArithmeticError("coordinate visibility no longer reconstructs the quotient-gain trace")
    return {
        "boundary": "posthoc: target coordinates and exact quotient labels are loaded only here",
        "target_count": len(targets),
        "chart_count": len(chart_oracle_rows),
        "declared_reduced_coordinate_height_bound": height_bound,
        "pairwise_visibility_confusion": confusion,
        "reduced_coordinate_visibility_matches_recorded_sources_exactly": exact_visibility,
        "coordinate_visibility_reconstructs_every_prefix_quotient_gain": prefix_matches,
        "reconstructed_final_quotient_rank": predicted_rank,
        "minimum_centered_half_height_predictor": centered_auc,
        "minimum_fiber_energy_predictor": energy_auc,
        "centered_height_strict_separation": {
            "maximum_over_productive_charts": decimal(max(productive_minima)),
            "minimum_over_nonproductive_charts": decimal(min(nonproductive_minima)),
            "gap_nonproductive_minus_productive": decimal(
                min(nonproductive_minima) - max(productive_minima)
            ),
        },
        "interpretation": (
            "On this sealed selected-chart ledger, a chart is productive exactly when "
            "one of the forty returned non-generic targets is sufficiently close to its "
            "half-lattice midpoint. This is a retrospective mechanism audit, not a "
            "prospective target-free classifier."
        ),
    }


def build_payload() -> dict[str, Any]:
    pari = Pari()
    pari.allocatemem(500_000_000)
    # Match the standalone GP default used by the sealed search runners.
    # hyperellred can choose a height-equivalent sign/translation at an LLL tie;
    # those harmless alternate representatives are explicitly counted below.
    pari.default("realprecision", 38)
    datasets, rank28, rank28_verify = detailed_datasets()
    rows, summaries, runtime = chart_census(pari, datasets)
    posthoc = rank28_posthoc_oracle(pari, rank28, rank28_verify, rows, runtime)
    compact_panel = compact_control_predictor_panel()

    event_datasets = [
        summary for summary in summaries if summary["first_independent_event_chart_count"] > 0
    ]
    predictor_sets = [
        summary["presearch_predictors_of_first_independent_event"]
        for summary in event_datasets
    ] + [
        case["presearch_predictors_of_first_independent_event"]
        for case in compact_panel["cases"]
    ]
    stable_predictors = []
    predictor_names = set.intersection(*(set(item) for item in predictor_sets))
    for name in sorted(predictor_names):
        available = [item[name] for item in predictor_sets if item[name] is not None]
        if len(available) != len(predictor_sets):
            continue
        orientations = {item["better_orientation"] for item in available}
        worst_auc = min(float(item["best_orientation_auc"]) for item in available)
        if len(orientations) == 1 and worst_auc >= 0.7:
            stable_predictors.append(name)

    return {
        "schema": "elliptic-curves.half-lattice-height-compression-analysis.v1",
        "status": "PASS_EXACT_IDENTITIES_AND_BOUNDED_RETROSPECTIVE_MECHANISM_AUDIT",
        "input_hashes": {
            str(path.relative_to(ROOT)): digest(path)
            for path in (
                RANK28_BLIND,
                RANK28_VERIFY,
                CURVE385_ITERATED,
                CURVE385_PRIMARY,
                CURVE398_ADAPTIVE,
                R17_MATRIX_BLIND,
                R17_MATRIX_VERIFY,
                RANK29_CONTROLS_BLIND,
                RANK29_CONTROLS_VERIFY,
                Path(__file__).resolve(),
            )
        },
        "theorem": {
            "chart": "C_Q: w^2=t^4-6*x_Q*t^2-8*y_Q*t-3*x_Q^2-4*A",
            "fiber_involution": "R -> Q-R",
            "exact_coordinate_identity": "t_Q(R)^2=x(R)+x(Q-R)+x(Q)",
            "exact_canonical_height_identity": (
                "hhat(R)+hhat(Q-R)=(hhat(Q)+hhat(2R-Q))/2"
            ),
            "centered_degree_two_height": (
                "2*hhat(R-Q/2)=hhat(2R-Q)/2"
            ),
            "half_lattice_empty_ball_identity": (
                "if Q_c is shortest in c modulo 2M, then "
                "min_{P in M} hhat(P-Q_c/2)=hhat(Q_c)/4"
            ),
            "invariants": "I=-48*A and J=-1728*B for every Q",
            "pgl2_reduction_effect": (
                "a reduced horizontal coordinate changes the Weil height by a "
                "chart-dependent bounded function but preserves the midpoint Q/2"
            ),
        },
        "chart_census": {
            "dataset_count": len(summaries),
            "chart_count": len(rows),
            "datasets": summaries,
            "rows": rows,
        },
        "compact_control_predictor_panel": compact_panel,
        "predictor_conclusion": {
            "positive_chart_order_count": len(predictor_sets),
            "common_target_free_scalar_predictors": sorted(predictor_names),
            "stable_target_free_scalar_predictors_at_auc_at_least_0_7": stable_predictors,
            "conclusion": (
                "No recorded target-free scalar is direction-stable with worst-case "
                "AUC at least 0.7 across the eleven positive chart orders. Depth is a "
                "coarse old-point exclusion selector, not a within-stratum discovery "
                "probability. Target-relative centered height is the mechanism visible "
                "in the sealed rank-28 outcome."
            ),
        },
        "rank28_posthoc_target_oracle": posthoc,
        "scheduler_consequence": {
            "primary_score": (
                "maximize a certified or enumerated old-point exclusion margin in the "
                "reduced coordinate, using actual current-lattice Voronoi/deep-hole "
                "centers and a chart-specific height-distortion bound"
            ),
            "secondary_rule": (
                "choose geometrically diverse centers; quotient Hamming weight is not "
                "a substitute for covering-radius coverage"
            ),
            "restart_rule": (
                "after every rank or finite-index enlargement, recompute the current "
                "lattice, its centers, reduced maps, and distortion margins"
            ),
            "curve385_primary_interpretation": (
                "the completed weight-one/two miss says that 3,116 additional coordinate "
                "charts exposed only the already known M29; it does not contradict the "
                "midpoint identity and supplies no rank upper bound"
            ),
        },
        "claim_boundary": [
            "The quartic formulas, invariant identities, rational maps, and coordinate visibility checks are exact.",
            "Canonical heights and AUC summaries are high-precision numerical diagnostics, not exact inequalities.",
            "The strict rank-28 separation is retrospective on a fixed sealed selected-chart ledger.",
            "No target-free success-probability model is claimed.",
            "No bounded miss implies point absence, saturation, a Selmer statement, or a rank upper bound.",
        ],
        "reproducing_command": (
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
            "elliptic-curves/cas/analyze_half_lattice_height_compression.sage"
        ),
    }


def deterministic_gzip(payload: dict[str, Any]) -> bytes:
    encoded = (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode()
    import io

    buffer = io.BytesIO()
    with gzip.GzipFile(filename="", mode="wb", fileobj=buffer, mtime=0) as handle:
        handle.write(encoded)
    return buffer.getvalue()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    payload = build_payload()
    if args.check:
        if not args.output.is_file():
            raise SystemExit(f"missing artifact: {args.output}")
        existing = load_json(args.output)
        if existing != payload:
            raise SystemExit("half-lattice height-compression artifact is stale")
        print(
            "HALFLATTICEHEIGHT|status=PASS|"
            f"charts={payload['chart_census']['chart_count']}|"
            "visibility=exact|prefix_gain=exact"
        )
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(deterministic_gzip(payload))
    print(
        "HALFLATTICEHEIGHT|status=WROTE|"
        f"charts={payload['chart_census']['chart_count']}|"
        f"output={args.output.relative_to(ROOT)}"
    )


if __name__ == "__main__":
    main()
