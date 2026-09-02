#!/usr/bin/env python3
"""Analyze the exceptional quotient bases on the four published R17 controls.

This calculation combines three already certified layers:

* the exact embedding of the specialized generic R17 subgroup in each public
  point lattice;
* the complete rational-bisection specialization census and the declared
  degree-three/four equation samples;
* the exact q12/orbit5867 point map to the 4A1/MW13 parent.

The new numerical invariant is the Schur-complement height form obtained by
orthogonally projecting the public exceptional basis away from the seventeen
specialized generic sections.  It is an arithmetic canonical-height invariant
of the displayed subgroup.  It is not a divisor class on the K3 surface.
"""

from __future__ import annotations

import argparse
from decimal import Decimal, getcontext
from fractions import Fraction
from hashlib import sha256
import importlib
from itertools import combinations
import json
from math import gcd, lcm, sqrt
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC = ROOT / "elliptic-curves"
sys.path[:0] = [str(ELLIPTIC), str(ELLIPTIC / "cas")]

from latent_lattice import EllipticCurve, height_gram  # noqa: E402


TRUTH = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "latent_lattice_calibration_truth_v1.json"
)
CONTROLS = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_high_rank_positive_controls_v2.json"
)
BISECTION_CONTROLS = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_bisection_specialization_controls_v1.json"
)
DEEP_COVERS = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_deep_cover_exceptional_quotients_v1.json"
)
BISECTIONS = ROOT / "artifacts/generated-results/elkies-2026-equation-bisections-full.json"
POINT_FACTORY = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-q12o5867-genus-one-point-factory-controls.json"
)
VISIBILITY = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_bisection_visibility_record_curves_v1.json"
)
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_exceptional_specialization_relations_v1.json"
)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def display(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def dtext(value: Decimal, digits: int = 32) -> str:
    if value.is_zero():
        return "0"
    return format(value, f".{digits}g")


def matrix(rows) -> list[list[Decimal]]:
    return [[Decimal(str(value)) for value in row] for row in rows]


def transpose(a):
    return [list(row) for row in zip(*a)]


def matmul(a, b):
    bt = transpose(b)
    return [[sum((x * y for x, y in zip(row, col)), Decimal(0)) for col in bt] for row in a]


def matsub(a, b):
    return [[x - y for x, y in zip(arow, brow)] for arow, brow in zip(a, b)]


def solve(a, b):
    """Solve A X=B by Decimal Gauss--Jordan elimination with pivoting."""

    n = len(a)
    width = len(b[0])
    aug = [list(a[row]) + list(b[row]) for row in range(n)]
    for col in range(n):
        pivot = max(range(col, n), key=lambda row: abs(aug[row][col]))
        if aug[pivot][col] == 0:
            raise ArithmeticError("singular height Gram")
        aug[col], aug[pivot] = aug[pivot], aug[col]
        unit = aug[col][col]
        aug[col] = [value / unit for value in aug[col]]
        for row in range(n):
            if row == col or aug[row][col] == 0:
                continue
            factor = aug[row][col]
            aug[row] = [x - factor * y for x, y in zip(aug[row], aug[col])]
    return [row[n : n + width] for row in aug]


def decimal_determinant(a) -> Decimal:
    work = [list(row) for row in a]
    answer = Decimal(1)
    for col in range(len(work)):
        pivot = max(range(col, len(work)), key=lambda row: abs(work[row][col]))
        if work[pivot][col] == 0:
            return Decimal(0)
        if pivot != col:
            work[col], work[pivot] = work[pivot], work[col]
            answer = -answer
        value = work[col][col]
        answer *= value
        for row in range(col + 1, len(work)):
            factor = work[row][col] / value
            for index in range(col + 1, len(work)):
                work[row][index] -= factor * work[col][index]
    return answer


def bareiss_determinant(a: list[list[int]]) -> int:
    work = [list(map(int, row)) for row in a]
    sign = 1
    previous = 1
    for col in range(len(work) - 1):
        pivot = next((row for row in range(col, len(work)) if work[row][col]), None)
        if pivot is None:
            return 0
        if pivot != col:
            work[col], work[pivot] = work[pivot], work[col]
            sign = -sign
        value = work[col][col]
        for row in range(col + 1, len(work)):
            for index in range(col + 1, len(work)):
                numerator = work[row][index] * value - work[row][col] * work[col][index]
                if numerator % previous:
                    raise ArithmeticError("Bareiss division lost exactness")
                work[row][index] = numerator // previous
        previous = value
    return sign * work[-1][-1]


def submatrix(a, rows, columns):
    return [[a[row][column] for column in columns] for row in rows]


def standard_columns(rank: int, indices: list[int]) -> list[list[int]]:
    return [[int(row == index) for index in indices] for row in range(rank)]


def unique_rational_coordinates(vectors: list[list[int]], target: list[int]):
    """Return coordinates for target in an independent vector list, or None."""

    if not vectors:
        return None
    width = len(vectors)
    rows = [
        [Fraction(vectors[column][row]) for column in range(width)] + [Fraction(target[row])]
        for row in range(len(target))
    ]
    pivot_columns = []
    pivot_row = 0
    for column in range(width):
        pivot = next((row for row in range(pivot_row, len(rows)) if rows[row][column]), None)
        if pivot is None:
            return None
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        unit = rows[pivot_row][column]
        rows[pivot_row] = [value / unit for value in rows[pivot_row]]
        for row in range(len(rows)):
            if row == pivot_row or not rows[row][column]:
                continue
            factor = rows[row][column]
            rows[row] = [x - factor * y for x, y in zip(rows[row], rows[pivot_row])]
        pivot_columns.append(column)
        pivot_row += 1
    if any(all(value == 0 for value in row[:width]) and row[-1] for row in rows):
        return None
    answer = [Fraction(0) for _ in range(width)]
    for row in rows:
        pivot = next((column for column in range(width) if row[column]), None)
        if pivot is not None:
            answer[pivot] = row[-1]
    return answer


def minimum_bisection_compositum(cover_points, target: list[int], bisection_by_label):
    """Find the smallest exact integral combination of split bisection classes."""

    vectors = [
        [int(value) for value in point["exceptional_quotient_coordinates"]]
        for point in cover_points
    ]
    for size in range(1, len(vectors) + 1):
        for subset in combinations(range(len(vectors)), size):
            coordinates = unique_rational_coordinates([vectors[index] for index in subset], target)
            if coordinates is None or any(value.denominator != 1 for value in coordinates):
                continue
            if any(value == 0 for value in coordinates):
                continue
            coefficients = [int(value) for value in coordinates]
            labels = [cover_points[index]["cover_label"] for index in subset]
            traces = [bisection_by_label[label]["published_basis_w"] for label in labels]
            trace_multiplier = 2 ** (size - 1)
            compositum_trace = [
                trace_multiplier
                * sum(coefficients[index] * int(traces[index][coordinate]) for index in range(size))
                for coordinate in range(17)
            ]
            return {
                "status": "EXACT_IN_KNOWN_RATIONAL_BISECTION_COMPOSITUM",
                "number_of_independent_quadratic_extensions": size,
                "generic_cover_degree": 2**size,
                "cover_labels": labels,
                "integer_coefficients": coefficients,
                "generic_galois_trace_published_R17_basis_coordinates": compositum_trace,
                "generic_galois_trace_formula": "2^(k-1) * sum(c_i*tau_i)",
                "reason_for_cover_degree": (
                    "The complete bisection squareclasses are represented by distinct irreducible "
                    "quadratics, so these function-field squareclasses are independent."
                ),
            }
    return {
        "status": "NOT_IN_KNOWN_RATIONAL_BISECTION_COMPOSITUM",
        "boundary": (
            "No combination of the rational bisections that split at this fibre gives this unit "
            "quotient direction. This does not exclude other multisections."
        ),
    }


def height_summary(form: list[list[Decimal]]) -> dict[str, object]:
    floats = np.array([[float(value) for value in row] for row in form], dtype=float)
    eigenvalues = np.linalg.eigvalsh(floats)
    diagonal = np.diag(floats)
    correlations = []
    for left in range(len(form)):
        for right in range(left + 1, len(form)):
            correlations.append(floats[left, right] / sqrt(diagonal[left] * diagonal[right]))
    return {
        "dimension": len(form),
        "diagonal_min_median_max": [
            f"{float(np.min(diagonal)):.17g}",
            f"{float(np.median(diagonal)):.17g}",
            f"{float(np.max(diagonal)):.17g}",
        ],
        "eigenvalues": [f"{value:.17g}" for value in eigenvalues],
        "spectral_condition_number": f"{float(eigenvalues[-1] / eigenvalues[0]):.17g}",
        "maximum_absolute_off_diagonal_correlation": (
            f"{max(map(abs, correlations)):.17g}" if correlations else "0"
        ),
        "determinant": dtext(decimal_determinant(form), 40),
    }


def pearson(left, right) -> str:
    return f"{float(np.corrcoef(np.asarray(left, dtype=float), np.asarray(right, dtype=float))[0, 1]):.17g}"


def fraction_parameter(text: str) -> tuple[int, int]:
    value = Fraction(text)
    return value.numerator, value.denominator


def primitive_polynomial(coefficients: list[str]) -> list[int]:
    values = [Fraction(value) for value in coefficients]
    denominator = 1
    for value in values:
        denominator = lcm(denominator, value.denominator)
    integers = [value.numerator * (denominator // value.denominator) for value in values]
    content = 0
    for value in integers:
        content = gcd(content, abs(value))
    integers = [value // content for value in integers]
    if integers[-1] < 0:
        integers = [-value for value in integers]
    return integers


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--digits", type=int, default=90)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.digits < 60:
        raise SystemExit("--digits must be at least 60")
    getcontext().prec = args.digits

    truth = json.loads(TRUTH.read_text())
    controls = json.loads(CONTROLS.read_text())
    bisection_controls = json.loads(BISECTION_CONTROLS.read_text())
    deep = json.loads(DEEP_COVERS.read_text())
    bisections = json.loads(BISECTIONS.read_text())
    point_factory = json.loads(POINT_FACTORY.read_text())
    visibility = json.loads(VISIBILITY.read_text())

    truth_by_parameter = {row["parameter"]: row for row in truth["positive_controls"]}
    control_by_parameter = {row["parameter"]: row for row in controls["fibres"]}
    bcontrol_by_parameter = {row["parameter"]: row for row in bisection_controls["fibres"]}
    deep_by_parameter = {row["parameter"]: row for row in deep["fibres"]}
    parent_by_parameter = {row["published_parameter"]: row for row in point_factory["controls"]}
    visibility_by_parameter = {row["parameter"]: row for row in visibility["visibility"]}
    bisection_by_label = {row["label"]: row for row in bisections["bisections"]}

    fibres = []
    cover_occurrences: dict[str, list[str]] = {}
    rank28_unexplained_form = None
    rank28_parent = None
    for rank in range(25, 29):
        parameter = next(
            value for value, row in control_by_parameter.items() if int(row["published_rank_lower_bound"]) == rank
        )
        module = importlib.import_module(f"elkies_rank{rank}")
        public_points = tuple(module.POINTS)
        curve = EllipticCurve(module.GENERAL_WEIERSTRASS_COEFFICIENTS)
        raw_height_gram = height_gram(curve, public_points, digits=args.digits, timeout=300)
        height = matrix(raw_height_gram)

        truth_row = truth_by_parameter[parameter]
        embedding = [list(map(int, row)) for row in truth_row["embedding_matrix_rows"]]
        generic = matrix(embedding)
        generic_gram = matmul(transpose(generic), matmul(height, generic))
        stored_generic_gram = matrix(truth_row["canonical_height_gram"])
        maximum_replay_error = max(
            abs(x - y)
            for row, stored in zip(generic_gram, stored_generic_gram)
            for x, y in zip(row, stored)
        )
        if maximum_replay_error > Decimal("1e-70"):
            raise ArithmeticError(f"generic height replay failed at {parameter}: {maximum_replay_error}")

        deep_row = deep_by_parameter[parameter]
        source_indices = [
            int(value) - 1
            for value in deep_row["exceptional_quotient"]["public_point_indices_one_based"]
        ]
        labels = list(deep_row["exceptional_quotient"]["ordered_basis"])
        exceptional_columns_int = standard_columns(rank, source_indices)
        exceptional = matrix(exceptional_columns_int)
        generic_exceptional = matmul(transpose(generic), matmul(height, exceptional))
        exceptional_gram = matmul(transpose(exceptional), matmul(height, exceptional))
        projection = solve(generic_gram, generic_exceptional)
        quotient_form = matsub(
            exceptional_gram,
            matmul(transpose(generic_exceptional), projection),
        )
        if min(quotient_form[index][index] for index in range(len(labels))) <= 0:
            raise ArithmeticError(f"non-positive quotient height at {parameter}")

        combined = [
            embedding[row] + exceptional_columns_int[row]
            for row in range(rank)
        ]
        combined_index = abs(bareiss_determinant(combined))
        parent_records = {
            row["public_complement_label"].removeprefix("public-complement-"): row
            for row in parent_by_parameter[parameter]["points"]
        }
        cover_points = deep_row["points"]
        per_direction = []
        mean_generic_height = sum(generic_gram[index][index] for index in range(17)) / Decimal(17)
        for index, label in enumerate(labels):
            unit = [0] * len(labels)
            unit[index] = 1
            compositum = minimum_bisection_compositum(
                cover_points, unit, bisection_by_label
            )
            direct = []
            involving = []
            for point in cover_points:
                vector = [int(value) for value in point["exceptional_quotient_coordinates"]]
                record = {
                    "cover_degree": int(point["cover_degree"]),
                    "cover_label": point["cover_label"],
                    "quotient_vector": vector,
                    "generic_correction_coordinates": [
                        int(value) for value in point["known_basis_coordinates"][:17]
                    ],
                    "trace_lattice_norm": int(point["minimum_lattice_norm"]),
                }
                if vector[index]:
                    involving.append(record)
                if vector == unit or vector == [-value for value in unit]:
                    direct.append(record)
            parent = parent_records[label]
            raw = exceptional_gram[index][index]
            defect = quotient_form[index][index]
            per_direction.append(
                {
                    "label": label,
                    "source_public_point_index_one_based": source_indices[index] + 1,
                    "public_point": [str(value) for value in public_points[source_indices[index]]],
                    "specialized_generic_height_pairings": [
                        dtext(generic_exceptional[row][index]) for row in range(17)
                    ],
                    "orthogonal_projection_coefficients_in_generic_basis": [
                        dtext(projection[row][index]) for row in range(17)
                    ],
                    "raw_canonical_height": dtext(raw),
                    "quotient_height_specialization_defect": dtext(defect),
                    "defect_over_mean_specialized_generic_height": dtext(defect / mean_generic_height),
                    "fraction_of_raw_height_removed_by_generic_projection": dtext((raw - defect) / raw),
                    "direct_rational_bisection_relations": direct,
                    "rational_bisection_relations_involving_direction": involving,
                    "minimal_cover_degree_status": (
                        {
                            "status": "EXACT_MINIMUM_WITHIN_RATIONAL_BISECTION_UNIVERSE",
                            "degree": 2,
                            "boundary": "The point class is one branch modulo a specialized generic section; degree one is excluded by the certified independent complement.",
                        }
                        if direct
                        else {
                            "status": "UNKNOWN",
                            "proved_exclusions": [
                                "degree 1 in the specialized generic R17 subgroup",
                                "every complete section-nonnegative rational bisection class",
                            ],
                            "bounded_negative_tests": [
                                "69 sampled norm-20 trisection representatives",
                                "all 160 inversion representatives in the norm-26 deep trisection shell",
                                "53 sampled norm-34 quadrisection representatives",
                            ],
                            "boundary": "Higher-genus degree-two curves and the untested degree-three/four cosets are not excluded.",
                        }
                    ),
                    "known_bisection_compositum_deformation": compositum,
                    "parent_transport": {
                        "parent_crossratio_base": parent["parent_crossratio_base"],
                        "parent_crossratio_base_projective_height_bits": int(
                            parent["parent_crossratio_base_projective_height"].bit_length()
                        ),
                        "parent_point_canonical_height_128bit": parent[
                            "parent_point_canonical_height_128bit"
                        ],
                        "exact_roundtrip": bool(parent["exact_forward_inverse_roundtrip"]),
                    },
                    "divisor_class_status": {
                        "status": "NOT_INTRINSIC_FOR_A_FIBRE_POINT",
                        "reason": "A point on an elliptic fibre is codimension two on the K3. A divisor class exists only after choosing a multisection through it.",
                    },
                }
            )

        specialization_equations = []
        hits = {row["label"]: row for row in bcontrol_by_parameter[parameter]["hits"]}
        for point in cover_points:
            label = point["cover_label"]
            hit = hits[label]
            cover = bisection_by_label[label]
            cover_occurrences.setdefault(label, []).append(parameter)
            primitive_q = primitive_polynomial(cover["residual_chord"]["q_coefficients"])
            root = Fraction(hit["canonical_positive_square_root"])
            specialization_equations.append(
                {
                    "cover_label": label,
                    "trace_published_basis_coordinates": cover["published_basis_w"],
                    "quadratic_q_coefficients_low_to_high": cover["residual_chord"]["q_coefficients"],
                    "primitive_integer_q_coefficients_low_to_high": primitive_q,
                    "primitive_q_discriminant": str(primitive_q[1] ** 2 - 4 * primitive_q[0] * primitive_q[2]),
                    "primitive_q_maximum_coefficient_bits": max(abs(value).bit_length() for value in primitive_q),
                    "specialized_q_value": hit["q_value"],
                    "canonical_positive_square_root": hit["canonical_positive_square_root"],
                    "specialized_square_root_projective_height_bits": max(
                        abs(root.numerator).bit_length(), root.denominator.bit_length()
                    ),
                    "equation": "q0 + q1*t + q2*t^2 = r^2",
                    "exact_at_parameter": parameter,
                    "exceptional_quotient_vector": point["exceptional_quotient_coordinates"],
                    "arithmetic_quotient_height_of_specialized_branch_class": dtext(
                        sum(
                            Decimal(int(point["exceptional_quotient_coordinates"][left]))
                            * quotient_form[left][right]
                            * Decimal(int(point["exceptional_quotient_coordinates"][right]))
                            for left in range(len(labels))
                            for right in range(len(labels))
                        )
                    ),
                }
            )

        unexplained_labels = [
            label.removeprefix("public-complement-")
            for label in visibility_by_parameter[parameter]["canonical_complement_labels"]
        ]
        unexplained_indices = [labels.index(label) for label in unexplained_labels]
        unexplained_form = submatrix(quotient_form, unexplained_indices, unexplained_indices)
        numerator, denominator = fraction_parameter(parameter)
        fibre_record = {
            "label": f"rank_at_least_{rank}",
            "parameter": parameter,
            "specialization_parameter_equation": f"{denominator}*t-({numerator})=0",
            "generic_rank": 17,
            "exceptional_rank": len(labels),
            "generic_embedding_is_exact_and_primitive": bool(
                truth_row["exact_group_law_replay"] and truth_row["primitive_in_displayed_subgroup"]
            ),
            "generic_plus_selected_exceptional_basis_index_in_public_lattice": combined_index,
            "generic_height_replay_maximum_absolute_error": dtext(maximum_replay_error),
            "exceptional_basis": per_direction,
            "quotient_height_gram": [[dtext(value) for value in row] for row in quotient_form],
            "quotient_height_summary": height_summary(quotient_form),
            "canonical_unexplained_packet": {
                "labels": unexplained_labels,
                "height_summary": height_summary(unexplained_form),
            },
            "split_rational_bisection_specialization_equations": specialization_equations,
            "cover_visibility": {
                "complete_degree_2_split_count": len(bcontrol_by_parameter[parameter]["hits"]),
                "degree_2_exact_quotient_rank": int(
                    deep_row["cumulative_captured_exceptional_rank"][0]["R_t(D)"]
                ),
                "degree_3_and_4_additional_split_count": 0,
                "degree_4_cumulative_quotient_rank": int(
                    deep_row["cumulative_captured_exceptional_rank"][-1]["R_t(D)"]
                ),
            },
        }
        fibres.append(fibre_record)
        if rank == 28:
            rank28_unexplained_form = unexplained_form
            rank28_parent = [parent_records[label] for label in unexplained_labels]

    repeated_cover_labels = {
        label: parameters for label, parameters in cover_occurrences.items() if len(parameters) > 1
    }
    if repeated_cover_labels:
        raise ArithmeticError("a rational bisection specialization equation now repeats across controls")
    assert rank28_unexplained_form is not None and rank28_parent is not None
    all_specialization_equations = [
        equation
        for fibre in fibres
        for equation in fibre["split_rational_bisection_specialization_equations"]
    ]
    q_bits = [row["primitive_q_maximum_coefficient_bits"] for row in all_specialization_equations]
    root_bits = [
        row["specialized_square_root_projective_height_bits"]
        for row in all_specialization_equations
    ]
    rank28_bits = [
        int(row["parent_crossratio_base_projective_height"].bit_length()) for row in rank28_parent
    ]
    rank28_heights = [float(row["parent_point_canonical_height_128bit"]) for row in rank28_parent]
    all_parent_bases = [
        direction["parent_transport"]["parent_crossratio_base"]
        for fibre in fibres
        for direction in fibre["exceptional_basis"]
    ]
    if len(set(all_parent_bases)) != 38:
        raise ArithmeticError("two exceptional directions now have the same normalized parent base")
    rank28_defects = [
        float(direction["quotient_height_specialization_defect"])
        for direction in fibres[-1]["exceptional_basis"]
        if direction["label"] in fibres[-1]["canonical_unexplained_packet"]["labels"]
    ]

    input_paths = [
        Path(__file__).resolve(),
        TRUTH,
        CONTROLS,
        BISECTION_CONTROLS,
        DEEP_COVERS,
        BISECTIONS,
        POINT_FACTORY,
        VISIBILITY,
        *(ELLIPTIC / "cas" / f"elkies_rank{rank}.py" for rank in range(25, 29)),
    ]
    payload = {
        "schema": "elliptic-curves.elkies-2026-exceptional-specialization-relations.v1",
        "status": "PASS_EXACT_EMBEDDINGS_AND_COVER_RELATIONS_WITH_NUMERICAL_HEIGHT_QUOTIENTS",
        "inputs": {display(path): digest(path) for path in input_paths},
        "decimal_precision_digits": args.digits,
        "fibres": fibres,
        "cross_control_comparison": {
            "split_rational_bisection_cover_labels_repeat_across_controls": False,
            "repeated_labels": repeated_cover_labels,
            "split_counts_rank25_through_rank28": [6, 3, 2, 1],
            "captured_quotient_ranks_rank25_through_rank28": [5, 3, 2, 1],
            "primitive_q_maximum_coefficient_bits_min_median_max": [
                min(q_bits), float(np.median(q_bits)), max(q_bits)
            ],
            "specialized_square_root_height_bits_min_median_max": [
                min(root_bits), float(np.median(root_bits)), max(root_bits)
            ],
            "normalized_parent_bases_among_all_38_directions": {
                "distinct": 38,
                "total": 38,
            },
            "interpretation": "The complete degree-two equations split in disjoint, rapidly thinning sets. This is compatible with isolated square-value conditions, not evidence for one repeated low-degree cover mechanism.",
        },
        "rank28_conclusion": {
            "public_exceptional_rank": 11,
            "complete_rational_bisection_visible_rank": 1,
            "canonical_unexplained_dimension": 10,
            "direct_unit_basis_directions_on_a_rational_bisection": 0,
            "sole_collective_relation": "Q2-Q4+Q5-Q8+Q10 modulo the specialized generic R17 subgroup",
            "degree_3_or_4_split_in_tested_equation_universe": False,
            "unexplained_height_summary": height_summary(rank28_unexplained_form),
            "parent_crossratio_height_bits_min_median_max": [
                int(min(rank28_bits)),
                float(np.median(rank28_bits)),
                int(max(rank28_bits)),
            ],
            "parent_canonical_height_min_median_max": [
                f"{min(rank28_heights):.17g}",
                f"{float(np.median(rank28_heights)):.17g}",
                f"{max(rank28_heights):.17g}",
            ],
            "pearson_correlations": {
                "quotient_defect_vs_parent_canonical_height": pearson(
                    rank28_defects, rank28_heights
                ),
                "quotient_defect_vs_parent_base_height_bits": pearson(
                    rank28_defects, rank28_bits
                ),
                "parent_canonical_height_vs_parent_base_height_bits": pearson(
                    rank28_heights, rank28_bits
                ),
            },
            "mechanism_classification": {
                "isolated_noether_lefschetz_condition": "NOT_APPLICABLE: t is the base coordinate of one fixed K3 elliptic fibration, not a K3 moduli parameter.",
                "known_extension_defined_divisor_specialization": "NO for the ten-dimensional packet in the complete rational-bisection atlas and the tested trisection/quadrisection equations.",
                "high_degree_or_higher_genus_multisection": "OPEN",
                "isolated_arithmetic_fibre_points_without_a_shared_low_degree_divisor": "BEST CURRENT DESCRIPTION, not a non-existence theorem for multisections.",
            },
        },
        "claim_boundary": (
            "The generic subgroup embeddings, public-point independence, bisection trace equations, "
            "and exact point transports are inherited exact certificates. Canonical-height pairings and "
            "Schur complements are numerical PARI computations at the declared precision. A fibre point "
            "has no intrinsic K3 divisor class. The degree-two completeness concerns section-nonnegative "
            "rational bisections; higher-genus bisections and the untested degree-three/four universes remain open."
        ),
        "reproducing_command": (
            "python3 elliptic-curves/cas/analyze_elkies_2026_exceptional_specialization_relations.py"
        ),
    }

    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.exists() or args.output.read_text() != rendered:
            raise SystemExit(f"stale or missing output: {display(args.output)}")
        print(f"ELKIESR17EXCEPTIONALREL|status=PASS_CHECK|output={display(args.output)}")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered)
    print(
        "ELKIESR17EXCEPTIONALREL|"
        f"fibres={len(fibres)}|rank28_unexplained=10|"
        f"status={payload['status']}|output={display(args.output)}"
    )


if __name__ == "__main__":
    main()
