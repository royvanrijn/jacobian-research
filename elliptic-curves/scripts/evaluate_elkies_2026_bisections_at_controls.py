#!/usr/bin/env python3
"""Evaluate the complete rootless-R17 bisection atlas at five exact fibres.

For every stored branch quadratic ``q_i(t)`` this replay tests whether
``q_i(t_0)`` is a nonzero rational square.  Each hit is specialized in both
square-root branches, transported to the pinned global minimal fibre, and
checked exactly.  The two points are also added on a short model and checked
to have the stored trace.

The finite-quotient calculation puts the seventeen generic points first,
followed by the deterministic public complement and the split-bisection
points.  It computes their images in the direct sum of exact quotients
``E(F_p)/2E(F_p)`` for every usable prime up to 1000.  Classes are reported
modulo the generic seventeen in the ordered public-complement basis, extended
by a named escape basis if a split point leaves the known span.

Failure to escape these finite quotients is not a Mordell--Weil dependence
proof.  The replay therefore uses a height/LLL calculation only to discover a
full relation block, verifies every relation by exact rational group addition,
and reports the resulting exact rank of the displayed generated subgroup
separately from the still-unknown rank of the full curve.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from math import isqrt
import json
from pathlib import Path
import re
import sys
from typing import Any, Sequence


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC_ROOT = ROOT / "elliptic-curves"
CAS = ELLIPTIC_ROOT / "cas"
sys.path[:0] = [str(ELLIPTIC_ROOT), str(CAS)]

from ecsearch.q12o5867_specialization import (  # noqa: E402
    evaluate_projective_specialization,
    global_minimal_model_with_change,
    load_q12o5867_data,
    short_certificate_model,
)
from elliptic_candidate_record import (  # noqa: E402
    build_finite_quotient_certificate,
    is_on_weierstrass_curve,
    matrix_rank_and_pivots_mod_prime,
    source_point_to_target,
    verify_finite_quotient_certificate,
)
from search_extra_points import gp_rational, run_gp  # noqa: E402
from elkies_rank25 import POINTS as RANK25_POINTS  # noqa: E402
from elkies_rank26 import POINTS as RANK26_POINTS  # noqa: E402
from elkies_rank27 import POINTS as RANK27_POINTS  # noqa: E402
from elkies_rank28 import POINTS as RANK28_POINTS  # noqa: E402
from verify_icarm_curve394_rank21 import (  # noqa: E402
    PUBLIC_COMPLEMENT_INDICES_ZERO_BASED as CURVE394_COMPLEMENT,
    PUBLIC_MODEL as CURVE394_MODEL,
    PUBLIC_POINTS as CURVE394_POINTS,
)


Q = Fraction
PROTOCOL = "ELKIES2026BISECTIONSPECIALIZE"
STATUS = "PASS_EXACT_ELKIES_2026_BISECTION_SPECIALIZATION_CONTROLS"

MODEL = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
SECTIONS = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json"
BISECTIONS = (
    ROOT / "artifacts/generated-results/elkies-2026-equation-bisections-full.json"
)
CONTROLS = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_high_rank_positive_controls_v2.json"
)
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_bisection_specialization_controls_v1.json"
)
REPRODUCING_COMMAND = (
    ".venv/bin/python "
    "elliptic-curves/scripts/evaluate_elkies_2026_bisections_at_controls.py"
)


PUBLIC_POINTS = {
    "-2/377": RANK25_POINTS,
    "-308/251": RANK26_POINTS,
    "2456/135": RANK27_POINTS,
    "-9529/5471": RANK28_POINTS,
    "3/8": CURVE394_POINTS,
}

PARAMETERS = (
    ("rank_at_least_25", -2, 377, 25),
    ("rank_at_least_26", -308, 251, 26),
    ("rank_at_least_27", 2456, 135, 27),
    ("rank_at_least_28", -9529, 5471, 28),
    ("icarm_curve394_rank_at_least_21", 3, 8, 21),
)

# These values make the complete census fail closed if the atlas, fibre
# normalization, public complement, or quotient convention changes.
EXPECTED = {
    "-2/377": {
        "labels": [
            "orbit-1cb25", "orbit-0cff7", "orbit-1ea09",
            "orbit-051a1", "orbit-0d4ca", "orbit-1d5bb",
        ],
        "classes": [
            "10000000", "10000010", "00000110",
            "00100000", "00000100", "00001000",
        ],
        "span": 5,
    },
    "-308/251": {
        "labels": ["orbit-0da89", "orbit-12c1b", "orbit-1ea54"],
        "classes": ["111010100", "100010111", "110010111"],
        "span": 3,
    },
    "2456/135": {
        "labels": ["orbit-195a4", "orbit-00edf"],
        "classes": ["0100000110", "1010000000"],
        "span": 2,
    },
    "-9529/5471": {
        "labels": ["orbit-15a68"],
        "classes": ["01011001010"],
        "span": 1,
    },
    "3/8": {
        "labels": [
            "orbit-05980", "orbit-04d17", "orbit-101f2", "orbit-07843",
            "orbit-090a3", "orbit-05443", "orbit-18fd5", "orbit-1f786",
            "orbit-02d31", "orbit-055ad", "orbit-0be21", "orbit-00ca6",
            "orbit-0976c", "orbit-045c2", "orbit-0fb68", "orbit-08e3a",
            "orbit-196a3", "orbit-01926", "orbit-0888a", "orbit-06faa",
            "orbit-10aaa", "orbit-06f04", "orbit-0eba4", "orbit-01e36",
            "orbit-126e6",
        ],
        "classes": [
            "1100", "0100", "0101", "0011", "1000",
            "1000", "0100", "0111", "0011", "1000",
            "0000", "0100", "1000", "1000", "1000",
            "0010", "1000", "1100", "0100", "1000",
            "0000", "0000", "0100", "1100", "0000",
        ],
        "span": 4,
    },
}


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def fraction_text(value: Fraction | int) -> str:
    value = Q(value)
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def point_record(point: tuple[Fraction, Fraction]) -> list[str]:
    return [fraction_text(point[0]), fraction_text(point[1])]


def evaluate_polynomial(coefficients: Sequence[str], value: Fraction) -> Fraction:
    answer = Q(0)
    for coefficient in reversed(coefficients):
        answer = answer * value + Q(coefficient)
    return answer


def nonzero_rational_square_root(value: Fraction) -> Fraction | None:
    value = Q(value)
    if value <= 0:
        return None
    numerator = isqrt(value.numerator)
    denominator = isqrt(value.denominator)
    if numerator**2 != value.numerator or denominator**2 != value.denominator:
        return None
    return Q(numerator, denominator)


def short_add(
    left: tuple[Fraction, Fraction] | None,
    right: tuple[Fraction, Fraction] | None,
    coefficient_a: Fraction,
) -> tuple[Fraction, Fraction] | None:
    if left is None:
        return right
    if right is None:
        return left
    x_left, y_left = left
    x_right, y_right = right
    if x_left == x_right:
        if y_left == -y_right:
            return None
        if y_left == 0:
            return None
        slope = (3 * x_left**2 + coefficient_a) / (2 * y_left)
    else:
        slope = (y_right - y_left) / (x_right - x_left)
    x_answer = slope**2 - x_left - x_right
    return x_answer, -y_left + slope * (x_left - x_answer)


def short_multiply(
    point: tuple[Fraction, Fraction] | None,
    scalar: int,
    coefficient_a: Fraction,
) -> tuple[Fraction, Fraction] | None:
    if scalar < 0:
        if point is not None:
            point = point[0], -point[1]
        return short_multiply(point, -scalar, coefficient_a)
    answer = None
    addend = point
    while scalar:
        if scalar & 1:
            answer = short_add(answer, addend, coefficient_a)
        addend = short_add(addend, addend, coefficient_a)
        scalar >>= 1
    return answer


def short_linear_combination(
    points: Sequence[tuple[Fraction, Fraction]],
    coefficients: Sequence[int],
    coefficient_a: Fraction,
) -> tuple[Fraction, Fraction] | None:
    if len(points) != len(coefficients):
        raise ValueError("point and coefficient counts differ")
    answer = None
    for point, coefficient in zip(points, coefficients):
        answer = short_add(
            answer, short_multiply(point, int(coefficient), coefficient_a), coefficient_a
        )
    return answer


def rational_matrix_rank(rows: Sequence[Sequence[int]]) -> int:
    if not rows:
        return 0
    matrix = [[Q(value) for value in row] for row in rows]
    column_count = len(matrix[0])
    rank = 0
    for column in range(column_count):
        pivot = next(
            (index for index in range(rank, len(matrix)) if matrix[index][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        pivot_value = matrix[rank][column]
        matrix[rank] = [value / pivot_value for value in matrix[rank]]
        for index in range(len(matrix)):
            if index == rank or not matrix[index][column]:
                continue
            multiple = matrix[index][column]
            matrix[index] = [
                left - multiple * right
                for left, right in zip(matrix[index], matrix[rank])
            ]
        rank += 1
    return rank


def discover_and_verify_relations(
    short_model: Sequence[Fraction],
    points: Sequence[tuple[Fraction, Fraction]],
    existing_count: int,
) -> dict[str, Any]:
    """Discover a relation lattice numerically, then prove it by exact addition."""

    model_text = ",".join(gp_rational(Q(value)) for value in short_model)
    point_text = ",".join(
        f"[{gp_rational(point[0])},{gp_rational(point[1])}]" for point in points
    )
    program = f"""default(realprecision,100);
E=ellinit([{model_text}]);
P=[{point_text}];
M=ellheightmatrix(E,P);
T=qflllgram(M);
print("HEIGHTRANK|",matrank(M));
for(j=1,#P,if(abs(T[,j]~*M*T[,j])<1e-50,print("REL|",Vec(T[,j]))));
quit
"""
    output, _elapsed = run_gp(program, timeout=60.0, stack_bytes=512_000_000)
    rank_match = re.search(r"^HEIGHTRANK\|(\d+)$", output, re.MULTILINE)
    if rank_match is None:
        raise ArithmeticError("PARI omitted the numerical height rank")
    discovered = []
    for payload in re.findall(r"^REL\|\[([^]]*)\]$", output, re.MULTILINE):
        vector = [int(value.strip()) for value in payload.split(",")]
        if len(vector) != len(points):
            raise ArithmeticError("PARI returned a relation of the wrong length")
        discovered.append(vector)
    new_count = len(points) - existing_count
    if len(discovered) < new_count:
        raise ArithmeticError("PARI did not expose enough candidate relations")

    selected = []
    new_rank = 0
    for vector in discovered:
        candidate = selected + [vector]
        candidate_rank = rational_matrix_rank(
            [row[existing_count:] for row in candidate]
        )
        if candidate_rank > new_rank:
            selected.append(vector)
            new_rank = candidate_rank
        if new_rank == new_count:
            break
    if new_rank != new_count:
        raise ArithmeticError("the candidate relations do not eliminate every new point")

    coefficient_a = Q(short_model[3])
    for vector in selected:
        if short_linear_combination(points, vector, coefficient_a) is not None:
            raise ArithmeticError("a numerical height relation failed exact group addition")
    return {
        "discovery_engine": "PARI/GP ellheightmatrix plus qflllgram",
        "real_precision_digits": 100,
        "zero_height_threshold": "1e-50",
        "numerical_height_matrix_rank": int(rank_match.group(1)),
        "relation_count": len(selected),
        "new_coefficient_block_rank_over_q": new_rank,
        "all_relations_verified_by_exact_group_addition": True,
        "relations": selected,
    }


def solve_columns_mod_two(
    rows: Sequence[Sequence[int]], basis_columns: Sequence[int], target_column: int
) -> list[int] | None:
    """Solve target as a unique combination of independent basis columns."""

    width = len(basis_columns)
    matrix = [
        [int(row[column]) & 1 for column in basis_columns]
        + [int(row[target_column]) & 1]
        for row in rows
    ]
    pivot_rows: dict[int, int] = {}
    next_row = 0
    for column in range(width):
        pivot = next(
            (
                row_index
                for row_index in range(next_row, len(matrix))
                if matrix[row_index][column]
            ),
            None,
        )
        if pivot is None:
            raise ArithmeticError("the declared quotient basis is not independent")
        matrix[next_row], matrix[pivot] = matrix[pivot], matrix[next_row]
        for row_index in range(len(matrix)):
            if row_index != next_row and matrix[row_index][column]:
                matrix[row_index] = [
                    left ^ right
                    for left, right in zip(matrix[row_index], matrix[next_row])
                ]
        pivot_rows[column] = next_row
        next_row += 1
    if any(not any(row[:width]) and row[width] for row in matrix):
        return None
    return [matrix[pivot_rows[column]][width] for column in range(width)]


def column_span_rank(columns: Sequence[Sequence[int]]) -> int:
    if not columns:
        return 0
    coordinate_count = len(columns[0])
    rows = [
        [int(column[coordinate]) for column in columns]
        for coordinate in range(coordinate_count)
    ]
    rank, _pivots = matrix_rank_and_pivots_mod_prime(rows, len(columns), 2)
    return rank


def signature_rows(certificate: dict[str, Any]) -> list[list[int]]:
    return [
        list(row)
        for signature in certificate["signatures"]
        for row in signature["rows"]
    ]


def rows_sha256(rows: Sequence[Sequence[int]]) -> str:
    digest = sha256()
    for row in rows:
        digest.update(("".join(str(int(value)) for value in row) + "\n").encode())
    return digest.hexdigest()


def public_complements(controls: dict[str, Any]) -> dict[str, tuple[int, ...]]:
    answer = {
        record["parameter"]: tuple(
            int(index) - 1
            for index in record["public_positive_control"][
                "selected_public_point_indices_one_based"
            ]
        )
        for record in controls["fibres"]
    }
    answer["3/8"] = tuple(int(index) for index in CURVE394_COMPLEMENT)
    return answer


def expected_minimal_models(controls: dict[str, Any]) -> dict[str, tuple[Fraction, ...]]:
    answer = {
        record["parameter"]: tuple(Q(value) for value in record["minimal_model"])
        for record in controls["fibres"]
    }
    answer["3/8"] = tuple(Q(value) for value in CURVE394_MODEL)
    return answer


def construct_hits(
    records: Sequence[dict[str, Any]], numerator: int, denominator: int
) -> list[dict[str, Any]]:
    parameter = Q(numerator, denominator)
    hits: list[dict[str, Any]] = []
    for record in records:
        q_value = evaluate_polynomial(
            record["residual_chord"]["q_coefficients"], parameter
        )
        square_root = nonzero_rational_square_root(q_value)
        if square_root is None:
            continue
        lifted = record["lifted_section"]
        x_zero = evaluate_polynomial(lifted["x0_coefficients"], parameter)
        x_one = evaluate_polynomial(lifted["x1_coefficients"], parameter)
        y_zero = evaluate_polynomial(lifted["y0_coefficients"], parameter)
        y_one = evaluate_polynomial(lifted["y1_coefficients"], parameter)
        scale_x = denominator**4
        scale_y = denominator**6
        positive_point = (
            Q(scale_x) * (x_zero + x_one * square_root),
            Q(scale_y) * (y_zero + y_one * square_root),
        )
        negative_point = (
            Q(scale_x) * (x_zero - x_one * square_root),
            Q(scale_y) * (y_zero - y_one * square_root),
        )
        trace = record["trace_section"]
        h_value = evaluate_polynomial(trace["h_coefficients"], parameter)
        if h_value == 0:
            raise ArithmeticError("a split hit specializes its trace to infinity")
        trace_point = (
            Q(scale_x)
            * evaluate_polynomial(trace["Nx_coefficients"], parameter)
            / h_value**2,
            Q(scale_y)
            * evaluate_polynomial(trace["Ny_coefficients"], parameter)
            / h_value**3,
        )
        hits.append(
            {
                "source_record": record,
                "q_value": q_value,
                "square_root": square_root,
                "positive_source_point": positive_point,
                "negative_source_point": negative_point,
                "trace_source_point": trace_point,
            }
        )
    return hits


def evaluate_fibre(
    *,
    data,
    records: Sequence[dict[str, Any]],
    complements: dict[str, tuple[int, ...]],
    minimal_models: dict[str, tuple[Fraction, ...]],
    label: str,
    numerator: int,
    denominator: int,
    known_rank: int,
) -> dict[str, Any]:
    parameter_text = fraction_text(Q(numerator, denominator))
    hits = construct_hits(records, numerator, denominator)
    expectation = EXPECTED[parameter_text]
    hit_labels = [hit["source_record"]["label"] for hit in hits]
    if hit_labels != expectation["labels"]:
        raise AssertionError(
            f"the split-bisection census changed at {parameter_text}: {hit_labels}"
        )

    specialization = evaluate_projective_specialization(data, numerator, denominator)
    for hit in hits:
        for point_name in (
            "positive_source_point", "negative_source_point", "trace_source_point"
        ):
            if not is_on_weierstrass_curve(specialization.model, hit[point_name]):
                raise ArithmeticError(
                    f"{hit['source_record']['label']} missed the source fibre"
                )

    minimal_model, minimal_change, minimalization = global_minimal_model_with_change(
        specialization.model
    )
    if minimal_model != minimal_models[parameter_text]:
        raise ArithmeticError(f"the pinned minimal model changed at {parameter_text}")
    generic_minimal = tuple(
        source_point_to_target(point, minimal_change)
        for point in specialization.points
    )
    public_points = tuple(PUBLIC_POINTS[parameter_text])
    complement_indices = complements[parameter_text]
    public_complement = tuple(public_points[index] for index in complement_indices)
    existing_minimal = generic_minimal + public_complement
    if len(existing_minimal) != known_rank:
        raise AssertionError("the known public basis has the wrong size")

    short_model, short_change = short_certificate_model(minimal_model)
    existing_short = tuple(
        source_point_to_target(point, short_change) for point in existing_minimal
    )
    positive_short = []
    negative_short = []
    trace_short = []
    for hit in hits:
        positive_minimal = source_point_to_target(
            hit["positive_source_point"], minimal_change
        )
        negative_minimal = source_point_to_target(
            hit["negative_source_point"], minimal_change
        )
        trace_minimal = source_point_to_target(
            hit["trace_source_point"], minimal_change
        )
        for point in (positive_minimal, negative_minimal, trace_minimal):
            if not is_on_weierstrass_curve(minimal_model, point):
                raise ArithmeticError("a transported point missed the minimal fibre")
        positive = source_point_to_target(positive_minimal, short_change)
        negative = source_point_to_target(negative_minimal, short_change)
        trace = source_point_to_target(trace_minimal, short_change)
        if short_add(positive, negative, Q(short_model[3])) != trace:
            raise ArithmeticError("the two split points do not have the stored trace")
        trace_vector = [int(value) for value in hit["source_record"]["published_basis_w"]]
        if (
            short_linear_combination(
                existing_short[:17], trace_vector, Q(short_model[3])
            )
            != trace
        ):
            raise ArithmeticError("the stored trace missed its published-basis word")
        positive_short.append(positive)
        negative_short.append(negative)
        trace_short.append(trace)
        hit["positive_minimal_point"] = positive_minimal
        hit["negative_minimal_point"] = negative_minimal

    all_short = existing_short + tuple(positive_short)
    joint = build_finite_quotient_certificate(
        short_model, all_short, relation_prime=2, prime_bound=1000
    )
    rows = signature_rows(joint)
    existing_rank, _ = matrix_rank_and_pivots_mod_prime(
        [row[:known_rank] for row in rows], known_rank, 2
    )
    if existing_rank != known_rank:
        raise ArithmeticError("the finite quotient lost a known public direction")
    joint_rank = int(joint["combined_rank_over_relation_field"])
    if joint_rank < known_rank:
        raise ArithmeticError("the joint finite quotient rank fell below the basis rank")

    quotient_basis_columns = list(range(known_rank))
    quotient_basis_labels = [
        f"public-complement-Q{index + 1}" for index in range(known_rank - 17)
    ]
    provisional_classes = []
    escape_labels = []
    for hit_index, hit in enumerate(hits):
        target_column = known_rank + hit_index
        coefficients = solve_columns_mod_two(
            rows, quotient_basis_columns, target_column
        )
        if coefficients is None:
            quotient_basis_columns.append(target_column)
            escape_label = f"escape-{hit['source_record']['label']}"
            quotient_basis_labels.append(escape_label)
            escape_labels.append(escape_label)
            coefficients = [0] * (len(quotient_basis_columns) - 1) + [1]
        provisional_classes.append(
            {
                "generic_correction": coefficients[:17],
                "quotient_coordinates": coefficients[17:],
            }
        )
    final_quotient_dimension = len(quotient_basis_labels)
    for item in provisional_classes:
        item["quotient_coordinates"].extend(
            [0] * (final_quotient_dimension - len(item["quotient_coordinates"]))
        )
    finite_gain = len(escape_labels)
    if joint_rank != known_rank + finite_gain:
        raise AssertionError("escape-basis dimension disagrees with joint rank")

    class_columns = [item["quotient_coordinates"] for item in provisional_classes]
    class_span = column_span_rank(class_columns)
    public_dimension = known_rank - 17
    public_class_columns = [column[:public_dimension] for column in class_columns]
    public_span = column_span_rank(public_class_columns)
    class_strings = [
        "".join(str(bit) for bit in column[:public_dimension])
        for column in class_columns
    ]
    if finite_gain != 0:
        raise AssertionError(
            f"a new finite-quotient escape appeared at {parameter_text}: {escape_labels}"
        )
    if class_strings != expectation["classes"] or public_span != expectation["span"]:
        raise AssertionError(f"the quotient-class atlas changed at {parameter_text}")

    independent_indices = tuple(int(index) for index in joint["pivot_columns_zero_based"])
    independent_points = tuple(all_short[index] for index in independent_indices)
    rank_certificate = build_finite_quotient_certificate(
        short_model, independent_points, relation_prime=2, prime_bound=1000
    )
    verify_finite_quotient_certificate(short_model, independent_points, rank_certificate)
    if (
        not rank_certificate["certified_independent"]
        or rank_certificate["certified_rank_lower_bound"] != joint_rank
    ):
        raise ArithmeticError("the independent joint subbasis was not certified")

    exact_relations = discover_and_verify_relations(
        short_model, all_short, known_rank
    )
    if (
        exact_relations["relation_count"] != len(hits)
        or exact_relations["new_coefficient_block_rank_over_q"] != len(hits)
        or exact_relations["numerical_height_matrix_rank"] != known_rank
    ):
        raise ArithmeticError("the exact relation closeout has the wrong rank")

    hit_records = []
    for hit, quotient_class in zip(hits, provisional_classes):
        source = hit["source_record"]
        support = [
            quotient_basis_labels[index]
            for index, bit in enumerate(quotient_class["quotient_coordinates"])
            if bit
        ]
        hit_records.append(
            {
                "label": source["label"],
                "lattice_orbit_mask": int(source["lattice_orbit_mask"]),
                "q_value": fraction_text(hit["q_value"]),
                "canonical_positive_square_root": fraction_text(hit["square_root"]),
                "positive_minimal_point": point_record(hit["positive_minimal_point"]),
                "negative_minimal_point": point_record(hit["negative_minimal_point"]),
                "exact_verification": {
                    "both_points_on_projective_source_fibre": True,
                    "both_points_on_global_minimal_fibre": True,
                    "sum_of_branches_equals_stored_trace": True,
                    "stored_trace_equals_published_basis_word": True,
                },
                "finite_quotient_class_modulo_generic_17": {
                    "basis": quotient_basis_labels,
                    "coordinates_over_f2": quotient_class["quotient_coordinates"],
                    "support": support,
                    "generic_correction_over_f2": quotient_class[
                        "generic_correction"
                    ],
                    "relation": (
                        "In the scanned direct sum of E(F_p)/2E(F_p), the point "
                        "equals the displayed generic correction plus quotient class."
                    ),
                },
            }
        )

    signature_summary = {
        "relation_prime": 2,
        "all_good_reduction_primes_scanned_through": 1000,
        "rank_increasing_primes": joint["certificate_primes"],
        "torsion_witness": joint["torsion_witness"],
        "rank_increasing_prime_data": [
            {
                "prime": signature["prime"],
                "group_order": signature["group_order"],
                "quotient_dimension": signature["quotient_dimension"],
            }
            for signature in joint["signatures"]
        ],
        "stacked_rows_sha256": rows_sha256(rows),
        "existing_public_basis_rank": existing_rank,
        "rank_after_adjoining_split_points": joint_rank,
        "gain_beyond_existing_public_basis": finite_gain,
        "pivot_columns_zero_based": list(independent_indices),
        "certified_independent_subbasis": {
            "point_count": rank_certificate["point_count"],
            "certificate_primes": rank_certificate["certificate_primes"],
            "torsion_witness": rank_certificate["torsion_witness"],
            "certified_rank_lower_bound": rank_certificate[
                "certified_rank_lower_bound"
            ],
        },
    }
    return {
        "label": label,
        "parameter": parameter_text,
        "projective_parameter": [numerator, denominator],
        "square_tests": len(records),
        "split_bisection_count": len(hits),
        "public_complement": {
            "dimension": public_dimension,
            "source_point_indices_one_based": [index + 1 for index in complement_indices],
            "ordered_basis_labels": quotient_basis_labels[:public_dimension],
        },
        "split_class_span": {
            "dimension_modulo_generic_17": class_span,
            "dimension_in_known_public_complement": public_span,
            "known_public_complement_dimension": public_dimension,
            "finite_quotient_escape_count": finite_gain,
            "escape_basis_labels": escape_labels,
            "distinct_known_public_classes": len(set(class_strings)),
        },
        "rank_result": {
            "existing_unconditional_rank_lower_bound": known_rank,
            "finite_quotient_signature_rank_after_adjoining": joint_rank,
            "unconditional_rank_lower_bound_after_adjoining": joint_rank,
            "exact_generated_subgroup_rank_after_adjoining": known_rank,
            "exact_generated_subgroup_rank_determined": True,
            "relation_basis": {
                "ordered_points": (
                    [f"generic-P{index + 1}" for index in range(17)]
                    + quotient_basis_labels[:public_dimension]
                    + [f"split-{hit['source_record']['label']}" for hit in hits]
                ),
                **exact_relations,
            },
            "boundary": (
                "The exact displayed subgroup rank is determined by independent "
                "known points and a full-rank exact relation block for all split "
                "points. This is not an upper bound for the full Mordell--Weil group."
            ),
        },
        "minimal_model": [fraction_text(value) for value in minimal_model],
        "minimalization": minimalization,
        "finite_quotient": signature_summary,
        "hits": hit_records,
    }


def build_payload() -> dict[str, Any]:
    sys.set_int_max_str_digits(0)
    data = load_q12o5867_data(MODEL, SECTIONS)
    atlas = json.loads(BISECTIONS.read_text())
    records = atlas["bisections"]
    if len(records) != 39120:
        raise AssertionError("the complete bisection atlas no longer has 39,120 rows")
    controls = json.loads(CONTROLS.read_text())
    complements = public_complements(controls)
    minimal_models = expected_minimal_models(controls)
    fibres = [
        evaluate_fibre(
            data=data,
            records=records,
            complements=complements,
            minimal_models=minimal_models,
            label=label,
            numerator=numerator,
            denominator=denominator,
            known_rank=known_rank,
        )
        for label, numerator, denominator, known_rank in PARAMETERS
    ]
    return {
        "schema": "elliptic-curves.elkies-2026-bisection-specialization-controls.v1",
        "status": STATUS,
        "claim": (
            "Complete exact specialization of all 39,120 stored bisections at the "
            "four rank-25--28 controls and ICARM curve 394."
        ),
        "claim_boundary": (
            "The square census, point identities, finite-quotient classes, class-span "
            "dimensions, exact group relations, and ranks of the displayed generated "
            "subgroups are exact. These subgroup ranks are not upper bounds for the "
            "full Mordell--Weil groups, so no exact curve-rank claim is made."
        ),
        "atlas_size": len(records),
        "total_square_tests": len(records) * len(PARAMETERS),
        "fibres": fibres,
        "summary": {
            "parameters": [fibre["parameter"] for fibre in fibres],
            "split_counts": [fibre["split_bisection_count"] for fibre in fibres],
            "known_complement_dimensions": [
                fibre["public_complement"]["dimension"] for fibre in fibres
            ],
            "split_class_span_dimensions": [
                fibre["split_class_span"]["dimension_in_known_public_complement"]
                for fibre in fibres
            ],
            "finite_quotient_escape_counts": [
                fibre["split_class_span"]["finite_quotient_escape_count"]
                for fibre in fibres
            ],
            "rank28_outcome": (
                "The sole split bisection spans one of eleven known exceptional "
                "finite-quotient directions; ten directions are invisible."
            ),
            "curve394_outcome": (
                "The 25 split bisections span all four known directions beyond R17."
            ),
        },
        "generation": {
            "command": REPRODUCING_COMMAND,
            "checker_sha256": file_sha256(Path(__file__)),
            "inputs": {
                relative(path): file_sha256(path)
                for path in (MODEL, SECTIONS, BISECTIONS, CONTROLS)
            },
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    payload = build_payload()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if not arguments.output.is_file() or arguments.output.read_text() != rendered:
            raise SystemExit(f"stale pinned certificate: rerun {REPRODUCING_COMMAND}")
        terminal = "PASS"
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered)
        terminal = "WROTE"
    summary = payload["summary"]
    print(
        f"{PROTOCOL}|hits={','.join(map(str, summary['split_counts']))}"
        f"|spans={','.join(map(str, summary['split_class_span_dimensions']))}"
        f"|escapes={','.join(map(str, summary['finite_quotient_escape_counts']))}"
        f"|artifact={arguments.output}|status={terminal}"
    )


if __name__ == "__main__":
    main()
