#!/usr/bin/env python3
"""Map split R17 cover points exactly into the known exceptional quotients.

The degree-two baseline is read from its existing exact relation certificate.
New degree-three points are transported from the published projective fibre to
the pinned minimal/short model.  PARI heights are used only to discover a full
relation block; every relation is then checked by exact rational group
addition.  Solving that block gives coordinates in the ordered known basis and
hence exact classes in L_t/M_t.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction
from hashlib import sha256
import json
from math import gcd
from pathlib import Path
import sys
from typing import Any, Sequence


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC_ROOT = ROOT / "elliptic-curves"
CAS = ELLIPTIC_ROOT / "cas"
sys.path[:0] = [str(ELLIPTIC_ROOT), str(CAS), str(ELLIPTIC_ROOT / "scripts")]

from ecsearch.q12o5867_specialization import (  # noqa: E402
    evaluate_projective_specialization,
    global_minimal_model_with_change,
    load_q12o5867_data,
    short_certificate_model,
)
from elliptic_candidate_record import (  # noqa: E402
    is_on_weierstrass_curve,
    source_point_to_target,
)
from elkies_rank25 import POINTS as RANK25_POINTS  # noqa: E402
from elkies_rank26 import POINTS as RANK26_POINTS  # noqa: E402
from elkies_rank27 import POINTS as RANK27_POINTS  # noqa: E402
from elkies_rank28 import POINTS as RANK28_POINTS  # noqa: E402
from evaluate_elkies_2026_bisections_at_controls import (  # noqa: E402
    discover_and_verify_relations,
    public_complements,
    short_linear_combination,
)


Q = Fraction
MODEL = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
SECTIONS = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json"
CONTROLS = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_high_rank_positive_controls_v2.json"
)
BISECTION_SPECIALIZATIONS = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_bisection_specialization_controls_v1.json"
)
BISECTION_ATLAS = ROOT / "artifacts/generated-results/elkies-2026-equation-bisections-full.json"
TRISECTIONS = ROOT / "artifacts/generated-results/elkies-k3-r17-deep-trisections-v1.json"
FRONTIER_TRISECTIONS = (
    ROOT / "artifacts/generated-results/elkies-k3-r17-sampled-frontier-trisections-v1.json"
)
QUADRISECTIONS = ROOT / "artifacts/generated-results/elkies-k3-r17-sampled-quadrisections-v1.json"
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_deep_cover_exceptional_quotients_v1.json"
)

PUBLIC_POINTS = {
    "-2/377": RANK25_POINTS,
    "-308/251": RANK26_POINTS,
    "2456/135": RANK27_POINTS,
    "-9529/5471": RANK28_POINTS,
}
PARAMETERS = (
    ("rank_at_least_28", -9529, 5471, 28),
    ("rank_at_least_27", 2456, 135, 27),
    ("rank_at_least_26", -308, 251, 26),
    ("rank_at_least_25", -2, 377, 25),
)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def qtext(value: Fraction | int) -> str:
    value = Q(value)
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def lcm_int(left: int, right: int) -> int:
    return abs(left * right) // gcd(left, right) if left and right else 0


def matrix_rank(rows: Sequence[Sequence[Fraction | int]]) -> int:
    if not rows:
        return 0
    work = [[Q(value) for value in row] for row in rows]
    row_count = len(work)
    column_count = len(work[0])
    rank = 0
    for column in range(column_count):
        pivot = next((row for row in range(rank, row_count) if work[row][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        divisor = work[rank][column]
        work[rank] = [value / divisor for value in work[rank]]
        for row in range(row_count):
            if row == rank or not work[row][column]:
                continue
            multiple = work[row][column]
            work[row] = [
                left - multiple * right for left, right in zip(work[row], work[rank])
            ]
        rank += 1
    return rank


def column_rank(columns: Sequence[Sequence[Fraction | int]]) -> int:
    if not columns:
        return 0
    return matrix_rank(
        [[column[row] for column in columns] for row in range(len(columns[0]))]
    )


def inverse(square: Sequence[Sequence[Fraction | int]]) -> list[list[Fraction]]:
    size = len(square)
    if any(len(row) != size for row in square):
        raise ValueError("matrix is not square")
    work = [
        [Q(value) for value in row]
        + [Q(int(row_index == column)) for column in range(size)]
        for row_index, row in enumerate(square)
    ]
    for column in range(size):
        pivot = next((row for row in range(column, size) if work[row][column]), None)
        if pivot is None:
            raise ArithmeticError("relation block is singular")
        work[column], work[pivot] = work[pivot], work[column]
        divisor = work[column][column]
        work[column] = [value / divisor for value in work[column]]
        for row in range(size):
            if row == column or not work[row][column]:
                continue
            multiple = work[row][column]
            work[row] = [
                left - multiple * right for left, right in zip(work[row], work[column])
            ]
    return [row[size:] for row in work]


def multiply(left, right):
    return [
        [
            sum((Q(left[row][k]) * Q(right[k][column]) for k in range(len(right))), Q(0))
            for column in range(len(right[0]))
        ]
        for row in range(len(left))
    ]


def coordinates_from_relations(
    relations: Sequence[Sequence[int]], existing_count: int, new_count: int
) -> list[list[Fraction]]:
    if len(relations) != new_count:
        raise ArithmeticError("relation count does not match the new point count")
    old = [[Q(value) for value in row[:existing_count]] for row in relations]
    new = [[Q(value) for value in row[existing_count:]] for row in relations]
    if any(len(row) != new_count for row in new):
        raise ArithmeticError("relation block width mismatch")
    solved = multiply(inverse(new), old)
    return [[-value for value in row] for row in solved]


def solve_independent_columns(columns, target):
    """Solve an exact full-column-rank system, or return None outside its span."""

    if not columns:
        return [] if not any(target) else None
    dimension = len(target)
    width = len(columns)
    augmented = [
        [Q(columns[column][row]) for column in range(width)] + [Q(target[row])]
        for row in range(dimension)
    ]
    pivot_row = 0
    pivot_for_column = {}
    for column in range(width):
        pivot = next(
            (row for row in range(pivot_row, dimension) if augmented[row][column]),
            None,
        )
        if pivot is None:
            raise ArithmeticError("declared independent columns are dependent")
        augmented[pivot_row], augmented[pivot] = augmented[pivot], augmented[pivot_row]
        divisor = augmented[pivot_row][column]
        augmented[pivot_row] = [value / divisor for value in augmented[pivot_row]]
        for row in range(dimension):
            if row == pivot_row or not augmented[row][column]:
                continue
            multiple = augmented[row][column]
            augmented[row] = [
                left - multiple * right
                for left, right in zip(augmented[row], augmented[pivot_row])
            ]
        pivot_for_column[column] = pivot_row
        pivot_row += 1
    if any(not any(row[:width]) and row[width] for row in augmented):
        return None
    return [augmented[pivot_for_column[column]][width] for column in range(width)]


def verify_coordinate_row(short_model, existing_short, point, coordinates):
    common = 1
    for value in coordinates:
        common = lcm_int(common, value.denominator)
    integer_coordinates = [int(common * value) for value in coordinates]
    relation_points = tuple(existing_short) + (point,)
    relation_coefficients = integer_coordinates + [-common]
    if short_linear_combination(
        relation_points, relation_coefficients, Q(short_model[3])
    ) is not None:
        raise ArithmeticError("a solved coordinate row failed exact group addition")
    return common


def baseline_classes(fibre, atlas_by_label):
    relation_basis = fibre["rank_result"]["relation_basis"]
    split_count = fibre["split_bisection_count"]
    existing_count = fibre["rank_result"]["existing_unconditional_rank_lower_bound"]
    coordinates = coordinates_from_relations(
        relation_basis["relations"], existing_count, split_count
    )
    exceptional_dimension = existing_count - 17
    records = []
    for hit, coordinate in zip(fibre["hits"], coordinates):
        source = atlas_by_label[hit["label"]]
        if any(value.denominator != 1 for value in coordinate):
            raise ArithmeticError("a baseline point is not integral in the known basis")
        records.append(
            {
                "point_label": f"d2-{hit['label']}",
                "cover_label": hit["label"],
                "cover_degree": 2,
                "minimum_lattice_norm": 10,
                "known_basis_coordinates": [qtext(value) for value in coordinate],
                "exceptional_quotient_coordinates": [
                    qtext(value) for value in coordinate[-exceptional_dimension:]
                ],
                "coordinate_denominator": 1,
                "equation_complexity": source["equation_complexity"],
                "exact_relation_reused_from_degree2_certificate": True,
            }
        )
    return records


def cover_points_for_parameter(payload, parameter, degree, minimum_norm):
    records = []
    for cover in payload["records"]:
        specialization = next(
            item for item in cover["specializations"] if item["parameter"] == parameter
        )
        for root_index, point in enumerate(specialization["projective_source_points"], start=1):
            records.append(
                {
                    "point_label": f"d{degree}-{cover['label']}-root{root_index}",
                    "cover_label": cover["label"],
                    "root_index": root_index,
                    "source_point": (Q(point[0]), Q(point[1])),
                    "cover_degree": degree,
                    "minimum_lattice_norm": minimum_norm,
                    "equation_complexity": cover["equation_complexity"],
                    "factor_degrees": specialization["factor_degrees"],
                }
            )
    return records


def complexity_key(record):
    complexity = record["equation_complexity"]
    return (
        int(record["cover_degree"]),
        int(record["minimum_lattice_norm"]),
        int(complexity["group_addition_upper_bound"]),
        int(complexity["support_count"]),
        int(complexity["dependency_count"]),
        int(complexity["coordinate_input_bits"]),
        int(complexity["maximum_absolute_coefficient"]),
        int(complexity["coefficient_l1"]),
        record["cover_label"],
        record["point_label"],
    )


def vector_from_record(record):
    return [Q(value) for value in record["exceptional_quotient_coordinates"]]


def matroid_structure(records):
    ordered = sorted(records, key=complexity_key)
    independent_vectors = []
    independent_labels = []
    circuits = []
    for record in ordered:
        value = vector_from_record(record)
        if not any(value):
            circuits.append(
                {"type": "loop", "point_labels": [record["point_label"]], "coefficients": ["1"]}
            )
            continue
        coefficients = solve_independent_columns(independent_vectors, value)
        if coefficients is None:
            independent_vectors.append(value)
            independent_labels.append(record["point_label"])
            continue
        support_labels = []
        support_coefficients = []
        for label, coefficient in zip(independent_labels, coefficients):
            if coefficient:
                support_labels.append(label)
                support_coefficients.append(qtext(coefficient))
        support_labels.append(record["point_label"])
        support_coefficients.append("-1")
        circuits.append(
            {
                "type": "fundamental_circuit",
                "point_labels": support_labels,
                "coefficients": support_coefficients,
            }
        )

    parallel = defaultdict(list)
    for record in ordered:
        value = vector_from_record(record)
        first = next((entry for entry in value if entry), None)
        if first is None:
            continue
        normalized = tuple(entry / first for entry in value)
        parallel[normalized].append((record["point_label"], qtext(first)))
    parallel_classes = [
        {
            "normalized_direction": [qtext(value) for value in direction],
            "members": [
                {"point_label": label, "scalar": scalar} for label, scalar in members
            ],
        }
        for direction, members in parallel.items()
        if len(members) > 1
    ]
    return {
        "rank": len(independent_vectors),
        "canonical_independent_point_labels": independent_labels,
        "parallel_duplicate_classes": parallel_classes,
        "fundamental_circuits": circuits,
        "boundary": (
            "The displayed fundamental circuits are relative to the deterministic "
            "complexity-ordered basis; they generate all linear dependencies but are "
            "not an enumeration of every circuit of the vector matroid."
        ),
    }


def exposure_profile(records, dimension):
    covers = defaultdict(list)
    representative = {}
    for record in records:
        covers[record["cover_label"]].append(vector_from_record(record))
        representative[record["cover_label"]] = record
    ordered_labels = sorted(covers, key=lambda label: complexity_key(representative[label]))
    accumulated = []
    rank_steps = []
    direction_exposure = [None] * dimension
    for label in ordered_labels:
        old_rank = column_rank(accumulated)
        accumulated.extend(covers[label])
        new_rank = column_rank(accumulated)
        if new_rank > old_rank:
            row = representative[label]
            rank_steps.append(
                {
                    "cover_label": label,
                    "exceptional_quotient_vectors": [
                        [qtext(value) for value in vector] for vector in covers[label]
                    ],
                    "rank_before": old_rank,
                    "rank_after": new_rank,
                    "minimum_lattice_norm": row["minimum_lattice_norm"],
                    "cover_degree": row["cover_degree"],
                    "equation_complexity": row["equation_complexity"],
                }
            )
        for index in range(dimension):
            if direction_exposure[index] is not None:
                continue
            target = [Q(int(row == index)) for row in range(dimension)]
            if column_rank(accumulated + [target]) == column_rank(accumulated):
                row = representative[label]
                direction_exposure[index] = {
                    "ordered_public_direction": f"Q{index + 1}",
                    "first_exposing_cover": label,
                    "minimum_lattice_norm": row["minimum_lattice_norm"],
                    "cover_degree": row["cover_degree"],
                    "equation_complexity": row["equation_complexity"],
                }
    return {
        "complexity_ordered_rank_increments": rank_steps,
        "minimum_exposure_of_ordered_public_directions": direction_exposure,
        "final_rank": column_rank(accumulated),
    }


def analyze_fibre(
    *, data, controls, bisection_fibre, atlas_by_label, frontier_trisections,
    trisections, quadrisections,
    label, numerator, denominator, known_rank
):
    parameter = qtext(Q(numerator, denominator))
    complements = public_complements(controls)
    specialization = evaluate_projective_specialization(data, numerator, denominator)
    minimal_model, minimal_change, _minimalization = global_minimal_model_with_change(
        specialization.model
    )
    generic_minimal = tuple(
        source_point_to_target(point, minimal_change) for point in specialization.points
    )
    public_points = tuple(PUBLIC_POINTS[parameter])
    complement_indices = complements[parameter]
    public_complement = tuple(public_points[index] for index in complement_indices)
    existing_minimal = generic_minimal + public_complement
    if len(existing_minimal) != known_rank:
        raise ArithmeticError("known basis has the wrong rank")
    short_model, short_change = short_certificate_model(minimal_model)
    existing_short = tuple(
        source_point_to_target(point, short_change) for point in existing_minimal
    )

    baseline = baseline_classes(bisection_fibre, atlas_by_label)
    frontier3 = cover_points_for_parameter(frontier_trisections, parameter, 3, 20)
    deep3 = frontier3 + cover_points_for_parameter(trisections, parameter, 3, 26)
    deep4 = cover_points_for_parameter(quadrisections, parameter, 4, 34)
    deep = deep3 + deep4
    deep_short = []
    for record in deep:
        if not is_on_weierstrass_curve(specialization.model, record["source_point"]):
            raise ArithmeticError("a stored trisection point misses the source fibre")
        minimal_point = source_point_to_target(record["source_point"], minimal_change)
        if not is_on_weierstrass_curve(minimal_model, minimal_point):
            raise ArithmeticError("a trisection point misses the minimal fibre")
        deep_short.append(source_point_to_target(minimal_point, short_change))

    if deep:
        relations = discover_and_verify_relations(
            short_model, existing_short + tuple(deep_short), known_rank
        )
        coordinates = coordinates_from_relations(
            relations["relations"], known_rank, len(deep)
        )
        for record, point, coordinate in zip(deep, deep_short, coordinates):
            common = verify_coordinate_row(short_model, existing_short, point, coordinate)
            record["known_basis_coordinates"] = [qtext(value) for value in coordinate]
            record["exceptional_quotient_coordinates"] = [
                qtext(value) for value in coordinate[17:]
            ]
            record["coordinate_denominator"] = common
            record["exact_relation_verified_by_rational_group_addition"] = True
            record.pop("source_point")
        relation_summary = {
            key: value for key, value in relations.items() if key != "relations"
        }
        relation_summary["relation_vectors_sha256"] = sha256(
            (json.dumps(relations["relations"], separators=(",", ":")) + "\n").encode()
        ).hexdigest()
    else:
        relation_summary = {
            "relation_count": 0,
            "all_relations_verified_by_exact_group_addition": True,
        }

    all_records = baseline + deep
    dimension = known_rank - 17
    baseline_rank = column_rank([vector_from_record(record) for record in baseline])
    degree3_records = baseline + deep3
    degree3_rank = column_rank([vector_from_record(record) for record in degree3_records])
    cumulative_rank = column_rank([vector_from_record(record) for record in all_records])
    structure = matroid_structure(all_records)
    exposure = exposure_profile(all_records, dimension)
    if exposure["final_rank"] != cumulative_rank or structure["rank"] != cumulative_rank:
        raise ArithmeticError("rank summaries disagree")
    return {
        "label": label,
        "parameter": parameter,
        "known_rank_lower_bound": known_rank,
        "exceptional_quotient": {
            "notation": "L_t/M_t",
            "dimension": dimension,
            "ordered_basis": [f"Q{index + 1}" for index in range(dimension)],
            "public_point_indices_one_based": [index + 1 for index in complement_indices],
        },
        "tested_cover_universe": {
            "degree_2_complete_bisection_atlas_size": 39120,
            "degree_3_complete_deep_translation_cosets": 320,
            "degree_3_tested_inversion_representatives": 160,
            "degree_3_frontier_sampled_cosets": 1025,
            "degree_3_frontier_norm_20_tested_inversion_representatives": 69,
            "degree_3_split_cover_count": len({record["cover_label"] for record in deep3}),
            "degree_3_rational_point_count": len(deep3),
            "degree_4_sampled_cosets": 1025,
            "degree_4_norm_34_tested_inversion_representatives": 53,
            "degree_4_split_cover_count": len({record["cover_label"] for record in deep4}),
            "degree_4_rational_point_count": len(deep4),
        },
        "cumulative_captured_exceptional_rank": [
            {"maximum_cover_degree": 2, "R_t(D)": baseline_rank},
            {"maximum_cover_degree": 3, "R_t(D)": degree3_rank},
            {"maximum_cover_degree": 4, "R_t(D)": cumulative_rank},
        ],
        "all_known_directions_explained": cumulative_rank == dimension,
        "degree_3_relation_discovery_and_exact_closeout": relation_summary,
        "points": sorted(all_records, key=complexity_key),
        "duplicate_and_circuit_structure": structure,
        "minimum_norm_and_equation_complexity_exposure": exposure,
    }


def build_payload():
    sys.set_int_max_str_digits(0)
    data = load_q12o5867_data(MODEL, SECTIONS)
    controls = json.loads(CONTROLS.read_text())
    bisection_specializations = json.loads(BISECTION_SPECIALIZATIONS.read_text())
    bisection_atlas = json.loads(BISECTION_ATLAS.read_text())
    trisections = json.loads(TRISECTIONS.read_text())
    frontier_trisections = json.loads(FRONTIER_TRISECTIONS.read_text())
    quadrisections = json.loads(QUADRISECTIONS.read_text())
    atlas_by_label = {record["label"]: record for record in bisection_atlas["bisections"]}
    bisection_by_parameter = {
        fibre["parameter"]: fibre for fibre in bisection_specializations["fibres"]
    }
    fibres = []
    for label, numerator, denominator, known_rank in PARAMETERS:
        parameter = qtext(Q(numerator, denominator))
        print(f"ELKIES2026DEEPCOVER|stage=quotient|parameter={parameter}", flush=True)
        fibres.append(
            analyze_fibre(
                data=data,
                controls=controls,
                bisection_fibre=bisection_by_parameter[parameter],
                atlas_by_label=atlas_by_label,
                trisections=trisections,
                frontier_trisections=frontier_trisections,
                quadrisections=quadrisections,
                label=label,
                numerator=numerator,
                denominator=denominator,
                known_rank=known_rank,
            )
        )
    rank28 = fibres[0]
    ranks = rank28["cumulative_captured_exceptional_rank"]
    stop_full = rank28["all_known_directions_explained"]
    stop_stalled = ranks[-1]["R_t(D)"] == ranks[-2]["R_t(D)"]
    return {
        "schema": "elliptic-curves.elkies-2026-deep-cover-exceptional-quotients.v1",
        "status": "PASS_EXACT_DEEP_COVER_EXCEPTIONAL_QUOTIENT_ANALYSIS",
        "claim": (
            "Exact degree-2 through degree-4 split-cover coordinates in the known "
            "exceptional quotients at the rank-25--28 R17 fibres."
        ),
        "claim_boundary": (
            "Ranks are ranks of displayed classes in the known quotient L_t/M_t. "
            "They neither prove that L_t is saturated nor give an upper bound for "
            "the full specialized Mordell--Weil group."
        ),
        "fibres": fibres,
        "degree_stopping_decision": {
            "rank28_all_eleven_known_directions_explained_by_degree_at_most_4": stop_full,
            "rank28_degree_4_sample_gained_rank": not stop_stalled,
            "degree_4_test_completed": True,
            "degree_5_or_6_test_required": not (stop_full or stop_stalled),
            "rule": (
                "Stop increasing degree once rank 11 is reached at the prioritized "
                "rank-28 fibre; otherwise continue until a tested degree gains no rank."
            ),
        },
        "generation": {
            "command": (
                ".venv/bin/python elliptic-curves/scripts/"
                "analyze_elkies_2026_deep_cover_quotients.py"
            ),
            "checker_sha256": digest(Path(__file__)),
            "inputs": {
                relative(path): digest(path)
                for path in (
                    MODEL, SECTIONS, CONTROLS, BISECTION_SPECIALIZATIONS,
                    BISECTION_ATLAS, FRONTIER_TRISECTIONS, TRISECTIONS,
                    QUADRISECTIONS,
                )
            },
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    payload = build_payload()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if not arguments.output.is_file() or arguments.output.read_text() != rendered:
            raise SystemExit("stale deep-cover quotient certificate")
        terminal = "PASS"
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered)
        terminal = "WROTE"
    print(
        "ELKIES2026DEEPCOVER|ranks={}|status={}".format(
            ",".join(
                str(fibre["cumulative_captured_exceptional_rank"][-1]["R_t(D)"])
                for fibre in payload["fibres"]
            ),
            terminal,
        )
    )


if __name__ == "__main__":
    main()
