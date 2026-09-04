#!/usr/bin/env python3
"""Build marked-U carrier-receptivity profiles from the exact ICARM controls."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"

ATLAS = GENERATED / "elkies-k3-r17-norm12-record-lineage-atlas-v1.json"
SWEEP = GENERATED / "elkies-k3-r17-norm12-icarm-database-sweep-v1.json"
CALIBRATION = GENERATED / "elkies-k3-r17-norm12-icarm-calibration-dataset-v1.json"
NATIVE = GENERATED / "elkies-k3-r17-norm12-native-icarm-quotient-audit-v1.json"
WGXLI = GENERATED / "elkies-k3-r17-norm12-wgxli-lineage-fibres-v1.json"
NORM8_CURVE12 = GENERATED / "elkies-k3-r17-norm12-curve12-norm8-incidence-v1.json"
NORM8_MORE = GENERATED / "elkies-k3-r17-norm12-icarm-norm8-incidence-v1.json"
DIVERSITY = GENERATED / "elkies-k3-r17-multisection-diversity-v1.json"
PAIR_GEOMETRY = GENERATED / "elkies-2026-bisection-pair-cover-geometry-full.json"
PAIR_RANKS = GENERATED / "elkies-2026-immediate-point-pair-rank-ledger.json"
ALTERNATE_LAB = (
    GENERATED
    / "elkies-k3-r17-norm12-11952-alternate-arithmetic-laboratory-cheapest-1024-v1.json"
)
ALTERNATE_BASE_RANKS = GENERATED / "elkies-k3-r17-norm12-11952-v4-base-rank-screen-v1.json"
ALTERNATE_TATE = GENERATED / "elkies-k3-r17-norm12-11952-product-tate-parity-v1.json"
ALTERNATE_CHART_SWEEP = (
    GENERATED / "elkies-k3-r17-norm12-alternate-chart-character-sweep-v1.json"
)
COMPLETE_CHARACTER_103B2 = (
    GENERATED / "elkies-k3-r17-norm12-103b2-complete-character-closure-v1.json"
)
RIGID_074D9 = (
    GENERATED / "elkies-k3-r17-074d9-cross-fibre-bisection-transfer-v1.json"
)
NORM8_074D9 = (
    GENERATED / "elkies-k3-r17-074d9-norm8-cross-fibre-transfer-v1.json"
)
SINGULAR_103B2 = (
    GENERATED / "elkies-k3-r17-norm12-103b2-singular-bisection-search-complete-v1.json"
)
SINGULAR_11952 = (
    GENERATED / "elkies-k3-r17-norm12-11952-singular-bisection-search-complete-v1.json"
)

OUTPUT = GENERATED / "elkies-k3-r17-carrier-receptivity-profiles-v1.json"
TABLE = GENERATED / "elkies-k3-r17-carrier-receptivity-profiles-v1.tsv"


def load(path: Path):
    return json.loads(path.read_text())


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def parameter_complexity(value: str) -> dict:
    q = Fraction(value)
    numerator = int(q.numerator)
    denominator = int(q.denominator)
    return {
        "affine_parameter": str(q),
        "numerator_bits": abs(numerator).bit_length(),
        "denominator_bits": denominator.bit_length(),
        "projective_height_bits": max(abs(numerator).bit_length(), denominator.bit_length()),
        "numerator_decimal_digits": len(str(abs(numerator))),
        "denominator_decimal_digits": len(str(denominator)),
    }


def evaluate_polynomial(coefficients: list[str], value: Fraction) -> Fraction:
    result = Fraction(0)
    for coefficient in reversed(coefficients):
        result = result * value + Fraction(coefficient)
    return result


def rational_square(value: Fraction) -> bool:
    if value < 0:
        return False
    numerator_root = math.isqrt(value.numerator)
    denominator_root = math.isqrt(value.denominator)
    return (
        numerator_root * numerator_root == value.numerator
        and denominator_root * denominator_root == value.denominator
    )


def coefficient_fingerprint(coefficients: list[str]) -> str:
    raw = json.dumps(coefficients, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def cover_audit(row: dict) -> dict:
    keys = [key for key in row if key.endswith("cover_audit")]
    if len(keys) != 1:
        raise ArithmeticError(f"curve {row['curve_id']} has {len(keys)} cover-audit blocks")
    return row[keys[0]]


def unknown_component(reason: str) -> dict:
    return {"status": "UNKNOWN", "reason": reason}


def build():
    atlas = load(ATLAS)
    sweep = load(SWEEP)
    calibration = load(CALIBRATION)
    native = load(NATIVE)
    wgxli = load(WGXLI)
    norm8_curve12 = load(NORM8_CURVE12)
    norm8_more = load(NORM8_MORE)
    diversity = load(DIVERSITY)
    pair_geometry = load(PAIR_GEOMETRY)
    pair_ranks = load(PAIR_RANKS)
    alternate_lab = load(ALTERNATE_LAB)
    alternate_base_ranks = load(ALTERNATE_BASE_RANKS)
    alternate_tate = load(ALTERNATE_TATE)
    alternate_chart_sweep = load(ALTERNATE_CHART_SWEEP)
    complete_character_103b2 = load(COMPLETE_CHARACTER_103B2)
    rigid_074d9 = load(RIGID_074D9)
    norm8_074d9 = load(NORM8_074D9)
    singular_103b2 = load(SINGULAR_103B2)
    singular_11952 = load(SINGULAR_11952)

    if (
        alternate_chart_sweep["status"]
        != "PASS_EXACT_ALTERNATE_CHART_SWEEP_NO_RANK19_OR_RANK20_CHARACTER_INPUT"
    ):
        raise ArithmeticError("the alternate marked-chart character sweep is not exact")
    if (
        complete_character_103b2["status"]
        != "PASS_EXACT_NO_COMPLETE_THREE_CHARACTER_CLOSURE"
        or complete_character_103b2["source_label"] != "norm12-orbit-103b2"
        or int(complete_character_103b2["character_count"]) != 39120
    ):
        raise ArithmeticError("the 103b2 complete branch-character layer is not exact")
    if (
        rigid_074d9["status"]
        != "PASS_EXACT_COMPLETE_074D9_CROSS_FIBRE_BISECTION_TRANSFER"
        or rigid_074d9["summary"]["target_ids"] != [351, 356, 376, 377, 385]
        or not rigid_074d9["cross_fibre_transfer"][
            "all_split_points_have_exact_displayed_quotient_coordinates"
        ]
    ):
        raise ArithmeticError("the 074d9 complete rigid-transfer layer is not exact")
    if (
        norm8_074d9["status"] != "PASS_EXACT_CANONICAL_NORM8_CROSS_FIBRE_TRANSFER"
        or not norm8_074d9["canonical_trace"]["selection_is_target_independent"]
        or int(norm8_074d9["canonical_trace"]["complete_minimum_norm_eight_class_count"])
        != 63925
        or norm8_074d9["different_family_control"]["status"]
        != "NOT_APPLICABLE_DIFFERENT_FIBRATION"
    ):
        raise ArithmeticError("the 074d9 norm-eight cross-fibre layer is not exact")

    charts = atlas["atlas"]["charts"]
    if len(charts) != 43:
        raise ArithmeticError("the exact lineage atlas no longer has 43 marked charts")
    chart_by_label = {row["label"]: row for row in charts}
    if len(chart_by_label) != 43:
        raise ArithmeticError("marked chart labels are not unique")

    family_by_chart = {}
    for family in atlas["atlas"]["pgl2_equivalence_classes"]:
        representative = family["representative"]
        for member in family["members"]:
            label = member["label"]
            if label in family_by_chart:
                raise ArithmeticError(f"duplicate PGL2 family assignment for {label}")
            family_by_chart[label] = representative
    if set(family_by_chart) != set(chart_by_label):
        raise ArithmeticError("PGL2 families do not partition the 43 marked charts")

    exact_rows = {
        int(row["curve_id"]): row
        for row in calibration["rows"]
        if row["displayed_exceptional_quotient_dimension"] is not None
    }
    unknown_rows = [
        row
        for row in calibration["rows"]
        if row["displayed_exceptional_quotient_dimension"] is None
    ]
    if len(exact_rows) != 12 or len(unknown_rows) != 57:
        raise ArithmeticError("expected twelve exact and 57 unknown quotient-labelled rows")

    native_by_curve = {int(row["curve_id"]): row for row in native["fibres"]}
    old_by_curve = {int(row["curve_id"]): row for row in wgxli["exceptional_quotients"]}
    if set(native_by_curve) & set(old_by_curve):
        raise ArithmeticError("native and prior exact quotient audits overlap")
    if set(native_by_curve) | set(old_by_curve) != set(exact_rows):
        raise ArithmeticError("exact quotient sources do not match the twelve calibration rows")

    old_chart = wgxli["representative"]["chart"]
    exact_chart_by_curve = {
        **{curve_id: row["native_chart"] for curve_id, row in native_by_curve.items()},
        **{curve_id: old_chart for curve_id in old_by_curve},
    }
    exact_targets_by_chart = {label: [] for label in chart_by_label}
    for curve_id, label in exact_chart_by_curve.items():
        exact_targets_by_chart[label].append(curve_id)
    for curve_ids in exact_targets_by_chart.values():
        curve_ids.sort()

    sweep_by_curve = {
        int(row["curve_id"]): row for row in sweep["rational_j_hits_and_twists"]
    }
    missing_sweep = set(exact_rows) - set(sweep_by_curve)
    if missing_sweep:
        raise ArithmeticError(f"exact curves missing from public atlas: {sorted(missing_sweep)}")
    specialization_by_chart = {label: [] for label in chart_by_label}
    specialization_lookup = {}
    for curve_id, exact_row in exact_rows.items():
        hit = sweep_by_curve[curve_id]
        for native_twist in hit["native_chart_twists"]:
            chart = native_twist["chart"]
            parameter = native_twist["native_parameter"]["affine_parameter"]
            marked_exact = exact_chart_by_curve[curve_id] == chart
            record = {
                "curve_id": curve_id,
                "snapshot_rank_lower_bound": int(exact_row["snapshot_rank_lower_bound"]),
                "known_curve_displayed_quotient_dimension": int(
                    exact_row["displayed_exceptional_quotient_dimension"]
                ),
                "marked_quotient_transport_status": (
                    "EXACT_CHART_SPECIFIC_SATURATED_TRANSPORT"
                    if marked_exact
                    else "UNKNOWN_NOT_TRANSPORTED_TO_THIS_MARKED_U"
                ),
                "marked_displayed_quotient_dimension": (
                    int(exact_row["displayed_exceptional_quotient_dimension"])
                    if marked_exact
                    else None
                ),
                "parameter_complexity": parameter_complexity(parameter),
                "QQ_isomorphism_status": native_twist["twist"]["status"],
            }
            specialization_by_chart[chart].append(record)
            specialization_lookup[(chart, curve_id)] = record
    for rows in specialization_by_chart.values():
        rows.sort(key=lambda row: row["curve_id"])

    fitted_sources = []
    fitted_sources.append(
        {
            "curve_id": int(norm8_curve12["curve"]["icarm_id"]),
            "native_chart": norm8_curve12["curve"]["native_chart"],
            "native_parameter": norm8_curve12["curve"]["native_parameter"],
            "directions": norm8_curve12["incidence_signature"]["directions"],
        }
    )
    for row in norm8_more["fibres"]:
        fitted_sources.append(
            {
                "curve_id": int(row["curve_id"]),
                "native_chart": row["native_chart"],
                "native_parameter": row["native_parameter"],
                "directions": row["incidence_signature"]["directions"],
            }
        )
    fitted_sources.sort(key=lambda row: row["curve_id"])
    if {row["curve_id"] for row in fitted_sources} != {12, 363, 364, 378, 395}:
        raise ArithmeticError("unexpected norm-eight fitted source set")

    fitted_by_chart = {label: [] for label in chart_by_label}
    fitted_graph_edges_by_chart = {label: [] for label in chart_by_label}
    fitted_matrix_by_chart = {label: [] for label in chart_by_label}
    for source in fitted_sources:
        chart = source["native_chart"]
        fitted_by_chart[chart].append(source)
        source_curve_id = source["curve_id"]
        for target_curve_id in exact_targets_by_chart[chart]:
            target_parameter = Fraction(
                specialization_lookup[(chart, target_curve_id)]["parameter_complexity"]
                ["affine_parameter"]
            )
            hit_directions = []
            for direction in source["directions"]:
                value = evaluate_polynomial(
                    direction["branch_quartic_coefficients_low_to_high"], target_parameter
                )
                if rational_square(value):
                    direction_label = direction["quotient_basis_direction"]
                    hit_directions.append(direction_label)
                    fitted_graph_edges_by_chart[chart].append(
                        {
                            "carrier": (
                                f"norm8-fitted:curve-{source_curve_id}:"
                                f"{direction_label}"
                            ),
                            "source_curve_id": source_curve_id,
                            "source_quotient_basis_direction": direction_label,
                            "target_curve_id": target_curve_id,
                            "diagonal": source_curve_id == target_curve_id,
                            "branch_quartic_sha256": coefficient_fingerprint(
                                direction["branch_quartic_coefficients_low_to_high"]
                            ),
                        }
                    )
            fitted_matrix_by_chart[chart].append(
                {
                    "source_curve_id": source_curve_id,
                    "target_curve_id": target_curve_id,
                    "source_direction_count": len(source["directions"]),
                    "split_source_direction_rank": len(hit_directions),
                    "split_source_directions": hit_directions,
                    "diagonal": source_curve_id == target_curve_id,
                }
            )

    for source in fitted_sources:
        diagonal = [
            row
            for row in fitted_matrix_by_chart[source["native_chart"]]
            if row["source_curve_id"] == source["curve_id"]
            and row["target_curve_id"] == source["curve_id"]
        ]
        if len(diagonal) != 1 or diagonal[0]["split_source_direction_rank"] != len(
            source["directions"]
        ):
            raise ArithmeticError(f"curve {source['curve_id']} lost its fitted diagonal")

    for source_experiment in norm8_074d9["source_experiments"]:
        chart = "norm12-orbit-074d9"
        source_curve_id = int(source_experiment["source_curve_id"])
        covers = source_experiment["covers"]
        if source_curve_id not in (356, 385) or len(covers) != 12:
            raise ArithmeticError("the 074d9 norm-eight source set changed")
        fitted_by_chart[chart].append(
            {
                "curve_id": source_curve_id,
                "native_chart": chart,
                "directions": covers,
            }
        )
        result_by_target = {
            int(row["curve_id"]): row
            for row in source_experiment["quotient_transfer_results"]
        }
        if set(result_by_target) != set(exact_targets_by_chart[chart]):
            raise ArithmeticError("the 074d9 norm-eight target set changed")
        for target_curve_id in exact_targets_by_chart[chart]:
            result = result_by_target[target_curve_id]
            records_by_cover = {
                row["cover_id"]: row for row in result.get("records", [])
            }
            split_directions = []
            for cover in covers:
                evaluation = next(
                    row
                    for row in cover["frozen_evaluations"]
                    if int(row["curve_id"]) == target_curve_id
                )
                if evaluation["status"] != "SPLIT_OVER_Q":
                    continue
                direction = cover["source_quotient_basis_direction"]
                split_directions.append(direction)
                exact_class = records_by_cover[cover["cover_id"]][
                    "exact_displayed_free_quotient_class"
                ]
                if exact_class["status"] != "PROVED_EXACT_IN_DISPLAYED_SUBGROUP":
                    raise ArithmeticError("a 074d9 norm-eight split lacks exact coordinates")
                fitted_graph_edges_by_chart[chart].append(
                    {
                        "carrier": f"norm8-fitted:curve-{source_curve_id}:{direction}",
                        "source_curve_id": source_curve_id,
                        "source_quotient_basis_direction": direction,
                        "target_curve_id": target_curve_id,
                        "diagonal": source_curve_id == target_curve_id,
                        "branch_quartic_sha256": coefficient_fingerprint(
                            cover["branch_quartic_coefficients_low_to_high"]
                        ),
                        "target_quotient_class": {
                            "basis": exact_class["basis"],
                            "coordinates": exact_class["coordinates"],
                        },
                    }
                )
            fitted_matrix_by_chart[chart].append(
                {
                    "source_curve_id": source_curve_id,
                    "target_curve_id": target_curve_id,
                    "source_direction_count": len(covers),
                    "split_source_direction_rank": int(result["class_span_rank"]),
                    "split_source_directions": split_directions,
                    "diagonal": source_curve_id == target_curve_id,
                }
            )
            if len(split_directions) != int(result["split_count"]):
                raise ArithmeticError("a 074d9 norm-eight split count changed")

    rigid_targets_by_chart = {label: [] for label in chart_by_label}
    rigid_edges_by_chart = {label: [] for label in chart_by_label}
    for curve_id, row in native_by_curve.items():
        chart = row["native_chart"]
        audit = cover_audit(row)
        rigid_targets_by_chart[chart].append(
            {
                "curve_id": curve_id,
                "inventory_size": int(audit["covers_evaluated"]),
                "rational_split_count": int(audit["rational_split_count"]),
                "transfer_rank_in_displayed_exceptional_quotient": int(
                    audit["exact_split_span_rank_in_exceptional_quotient"]
                ),
                "transfer_span_primitive": bool(audit["split_span_is_primitive"]),
            }
        )
        for split in audit["splits"]:
            rigid_edges_by_chart[chart].append(
                {
                    "carrier": f"rigid:{split['label']}",
                    "target_curve_id": curve_id,
                    "priority_rank": int(split["priority_rank"]),
                    "target_quotient_class": {
                        "basis": "deterministic Smith quotient basis",
                        "coordinates": split[
                            "plus_quotient_vector_in_deterministic_smith_basis"
                        ],
                    },
                }
            )
    for fibre in rigid_074d9["fibres"]:
        curve_id = int(fibre["curve_id"])
        if curve_id not in exact_targets_by_chart["norm12-orbit-074d9"]:
            raise ArithmeticError(f"unexpected 074d9 rigid target {curve_id}")
        rigid_targets_by_chart["norm12-orbit-074d9"].append(
            {
                "curve_id": curve_id,
                "inventory_size": int(rigid_074d9["lattice_transport"]["record_count"]),
                "rational_split_count": int(fibre["split_count"]),
                "transfer_rank_in_displayed_exceptional_quotient": int(
                    fibre["class_span_rank"]
                ),
                "transfer_span_primitive": None,
                "transfer_span_primitivity_status": (
                    "UNKNOWN_NOT_COMPUTED_BY_THE_074D9_CERTIFICATE"
                ),
            }
        )
        for split in fibre["records"]:
            exact_class = split["exact_displayed_free_quotient_class"]
            if exact_class["status"] != "PROVED_EXACT_IN_DISPLAYED_SUBGROUP":
                raise ArithmeticError("a 074d9 split lacks an exact displayed quotient class")
            rigid_edges_by_chart["norm12-orbit-074d9"].append(
                {
                    "carrier": f"rigid:{split['label']}",
                    "target_curve_id": curve_id,
                    "priority_rank": int(split["priority_rank"]),
                    "target_quotient_class": {
                        "basis": exact_class["basis"],
                        "coordinates": exact_class["coordinates"],
                    },
                }
            )
    for rows in rigid_targets_by_chart.values():
        rows.sort(key=lambda row: row["curve_id"])
    for chart, edges in rigid_edges_by_chart.items():
        carrier_targets = {}
        for edge in edges:
            carrier_targets.setdefault(edge["carrier"], set()).add(edge["target_curve_id"])
        shared = {carrier: targets for carrier, targets in carrier_targets.items() if len(targets) > 1}
        if shared:
            raise ArithmeticError(f"rigid carrier shared across exact targets in {chart}: {shared}")

    r17_histogram = {
        str(key): int(value)
        for key, value in diversity["degree_two"]["minimum_norm_histogram"].items()
    }
    alternate_histogram = {
        str(key): int(value)
        for key, value in alternate_tate["invariant_trace_parity"][
            "minimum_norm_histogram"
        ].items()
    }
    if sum(r17_histogram.values()) != 2**17 or sum(alternate_histogram.values()) != 2**17:
        raise ArithmeticError("minimum-norm histogram does not cover M/2M")
    if int(singular_103b2["minimum_translation_class_count"]) != r17_histogram["8"]:
        raise ArithmeticError("103b2 marked norm-eight count disagrees with the frame census")
    if int(singular_11952["minimum_translation_class_count"]) != alternate_histogram["8"]:
        raise ArithmeticError("11952 marked norm-eight count disagrees with the frame census")

    inherited_character_by_chart = {
        row["label"]: row for row in alternate_chart_sweep["inherited_layer"]["records"]
    }
    alternate_complete_character_by_chart = {
        row["label"]: row
        for row in alternate_chart_sweep["complete_smooth_layer"]["records"]
    }
    obstruction_by_chart = {
        row["label"]: row
        for row in alternate_chart_sweep["remaining_full_atlas_marking_obstructions"]
    }
    alternate_charts = {
        row["label"] for row in charts if row["frame_class"] == "alternate-Q80"
    }
    if set(inherited_character_by_chart) != alternate_charts:
        raise ArithmeticError("the inherited character layer does not cover ten alternate charts")
    if set(alternate_complete_character_by_chart) | set(obstruction_by_chart) != alternate_charts:
        raise ArithmeticError("complete and obstructed alternate character layers do not partition")
    if set(alternate_complete_character_by_chart) & set(obstruction_by_chart):
        raise ArithmeticError("an alternate chart is both complete and marking-obstructed")
    complete_character_by_chart = dict(alternate_complete_character_by_chart)
    complete_character_by_chart["norm12-orbit-103b2"] = {
        "label": "norm12-orbit-103b2",
        "cover_count": int(complete_character_103b2["character_count"]),
        "distinct_character_count": int(
            complete_character_103b2["distinct_character_count"]
        ),
        "collision_count": (
            int(complete_character_103b2["character_count"])
            - int(complete_character_103b2["distinct_character_count"])
        ),
        "internal_relation_count": int(
            complete_character_103b2["three_character_relation_count"]
        ),
        "committed_catalog_match_count": int(
            complete_character_103b2["older_published_base_catalog"]
            ["formal_product_match_count"]
        ),
        "evidence": relative(COMPLETE_CHARACTER_103B2),
    }

    pair_rank_summary = pair_ranks["summary"]
    exact_rank_one_pairs = [
        row
        for row in alternate_base_ranks["results"]
        if row.get("rank_lower_bound") == 1 and row.get("rank_upper_bound") == 1
    ]
    if len(exact_rank_one_pairs) != 17:
        raise ArithmeticError("alternate-Q80 exact rank-one base count changed")

    low_genus_by_chart = {
        "norm12-orbit-074d9": {
            "status": "EXACT_COMPLETE_RIGID_CROSS_FIBRE_ATLAS",
            "individual_quadratic_bases": {
                "base_curve_genus": 0,
                "complete_carrier_count": int(
                    rigid_074d9["lattice_transport"]["record_count"]
                ),
                "rational_point_incidence_count_on_current_targets": sum(
                    int(row["split_count"]) for row in rigid_074d9["fibres"]
                ),
                "incidences_by_target": [
                    {
                        "curve_id": int(row["curve_id"]),
                        "rational_carrier_count": int(row["split_count"]),
                        "surface_generic_MW_rank_lower_bound": 18,
                    }
                    for row in rigid_074d9["fibres"]
                ],
            },
            "rank_29_record_corridor": {
                "target_curve_ids": [356, 385],
                "common_rational_carrier_count": int(
                    rigid_074d9["cross_fibre_transfer"][
                        "record_anchored_rank_at_least_18_cover_count"
                    ]
                ),
                "rank_at_least_19_two_character_compositum_available": bool(
                    rigid_074d9["cross_fibre_transfer"][
                        "rank_at_least_19_two_character_compositum_available"
                    ]
                ),
            },
            "boundary": rigid_074d9["claim_boundary"],
        },
        "norm12-orbit-103b2": {
            "status": "EXACT_COMPLETE_NATIVE_RIGID_ATLAS",
            "individual_quadratic_bases": {
                "base_curve_genus": 0,
                "cover_count": int(
                    pair_geometry["complete_conic_classification"]["record_count"]
                ),
                "all_Q_rational": bool(
                    pair_geometry["complete_conic_classification"]["all_conics_Q_rational"]
                ),
                "surface_generic_MW_rank_lower_bound": 18,
            },
            "paired_V4_bases": {
                "base_curve_genus": int(
                    pair_geometry["all_distinct_pairs"]["fiber_product_genus"]
                ),
                "pair_count": int(pair_geometry["all_distinct_pairs"]["pair_count"]),
                "surface_generic_MW_rank_lower_bound": int(
                    pair_geometry["all_distinct_pairs"]["generic_mw_rank_lower_bound"]
                ),
                "immediate_Q_point_pair_count": int(
                    pair_geometry["all_distinct_pairs"]["immediate_Q_point_subfamily"]
                    ["distinct_pair_count"]
                ),
                "base_Jacobian_rank_lower_bound_distribution": {
                    str(key): int(value)
                    for key, value in pair_rank_summary[
                        "certified_rank_lower_bound_counts"
                    ].items()
                },
                "maximum_certified_base_Jacobian_rank_lower_bound": int(
                    pair_rank_summary["maximum_certified_rank_lower_bound"]
                ),
            },
            "boundary": (
                "Base-Jacobian values are certified lower bounds; zero means no point "
                "was found in the declared pass, not rank zero."
            ),
        },
        "norm12-orbit-11952": {
            "status": "EXACT_BOUNDED_1143_CLASS_LABORATORY",
            "individual_quadratic_bases": {
                "base_curve_genus": 0,
                "cover_count": int(
                    alternate_lab["individual_cover_arithmetic"]["cover_count"]
                ),
                "all_Q_rational": bool(
                    alternate_lab["individual_cover_arithmetic"]["all_conics_Q_rational"]
                ),
                "surface_generic_MW_rank_lower_bound": int(
                    alternate_lab["individual_cover_arithmetic"]
                    ["generic_rank_lower_bound_on_each_Q_rational_cover"]
                ),
            },
            "paired_V4_bases": {
                "base_curve_genus": int(
                    alternate_lab["pair_and_product_search"]
                    ["minimum_connected_pair_fibre_product_genus"]
                ),
                "pair_count": int(
                    alternate_lab["pair_and_product_search"]
                    ["disjoint_branch_genus_one_pair_count"]
                ),
                "surface_generic_MW_rank_lower_bound": int(
                    alternate_lab["pair_and_product_search"]
                    ["pair_cover_generic_rank_lower_bound"]
                ),
                "bounded_base_rank_screen_size": int(
                    alternate_base_ranks["limits"]["shortlist_prefix"]
                ),
                "completed_base_rank_screens": int(
                    alternate_base_ranks["summary"]["completed"]
                ),
                "exact_rank_one_base_count": len(exact_rank_one_pairs),
                "maximum_certified_base_Jacobian_rank_lower_bound_in_screen": 1,
            },
            "genus_five_three_character_bases": {
                "base_curve_genus": int(
                    alternate_lab["independent_character_triples"]["common_cover_genus"]
                ),
                "count": int(alternate_lab["independent_character_triples"]["count"]),
                "surface_generic_MW_rank_lower_bound": int(
                    alternate_lab["independent_character_triples"]
                    ["generic_rank_lower_bound"]
                ),
            },
            "boundary": alternate_lab["proof_boundary"],
        },
    }

    frame_histogram = {
        "published-R17": r17_histogram,
        "alternate-Q80": alternate_histogram,
    }
    marked_minimum_status = {
        "norm12-orbit-074d9": {
            "status": "EXACT_MARKED_NORM_EIGHT_AND_NORM_TEN_COUNTS_ONLY",
            "minimum_norm_histogram": None,
            "marked_norm_eight_class_count": int(
                norm8_074d9["canonical_trace"]["complete_minimum_norm_eight_class_count"]
            ),
            "marked_rational_norm_ten_class_count": int(
                rigid_074d9["lattice_transport"]["record_count"]
            ),
            "marked_transport_evidence": [
                relative(NORM8_074D9),
                relative(RIGID_074D9),
            ],
        },
        "norm12-orbit-103b2": {
            "status": "EXACT_MARKED_MINIMUM_CLASSES",
            "minimum_norm_histogram": r17_histogram,
            "marked_norm_eight_class_count": r17_histogram["8"],
            "marked_rational_norm_ten_class_count": r17_histogram["10"],
            "marked_transport_evidence": relative(SINGULAR_103B2),
        },
        "norm12-orbit-11952": {
            "status": "EXACT_MARKED_MINIMUM_CLASSES_AND_ZERO_CLASS_REDUCTION",
            "minimum_norm_histogram": alternate_histogram,
            "marked_norm_eight_class_count": alternate_histogram["8"],
            "marked_rational_norm_ten_class_count": alternate_histogram["10"],
            "marked_transport_evidence": relative(SINGULAR_11952),
            "deep_norm_twelve_trace_classes": int(
                alternate_tate["invariant_trace_parity"]["deep_norm12_class_count"]
            ),
            "product_targets_with_zero_Tate_class_excluded_under_stated_gates": int(
                alternate_tate["interaction_with_completed_inversion"][
                    "product_target_count"
                ]
            ),
            "residual_norm_twelve_target_trace_cases": int(
                alternate_tate["interaction_with_completed_inversion"][
                    "residual_trace_target_pairs"
                ]
            ),
        },
        "norm12-orbit-08f72": {
            "status": "EXACT_MARKED_NORM_EIGHT_AND_NORM_TEN_COUNTS_ONLY",
            "minimum_norm_histogram": None,
            "marked_norm_eight_class_count": 63917,
            "marked_rational_norm_ten_class_count": 39147,
            "marked_transport_evidence": [relative(NORM8_MORE), relative(NATIVE)],
        },
    }
    for label, character_record in alternate_complete_character_by_chart.items():
        if label in marked_minimum_status:
            continue
        if int(character_record["cover_count"]) != alternate_histogram["10"]:
            raise ArithmeticError(f"{label} complete norm-ten count disagrees with Q80")
        marked_minimum_status[label] = {
            "status": "EXACT_MARKED_NORM_TEN_COUNT_ONLY",
            "minimum_norm_histogram": None,
            "marked_norm_eight_class_count": None,
            "marked_rational_norm_ten_class_count": int(character_record["cover_count"]),
            "marked_transport_evidence": relative(ALTERNATE_CHART_SWEEP),
        }

    profiles = []
    for chart in charts:
        label = chart["label"]
        frame_class = chart["frame_class"]
        exact_target_ids = exact_targets_by_chart[label]
        rigid_targets = rigid_targets_by_chart[label]
        rigid_edges = rigid_edges_by_chart[label]
        fitted_edges = fitted_graph_edges_by_chart[label]
        fitted_matrix = fitted_matrix_by_chart[label]
        off_diagonal = [row for row in fitted_matrix if not row["diagonal"]]
        diagonal = [row for row in fitted_matrix if row["diagonal"]]
        inherited_character = inherited_character_by_chart.get(label)
        complete_character = complete_character_by_chart.get(label)
        marking_obstruction = obstruction_by_chart.get(label)

        if rigid_targets:
            rigid_component = {
                "status": "EXACT_COMPLETE_FIXED_NATIVE_INVENTORY_ON_LABELLED_TARGETS",
                "targets": rigid_targets,
                "rank_vector_in_curve_id_order": [
                    {
                        "curve_id": row["curve_id"],
                        "rank": row[
                            "transfer_rank_in_displayed_exceptional_quotient"
                        ],
                    }
                    for row in rigid_targets
                ],
                "boundary": (
                    "Ranks are spans inside the displayed exceptional quotients, not "
                    "upper bounds for the full specialized Mordell-Weil groups."
                ),
            }
        elif exact_target_ids:
            rigid_component = unknown_component(
                "Exact quotient labels exist, but the fixed native carrier inventory "
                "has not been transported and audited on this marked U."
            )
            rigid_component["target_curve_ids"] = exact_target_ids
        else:
            rigid_component = unknown_component(
                "No chart-specific exact quotient target is currently available."
            )

        if off_diagonal:
            off_diagonal_component = {
                "status": "EXACT_FOR_CURRENT_FITTED_NORM_EIGHT_CARRIERS_AND_TARGETS",
                "rank_definition": (
                    "For a source fibre, count the independent source quotient-basis "
                    "directions whose fitted branch quartics split at the target. A "
                    "zero split count proves transfer rank zero; positive values are "
                    "source-carrier ranks and are not target quotient ranks without a "
                    "separate target point transport."
                ),
                "matrix": fitted_matrix,
                "diagonal_fitted_rank_vector": [
                    {
                        "curve_id": row["source_curve_id"],
                        "rank": row["split_source_direction_rank"],
                    }
                    for row in diagonal
                ],
                "off_diagonal_direction_tests": sum(
                    row["source_direction_count"] for row in off_diagonal
                ),
                "off_diagonal_split_direction_count": sum(
                    row["split_source_direction_rank"] for row in off_diagonal
                ),
                "all_off_diagonal_transfer_ranks_zero": all(
                    row["split_source_direction_rank"] == 0 for row in off_diagonal
                ),
                "boundary": (
                    "The source pencils were fitted using their diagonal exceptional "
                    "points. This is an exact cross-incidence audit, not a held-out "
                    "discovery experiment or an exhaustion of the norm-eight frontier."
                ),
            }
            if label == "norm12-orbit-074d9":
                off_diagonal_component["canonical_trace_selection"] = {
                    "selection_is_target_independent": bool(
                        norm8_074d9["canonical_trace"]["selection_is_target_independent"]
                    ),
                    "selected_priority_rank": int(
                        norm8_074d9["canonical_trace"]["selected_priority_rank"]
                    ),
                    "complete_marked_norm_eight_class_count": int(
                        norm8_074d9["canonical_trace"]
                        ["complete_minimum_norm_eight_class_count"]
                    ),
                }
                off_diagonal_component["different_fibration_control"] = norm8_074d9[
                    "different_family_control"
                ]
                off_diagonal_component["evidence"] = relative(NORM8_074D9)
        elif fitted_by_chart[label]:
            off_diagonal_component = unknown_component(
                "Only one fitted source/target is present, so there is no off-diagonal pair."
            )
        else:
            off_diagonal_component = unknown_component(
                "No target-fitted norm-eight carriers have been compiled on this marked U."
            )

        graph_component = {
            "status": (
                "EXACT_FOR_STORED_RIGID_AND_FITTED_CARRIER_LAYERS"
                if rigid_edges or fitted_edges
                else (
                    "EXACT_FOR_STORED_NATIVE_CHARACTER_LAYER"
                    if inherited_character is not None
                    else "UNKNOWN_NO_STORED_CARRIER_INCIDENCE_LAYER"
                )
            ),
            "target_curve_ids": exact_target_ids,
            "rigid_cover_layer": {
                "carrier_nodes_with_incidence": len({edge["carrier"] for edge in rigid_edges}),
                "edges": rigid_edges,
                "edge_count": len(rigid_edges),
                "carriers_incident_to_multiple_current_targets": 0,
            },
            "fitted_norm_eight_layer": {
                "carrier_nodes": len(
                    {
                        edge["carrier"]
                        for edge in fitted_graph_edges_by_chart[label]
                    }
                ),
                "edges": fitted_edges,
                "edge_count": len(fitted_edges),
                "diagonal_edge_count": sum(edge["diagonal"] for edge in fitted_edges),
                "off_diagonal_edge_count": sum(
                    not edge["diagonal"] for edge in fitted_edges
                ),
            },
            "native_branch_character_layer": {
                "status": (
                    "EXACT_COMPLETE_SMOOTH_LAYER"
                    if complete_character is not None
                    else (
                        "EXACT_INHERITED_LAYER_FULL_SMOOTH_LAYER_MARKING_OBSTRUCTED"
                        if marking_obstruction is not None
                        else (
                            "EXACT_INHERITED_LAYER_ONLY"
                            if inherited_character is not None
                            else "UNKNOWN"
                        )
                    )
                ),
                "inherited_height_four": inherited_character,
                "complete_smooth_norm_ten": complete_character,
                "full_smooth_layer_marking_obstruction": marking_obstruction,
                "evidence": (
                    relative(ALTERNATE_CHART_SWEEP)
                    if inherited_character is not None
                    else None
                ),
            },
            "boundary": (
                "This bipartite graph has carrier nodes and exact quotient-labelled "
                "specialization nodes. It is scoped to the stored fixed inventories "
                "and 75 fitted norm-eight carriers. The separate native-character "
                "layer records within-chart squareclass collisions and relations, not "
                "specialization edges. Missing edges outside those scopes are UNKNOWN, "
                "not absent."
            ),
        }

        minimum_data = marked_minimum_status.get(label)
        tate_component = {
            "status": (
                minimum_data["status"]
                if minimum_data is not None
                else "FRAME_REFERENCE_ONLY_MARKED_MINIMUM_TRANSPORT_UNKNOWN"
            ),
            "abstract_frame_minimum_norm_histogram": frame_histogram[frame_class],
            "marked_minimum_classes": minimum_data,
            "Tate_quotient": {
                "status": "UNKNOWN",
                "reason": (
                    "No full anti-invariant Mordell-Weil lattice and integral "
                    "base-change glue have been computed for this marked-U target."
                ),
            },
            "boundary": (
                "Frame histograms are exact lattice references but do not become "
                "chart-specific carrier data without the displayed marked transport."
            ),
        }
        if label == "norm12-orbit-11952":
            tate_component["Tate_quotient"]["zero_class_result"] = (
                "For the seventeen separately selected product characters, the zero "
                "class is excluded for height-eight sections under the direct-polynomial "
                "and local-component gates; quotient dimensions and nonzero classes remain UNKNOWN."
            )

        specializations = specialization_by_chart[label]
        exact_marked_specializations = [
            row
            for row in specializations
            if row["marked_quotient_transport_status"]
            == "EXACT_CHART_SPECIFIC_SATURATED_TRANSPORT"
        ]
        parameter_bits = [
            row["parameter_complexity"]["projective_height_bits"]
            for row in specializations
        ]
        equation_parameter_component = {
            "status": "EXACT_EQUATION_COMPLEXITY_AND_PUBLIC_SPECIALIZATION_PARAMETERS",
            "equation_complexity": chart["equation_complexity"],
            "current_exact_curve_matches_in_this_PGL2_chart": specializations,
            "chart_specific_exact_quotient_transport_count": len(
                exact_marked_specializations
            ),
            "matched_parameter_projective_height_bits_range": (
                [min(parameter_bits), max(parameter_bits)] if parameter_bits else None
            ),
            "boundary": (
                "A QQ-isomorphic public specialization supplies an exact parameter, "
                "but its quotient and carrier ranks remain UNKNOWN on this marked U "
                "until the generic saturated section transport is compiled."
            ),
        }

        profiles.append(
            {
                "marked_U": {
                    "label": label,
                    "identity": "primitive ordered pair <F,O> in the pinned NS marking",
                    "trace_vector_in_pinned_rank17_basis": chart["trace_vector"],
                    "frame_class": frame_class,
                    "PGL2_family_representative": family_by_chart[label],
                },
                "exact_quotient_labelled_curve_ids": exact_target_ids,
                "carrier_receptivity_profile": {
                    "rigid_cover_transfer_ranks": rigid_component,
                    "off_diagonal_norm_eight_transfer_ranks": off_diagonal_component,
                    "branch_incidence_graph": graph_component,
                    "low_genus_base_ranks": low_genus_by_chart.get(
                        label,
                        unknown_component(
                            "No chart-specific low-genus carrier-base rank audit is stored."
                        ),
                    ),
                    "Tate_quotient_and_minimum_classes": tate_component,
                    "equation_and_parameter_complexity": equation_parameter_component,
                },
            }
        )

    profiles.sort(
        key=lambda row: row["carrier_receptivity_profile"]
        ["equation_and_parameter_complexity"]["equation_complexity"][
            "deep_equation_rank"
        ]
    )

    highest_unknown_rank = max(int(row["snapshot_rank_lower_bound"]) for row in unknown_rows)
    selected_unknown_rows = sorted(
        [
            row
            for row in unknown_rows
            if int(row["snapshot_rank_lower_bound"]) == highest_unknown_rank
        ],
        key=lambda row: int(row["curve_id"]),
    )
    if highest_unknown_rank != 28 or [int(row["curve_id"]) for row in selected_unknown_rows] != [
        11,
        391,
        423,
    ]:
        raise ArithmeticError("the declared highest-rank transport tranche changed")

    transport_tranche = []
    for row in selected_unknown_rows:
        curve_id = int(row["curve_id"])
        marked_u = row["family"]
        hit = sweep_by_curve[curve_id]
        choices = [entry for entry in hit["native_chart_twists"] if entry["chart"] == marked_u]
        if len(choices) != 1:
            raise ArithmeticError(f"curve {curve_id} has no unique representative-chart parameter")
        parameter = choices[0]["native_parameter"]["affine_parameter"]
        transport_tranche.append(
            {
                "curve_id": curve_id,
                "marked_U": marked_u,
                "PGL2_family_representative": hit["representative"],
                "snapshot_rank_lower_bound": int(row["snapshot_rank_lower_bound"]),
                "rank_jump_lower_bound_over_generic_17": int(
                    row["rank_jump_lower_bound_over_generic_17"]
                ),
                "transport_status": "UNKNOWN_SELECTED_NOT_YET_COMPILED",
                "parameter_complexity": parameter_complexity(parameter),
                "marked_U_equation_complexity": chart_by_label[marked_u][
                    "equation_complexity"
                ],
                "required_completion": (
                    "compile the saturated chart-specific generic section transport, "
                    "then certify the displayed quotient and fixed native carrier span"
                ),
            }
        )

    selected_by_marked_u = {label: [] for label in chart_by_label}
    for candidate in transport_tranche:
        selected_by_marked_u[candidate["marked_U"]].append(candidate)
    for profile in profiles:
        label = profile["marked_U"]["label"]
        complexity = profile["carrier_receptivity_profile"][
            "equation_and_parameter_complexity"
        ]
        complexity["selected_initial_transport_candidates"] = selected_by_marked_u[
            label
        ]
        complexity["selected_initial_transport_count"] = len(
            selected_by_marked_u[label]
        )

    table_rows = []
    for profile in profiles:
        marked = profile["marked_U"]
        receptivity = profile["carrier_receptivity_profile"]
        rigid = receptivity["rigid_cover_transfer_ranks"]
        offdiag = receptivity["off_diagonal_norm_eight_transfer_ranks"]
        graph = receptivity["branch_incidence_graph"]
        native_characters = graph["native_branch_character_layer"]
        low_genus = receptivity["low_genus_base_ranks"]
        tate = receptivity["Tate_quotient_and_minimum_classes"]
        complexity = receptivity["equation_and_parameter_complexity"]
        equation = complexity["equation_complexity"]
        marked_minima = tate["marked_minimum_classes"] or {}
        table_rows.append(
            {
                "marked_U": marked["label"],
                "frame_class": marked["frame_class"],
                "PGL2_family": marked["PGL2_family_representative"],
                "deep_equation_rank": equation["deep_equation_rank"],
                "equation_support_count": equation["support_count"],
                "equation_coordinate_input_bits": equation["coordinate_input_bits"],
                "exact_quotient_target_count": len(
                    profile["exact_quotient_labelled_curve_ids"]
                ),
                "exact_quotient_curve_ids": ",".join(
                    map(str, profile["exact_quotient_labelled_curve_ids"])
                ),
                "rigid_transfer_status": rigid["status"],
                "rigid_transfer_rank_vector": ",".join(
                    f"{row['curve_id']}:{row['rank']}"
                    for row in rigid.get("rank_vector_in_curve_id_order", [])
                ),
                "norm8_off_diagonal_status": offdiag["status"],
                "norm8_off_diagonal_tests": offdiag.get(
                    "off_diagonal_direction_tests"
                ),
                "norm8_off_diagonal_hits": offdiag.get(
                    "off_diagonal_split_direction_count"
                ),
                "rigid_branch_graph_edges": graph["rigid_cover_layer"]["edge_count"],
                "fitted_norm8_branch_graph_edges": graph["fitted_norm_eight_layer"]
                ["edge_count"],
                "native_branch_character_status": native_characters["status"],
                "inherited_branch_character_count": (
                    native_characters["inherited_height_four"] or {}
                ).get("distinct_character_count"),
                "complete_smooth_norm10_character_count": (
                    native_characters["complete_smooth_norm_ten"] or {}
                ).get("distinct_character_count"),
                "full_smooth_marking_saturation_index": (
                    native_characters["full_smooth_layer_marking_obstruction"] or {}
                ).get("complete_old_degree_one_lattice_saturation_index"),
                "low_genus_base_rank_status": low_genus["status"],
                "maximum_certified_base_Jacobian_rank_lower_bound": low_genus.get(
                    "paired_V4_bases", {}
                ).get("maximum_certified_base_Jacobian_rank_lower_bound")
                or low_genus.get("paired_V4_bases", {}).get(
                    "maximum_certified_base_Jacobian_rank_lower_bound_in_screen"
                ),
                "marked_minimum_class_status": tate["status"],
                "marked_norm8_class_count": marked_minima.get(
                    "marked_norm_eight_class_count"
                ),
                "Tate_quotient_status": tate["Tate_quotient"]["status"],
                "matched_exact_curve_parameter_count": len(
                    complexity["current_exact_curve_matches_in_this_PGL2_chart"]
                ),
                "marked_quotient_transport_count": complexity[
                    "chart_specific_exact_quotient_transport_count"
                ],
                "matched_parameter_height_bits_range": (
                    ""
                    if complexity["matched_parameter_projective_height_bits_range"]
                    is None
                    else "-".join(
                        map(
                            str,
                            complexity[
                                "matched_parameter_projective_height_bits_range"
                            ],
                        )
                    )
                ),
                "selected_initial_transport_count": complexity[
                    "selected_initial_transport_count"
                ],
            }
        )

    fields = list(table_rows[0])
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(table_rows)
    table_text = stream.getvalue()

    exact_profile_labels = [
        row["marked_U"]["label"]
        for row in profiles
        if row["exact_quotient_labelled_curve_ids"]
    ]
    total_off_diagonal_tests = sum(
        row["carrier_receptivity_profile"]["off_diagonal_norm_eight_transfer_ranks"].get(
            "off_diagonal_direction_tests", 0
        )
        for row in profiles
    )
    total_off_diagonal_hits = sum(
        row["carrier_receptivity_profile"]["off_diagonal_norm_eight_transfer_ranks"].get(
            "off_diagonal_split_direction_count", 0
        )
        for row in profiles
    )
    total_fitted_edges = sum(
        row["carrier_receptivity_profile"]["branch_incidence_graph"]
        ["fitted_norm_eight_layer"]["edge_count"]
        for row in profiles
    )
    total_rigid_edges = sum(
        row["carrier_receptivity_profile"]["branch_incidence_graph"]
        ["rigid_cover_layer"]["edge_count"]
        for row in profiles
    )
    exact_rigid_profile_count = sum(
        row["carrier_receptivity_profile"]["rigid_cover_transfer_ranks"]["status"].startswith(
            "EXACT_"
        )
        for row in profiles
    )
    exact_off_diagonal_profile_count = sum(
        row["carrier_receptivity_profile"][
            "off_diagonal_norm_eight_transfer_ranks"
        ]["status"].startswith("EXACT_")
        for row in profiles
    )
    exact_branch_profile_count = sum(
        not row["carrier_receptivity_profile"]["branch_incidence_graph"][
            "status"
        ].startswith("UNKNOWN")
        for row in profiles
    )
    exact_low_genus_profile_count = sum(
        row["carrier_receptivity_profile"]["low_genus_base_ranks"]["status"].startswith(
            "EXACT_"
        )
        for row in profiles
    )
    marked_minimum_profile_count = sum(
        row["carrier_receptivity_profile"]["Tate_quotient_and_minimum_classes"][
            "marked_minimum_classes"
        ]
        is not None
        for row in profiles
    )
    parameter_match_profile_count = sum(
        bool(
            row["carrier_receptivity_profile"]["equation_and_parameter_complexity"]
            ["current_exact_curve_matches_in_this_PGL2_chart"]
        )
        for row in profiles
    )
    parameter_match_count = sum(
        len(
            row["carrier_receptivity_profile"]["equation_and_parameter_complexity"]
            ["current_exact_curve_matches_in_this_PGL2_chart"]
        )
        for row in profiles
    )
    if total_off_diagonal_tests != 175 or total_off_diagonal_hits != 0:
        raise ArithmeticError("the first norm-eight off-diagonal audit changed")
    if total_fitted_edges != 75 or total_rigid_edges != 70:
        raise ArithmeticError("the first branch-incidence graph counts changed")

    input_paths = (
        ATLAS,
        SWEEP,
        CALIBRATION,
        NATIVE,
        WGXLI,
        NORM8_CURVE12,
        NORM8_MORE,
        DIVERSITY,
        PAIR_GEOMETRY,
        PAIR_RANKS,
        ALTERNATE_LAB,
        ALTERNATE_BASE_RANKS,
        ALTERNATE_TATE,
        ALTERNATE_CHART_SWEEP,
        COMPLETE_CHARACTER_103B2,
        RIGID_074D9,
        NORM8_074D9,
        SINGULAR_103B2,
        SINGULAR_11952,
    )
    payload = {
        "schema": "elkies-k3.r17-carrier-receptivity-profiles.v1",
        "status": "PASS_FIRST_FAIL_CLOSED_MARKED_U_CARRIER_RECEPTIVITY_PROFILE",
        "definition": {
            "unit": "a primitive ordered marked U=<F,O> in the pinned Neron-Severi marking",
            "coordinates": [
                "rigid-cover transfer ranks",
                "off-diagonal norm-eight transfer ranks",
                "branch-incidence graph",
                "low-genus base ranks",
                "Tate quotient/minimum classes",
                "equation and parameter complexity",
            ],
            "non_invariance_warning": (
                "The profile is not an invariant of the abstract Mordell-Weil lattice. "
                "No value is copied between marked U embeddings merely because their "
                "frame lattices or PGL2 family labels agree."
            ),
        },
        "summary": {
            "marked_U_profiles": len(profiles),
            "PGL2_families": len(set(family_by_chart.values())),
            "frame_class_counts": atlas["atlas"]["frame_class_counts"],
            "exact_quotient_labelled_fibres": len(exact_rows),
            "marked_U_profiles_with_exact_quotient_targets": len(exact_profile_labels),
            "marked_U_labels_with_exact_quotient_targets": sorted(exact_profile_labels),
            "marked_U_profiles_with_exact_rigid_transfer_data": exact_rigid_profile_count,
            "marked_U_profiles_with_exact_off_diagonal_norm_eight_data": (
                exact_off_diagonal_profile_count
            ),
            "marked_U_profiles_with_any_exact_branch_layer": exact_branch_profile_count,
            "marked_U_profiles_with_exact_low_genus_base_data": (
                exact_low_genus_profile_count
            ),
            "marked_U_profiles_with_marked_minimum_class_data": (
                marked_minimum_profile_count
            ),
            "marked_U_profiles_with_exact_Tate_quotient": 0,
            "marked_U_profiles_with_exact_curve_parameter_matches": (
                parameter_match_profile_count
            ),
            "exact_curve_chart_parameter_matches": parameter_match_count,
            "complete_fixed_inventory_rigid_transfer_targets": sum(
                len(rows) for rows in rigid_targets_by_chart.values()
            ),
            "fitted_norm_eight_carriers": total_fitted_edges,
            "norm_eight_off_diagonal_direction_tests": total_off_diagonal_tests,
            "norm_eight_off_diagonal_split_directions": total_off_diagonal_hits,
            "rigid_branch_incidence_edges": total_rigid_edges,
            "fitted_norm_eight_branch_incidence_edges": total_fitted_edges,
            "alternate_marked_U_with_complete_smooth_norm_ten_character_layer": len(
                alternate_complete_character_by_chart
            ),
            "alternate_marked_U_with_inherited_character_layer": len(
                inherited_character_by_chart
            ),
            "alternate_marked_U_with_full_smooth_layer_marking_obstruction": len(
                obstruction_by_chart
            ),
            "unknown_quotient_rows_before_initial_tranche": len(unknown_rows),
            "selected_initial_transport_rows": len(transport_tranche),
            "deferred_unknown_quotient_rows": len(unknown_rows) - len(transport_tranche),
        },
        "profiles": profiles,
        "initial_transport_tranche": {
            "selection_rule": (
                "Take every currently UNKNOWN quotient row at the maximum snapshot rank "
                "lower bound. Do not use family hit counts or inferred historical exposure."
            ),
            "maximum_unknown_snapshot_rank_lower_bound": highest_unknown_rank,
            "selected": transport_tranche,
            "selected_count": len(transport_tranche),
            "deferred_count": len(unknown_rows) - len(transport_tranche),
            "metadata_completion_is_not_the_objective": True,
            "completion_gate": {
                "required": [
                    "an exact saturated rank-17 generic section basis on the selected marked U",
                    "exact specialization of that basis at the recorded rational parameter",
                    "exact independence of the displayed public subgroup",
                    "integral relation recovery followed by exact elliptic-curve group-law replay",
                    "Smith proof of the displayed-subgroup quotient and its primitive generic image",
                ],
                "forbidden_inference": (
                    "A displayed-subgroup quotient is not the quotient of the full "
                    "Mordell-Weil group and supplies no rank upper bound."
                ),
            },
        },
        "identifiability_boundary": {
            "historical_search_exposure": "UNKNOWN",
            "public_hit_counts_used_as_denominators": False,
            "warning": (
                "Successful public fibres have no known denominator of failed historical "
                "trials. Hit counts are omitted from the profile and transport priority."
            ),
        },
        "claim_boundary": {
            "exact": [
                "the 43 marked-U identities, PGL2 family partition, and equation complexities",
                "the twelve chart-specific displayed quotient labels on four marked U embeddings",
                "the twelve complete fixed-inventory rigid transfer ranks",
                "the 175 off-diagonal fitted norm-eight square tests and their zero hits",
                "the 70 rigid and 75 fitted norm-eight specialization-incidence edges",
                "the inherited branch-character layer on ten alternate marked U embeddings",
                "the complete smooth norm-ten branch-character layer on six alternate marked U embeddings",
                "the complete smooth norm-ten branch-character layer on marked U=103b2",
                "the stored low-genus base lower bounds on marked U=103b2 and marked U=11952",
                "the displayed marked minimum-class data and fail-closed Tate status",
                "the rank-28 three-row initial transport tranche",
            ],
            "unknown_or_not_inferred": [
                "the full Mordell-Weil groups of the specialized public fibres",
                "positive target quotient ranks for off-diagonal fitted carriers",
                "carrier data on a PGL2-equivalent chart without a saturated marked transport",
                "every Tate quotient dimension and nonzero minimum Tate class",
                "the 54 lower-rank deferred quotient transports",
                "historical search exposure and causal family quality",
            ],
        },
        "table": relative(TABLE),
        "table_sha256": hashlib.sha256(table_text.encode()).hexdigest(),
        "inputs": {relative(path): digest(path) for path in input_paths},
        "reproducing_command": (
            ".venv/bin/python elkies-k3/scripts/"
            "build_r17_carrier_receptivity_profiles.py"
        ),
    }
    return payload, table_text


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--table-output", type=Path, default=TABLE)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    output = args.output.resolve()
    table_output = args.table_output.resolve()
    payload, table_text = build()
    payload["table"] = relative(table_output)
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not output.exists() or output.read_text() != serialized:
            raise ArithmeticError("stored carrier-receptivity JSON differs from replay")
        if not table_output.exists() or table_output.read_text() != table_text:
            raise ArithmeticError("stored carrier-receptivity table differs from replay")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        table_output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(serialized)
        table_output.write_text(table_text)
    print(
        "R17CARRIERRECEPTIVITY|marked_U=43|exact_quotients=12|"
        "offdiag_tests=175|offdiag_hits=0|initial_transports=3|"
        f"status=PASS|output={relative(output)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
