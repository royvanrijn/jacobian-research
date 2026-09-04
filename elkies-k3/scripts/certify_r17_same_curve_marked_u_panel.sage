#!/usr/bin/env sage-python
"""Certify a same-curve/different-marked-U carrier-affordance panel.

The five alternate-Q80 control curves occur on every chart in their rational
PGL2 class.  This replay holds the displayed public elliptic-curve subgroup
fixed, specializes a saturated generic MW17 basis on each marked U, and
exhausts the same-size 39,147-cover inventory on every chart.  Canonical
heights are used only to propose integral relations; every retained relation,
Smith calculation, and carrier split is checked exactly.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from fractions import Fraction
from pathlib import Path
import runpy
import sys

from sage.all import PolynomialRing, QQ, ZZ, matrix
from sage.env import SAGE_VERSION


sys.set_int_max_str_digits(0)
ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"

HELPER = ROOT / "elkies-k3/scripts/certify_r17_norm12_native_icarm_quotient_audit.sage"
ATLAS = GENERATED / "elkies-k3-r17-norm12-record-lineage-atlas-v1.json"
CALIBRATION = GENERATED / "elkies-k3-r17-norm12-icarm-calibration-dataset-v1.json"
SWEEP = GENERATED / "elkies-k3-r17-norm12-icarm-database-sweep-v1.json"
PUBLIC = GENERATED / "elkies-k3-r17-norm12-icarm-public-fibres-v1.json"
NATIVE = GENERATED / "elkies-k3-r17-norm12-native-icarm-quotient-audit-v1.json"
CHART_SWEEP = GENERATED / "elkies-k3-r17-norm12-alternate-chart-character-sweep-v1.json"
OUTPUT = GENERATED / "elkies-k3-r17-same-curve-marked-u-panel-v1.json"
TABLE = GENERATED / "elkies-k3-r17-same-curve-marked-u-panel-v1.tsv"

PREFIX_LIMITS = (1000, 5000, 10000, 20000, 39147)
FAMILIES = (
    {
        "representative": "norm12-orbit-11952",
        "curve_ids": (12, 395),
        "charts": ("11952", "08ab4", "091e4", "10f72", "1183a", "098fc"),
    },
    {
        "representative": "norm12-orbit-08f72",
        "curve_ids": (363, 364, 378),
        "charts": ("08f72", "135b7", "09952", "0ae21"),
    },
)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def parameter_height_bits(value: str) -> int:
    value = Fraction(value)
    return max(abs(value.numerator).bit_length(), value.denominator.bit_length())


def direct_path(short_label: str) -> Path:
    saturated = GENERATED / (
        f"elkies-k3-r17-norm12-orbit{short_label}-direct-fibration-saturated-v1.json"
    )
    ordinary = GENERATED / (
        f"elkies-k3-r17-norm12-orbit{short_label}-direct-fibration-v1.json"
    )
    return saturated if saturated.exists() else ordinary


def cover_path(short_label: str) -> Path:
    return GENERATED / (
        f"elkies-k3-r17-norm12-{short_label}-alternate-bisections-full-v1.json"
    )


def exact_span_intersection(left, right) -> int:
    return int(left.ncols() + right.ncols() - left.augment(right).rank())


def generic_matrix(cell: dict):
    return matrix(
        QQ,
        cell["specialized_generic_subgroup"][
            "coordinate_matrix_rows_in_ordered_public_points"
        ],
    )


def visible_matrix(cell: dict):
    generic = generic_matrix(cell)
    splits = cell["complete_carrier_audit"]["splits"]
    if not splits:
        return generic
    split_columns = matrix(
        QQ,
        [record["plus_point_coordinates_in_ordered_public_points"] for record in splits],
    ).transpose()
    answer = generic.augment(split_columns)
    return answer.matrix_from_columns(answer.pivots())


def generator_matrix_integer(cell: dict, include_carriers: bool):
    answer = matrix(
        ZZ,
        cell["specialized_generic_subgroup"][
            "coordinate_matrix_rows_in_ordered_public_points"
        ],
    )
    if include_carriers:
        splits = cell["complete_carrier_audit"]["splits"]
        if splits:
            answer = answer.augment(
                matrix(
                    ZZ,
                    [
                        record["plus_point_coordinates_in_ordered_public_points"]
                        for record in splits
                    ],
                ).transpose()
            )
    return answer


def row_lattice_hnf_of_columns(generators):
    return generators.transpose().hermite_form(include_zero_rows=False)


def pairwise_intersections(cells: list[dict], matrix_builder) -> list[dict]:
    answer = []
    for left_index, left in enumerate(cells):
        left_matrix = matrix_builder(left)
        for right in cells[left_index + 1 :]:
            right_matrix = matrix_builder(right)
            answer.append(
                {
                    "left_marked_U": left["marked_U"],
                    "right_marked_U": right["marked_U"],
                    "left_dimension": int(left_matrix.ncols()),
                    "right_dimension": int(right_matrix.ncols()),
                    "rational_intersection_dimension": exact_span_intersection(
                        left_matrix, right_matrix
                    ),
                }
            )
    return answer


def build():
    helper = runpy.run_path(str(HELPER), run_name="same_curve_marked_u_helper")
    special_fibre = helper["special_fibre"]

    atlas = json.loads(ATLAS.read_text())
    calibration = json.loads(CALIBRATION.read_text())
    sweep = json.loads(SWEEP.read_text())
    public = json.loads(PUBLIC.read_text())
    native = json.loads(NATIVE.read_text())
    chart_sweep = json.loads(CHART_SWEEP.read_text())

    if (
        chart_sweep["status"]
        != "PASS_EXACT_ALTERNATE_CHART_SWEEP_NO_RANK19_OR_RANK20_CHARACTER_INPUT"
    ):
        raise ArithmeticError("the ten-chart complete carrier layer is not exact")
    complete_counts = {
        row["label"]: int(row["cover_count"])
        for row in chart_sweep["complete_smooth_layer"]["records"]
    }
    expected_charts = {
        f"norm12-orbit-{short}"
        for family in FAMILIES
        for short in family["charts"]
    }
    if set(complete_counts) != expected_charts or set(complete_counts.values()) != {39147}:
        raise ArithmeticError("the balanced ten-chart 39,147-cover layer changed")

    chart_by_label = {row["label"]: row for row in atlas["atlas"]["charts"]}
    exact_curve_ids = {
        int(row["curve_id"])
        for row in calibration["rows"]
        if row["displayed_exceptional_quotient_dimension"] is not None
    }
    parameter_cells = [
        (int(hit["curve_id"]), record["chart"])
        for hit in sweep["rational_j_hits_and_twists"]
        if int(hit["curve_id"]) in exact_curve_ids
        for record in hit["native_chart_twists"]
    ]
    if len(parameter_cells) != 112 or len({chart for _, chart in parameter_cells}) != 43:
        raise ArithmeticError("the declared 112-cell, 43-chart parameter panel changed")

    hit_by_curve = {
        int(record["curve_id"]): record
        for record in sweep["rational_j_hits_and_twists"]
    }
    public_by_curve = {int(record["id"]): record for record in public["records"]}
    native_by_cell = {
        (int(record["curve_id"]), record["native_chart"]): record
        for record in native["fibres"]
    }
    ring = PolynomialRing(QQ, "u")

    cells = []
    inputs = {
        relative(path): digest(path)
        for path in (HELPER, ATLAS, CALIBRATION, SWEEP, PUBLIC, NATIVE, CHART_SWEEP)
    }
    for family in FAMILIES:
        representative = family["representative"]
        for short_label in family["charts"]:
            marked_u = f"norm12-orbit-{short_label}"
            d_path = direct_path(short_label)
            c_path = cover_path(short_label)
            direct = json.loads(d_path.read_text())
            covers = json.loads(c_path.read_text())
            if direct["sections"]["status"] != "PASS_EXACT_SATURATED_RANK17_BASIS":
                raise ArithmeticError(f"{marked_u} does not have a saturated MW17 basis")
            if (
                len(covers["bisections"]) != 39147
                or int(covers["construction"]["record_count"]) != 39147
            ):
                raise ArithmeticError(f"{marked_u} does not have the complete carrier layer")
            inputs[relative(d_path)] = digest(d_path)
            inputs[relative(c_path)] = digest(c_path)

            orbit_mask_by_label = {
                record["label"]: int(record["lattice_orbit_mask"])
                for record in covers["bisections"]
            }
            for curve_id in family["curve_ids"]:
                print(
                    f"R17SAMECURVE|curve={curve_id}|marked_U={marked_u}|"
                    "stage=exact_transport",
                    flush=True,
                )
                config = {
                    "source_chart": marked_u,
                    "representative": representative,
                }
                record = special_fibre(
                    config,
                    hit_by_curve[curve_id],
                    public_by_curve[curve_id],
                    direct,
                    covers,
                    ring,
                )
                carrier = record.pop("alternate_q80_cover_audit")
                record.pop("native_chart")
                for split in carrier["splits"]:
                    split["lattice_orbit_mask"] = orbit_mask_by_label[split["label"]]
                priorities = sorted(int(split["priority_rank"]) for split in carrier["splits"])
                carrier["split_priority_ranks"] = priorities
                carrier["first_split_priority_rank"] = priorities[0] if priorities else None
                carrier["prefix_split_counts"] = {
                    str(limit): sum(priority <= limit for priority in priorities)
                    for limit in PREFIX_LIMITS
                }

                equation = chart_by_label[marked_u]["equation_complexity"]
                parameter = record["native_parameter"]
                record["marked_U"] = marked_u
                record["native_anchor_for_curve"] = marked_u == representative
                record["parameter_height_bits"] = parameter_height_bits(parameter)
                record["marked_U_equation_complexity"] = equation
                record["complete_carrier_audit"] = carrier

                native_record = native_by_cell.get((curve_id, marked_u))
                if native_record is not None:
                    if (
                        record["specialized_generic_subgroup"]
                        != native_record["specialized_generic_subgroup"]
                        or record["displayed_exceptional_quotient"]
                        != native_record["displayed_exceptional_quotient"]
                    ):
                        raise ArithmeticError(
                            f"native transport regression at curve {curve_id}, {marked_u}"
                        )
                    native_carrier = native_record["alternate_q80_cover_audit"]
                    for key in (
                        "split_indicator_sha256_in_inventory_order",
                        "rational_split_count",
                        "exact_split_span_rank_in_exceptional_quotient",
                        "split_span_nonzero_smith_invariant_factors",
                    ):
                        if carrier[key] != native_carrier[key]:
                            raise ArithmeticError(
                                f"native carrier regression at curve {curve_id}, {marked_u}: {key}"
                            )
                cells.append(record)

            del direct
            del covers

    if len(cells) != 24 or len({cell["marked_U"] for cell in cells}) != 10:
        raise ArithmeticError("the balanced exact panel is no longer 24 cells on ten charts")

    curve_summaries = []
    table_rows = []
    for family in FAMILIES:
        anchor = family["representative"]
        for curve_id in family["curve_ids"]:
            curve_cells = [cell for cell in cells if cell["curve_id"] == curve_id]
            curve_cells.sort(key=lambda cell: cell["marked_U"])
            anchor_cell = next(cell for cell in curve_cells if cell["marked_U"] == anchor)
            anchor_generic = generic_matrix(anchor_cell)
            anchor_visible = visible_matrix(anchor_cell)
            anchor_generic_hnf = row_lattice_hnf_of_columns(
                generator_matrix_integer(anchor_cell, include_carriers=False)
            )
            anchor_visible_hnf = row_lattice_hnf_of_columns(
                generator_matrix_integer(anchor_cell, include_carriers=True)
            )
            split_counts = {
                cell["complete_carrier_audit"]["rational_split_count"]
                for cell in curve_cells
            }
            span_ranks = {
                cell["complete_carrier_audit"][
                    "exact_split_span_rank_in_exceptional_quotient"
                ]
                for cell in curve_cells
            }
            if len(split_counts) != 1 or len(span_ranks) != 1:
                raise ArithmeticError(
                    f"complete carrier totals changed across presentations of curve {curve_id}"
                )
            first_hits = [
                cell["complete_carrier_audit"]["first_split_priority_rank"]
                for cell in curve_cells
                if cell["complete_carrier_audit"]["first_split_priority_rank"] is not None
            ]
            generic_pairs = pairwise_intersections(curve_cells, generic_matrix)
            visible_pairs = pairwise_intersections(curve_cells, visible_matrix)
            basis_changes = []
            for cell in curve_cells:
                cell_generic = generic_matrix(cell)
                change = anchor_generic.solve_right(cell_generic)
                if not all(value in ZZ for value in change.list()) or abs(change.det()) != 1:
                    raise ArithmeticError(
                        f"curve {curve_id} generic subgroup changed at {cell['marked_U']}"
                    )
                if (
                    row_lattice_hnf_of_columns(
                        generator_matrix_integer(cell, include_carriers=False)
                    )
                    != anchor_generic_hnf
                ):
                    raise ArithmeticError(
                        f"curve {curve_id} generic lattice HNF changed at {cell['marked_U']}"
                    )
                if (
                    row_lattice_hnf_of_columns(
                        generator_matrix_integer(cell, include_carriers=True)
                    )
                    != anchor_visible_hnf
                ):
                    raise ArithmeticError(
                        f"curve {curve_id} visible extension lattice changed at {cell['marked_U']}"
                    )
                change_record = {
                    "marked_U": cell["marked_U"],
                    "determinant": int(change.det()),
                    "maximum_absolute_entry": max(
                        abs(int(value)) for value in change.list()
                    ),
                    "matrix_columns_are_marked_basis_in_native_anchor_basis": [
                        [int(value) for value in row] for row in change.rows()
                    ],
                }
                basis_changes.append(change_record)
                cell["generic_basis_change_from_native_anchor"] = change_record
            quotient_rank = curve_cells[0]["displayed_exceptional_quotient"]["free_rank"]
            curve_summary = {
                "curve_id": curve_id,
                "presentation_count": len(curve_cells),
                "marked_U_labels": [cell["marked_U"] for cell in curve_cells],
                "fixed_snapshot_rank_lower_bound": curve_cells[0][
                    "snapshot_rank_lower_bound"
                ],
                "fixed_displayed_public_subgroup_rank": curve_cells[0][
                    "displayed_point_count"
                ],
                "generic_rank_in_every_presentation": 17,
                "saturated_generic_sublattice_constant_across_presentations": True,
                "generic_basis_changes_from_native_anchor": basis_changes,
                "displayed_quotient_rank_in_every_presentation": quotient_rank,
                "complete_carrier_inventory_size_in_every_presentation": 39147,
                "complete_rational_split_count_in_every_presentation": next(
                    iter(split_counts)
                ),
                "complete_split_span_rank_in_every_presentation": next(iter(span_ranks)),
                "carrier_visible_extension_sublattice_constant_across_presentations": True,
                "first_split_priority_rank_range": (
                    [min(first_hits), max(first_hits)] if first_hits else None
                ),
                "first_split_priority_spread_ratio": (
                    str(Fraction(max(first_hits), min(first_hits))) if first_hits else None
                ),
                "generic_subspace_pairwise_intersections": generic_pairs,
                "generic_subspace_intersection_dimension_range": [
                    min(row["rational_intersection_dimension"] for row in generic_pairs),
                    max(row["rational_intersection_dimension"] for row in generic_pairs),
                ],
                "visible_extension_subspace_pairwise_intersections": visible_pairs,
                "visible_extension_intersection_dimension_range": [
                    min(row["rational_intersection_dimension"] for row in visible_pairs),
                    max(row["rational_intersection_dimension"] for row in visible_pairs),
                ],
            }
            curve_summaries.append(curve_summary)

            for cell in curve_cells:
                carrier = cell["complete_carrier_audit"]
                generic_intersection = exact_span_intersection(
                    anchor_generic, generic_matrix(cell)
                )
                visible_intersection = exact_span_intersection(
                    anchor_visible, visible_matrix(cell)
                )
                table_rows.append(
                    {
                        "curve_id": curve_id,
                        "snapshot_rank_lower_bound": cell["snapshot_rank_lower_bound"],
                        "marked_U": cell["marked_U"],
                        "native_anchor": int(cell["native_anchor_for_curve"]),
                        "parameter_height_bits": cell["parameter_height_bits"],
                        "deep_equation_rank": cell["marked_U_equation_complexity"][
                            "deep_equation_rank"
                        ],
                        "generic_intersection_with_native_anchor": generic_intersection,
                        "displayed_quotient_rank": cell["displayed_exceptional_quotient"][
                            "free_rank"
                        ],
                        "complete_split_count": carrier["rational_split_count"],
                        "complete_split_span_rank": carrier[
                            "exact_split_span_rank_in_exceptional_quotient"
                        ],
                        "visible_extension_intersection_with_native_anchor": visible_intersection,
                        "first_split_priority_rank": carrier[
                            "first_split_priority_rank"
                        ],
                        **{
                            f"splits_in_first_{limit}": carrier["prefix_split_counts"][
                                str(limit)
                            ]
                            for limit in PREFIX_LIMITS
                        },
                    }
                )

    table_stream = io.StringIO()
    writer = csv.DictWriter(
        table_stream,
        fieldnames=list(table_rows[0]),
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(table_rows)
    table_text = table_stream.getvalue()

    nonzero_first_hit_summaries = [
        row for row in curve_summaries if row["first_split_priority_rank_range"] is not None
    ]
    payload = {
        "schema": "elkies-k3.r17-same-curve-marked-u-panel.v1",
        "status": "PASS_EXACT_SAME_CURVE_DIFFERENT_MARKED_U_PANEL",
        "experimental_unit": (
            "one fixed public elliptic curve with its ordered displayed independent "
            "point subgroup, observed under several marked U presentations"
        ),
        "treatment": (
            "the primitive ordered marked U=<F,O>, its chart-specific generic MW17 "
            "basis, and its complete carrier geometry and deterministic cost order"
        ),
        "fixed_within_curve": [
            "the QQ-isomorphism class of the specialized elliptic curve",
            "the ordered displayed independent public-point subgroup",
            "the snapshot rank lower bound",
        ],
        "summary": {
            "available_exact_parameter_panel_cells": 112,
            "available_exact_parameter_panel_marked_U_count": 43,
            "balanced_complete_inventory_subpanel_cells": len(cells),
            "balanced_complete_inventory_subpanel_marked_U_count": 10,
            "balanced_complete_inventory_subpanel_curve_count": 5,
            "new_non_native_exact_transports": sum(
                not cell["native_anchor_for_curve"] for cell in cells
            ),
            "native_anchor_transports_replayed": sum(
                cell["native_anchor_for_curve"] for cell in cells
            ),
            "complete_carriers_per_cell": 39147,
            "all_within_curve_complete_split_counts_constant": True,
            "all_within_curve_complete_split_span_ranks_constant": True,
            "all_within_curve_saturated_generic_sublattices_constant": True,
            "all_within_curve_carrier_visible_extension_sublattices_constant": True,
            "largest_first_split_priority_spread": max(
                row["first_split_priority_rank_range"][1]
                / row["first_split_priority_rank_range"][0]
                for row in nonzero_first_hit_summaries
            ),
            "largest_first_split_priority_spread_curve_id": max(
                nonzero_first_hit_summaries,
                key=lambda row: row["first_split_priority_rank_range"][1]
                / row["first_split_priority_rank_range"][0],
            )["curve_id"],
        },
        "curve_summaries": curve_summaries,
        "cells": cells,
        "interpretation": {
            "exact_experimental_result": (
                "For each fixed curve, the complete inventory split count and its exact "
                "quotient-span rank are constant across the tested presentations. More "
                "strongly, the saturated generic subgroup and the carrier-visible extension "
                "sublattice are identical inside the fixed displayed public group. The "
                "deterministic first-hit cost position nevertheless moves by more than an "
                "order of magnitude."
            ),
            "search_affordance_conclusion": (
                "Variation in prefix visibility, first-hit priority, equation cost, or "
                "carrier labels is presentation-dependent search exposure. It is not "
                "evidence that the fixed specialized curve acquired rank."
            ),
            "mechanism_boundary": (
                "Only a quantity stable under the within-curve presentation changes is "
                "eligible here as an intrinsic rank-mechanism candidate. Constancy in "
                "this five-curve panel is evidence for this panel, not a general theorem."
            ),
        },
        "claim_boundary": {
            "proved": [
                "24 exact saturated generic-subgroup transports on five fixed curves and ten marked U charts",
                "24 exhaustive 39,147-cover split censuses and exact displayed-quotient span calculations",
                "the exact within-curve generic and visible-extension rational-subspace intersections",
                "equality of the saturated generic and carrier-visible extension sublattices within each fixed curve",
                "the exact chart-cost positions of every rationally splitting carrier",
            ],
            "not_proved": [
                "that any displayed public subgroup is the full Mordell-Weil group",
                "an exact Mordell-Weil rank or upper bound for any specialized curve",
                "presentation invariance beyond the five curves and ten complete inventories tested",
                "a causal model of historical search exposure or a population-level rank predictor",
                "the remaining 88 exact parameter matches without a balanced complete carrier inventory",
            ],
        },
        "table": relative(TABLE),
        "table_sha256": hashlib.sha256(table_text.encode()).hexdigest(),
        "inputs": inputs,
        "software_assumptions": {
            "sage_version": SAGE_VERSION,
            "required_features": [
                "exact QQ elliptic-curve group law",
                "exact Smith normal form",
                "exact rational matrix rank and intersections",
                "canonical heights used only for candidate relation recovery",
            ],
        },
        "reproducing_command": (
            "PYTHONPATH=elliptic-curves/cas sage -python elkies-k3/scripts/"
            "certify_r17_same_curve_marked_u_panel.sage"
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
            raise ArithmeticError("stored same-curve marked-U panel differs from replay")
        if not table_output.exists() or table_output.read_text() != table_text:
            raise ArithmeticError("stored same-curve marked-U table differs from replay")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        table_output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(serialized)
        table_output.write_text(table_text)
    print(
        "R17SAMECURVE|cells=24|charts=10|curves=5|new_transports=19|"
        f"status=PASS|output={relative(output)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
