#!/usr/bin/env python3
"""Assemble the exact/unknown-labelled 69-fibre R17 calibration table."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-icarm-local-fingerprints-v1.json"
NATIVE = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-native-icarm-quotient-audit-v1.json"
WGXLI = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-wgxli-lineage-fibres-v1.json"
NORM8 = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-curve12-norm8-incidence-v1.json"
NORM8_MORE = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-icarm-norm8-incidence-v1.json"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-icarm-calibration-dataset-v1.json"
TABLE = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-icarm-calibration-dataset-v1.tsv"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def build():
    local = json.loads(LOCAL.read_text())
    native = json.loads(NATIVE.read_text())
    wgxli = json.loads(WGXLI.read_text())
    norm8 = json.loads(NORM8.read_text())
    norm8_more = json.loads(NORM8_MORE.read_text())
    if len(local["fibres"]) != 69:
        raise ArithmeticError("the local feature table no longer has 69 rows")
    native_by_id = {int(row["curve_id"]): row for row in native["fibres"]}
    wgxli_by_id = {
        int(row["curve_id"]): row for row in wgxli["exceptional_quotients"]
    }
    overlap = set(native_by_id) & set(wgxli_by_id)
    if overlap:
        raise ArithmeticError(f"duplicate exact quotient audits: {sorted(overlap)}")
    exact_ids = set(native_by_id) | set(wgxli_by_id)
    if len(exact_ids) != 12:
        raise ArithmeticError("expected twelve exact displayed-quotient audits")
    hit_count = {
        row["representative"]: int(row["recognized_public_hit_count"])
        for row in local["class_summaries"]
    }
    norm8_by_id = {
        int(norm8["curve"]["icarm_id"]): int(
            norm8["incidence_signature"]["successful_directions"]
        ),
        **{
            int(row["curve_id"]): int(
                row["incidence_signature"]["successful_directions"]
            )
            for row in norm8_more["fibres"]
        },
    }
    if set(norm8_by_id) != {12, 363, 364, 378, 395}:
        raise ArithmeticError("expected five exact norm-eight incidence audits")

    rows = []
    for source in local["fibres"]:
        curve_id = int(source["curve_id"])
        native_row = native_by_id.get(curve_id)
        old_row = wgxli_by_id.get(curve_id)
        if native_row is not None:
            quotient_dimension = int(
                native_row["displayed_exceptional_quotient"]["free_rank"]
            )
            cover = native_row["alternate_q80_cover_audit"]
            quotient_source = "new native exact audit"
            fixed_cover_count = int(cover["covers_evaluated"])
            fixed_split_count = int(cover["rational_split_count"])
            fixed_visibility_rank = int(
                cover["exact_split_span_rank_in_exceptional_quotient"]
            )
            fixed_visibility_primitive = bool(cover["split_span_is_primitive"])
        elif old_row is not None:
            quotient_dimension = int(old_row["free_rank"])
            quotient_source = "prior exact 074d9/wgxli audit"
            fixed_cover_count = None
            fixed_split_count = None
            fixed_visibility_rank = None
            fixed_visibility_primitive = None
        else:
            quotient_dimension = None
            quotient_source = "UNKNOWN: saturated chart-specific section transport not yet compiled"
            fixed_cover_count = None
            fixed_split_count = None
            fixed_visibility_rank = None
            fixed_visibility_primitive = None

        representative = source["representative"]
        row = {
            "curve_id": curve_id,
            "parameter": source["parameter"],
            "family": representative,
            "family_frame_class": source["family_frame_class"],
            "leave_one_family_out_group": representative,
            "snapshot_rank_lower_bound": int(source["snapshot_rank_lower_bound"]),
            "rank_jump_lower_bound_over_generic_17": int(source["snapshot_rank_lower_bound"]) - 17,
            "recognized_public_hit_count_in_family": hit_count[representative],
            "historical_search_exposure": None,
            "search_exposure_status": (
                "UNKNOWN: public hit count is a selection/effort proxy, not a denominator "
                "of searched parameters"
            ),
            "displayed_exceptional_quotient_dimension": quotient_dimension,
            "quotient_audit_source": quotient_source,
            "fixed_native_cover_inventory_size": fixed_cover_count,
            "fixed_native_rational_split_count": fixed_split_count,
            "fixed_native_geometric_visibility_rank": fixed_visibility_rank,
            "fixed_native_visibility_span_primitive": fixed_visibility_primitive,
            "norm8_fitted_genus_one_incidence_directions": norm8_by_id.get(curve_id),
            "norm8_incidence_status": (
                "exact fitted positive control"
                if curve_id in norm8_by_id
                else "not run in this audit"
            ),
            "local_prime_block_scores_units_1e12": source["block_score_units_1e12"],
            "local_total_nagao_score_units_1e12": int(source["total_score_units_1e12"]),
            "local_good_prime_count": int(source["good_prime_count"]),
            "local_bad_reduction_prime_count": int(source["bad_reduction_prime_count"]),
            "bad_prime_count": int(source["bad_prime_count"]),
            "conductor_decimal_digits": int(source["conductor_decimal_digits"]),
        }
        if quotient_dimension is not None and quotient_dimension != row["rank_jump_lower_bound_over_generic_17"]:
            raise ArithmeticError(f"curve {curve_id} quotient/rank lower-bound mismatch")
        if (
            row["norm8_fitted_genus_one_incidence_directions"] is not None
            and row["norm8_fitted_genus_one_incidence_directions"] != quotient_dimension
        ):
            raise ArithmeticError(f"curve {curve_id} norm-eight/quotient mismatch")
        rows.append(row)

    fields = list(rows[0])
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for row in rows:
        serialized = dict(row)
        serialized["local_prime_block_scores_units_1e12"] = ",".join(
            map(str, row["local_prime_block_scores_units_1e12"])
        )
        writer.writerow(serialized)
    table_text = stream.getvalue()

    payload = {
        "schema": "elkies-k3.r17-norm12-icarm-calibration-dataset.v1",
        "status": "PASS_FAIL_CLOSED_69_FIBRE_CALIBRATION_DATASET",
        "summary": {
            "rows": len(rows),
            "families": 6,
            "rows_with_exact_displayed_exceptional_quotient": sum(
                row["displayed_exceptional_quotient_dimension"] is not None for row in rows
            ),
            "rows_with_complete_fixed_native_cover_visibility": sum(
                row["fixed_native_cover_inventory_size"] is not None for row in rows
            ),
            "rows_with_exact_norm8_fitted_incidence": len(norm8_by_id),
            "rows_with_local_fingerprints": len(rows),
            "unknown_quotient_rows_are_literal_nulls": True,
        },
        "rows": rows,
        "table": relative(TABLE),
        "table_sha256": hashlib.sha256(table_text.encode()).hexdigest(),
        "recommended_evaluation_split": (
            "Group by leave_one_family_out_group. Never randomly split fibres from the "
            "same PGL2 family between train and test."
        ),
        "identifiability_warning": (
            "The endpoint contains successful public fibres but no count or parameter "
            "distribution for failed historical trials. Family quality and search effort "
            "are therefore not separately identifiable from hit counts alone."
        ),
        "claim_boundary": {
            "exact": [
                "all local-symbol features",
                "the twelve non-null displayed exceptional quotient dimensions",
                "the seven non-null fixed-cover visibility records",
                "the five non-null norm-eight fitted incidence counts",
            ],
            "heuristic_or_unknown": [
                "Nagao scores as predictors",
                "historical search exposure",
                "family-quality attribution",
                "the 57 null displayed quotient dimensions",
            ],
        },
        "inputs": {
            relative(path): digest(path)
            for path in (LOCAL, NATIVE, WGXLI, NORM8, NORM8_MORE)
        },
        "reproducing_command": (
            ".venv/bin/python elkies-k3/scripts/"
            "build_r17_norm12_icarm_calibration_dataset.py"
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
            raise ArithmeticError("stored calibration JSON differs from replay")
        if not table_output.exists() or table_output.read_text() != table_text:
            raise ArithmeticError("stored calibration table differs from replay")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        table_output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(serialized)
        table_output.write_text(table_text)
    print(
        "R17ICARMCALIBRATION|rows=69|exact_quotients=12|native_visibility=7|"
        f"norm8=5|local=69|status=PASS|output={relative(output)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
