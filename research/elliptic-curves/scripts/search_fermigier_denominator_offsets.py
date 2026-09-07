#!/usr/bin/env python3
"""Denominator-aware exact quartic search on mixed-small-prime survivors."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import platform
import sys


PROGRAM_ROOT = Path(__file__).resolve().parents[1]
ROOT = PROGRAM_ROOT.parent
sys.path.insert(0, str(PROGRAM_ROOT))

from ecsearch.fermigier_offset_search import (  # noqa: E402
    denominator_offset_points,
    point_stream_digest,
)
from ecsearch.fermigier_rank import (  # noqa: E402
    section_and_point_cloud_differences,
    specialize_fermigier_rank_sections,
)
from ecsearch.rank_certification import select_independent_subset  # noqa: E402


DEFAULT_INPUT = ROOT / "artifacts/generated-results/elliptic-curves/fermigier_mixed_small_prime_crt_gauss_v1.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elliptic-curves/fermigier_denominator_offsets_v1.json"
EXPECTED_INPUT_SHA256 = "e5c50c236c86bcccd55b5ad45a578202b086ff34d5b9cd068c06f3c541ef0690"
MAXIMUM_DENOMINATOR = 64
MAXIMUM_ABS_NUMERATOR = 5_000
DEEP_MINIMUM_DENOMINATOR = 65
DEEP_MAXIMUM_DENOMINATOR = 256
DEEP_MAXIMUM_ABS_NUMERATOR = 20_000
EXPECTED_CANDIDATES = (
    "119/2", "595/2", "1155/2", "91/5", "1015", "245/3",
    "1785/4", "350", "385/12", "448/5", "1491/4", "539",
)


def rational_text(value: Q) -> str:
    return str(value.numerator) if value.denominator == 1 else str(value)


def stable_digest(artifact: dict) -> str:
    stable = {key: value for key, value in artifact.items() if key not in {"generated_at_utc", "result_sha256"}}
    return hashlib.sha256(json.dumps(stable, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def best_certificate(specialization, quartic_points):
    cloud = section_and_point_cloud_differences(specialization, quartic_points)
    attempts = []
    for relation_prime in (5, 3, 7):
        selected, certificate = select_independent_subset(
            specialization.canonical_model,
            cloud,
            relation_prime=relation_prime,
            maximum_reduction_prime=2_000,
        )
        attempts.append((len(selected), selected, certificate))
        if len(selected) == len(cloud):
            break
    rank, selected, certificate = max(
        attempts, key=lambda item: (item[0], -item[2].relation_prime)
    )
    return cloud, rank, selected, certificate


def run(input_path: Path, *, maximum_denominator: int, maximum_abs_numerator: int) -> dict:
    raw = input_path.read_bytes()
    observed_sha = hashlib.sha256(raw).hexdigest()
    if input_path == DEFAULT_INPUT and observed_sha != EXPECTED_INPUT_SHA256:
        raise AssertionError("the pinned mixed-small-prime input changed")
    source = json.loads(raw)
    source_rows = {
        row["adapter_u"]: row
        for row in source["selected"]
        if row["point_screen"]["certified_rank_lower_bound"] == 13
    }
    if tuple(source_rows) != EXPECTED_CANDIDATES:
        raise AssertionError("the pinned rank-13 candidate order changed")
    records = []
    for parameter_text, source_row in source_rows.items():
        parameter = Q(parameter_text)
        specialization = specialize_fermigier_rank_sections(parameter)
        new_points = denominator_offset_points(
            parameter,
            maximum_denominator=maximum_denominator,
            maximum_abs_numerator=maximum_abs_numerator,
        )
        prior_points = tuple(
            (Q(x_value), Q(y_value))
            for x_value, y_value in source_row["point_screen"]["extra_quartic_points"]
        )
        prior_x = {point[0] for point in specialization.quartic_points}
        prior_x.update(point[0] for point in prior_points)
        genuinely_new = tuple(point for point in new_points if point.x not in prior_x)
        combined = prior_points + tuple((point.x, point.raw_y) for point in genuinely_new)
        cloud, rank, selected, certificate = best_certificate(specialization, combined)
        records.append(
            {
                "adapter_u": parameter_text,
                "source_integer_offset_extra_points": len(prior_points),
                "denominator_offset_points_before_prior_exclusion": len(new_points),
                "new_denominator_offset_point_count": len(genuinely_new),
                "new_point_stream_sha256": point_stream_digest(genuinely_new),
                "new_points": [
                    {
                        "sign": point.sign,
                        "offset_numerator": point.offset_numerator,
                        "offset_denominator": point.offset_denominator,
                        "x": rational_text(point.x),
                        "raw_y": rational_text(point.raw_y),
                    }
                    for point in genuinely_new
                ],
                "deduplicated_jacobian_difference_count": len(cloud),
                "selected_indices": list(selected),
                "certified_rank_lower_bound": rank,
                "finite_reduction_certificate": certificate.to_json_object(),
            }
        )
    first_stage_maximum = max(row["certified_rank_lower_bound"] for row in records)
    for row in records:
        if row["certified_rank_lower_bound"] != first_stage_maximum or first_stage_maximum <= 13:
            row["deep_followup"] = {"status": "not_selected"}
            continue
        parameter_text = row["adapter_u"]
        parameter = Q(parameter_text)
        source_row = source_rows[parameter_text]
        specialization = specialize_fermigier_rank_sections(parameter)
        deep_points = denominator_offset_points(
            parameter,
            minimum_denominator=DEEP_MINIMUM_DENOMINATOR,
            maximum_denominator=DEEP_MAXIMUM_DENOMINATOR,
            maximum_abs_numerator=DEEP_MAXIMUM_ABS_NUMERATOR,
        )
        prior_points = tuple(
            (Q(x_value), Q(y_value))
            for x_value, y_value in source_row["point_screen"]["extra_quartic_points"]
        )
        first_stage_points = tuple(
            (Q(point["x"]), Q(point["raw_y"])) for point in row["new_points"]
        )
        prior_x = {point[0] for point in specialization.quartic_points}
        prior_x.update(point[0] for point in prior_points)
        prior_x.update(point[0] for point in first_stage_points)
        genuinely_new = tuple(point for point in deep_points if point.x not in prior_x)
        combined = prior_points + first_stage_points + tuple(
            (point.x, point.raw_y) for point in genuinely_new
        )
        cloud, rank, selected, certificate = best_certificate(specialization, combined)
        row["deep_followup"] = {
            "status": "completed",
            "minimum_denominator": DEEP_MINIMUM_DENOMINATOR,
            "maximum_denominator": DEEP_MAXIMUM_DENOMINATOR,
            "maximum_abs_offset_numerator": DEEP_MAXIMUM_ABS_NUMERATOR,
            "new_point_count": len(genuinely_new),
            "new_point_stream_sha256": point_stream_digest(genuinely_new),
            "new_points": [
                {
                    "sign": point.sign,
                    "offset_numerator": point.offset_numerator,
                    "offset_denominator": point.offset_denominator,
                    "x": rational_text(point.x),
                    "raw_y": rational_text(point.raw_y),
                }
                for point in genuinely_new
            ],
            "deduplicated_jacobian_difference_count": len(cloud),
            "selected_indices": list(selected),
            "certified_rank_lower_bound": rank,
            "finite_reduction_certificate": certificate.to_json_object(),
        }
    final_ranks = [
        row["deep_followup"].get(
            "certified_rank_lower_bound", row["certified_rank_lower_bound"]
        )
        for row in records
    ]
    artifact = {
        "schema": "elliptic-curves.fermigier-denominator-offsets.v1",
        "status": "complete bounded exact denominator-offset search",
        "claim_level": "finite-reduction rank lower bounds only",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "input": {
            "path": str(input_path.relative_to(ROOT)),
            "sha256": observed_sha,
            "source_result_sha256": source["result_sha256"],
            "selection": "all twelve source rows with certified rank lower bound 13",
        },
        "search": {
            "charts": "x=+/-2u+n/d",
            "minimum_denominator": 2,
            "maximum_denominator": maximum_denominator,
            "maximum_abs_offset_numerator": maximum_abs_numerator,
            "offsets_reduced": True,
            "square_test": "exact integer clearing B^6*d^4 plus modular square filters and isqrt",
            "deep_followup": {
                "selection": "all first-stage maximum-rank fibers when that maximum exceeds 13",
                "minimum_denominator": DEEP_MINIMUM_DENOMINATOR,
                "maximum_denominator": DEEP_MAXIMUM_DENOMINATOR,
                "maximum_abs_offset_numerator": DEEP_MAXIMUM_ABS_NUMERATOR,
            },
        },
        "candidates": records,
        "outcome": {
            "candidates_searched": len(records),
            "candidates_with_new_denominator_points": sum(bool(row["new_denominator_offset_point_count"]) for row in records),
            "new_denominator_offset_points": sum(row["new_denominator_offset_point_count"] for row in records),
            "deep_followup_candidates": sum(row["deep_followup"]["status"] == "completed" for row in records),
            "deep_followup_new_points": sum(row["deep_followup"].get("new_point_count", 0) for row in records),
            "maximum_certified_rank_lower_bound": max(final_ranks),
            "target_met": any(rank >= 21 for rank in final_ranks),
            "boundary": "This is a bounded rational-offset slice, not a complete quartic point search or a rank upper bound.",
        },
        "software": {"python": platform.python_version()},
        "reproducing_command": (
            "PYTHONPATH=elliptic-curves python3 "
            "elliptic-curves/scripts/search_fermigier_denominator_offsets.py"
        ),
    }
    artifact["result_sha256"] = stable_digest(artifact)
    return artifact


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--maximum-denominator", type=int, default=MAXIMUM_DENOMINATOR)
    parser.add_argument("--maximum-abs-numerator", type=int, default=MAXIMUM_ABS_NUMERATOR)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if not 2 <= args.maximum_denominator <= MAXIMUM_DENOMINATOR:
        raise SystemExit("maximum denominator outside the pinned range")
    if not 1 <= args.maximum_abs_numerator <= MAXIMUM_ABS_NUMERATOR:
        raise SystemExit("maximum numerator outside the pinned range")
    artifact = run(
        args.input,
        maximum_denominator=args.maximum_denominator,
        maximum_abs_numerator=args.maximum_abs_numerator,
    )
    if args.check:
        expected = json.loads(args.output.read_text())
        if artifact["result_sha256"] != expected["result_sha256"]:
            raise AssertionError("the denominator-offset result digest changed")
        print(f"PASS {artifact['result_sha256']}")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(json.dumps(artifact["outcome"], sort_keys=True))
    print(f"result_sha256={artifact['result_sha256']}")


if __name__ == "__main__":
    main()
