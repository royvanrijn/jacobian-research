#!/usr/bin/env python3
"""Second rank-aware/diversity tranche from the frozen max-root-300 panel.

The complete max-root-300 artifact contains exact visible-point and mod-3
signatures for all T=1,...,8 fibers of all 1,291 new nonsingular families.
Its first 64-family follow-up found one exact rank-15 specialization.  This
standalone continuation removes those 64 families before selection, then
applies the same outcome-blind 34 global + 3 per diameter-decile rule to the
remaining exact panel records.  Every selected fiber receives conductor-first
PARI replay and one H=5000 point/height screen.  Stable rank at least 15 is
immediately tested by exact finite reductions.
"""

from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
from decimal import Decimal
import hashlib
import json
from pathlib import Path
import platform
import shlex
import sys
import time
from typing import Any, Sequence

from search_mestre_root_tuple_scale import (
    CappedProcessTimeout,
    TARGET_LOG_CONDUCTOR,
    capped_minimal_curve_data,
    point_record,
    sha256_file,
)
from search_mestre_root_tuple_scale_max100 import search_h5000, stable_json_digest
from search_mestre_root_tuple_scale_max200 import (
    fiber_rank_key,
    mod3_independence_certificate,
    visible_points_and_coefficients,
)


ROOT = Path(__file__).resolve().parents[2]
SOURCE_ARTIFACT = (
    ROOT / "artifacts/generated-results/elliptic_mestre_root_tuple_scale_max300.json"
)
EXPECTED_SOURCE_ARTIFACT_SHA256 = (
    "e4a7be774ac0cae3c636c70bde7490e7ced7e313971dfba3c9017e48d730fca7"
)
SOURCE_DRIVER = ROOT / "elliptic-curves/cas/search_mestre_root_tuple_scale_max300.py"
EXPECTED_SOURCE_DRIVER_SHA256 = (
    "922cd33621e882dbb5483b041c547f568d3f4fdfea2bffc8cab0d3741a3445b4"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/"
    "elliptic_mestre_root_tuple_scale_max300_second_tranche.json"
)
SELECTED_COUNT = 64
GLOBAL_KEEP = 34
DIVERSITY_PER_DECILE = 3
STACK_BYTES = 256_000_000
GAIN_TRIGGER = 15
CERTIFICATE_PRIME_BOUND = 499


def best_fiber_records(
    family_records: Sequence[dict[str, Any]],
) -> list[dict[str, Any]]:
    answer = []
    for family in family_records:
        admissible = [record for record in family["fibers"] if record["admissible"]]
        if not admissible:
            continue
        best = dict(sorted(admissible, key=fiber_rank_key)[0])
        best.update(
            {
                "roots": family["roots"],
                "diameter": family["diameter"],
                "diameter_decile": family["diameter_decile"],
            }
        )
        answer.append(best)
    return answer


def select_second_tranche(source: dict[str, Any]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    family_records = source["complete_panel_screen"]["family_records"]
    first = source["rank_aware_diversity_selection"]["selected_records"]
    first_roots = {tuple(record["roots"]) for record in first}
    candidates = [
        record
        for record in best_fiber_records(family_records)
        if tuple(record["roots"]) not in first_roots
    ]
    if len(candidates) != 1_227:
        raise AssertionError("the max300 post-first-tranche family count changed")
    selected: dict[tuple[int, ...], dict[str, Any]] = {}
    for record in sorted(candidates, key=fiber_rank_key)[:GLOBAL_KEEP]:
        item = dict(record)
        item["selection_stratum"] = "second-tranche top-34 exact panel leader"
        selected[tuple(item["roots"])] = item
    diversity_counts = {}
    for lower in range(201, 301, 10):
        upper = lower + 9
        pool = [
            record
            for record in candidates
            if lower <= record["diameter"] <= upper
            and tuple(record["roots"]) not in selected
        ]
        chosen = sorted(pool, key=fiber_rank_key)[:DIVERSITY_PER_DECILE]
        if len(chosen) != DIVERSITY_PER_DECILE:
            raise AssertionError("a second-tranche diameter decile is undersized")
        diversity_counts[f"{lower}-{upper}"] = len(chosen)
        for record in chosen:
            item = dict(record)
            item["selection_stratum"] = f"second-tranche diversity-{lower}-{upper}"
            selected[tuple(item["roots"])] = item
    result = sorted(selected.values(), key=lambda record: record["identifier"])
    if len(result) != SELECTED_COUNT:
        raise AssertionError("the second max300 tranche must have 64 families")
    identifier_digest = hashlib.sha256(
        "\n".join(record["identifier"] for record in result).encode()
    ).hexdigest()
    return result, {
        "source_family_count": len(family_records),
        "first_tranche_excluded_family_count": len(first_roots),
        "remaining_best-fiber_count": len(candidates),
        "global_keep": GLOBAL_KEEP,
        "diversity_per_decile": DIVERSITY_PER_DECILE,
        "diversity_decile_counts": diversity_counts,
        "selected_count": len(result),
        "selected_identifier_sha256": identifier_digest,
        "selection_uses_no_second_tranche_conductor_or_point_outcomes": True,
        "records": result,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--conductor-timeout", type=float, default=8.0)
    parser.add_argument("--point-timeout", type=float, default=12.0)
    parser.add_argument("--height-timeout", type=float, default=12.0)
    args = parser.parse_args()
    if not all(
        0 < value <= 30
        for value in (args.conductor_timeout, args.point_timeout, args.height_timeout)
    ):
        raise SystemExit("all caps must lie in (0,30]")
    if sha256_file(SOURCE_ARTIFACT) != EXPECTED_SOURCE_ARTIFACT_SHA256:
        raise AssertionError("the complete max300 source artifact changed")
    if sha256_file(SOURCE_DRIVER) != EXPECTED_SOURCE_DRIVER_SHA256:
        raise AssertionError("the complete max300 source driver changed")
    source = json.loads(SOURCE_ARTIFACT.read_text(encoding="utf-8"))
    selected, selection = select_second_tranche(source)
    started = time.monotonic()
    records: list[dict[str, Any]] = []
    runtime: dict[str, tuple[tuple[int, ...], int]] = {}
    conductor_started = time.monotonic()
    for selected_record in selected:
        identifier = selected_record["identifier"]
        roots = tuple(map(int, selected_record["roots"]))
        parameter = int(selected_record["parameter"])
        record: dict[str, Any] = {
            "identifier": identifier,
            "roots": list(roots),
            "parameter": parameter,
            "selection_stratum": selected_record["selection_stratum"],
            "panel_visible_rank_lower_bound": selected_record[
                "mod3_finite_reduction_certificate"
            ]["certified_algebraic_rank_lower_bound"],
        }
        try:
            _, coefficients, visible = visible_points_and_coefficients(roots, parameter)
            record["exact_visible_points"] = [point_record(point) for point in visible]
            conductor = capped_minimal_curve_data(
                coefficients,
                timeout=args.conductor_timeout,
                stack_bytes=STACK_BYTES,
            )
            record["conductor_phase"] = {
                "status": "completed exact PARI minimal-model/conductor computation",
                **conductor,
                "below_strict_log_conductor_target_numerically": (
                    Decimal(conductor["log_conductor"]) < TARGET_LOG_CONDUCTOR
                ),
            }
            runtime[identifier] = roots, parameter
        except CappedProcessTimeout:
            record["conductor_phase"] = {
                "status": "timeout-no-retry",
                "timeout_seconds": args.conductor_timeout,
            }
        except Exception as error:
            record["conductor_phase"] = {
                "status": "error-no-retry",
                "error": str(error)[:1000],
            }
        records.append(record)
    conductor_wall = time.monotonic() - conductor_started

    point_started = time.monotonic()
    target_hits = []
    for position, record in enumerate(records, start=1):
        identifier = record["identifier"]
        if identifier not in runtime:
            record["point_triage"] = {"status": "not attempted after incomplete conductor"}
            record["exact_gain_attempt"] = {"status": "not attempted"}
            continue
        roots, parameter = runtime[identifier]
        try:
            triage, subset = search_h5000(
                roots,
                parameter,
                point_timeout=args.point_timeout,
                height_timeout=args.height_timeout,
            )
            record["point_triage"] = triage
            stable_rank = int(triage["stable_numerical_rank"])
            if stable_rank >= GAIN_TRIGGER and subset is not None:
                print(f"EARLY_SIGNAL {identifier} stable_rank={stable_rank}", flush=True)
                _, coefficients, _ = visible_points_and_coefficients(roots, parameter)
                certificate = mod3_independence_certificate(
                    coefficients, subset, prime_bound=CERTIFICATE_PRIME_BOUND
                )
                record["exact_gain_attempt"] = {"mod3": certificate}
                certified = certificate["certified_algebraic_rank_lower_bound"]
                print(f"EXACT_SIGNAL {identifier} certified_rank={certified}", flush=True)
                conductor = record["conductor_phase"]
                if certified >= 30 or (
                    certified >= 21
                    and conductor["below_strict_log_conductor_target_numerically"]
                ):
                    target_hits.append(
                        {
                            "identifier": identifier,
                            "certified_algebraic_rank_lower_bound": certified,
                            "conductor": conductor["conductor"],
                            "log_conductor": conductor["log_conductor"],
                        }
                    )
            else:
                record["exact_gain_attempt"] = {
                    "status": "not triggered",
                    "trigger_stable_numerical_rank": GAIN_TRIGGER,
                }
        except CappedProcessTimeout:
            record["point_triage"] = {
                "status": "timeout-no-retry",
                "point_timeout_seconds": args.point_timeout,
                "height_timeout_seconds": args.height_timeout,
            }
            record["exact_gain_attempt"] = {"status": "not attempted"}
        except Exception as error:
            record["point_triage"] = {
                "status": "error-no-retry",
                "error": str(error)[:1000],
            }
            record["exact_gain_attempt"] = {"status": "not attempted"}
        if position % 16 == 0:
            print(f"H5000 {position}/{len(records)}", flush=True)
    point_wall = time.monotonic() - point_started

    completed = [
        record
        for record in records
        if record["point_triage"]["status"].startswith("completed")
    ]
    rank_histogram = Counter(
        str(record["point_triage"]["stable_numerical_rank"])
        for record in completed
    )
    script = Path(__file__).resolve()
    artifact: dict[str, Any] = {
        "schema_version": 1,
        "status": "complete bounded second max-root-300 leader tranche",
        "target": {
            "rank_at_least": 21,
            "strict_log_conductor_upper_bound": "182.72",
            "alternative_rank_at_least": 30,
            "hits": target_hits,
        },
        "input": {
            "source_artifact": str(SOURCE_ARTIFACT.relative_to(ROOT)),
            "source_artifact_sha256": EXPECTED_SOURCE_ARTIFACT_SHA256,
            "source_driver_sha256": EXPECTED_SOURCE_DRIVER_SHA256,
        },
        "selection": selection,
        "followup": {
            "protocol": {
                "conductor_first": True,
                "one_H5000_search_per_completed_conductor": True,
                "no_retries_or_adaptive_broadening": True,
                "exact_gain_trigger": GAIN_TRIGGER,
            },
            "population": {
                "selected": len(records),
                "conductor_completed": sum(
                    record["conductor_phase"]["status"].startswith("completed")
                    for record in records
                ),
                "conductor_timeouts": sum(
                    record["conductor_phase"]["status"] == "timeout-no-retry"
                    for record in records
                ),
                "subtarget_conductors": sum(
                    record["conductor_phase"].get(
                        "below_strict_log_conductor_target_numerically", False
                    )
                    for record in records
                ),
                "point_search_completed": len(completed),
                "stable_numerical_rank_histogram": dict(sorted(rank_histogram.items())),
                "maximum_stable_numerical_rank": max(
                    (record["point_triage"]["stable_numerical_rank"] for record in completed),
                    default=None,
                ),
                "exact_gain_attempts": sum(
                    "mod3" in record["exact_gain_attempt"] for record in records
                ),
            },
            "records": records,
        },
        "parameters": {
            "conductor_timeout_seconds": args.conductor_timeout,
            "point_timeout_seconds": args.point_timeout,
            "height_timeout_seconds": args.height_timeout,
            "stack_bytes": STACK_BYTES,
        },
        "timings": {
            "conductor_wall_seconds": conductor_wall,
            "point_wall_seconds": point_wall,
            "total_wall_seconds": time.monotonic() - started,
        },
        "provenance": {
            "script": str(script.relative_to(ROOT)),
            "script_sha256": sha256_file(script),
            "command": " ".join(
                shlex.quote(part) for part in [sys.executable, *sys.argv]
            ),
            "python": platform.python_version(),
            "owned_processes_remaining": 0,
        },
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    artifact["result_sha256"] = stable_json_digest(
        {
            "source": EXPECTED_SOURCE_ARTIFACT_SHA256,
            "selection": selection["selected_identifier_sha256"],
            "records": [
                [
                    record["identifier"],
                    record["conductor_phase"]["status"],
                    record["conductor_phase"].get("conductor"),
                    record["point_triage"]["status"],
                    record["point_triage"].get("stable_numerical_rank"),
                    record["point_triage"].get("pool_point_sha256"),
                ]
                for record in records
            ],
            "target": artifact["target"],
        }
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {args.output}")
    print(
        f"completed={len(completed)}/{len(records)}; "
        f"max_rank={artifact['followup']['population']['maximum_stable_numerical_rank']}; "
        f"hits={len(target_hits)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
