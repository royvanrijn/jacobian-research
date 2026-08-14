#!/usr/bin/env python3
"""Escalate the three rank-14 leaders from max300 tranche two.

The frozen second-tranche artifact fixes exactly three H=5000 stable-rank-14
fibers, all well below the conductor target.  Replay each at H=50000,
H=250000, and H=1000000 with exact quartic membership/mapping and independent
72/120-digit height matrices.  Any stable rank at least 15 receives an exact
mod-3 finite-reduction certificate immediately.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import shlex
import sys
import time
from typing import Any

from mestre_root_tuples import SixRootMestreConstruction
from search_mestre_root_tuple_scale import (
    CappedProcessTimeout,
    bounded_quartic_points,
    canonical_signless_points,
    height_matrix_replay,
    numerical_subset,
    point_digest,
    point_record,
    primitive_visible_points,
    quartic_point_to_jacobian,
    quartic_value,
    sha256_file,
)
from search_mestre_root_tuple_scale_max100 import stable_json_digest
from search_mestre_root_tuple_scale_max200 import (
    mod3_independence_certificate,
    visible_points_and_coefficients,
)


Q = Fraction
ROOT = Path(__file__).resolve().parents[2]
SOURCE_SCRIPT = (
    ROOT
    / "elliptic-curves/cas/"
    "search_mestre_root_tuple_scale_max300_second_tranche.py"
)
EXPECTED_SOURCE_SCRIPT_SHA256 = (
    "a03456a4b844653432b9435cd9193cf721920ac3b52ada32e771bfd084009596"
)
SOURCE_ARTIFACT = (
    ROOT
    / "artifacts/generated-results/"
    "elliptic_mestre_root_tuple_scale_max300_second_tranche.json"
)
EXPECTED_SOURCE_ARTIFACT_SHA256 = (
    "94662ac54a1666f3ee238a7849db3b07cf3a433308b940d682b918a8cd58b796"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/"
    "elliptic_mestre_root_tuple_scale_max300_second_tranche_extension.json"
)
HEIGHT_STAGES = (50_000, 250_000, 1_000_000)
EXPECTED_LEADER_IDS = (
    "r0_17_167_204_209_255_t4",
    "r0_27_116_143_189_225_t2",
    "r0_28_161_205_228_281_t1",
)
STACK_BYTES = 256_000_000
MAPPING_CAP = 512
CERTIFICATE_PRIME_BOUND = 499


def search_height(
    roots: tuple[int, ...],
    parameter: int,
    height_bound: int,
    *,
    point_timeout: float,
    height_timeout: float,
) -> tuple[dict[str, Any], tuple[tuple[Fraction, Fraction], ...]]:
    construction = SixRootMestreConstruction(tuple(Q(root) for root in roots))
    parameter_q = Q(parameter)
    coefficients = construction.primitive_jacobian_coefficients(parameter_q)
    quartic_coefficients = construction.primitive_quartic_coefficients(parameter_q)
    visible_quartic = primitive_visible_points(construction, parameter_q)
    visible_jacobian = tuple(
        quartic_point_to_jacobian(construction, parameter_q, point)
        for point in visible_quartic
    )
    raw = bounded_quartic_points(
        quartic_coefficients,
        height_bound=height_bound,
        timeout=point_timeout,
        stack_bytes=STACK_BYTES,
    )
    signless = canonical_signless_points(raw)
    retained = signless[:MAPPING_CAP]
    if len(signless) > len(retained):
        raise AssertionError("the declared mapping cap truncated a leader pool")
    if any(
        point[1] ** 2 != quartic_value(quartic_coefficients, point[0])
        for point in retained
    ):
        raise AssertionError("a bounded point lies off the exact quartic")
    searched = tuple(
        quartic_point_to_jacobian(construction, parameter_q, point)
        for point in retained
    )
    pool_by_x = {point[0]: point for point in visible_jacobian}
    for point in searched:
        pool_by_x.setdefault(point[0], point)
    pool = tuple(pool_by_x.values())
    heights = height_matrix_replay(
        coefficients,
        pool,
        precisions=(72, 120),
        timeout=height_timeout,
        stack_bytes=STACK_BYTES,
    )
    subset = numerical_subset(pool, heights)
    return (
        {
            "status": "completed exact bounded point checks and numerical height triage",
            "height_bound": height_bound,
            "signed_points_returned": len(raw),
            "distinct_nonzero_ordinate_abscissas": len(signless),
            "visible_jacobian_point_count": len(visible_jacobian),
            "pool_point_count_modulo_inverse": len(pool),
            "pool_point_sha256": point_digest(pool),
            "height_matrix_runs": list(heights),
            "stable_numerical_rank": int(heights[-1]["numerical_rank"]),
            "numerical_subset": [point_record(point) for point in subset],
            "mapping_cap": MAPPING_CAP,
            "mapping_truncated": False,
            "numerical_rank_is_not_an_independence_certificate": True,
        },
        subset,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--point-timeout", type=float, default=90.0)
    parser.add_argument("--height-timeout", type=float, default=30.0)
    args = parser.parse_args()
    if not 0 < args.point_timeout <= 120 or not 0 < args.height_timeout <= 60:
        raise SystemExit("invalid point/height caps")
    if sha256_file(SOURCE_SCRIPT) != EXPECTED_SOURCE_SCRIPT_SHA256:
        raise AssertionError("the second-tranche source script changed")
    if sha256_file(SOURCE_ARTIFACT) != EXPECTED_SOURCE_ARTIFACT_SHA256:
        raise AssertionError("the second-tranche source artifact changed")
    source = json.loads(SOURCE_ARTIFACT.read_text(encoding="utf-8"))
    leaders = [
        record
        for record in source["followup"]["records"]
        if record["point_triage"].get("stable_numerical_rank") == 14
    ]
    if tuple(record["identifier"] for record in leaders) != EXPECTED_LEADER_IDS:
        raise AssertionError("the exact rank-14 leader set changed")

    started = time.monotonic()
    records = []
    target_hits = []
    for source_record in leaders:
        roots = tuple(map(int, source_record["roots"]))
        parameter = int(source_record["parameter"])
        record: dict[str, Any] = {
            "identifier": source_record["identifier"],
            "roots": list(roots),
            "parameter": parameter,
            "conductor": source_record["conductor_phase"]["conductor"],
            "log_conductor": source_record["conductor_phase"]["log_conductor"],
            "root_number": source_record["conductor_phase"]["root_number"],
            "source_H5000": source_record["point_triage"],
            "stages": [],
            "exact_gain_certificate": {"status": "not triggered"},
        }
        best_subset: tuple[tuple[Fraction, Fraction], ...] = ()
        best_rank = 14
        for height_bound in HEIGHT_STAGES:
            try:
                stage, subset = search_height(
                    roots,
                    parameter,
                    height_bound,
                    point_timeout=args.point_timeout,
                    height_timeout=args.height_timeout,
                )
                record["stages"].append(stage)
                if stage["stable_numerical_rank"] > best_rank:
                    best_rank = stage["stable_numerical_rank"]
                    best_subset = subset
                    print(
                        f"EARLY_SIGNAL {record['identifier']} H={height_bound} rank={best_rank}",
                        flush=True,
                    )
            except CappedProcessTimeout:
                record["stages"].append(
                    {
                        "status": "timeout-no-retry",
                        "height_bound": height_bound,
                        "point_timeout_seconds": args.point_timeout,
                        "height_timeout_seconds": args.height_timeout,
                    }
                )
        record["maximum_stable_numerical_rank"] = best_rank
        if best_rank >= 15 and best_subset:
            _, coefficients, _ = visible_points_and_coefficients(roots, parameter)
            certificate = mod3_independence_certificate(
                coefficients, best_subset, prime_bound=CERTIFICATE_PRIME_BOUND
            )
            record["exact_gain_certificate"] = certificate
            certified = certificate["certified_algebraic_rank_lower_bound"]
            print(f"EXACT_SIGNAL {record['identifier']} rank={certified}", flush=True)
            if certified >= 21 and Q(record["log_conductor"]) < Q(4568, 25):
                target_hits.append(
                    {
                        "identifier": record["identifier"],
                        "certified_algebraic_rank_lower_bound": certified,
                        "conductor": record["conductor"],
                        "log_conductor": record["log_conductor"],
                    }
                )
        records.append(record)

    script = Path(__file__).resolve()
    histogram = Counter(str(record["maximum_stable_numerical_rank"]) for record in records)
    artifact: dict[str, Any] = {
        "schema_version": 1,
        "status": "complete bounded escalation of three max300 tranche-two leaders",
        "target": {
            "rank_at_least": 21,
            "strict_log_conductor_upper_bound": "182.72",
            "hits": target_hits,
        },
        "input": {
            "source_script": str(SOURCE_SCRIPT.relative_to(ROOT)),
            "source_script_sha256": EXPECTED_SOURCE_SCRIPT_SHA256,
            "source_artifact": str(SOURCE_ARTIFACT.relative_to(ROOT)),
            "source_artifact_sha256": EXPECTED_SOURCE_ARTIFACT_SHA256,
            "leader_identifiers": list(EXPECTED_LEADER_IDS),
        },
        "protocol": {
            "height_stages": list(HEIGHT_STAGES),
            "all_three_leaders_receive_every_stage_once": True,
            "no_retries_or_adaptive_broadening": True,
            "exact_mod3_trigger": 15,
        },
        "result": {
            "records": records,
            "maximum_stable_numerical_rank": max(
                record["maximum_stable_numerical_rank"] for record in records
            ),
            "rank_histogram": dict(sorted(histogram.items())),
            "exact_gain_certificate_count": sum(
                "certified_algebraic_rank_lower_bound" in record["exact_gain_certificate"]
                for record in records
            ),
            "target_hit": bool(target_hits),
        },
        "timing": {"wall_seconds": time.monotonic() - started},
        "provenance": {
            "script": str(script.relative_to(ROOT)),
            "script_sha256": sha256_file(script),
            "command": " ".join(
                shlex.quote(part) for part in [sys.executable, *sys.argv]
            ),
            "python": platform.python_version(),
            "owned_processes_remaining": 0,
        },
    }
    artifact["result_sha256"] = stable_json_digest(
        {
            "input": artifact["input"],
            "records": [
                [
                    record["identifier"],
                    [
                        [stage["height_bound"], stage["status"], stage.get("stable_numerical_rank")]
                        for stage in record["stages"]
                    ],
                    record["maximum_stable_numerical_rank"],
                    record["exact_gain_certificate"].get(
                        "certified_algebraic_rank_lower_bound"
                    ),
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
        f"max_rank={artifact['result']['maximum_stable_numerical_rank']}; "
        f"certificates={artifact['result']['exact_gain_certificate_count']}; "
        f"hits={len(target_hits)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
