#!/usr/bin/env python3
"""Pin the completed b=50001..500000 rank-30 denominator-sieve run.

The reusable sieve driver records its default command string even when called
with nondefault bounds.  This metadata-only normalizer preserves the exact raw
execution artifact, validates every declared boundary/count, and writes a
canonical artifact with the full two-command replay.  It never reruns or
extends the bounded search.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RAW = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "elliptic_elkies_klagsbrun_rank30_denominator_sieve_b500000_raw.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "elliptic_elkies_klagsbrun_rank30_denominator_sieve_b500000.json"
)
ENGINE = (
    ROOT
    / "elliptic-curves"
    / "cas"
    / "search_elkies_klagsbrun_rank30_denominator_sieve.py"
)
EXPECTED_ENGINE_SHA256 = (
    "73d46aa761e912feee18667a7ffea240db9a211ce9401ec73b8722fe54f09878"
)
EXPECTED_RAW_SHA256 = (
    "39b2ce7834d73cfd04f6dca82b40e64a40d9595069f7e56865d2015dee216457"
)
EXPECTED_PROCESSED_COUNT = 257_729_687_722
EXPECTED_SURVIVOR_COUNT = 858_583
EXPECTED_SURVIVOR_SHA256 = (
    "bfda86a2212043454ef884f7b9a6a198ef1645e499fd9c68d5c2e8874d8f19f9"
)


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def validate(data: dict[str, Any]) -> None:
    if data["parameters"]["denominator_interval"] != [50_001, 500_000]:
        raise AssertionError("the denominator interval changed")
    if data["parameters"]["nonzero_offset_interval"] != [-16_384, 16_384]:
        raise AssertionError("the offset interval changed")
    result = data["search_result"]
    if not result["search_complete"] or result["wall_cap_reached"]:
        raise AssertionError("the bounded search did not complete")
    if result["completed_denominator_interval"] != [50_001, 500_000]:
        raise AssertionError("the completed interval changed")
    if result["declared_primitive_candidate_count"] != EXPECTED_PROCESSED_COUNT:
        raise AssertionError("the declared primitive count changed")
    if result["processed_primitive_candidate_count"] != EXPECTED_PROCESSED_COUNT:
        raise AssertionError("the processed primitive count changed")
    if result["modular_survivor_count_after_primitivity"] != EXPECTED_SURVIVOR_COUNT:
        raise AssertionError("the modular survivor count changed")
    if result["modular_survivor_manifest_sha256"] != EXPECTED_SURVIVOR_SHA256:
        raise AssertionError("the survivor manifest changed")
    if result["exact_nonsquare_count_after_sieve"] != EXPECTED_SURVIVOR_COUNT:
        raise AssertionError("not every survivor was replayed as a nonsquare")
    if result["exact_square_abscissa_count"] != 0:
        raise AssertionError("the run unexpectedly found a square abscissa")
    if result["certified_independent_30th_point_count"] != 0:
        raise AssertionError("the run unexpectedly certified rank 30")
    if result["rank30_target_hit"]:
        raise AssertionError("the negative run is marked as a target hit")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw", type=Path, default=DEFAULT_RAW)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--capture-current-output-as-raw",
        action="store_true",
        help="copy the exact pre-normalization output to --raw before validation",
    )
    args = parser.parse_args()
    if sha256_file(ENGINE) != EXPECTED_ENGINE_SHA256:
        raise AssertionError("the bounded sieve engine changed")
    if args.capture_current_output_as_raw:
        current = args.output.read_bytes()
        if sha256_bytes(current) != EXPECTED_RAW_SHA256:
            raise AssertionError("the current output is not the pinned raw execution")
        args.raw.parent.mkdir(parents=True, exist_ok=True)
        args.raw.write_bytes(current)
    raw_bytes = args.raw.read_bytes()
    if sha256_bytes(raw_bytes) != EXPECTED_RAW_SHA256:
        raise AssertionError("the raw bounded execution artifact changed")
    data = json.loads(raw_bytes)
    validate(data)
    normalizer = Path(__file__).resolve()
    data["reproduction"] = {
        "bounded_execution_command": (
            "PYTHONPATH=elliptic-curves/cas .venv/bin/python -u "
            "elliptic-curves/cas/"
            "search_elkies_klagsbrun_rank30_denominator_sieve.py "
            "--denominator-min 50001 --denominator-max 500000 "
            "--offset-radius 16384 --wall-cap-seconds 900 "
            "--progress-every 10000 --output "
            "artifacts/generated-results/"
            "elliptic_elkies_klagsbrun_rank30_denominator_sieve_b500000_raw.json"
        ),
        "normalization_command": (
            "python3 elliptic-curves/cas/"
            "normalize_elkies_klagsbrun_rank30_denominator_sieve_b500000.py"
        ),
        "bounded_execution_engine": str(ENGINE.relative_to(ROOT)),
        "bounded_execution_engine_sha256": EXPECTED_ENGINE_SHA256,
        "raw_execution_artifact": str(args.raw.relative_to(ROOT)),
        "raw_execution_artifact_sha256": EXPECTED_RAW_SHA256,
        "normalizer": str(normalizer.relative_to(ROOT)),
        "normalizer_sha256": sha256_file(normalizer),
        "metadata_normalization_reran_bounded_search": False,
        "python": platform.python_version(),
    }
    data["claim_scope"]["new_region"] = (
        "For each of the 29 public x-coordinates, every primitive offset "
        "x=x_i+k/b^2 with 50001<=b<=500000 and 0<|k|<=16384; this is "
        "disjoint from the completed b<=50000 denominator-sieve region."
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"validated {EXPECTED_PROCESSED_COUNT} primitive abscissas")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
