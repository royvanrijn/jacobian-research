#!/usr/bin/env python3
"""Pin the completed b=50001..500000 companion-center rank-30 sieve.

The reusable sieve stores its default command even under nondefault bounds.
This metadata-only normalizer preserves the raw execution, verifies its exact
population and survivor manifest, and records the full bounded replay.  It
does not rerun or broaden the search.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts" / "generated-results"
DEFAULT_RAW = GENERATED / "elliptic_elkies_klagsbrun_rank30_companion_center_sieve_b500000_raw.json"
DEFAULT_OUTPUT = GENERATED / "elliptic_elkies_klagsbrun_rank30_companion_center_sieve_b500000.json"
ENGINE = ROOT / "elliptic-curves" / "cas" / "search_elkies_klagsbrun_rank30_companion_center_sieve.py"
EXPECTED_ENGINE_SHA256 = "89f63d500bdea8b3b6602d4321b2395f76740c0123e6583cb9adce8baf91e680"
EXPECTED_RAW_SHA256 = "1d1acf0356c13494d023eaf996787e32a684edb064a7a6865855aeb0e0a5a987"
EXPECTED_CENTER_COUNT = 32
EXPECTED_CENTER_SHA256 = "f37a7ab9e8c1607f9fec18d8a9fd6eac095b5526d8c9c391f838d70f33f815eb"
EXPECTED_PROCESSED_COUNT = 273_248_070_766
EXPECTED_SURVIVOR_COUNT = 1_222_922
EXPECTED_SURVIVOR_SHA256 = "69e4359560792639a84090faf2f5d790e557102e3550171c7dda3b8c0f729381"


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def validate(data: dict[str, Any]) -> None:
    if data["center_manifest"]["count"] != EXPECTED_CENTER_COUNT:
        raise AssertionError("the companion-center count changed")
    if data["center_manifest"]["sha256"] != EXPECTED_CENTER_SHA256:
        raise AssertionError("the companion-center manifest changed")
    if data["parameters"]["denominator_interval"] != [50_001, 500_000]:
        raise AssertionError("the denominator interval changed")
    if data["parameters"]["nonzero_offset_interval"] != [-16_384, 16_384]:
        raise AssertionError("the offset interval changed")
    result = data["search_result"]
    if not result["one_pass_no_retry"]:
        raise AssertionError("the execution was not one-pass")
    if not result["search_complete"] or result["wall_cap_reached"]:
        raise AssertionError("the bounded search did not complete")
    if result["completed_denominator_interval"] != [50_001, 500_000]:
        raise AssertionError("the completed interval changed")
    if result["declared_primitive_candidate_count"] != EXPECTED_PROCESSED_COUNT:
        raise AssertionError("the declared primitive count changed")
    if result["processed_primitive_candidate_count"] != EXPECTED_PROCESSED_COUNT:
        raise AssertionError("the processed primitive count changed")
    if result["modular_survivor_count_after_primitivity"] != EXPECTED_SURVIVOR_COUNT:
        raise AssertionError("the survivor count changed")
    if result["modular_survivor_manifest_sha256"] != EXPECTED_SURVIVOR_SHA256:
        raise AssertionError("the survivor manifest changed")
    if result["exact_nonsquare_count_after_sieve"] != EXPECTED_SURVIVOR_COUNT:
        raise AssertionError("not every survivor was replayed as a nonsquare")
    if result["exact_square_abscissa_count"] != 0:
        raise AssertionError("the run unexpectedly found a point")
    if result["rank30_target_hit"]:
        raise AssertionError("the negative run is marked as a hit")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw", type=Path, default=DEFAULT_RAW)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--capture-current-output-as-raw", action="store_true")
    args = parser.parse_args()
    if sha256_file(ENGINE) != EXPECTED_ENGINE_SHA256:
        raise AssertionError("the companion-center sieve engine changed")
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
            "elliptic-curves/cas/search_elkies_klagsbrun_rank30_companion_center_sieve.py "
            "--denominator-min 50001 --denominator-max 500000 "
            "--offset-radius 16384 --wall-cap-seconds 900 "
            "--progress-every 10000 --output artifacts/generated-results/"
            "elliptic_elkies_klagsbrun_rank30_companion_center_sieve_b500000_raw.json"
        ),
        "normalization_command": (
            "python3 elliptic-curves/cas/"
            "normalize_elkies_klagsbrun_rank30_companion_center_sieve_b500000.py"
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
        "Around each of 32 exact nonpublic subgroup-companion centers, every "
        "primitive x=x_i+k/b^2 with 50001<=b<=500000 and "
        "0<|k|<=16384; disjoint from the completed companion b<=50000 boxes."
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"validated {EXPECTED_PROCESSED_COUNT} primitive abscissas")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
