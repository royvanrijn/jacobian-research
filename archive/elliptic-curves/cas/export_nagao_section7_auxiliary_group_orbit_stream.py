#!/usr/bin/env python3
"""Export the complete canonical T stream from the section-7 group orbits.

The main search artifact stores the top conductor proxies and a digest.  This
companion artifact stores all 4,251 canonical rational parameters so later
search lanes can perform exact population exclusions without repeating the
auxiliary group-law construction.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import platform
import shlex
import sys

from search_nagao_section7_auxiliary_group_orbits import (
    DEFAULT_DEV_GP,
    EXPECTED_DEV_GP_SHA256,
    EXPECTED_INPUT_SHA256,
    INPUT_ARTIFACT_RELATIVE,
    NAIVE_HEIGHT_EXCLUSION,
    ORBIT_SPECIFICATIONS,
    build_gp_orbit_program,
    file_sha256,
    h200000_parameters,
    load_inputs,
    projective_height,
    rational_string,
    run_gp_orbit,
    stream_sha256,
    validate_and_deduplicate_orbits,
)


OUTPUT_RELATIVE = Path(
    "artifacts/generated-results/elliptic_nagao_section7_auxiliary_group_orbit_stream.json"
)
EXPECTED_COUNT = 4251
EXPECTED_STREAM_SHA256 = (
    "b3209a0bf363b66367bbd1765b7e869f9e95260a5002f8911870313a81242b35"
)
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/export_nagao_section7_auxiliary_group_orbit_stream.py "
    "--gp /private/tmp/pari-map-src.33iJSU/pari/Odarwin-aarch64/gp-dyn"
)


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gp", type=Path, default=DEFAULT_DEV_GP)
    parser.add_argument("--orbit-timeout", type=float, default=60.0)
    parser.add_argument("--gp-stack-bytes", type=int, default=1_000_000_000)
    parser.add_argument("--output", type=Path, default=root / OUTPUT_RELATIVE)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if not 0 < args.orbit_timeout <= 60:
        raise SystemExit("--orbit-timeout must be in (0,60]")
    if not 256_000_000 <= args.gp_stack_bytes <= 1_000_000_000:
        raise SystemExit("--gp-stack-bytes must be in [256MB,1GB]")
    if not args.gp.is_file() or file_sha256(args.gp) != EXPECTED_DEV_GP_SHA256:
        raise SystemExit("the pinned development GP binary is absent or changed")

    root = Path(__file__).resolve().parents[2]
    input_artifact, priority = load_inputs(root)
    raw_orbits = {}
    orbit_records = []
    for specification in ORBIT_SPECIFICATIONS:
        item = priority[specification.priority_index_zero_based]
        known = h200000_parameters(
            input_artifact, specification.priority_index_zero_based
        )
        record, orbit = run_gp_orbit(
            args.gp,
            build_gp_orbit_program(specification, item, known),
            timeout=args.orbit_timeout,
            stack_bytes=args.gp_stack_bytes,
        )
        raw_orbits[specification.label] = orbit
        orbit_records.append(
            {
                "slice": specification.label,
                "status": record["status"],
                "coefficient_vector_count": record["coefficient_vector_count"],
                "raw_pullback_count": record["raw_pullback_count"],
                "wall_seconds": record["wall_seconds"],
                "timeout_seconds": record["timeout_seconds"],
                "one_attempt_no_retry": True,
            }
        )

    candidates, exclusions = validate_and_deduplicate_orbits(
        input_artifact, priority, raw_orbits
    )
    ordered = sorted(
        candidates,
        key=lambda value: (
            projective_height(value),
            value.numerator,
            value.denominator,
        ),
    )
    digest = stream_sha256(ordered)
    if len(ordered) != EXPECTED_COUNT or digest != EXPECTED_STREAM_SHA256:
        raise AssertionError("the canonical group-orbit stream changed")
    if any(projective_height(value) <= NAIVE_HEIGHT_EXCLUSION for value in ordered):
        raise AssertionError("the exported stream overlaps the naive H=200000 box")

    output = {
        "schema_version": 1,
        "status": "exact_canonical_parameter_stream_complete",
        "source_search": (
            "artifacts/generated-results/"
            "elliptic_nagao_section7_auxiliary_group_orbits.json"
        ),
        "input_accidental_slice_artifact": str(INPUT_ARTIFACT_RELATIVE),
        "input_accidental_slice_artifact_sha256": EXPECTED_INPUT_SHA256,
        "canonicalization": "T -> |T|",
        "naive_projective_height_exclusion": NAIVE_HEIGHT_EXCLUSION,
        "parameter_count": len(ordered),
        "parameter_stream_sha256": digest,
        "parameters": [rational_string(value) for value in ordered],
        "exclusion_counts": exclusions,
        "orbit_records": orbit_records,
        "reproduction": {
            "command": REPRODUCING_COMMAND,
            "actual_command": " ".join(
                shlex.quote(part) for part in [sys.executable, *sys.argv]
            ),
            "script_sha256": file_sha256(Path(__file__).resolve()),
            "group_orbit_module_sha256": file_sha256(
                Path(__file__).with_name(
                    "search_nagao_section7_auxiliary_group_orbits.py"
                )
            ),
        },
        "software": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "development_gp_path": str(args.gp),
            "development_gp_sha256": EXPECTED_DEV_GP_SHA256,
        },
        "bounded_process_policy": {
            "orbit_timeout_seconds": args.orbit_timeout,
            "orbit_stack_bytes": args.gp_stack_bytes,
            "one_attempt_no_retry_per_slice": True,
            "foreground_only": True,
            "no_detached_processes": True,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        f"wrote {args.output}: parameters={len(ordered)} sha256={digest}",
        flush=True,
    )


if __name__ == "__main__":
    main()
