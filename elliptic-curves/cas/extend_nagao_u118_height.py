#!/usr/bin/env python3
"""Single capped H=10^6 point/height extension for Nagao u=118.

There is exactly one ``hyperellratpoints`` attempt, capped at 60 seconds.  All
nonvisible points are mapped and checked exactly if it completes, and the
height matrix is replayed at two precisions.  This script performs no
``ellrank`` or saturation computation and makes no retry.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform
import shlex
import subprocess
import sys

from nagao_1994 import PRIMARY_SOURCE
from pari_bridge import pari_version
from search_extra_points import parse_precisions
from triage_nagao_rank13_finalists import exact_candidate_triage
from triage_nagao_rank13_local_candidates import load_local_integer_candidates


PARAMETER_U = 118


def validate_checkpoint(path: Path, *, expected_parameter_t: str) -> None:
    """Pin the lower-height result that justified this one-shot extension."""

    data = json.loads(path.read_text())
    record = next(
        (
            candidate
            for candidate in data.get("candidates", ())
            if int(candidate.get("parameter_u", -1)) == PARAMETER_U
        ),
        None,
    )
    if record is None:
        raise ValueError("the checkpoint artifact has no u=118 record")
    if record.get("parameter_t") != expected_parameter_t:
        raise ValueError("the checkpoint u=118 parameter changed")
    if int(record.get("frontier_stable_pool_numerical_rank", -1)) != 15:
        raise ValueError("the checkpoint u=118 numerical-rank frontier changed")


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    generated = root / "artifacts" / "generated-results"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=generated / "elliptic_nagao_rank13_local_crt.json",
    )
    parser.add_argument(
        "--checkpoint-input",
        type=Path,
        default=generated / "elliptic_nagao_rank13_local_candidate_triage.json",
    )
    parser.add_argument("--height-bound", type=int, default=1_000_000)
    parser.add_argument("--search-timeout", type=float, default=60.0)
    parser.add_argument("--height-timeout", type=float, default=30.0)
    parser.add_argument("--precisions", type=parse_precisions, default=(72, 120))
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=generated / "elliptic_nagao_u118_height_1000000.json",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.height_bound != 1_000_000:
        raise SystemExit("this one-shot extension pins --height-bound=1000000")
    if args.search_timeout <= 0 or args.search_timeout > 60:
        raise SystemExit("--search-timeout must be in (0,60]")
    if args.height_timeout <= 0:
        raise SystemExit("--height-timeout must be positive")
    if args.stack_bytes < 8_000_000:
        raise SystemExit("the PARI stack bound is too small")

    candidate = load_local_integer_candidates(args.input, (PARAMETER_U,))[0]
    validate_checkpoint(
        args.checkpoint_input, expected_parameter_t=str(candidate.parameter_t)
    )
    try:
        result, _ = exact_candidate_triage(
            candidate,
            height_bound=args.height_bound,
            precisions=args.precisions,
            search_timeout=args.search_timeout,
            height_timeout=args.height_timeout,
            stack_bytes=args.stack_bytes,
        )
        extension = {"status": "completed", "result": result}
        print(
            f"u=118 H={args.height_bound} "
            f"rank={result['frontier_stable_pool_numerical_rank']} "
            f"new={result['bounded_search']['new_distinct_jacobian_images']}",
            flush=True,
        )
    except subprocess.TimeoutExpired:
        extension = {
            "status": "timeout",
            "height_bound": args.height_bound,
            "timeout_seconds": args.search_timeout,
            "interpretation": (
                "the sole declared H=1000000 enumeration attempt timed out; "
                "no bounded-search conclusion is drawn"
            ),
        }
        print("u=118 H=1000000 search timed out", flush=True)
    except RuntimeError as error:
        extension = {
            "status": "pari_error",
            "height_bound": args.height_bound,
            "timeout_seconds": args.search_timeout,
            "error": str(error)[:1000],
        }
        print("u=118 H=1000000 search failed", flush=True)

    script_path = Path(__file__).resolve()
    engine_path = script_path.with_name("triage_nagao_rank13_finalists.py")
    command = " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv])
    artifact = {
        "schema_version": 1,
        "status": "single capped point/height extension; no descent or saturation",
        "primary_source": PRIMARY_SOURCE,
        "candidate": {
            "parameter_u": PARAMETER_U,
            "parameter_t": str(candidate.parameter_t),
            "log_conductor": candidate.log_conductor,
            "root_number": candidate.root_number,
        },
        "inputs": [
            {
                "path": str(path),
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
            for path in (args.input, args.checkpoint_input)
        ],
        "height_1000000_extension": extension,
        "target_status": {
            "rank21_log_conductor_target_certified": False,
            "rank30_target_certified": False,
        },
        "interpretation": (
            "exact point membership and bounded enumeration are exact; any "
            "reported height rank is numerical evidence only"
        ),
        "software": {
            "python": platform.python_version(),
            "python_implementation": platform.python_implementation(),
            "pari_gp": pari_version(),
        },
        "parameters": {
            "height_bound": args.height_bound,
            "height_precisions": list(args.precisions),
            "search_timeout_seconds": args.search_timeout,
            "height_timeout_seconds": args.height_timeout,
            "pari_stack_bytes": args.stack_bytes,
            "output": str(args.output),
        },
        "reproducing_command": command,
        "script_sha256": hashlib.sha256(script_path.read_bytes()).hexdigest(),
        "triage_engine_sha256": hashlib.sha256(engine_path.read_bytes()).hexdigest(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
