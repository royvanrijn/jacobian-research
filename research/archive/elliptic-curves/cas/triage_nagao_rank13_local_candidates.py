#!/usr/bin/env python3
"""H=50,000 exact-section triage for local-CRT Nagao candidates.

This is an incremental companion to ``triage_nagao_rank13_finalists.py``.  It
consumes named *integer* ``u`` candidates from the independently generated
local-CRT artifact, so new conductor-engineered candidates can be checkpointed
without repeating the expensive million-height finalist searches.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import shlex
import sys
from typing import Any, Sequence

from nagao_1994 import PRIMARY_SOURCE, rank13_base_parameter
from pari_bridge import pari_version
from search_extra_points import parse_precisions
from triage_nagao_rank13_finalists import (
    Finalist,
    exact_candidate_triage,
    parse_u_values,
)


Q = Fraction
DEFAULT_U = (118, 316)
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/triage_nagao_rank13_local_candidates.py"
)


def load_local_integer_candidates(
    path: Path, parameter_u_values: Sequence[int]
) -> tuple[Finalist, ...]:
    data = json.loads(path.read_text())
    records = data.get("final_conductor_records")
    if not isinstance(records, list):
        raise ValueError(f"{path} has no final_conductor_records list")
    by_u: dict[int, dict[str, Any]] = {}
    for record in records:
        try:
            parameter_u = Q(str(record["parameter_u"]))
        except (KeyError, ValueError, ZeroDivisionError):
            continue
        if parameter_u.denominator == 1:
            by_u[int(parameter_u)] = record
    answer = []
    for parameter_u in parameter_u_values:
        if parameter_u not in by_u:
            raise ValueError(f"integer u={parameter_u} is absent from {path.name}")
        record = by_u[parameter_u]
        parameter_t = Q(record["parameter_t"])
        if parameter_t != rank13_base_parameter(Q(parameter_u)):
            raise AssertionError(f"u={parameter_u} has an inconsistent T value")
        curve = record["curve"]
        answer.append(
            Finalist(
                parameter_u=parameter_u,
                parameter_t=parameter_t,
                score=str(record["score"]),
                last_numerical_prime=int(record["last_numerical_prime"]),
                log_conductor=str(curve["log_conductor"]),
                root_number=int(curve["root_number"]),
                source_artifact=path.name,
            )
        )
    return tuple(answer)


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=root
        / "artifacts"
        / "generated-results"
        / "elliptic_nagao_rank13_local_crt.json",
    )
    parser.add_argument("--u", type=parse_u_values, default=DEFAULT_U)
    parser.add_argument("--height-bound", type=int, default=50_000)
    parser.add_argument("--precisions", type=parse_precisions, default=(72, 120))
    parser.add_argument("--search-timeout", type=float, default=30.0)
    parser.add_argument("--height-timeout", type=float, default=30.0)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=root
        / "artifacts"
        / "generated-results"
        / "elliptic_nagao_rank13_local_candidate_triage.json",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if not args.u:
        raise SystemExit("--u must contain at least one integer parameter")
    if args.height_bound <= 0:
        raise SystemExit("--height-bound must be positive")
    if min(args.search_timeout, args.height_timeout) <= 0:
        raise SystemExit("timeouts must be positive")
    if args.stack_bytes < 8_000_000:
        raise SystemExit("the PARI stack bound is too small")

    candidates = load_local_integer_candidates(args.input, args.u)
    results = []
    for candidate in candidates:
        result, _ = exact_candidate_triage(
            candidate,
            height_bound=args.height_bound,
            precisions=args.precisions,
            search_timeout=args.search_timeout,
            height_timeout=args.height_timeout,
            stack_bytes=args.stack_bytes,
        )
        results.append(result)
        print(
            f"u={candidate.parameter_u} logN={candidate.log_conductor} "
            f"pool={result['frontier_stable_pool_numerical_rank']} "
            f"new={result['bounded_search']['new_distinct_jacobian_images']}",
            flush=True,
        )

    results.sort(
        key=lambda record: (
            -int(record["frontier_stable_pool_numerical_rank"]),
            Fraction(str(record["parameter_u"])),
        )
    )
    maximum = max(
        int(record["frontier_stable_pool_numerical_rank"]) for record in results
    )
    script_path = Path(__file__).resolve()
    engine_path = script_path.with_name("triage_nagao_rank13_finalists.py")
    command = " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv])
    artifact = {
        "schema_version": 1,
        "status": (
            "incremental bounded exact-point and numerical-height triage; no "
            "rank target hit is certified"
        ),
        "primary_source": PRIMARY_SOURCE,
        "input": {
            "path": str(args.input),
            "sha256": hashlib.sha256(args.input.read_bytes()).hexdigest(),
        },
        "method": {
            "parameters_u": list(args.u),
            "exact_affine_sections_per_candidate": 13,
            "split_infinity_covariant_limit_added": True,
            "quartic_naive_height_bound": args.height_bound,
            "height_decimal_precisions": list(args.precisions),
        },
        "summary": {
            "candidate_count": len(results),
            "maximum_stable_pool_numerical_rank": maximum,
            "candidates_at_maximum_numerical_rank": [
                int(record["parameter_u"])
                for record in results
                if int(record["frontier_stable_pool_numerical_rank"]) == maximum
            ],
            "interpretation": (
                "exact point membership plus precision-stable numerical height "
                "rank, not an exact independence or rank certificate"
            ),
        },
        "candidates": results,
        "software": {
            "python": platform.python_version(),
            "python_implementation": platform.python_implementation(),
            "pari_gp": pari_version(),
        },
        "parameters": {
            "search_timeout_seconds_per_candidate": args.search_timeout,
            "height_timeout_seconds_per_replay": args.height_timeout,
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
