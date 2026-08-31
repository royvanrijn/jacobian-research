#!/usr/bin/env python3
"""Specialize and certify the q12o5867 rootless rank-17 basis."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


REPOSITORY = Path(__file__).resolve().parents[2]
ELLIPTIC_ROOT = REPOSITORY / "elliptic-curves"
CAS = ELLIPTIC_ROOT / "cas"
sys.path.insert(0, str(ELLIPTIC_ROOT))
sys.path.insert(0, str(CAS))

from ecsearch.q12o5867_specialization import (  # noqa: E402
    build_specialization_record,
    load_q12o5867_data,
    normalize_projective_parameter,
)


DEFAULT_MODEL = REPOSITORY / "artifacts/local/elkies-k3/q12o5867-smooth-rr-qq.json"
DEFAULT_SECTIONS = (
    REPOSITORY
    / "artifacts/local/elkies-k3/q12o5867-rootless-selected-basis-qq.json"
)


def parse_relation_primes(text: str) -> tuple[int, ...]:
    try:
        values = tuple(int(value) for value in text.split(",") if value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("relation primes must be comma-separated integers") from error
    if not values:
        raise argparse.ArgumentTypeError("at least one relation prime is required")
    return values


def default_output(a: int, b: int) -> Path:
    a, b = normalize_projective_parameter(a, b)
    label = f"{a}_{b}".replace("-", "m")
    return (
        REPOSITORY
        / "artifacts/local/elliptic-curves/q12o5867-specializations"
        / f"q12o5867-specialization-{label}.json"
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Homogeneously specialize q12o5867 at the projective parameter "
            "(a:b), minimize exactly with PARI, and certify its 17 sections."
        )
    )
    parser.add_argument("--a", type=int, required=True, help="projective numerator")
    parser.add_argument("--b", type=int, required=True, help="projective denominator")
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--sections", type=Path, default=DEFAULT_SECTIONS)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--relation-primes", type=parse_relation_primes, default=(2, 3, 5))
    parser.add_argument("--reduction-prime-bound", type=int, default=500)
    parser.add_argument("--gp-timeout", type=float, default=300.0)
    parser.add_argument("--gp-stack-bytes", type=int, default=512_000_000)
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="replace the exact requested output path if it already exists",
    )
    args = parser.parse_args()
    if args.reduction_prime_bound < 3:
        parser.error("--reduction-prime-bound must be at least 3")
    if args.gp_timeout <= 0:
        parser.error("--gp-timeout must be positive")
    if args.gp_stack_bytes < 8_000_000:
        parser.error("--gp-stack-bytes must be at least 8000000")

    # The source artifacts contain decimal integers beyond Python's default
    # defensive conversion limit; this computation intentionally parses them.
    sys.set_int_max_str_digits(0)
    data = load_q12o5867_data(args.model.resolve(), args.sections.resolve())
    record = build_specialization_record(
        data,
        args.a,
        args.b,
        relation_primes=args.relation_primes,
        reduction_prime_bound=args.reduction_prime_bound,
        gp_timeout=args.gp_timeout,
        gp_stack_bytes=args.gp_stack_bytes,
    )
    output = (args.output or default_output(args.a, args.b)).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    mode = "w" if args.overwrite else "x"
    with output.open(mode) as handle:
        json.dump(record, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(f"status={record['status']}")
    print(f"parameter={record['parameter']['normalized_projective']}")
    print(f"output={output}")
    certificate = record.get("finite_quotient_independence")
    if certificate is not None:
        print(f"certified_independent={certificate['certified_independent']}")
        print(f"certified_rank_lower_bound={certificate['certified_rank_lower_bound']}")
        print(f"relation_prime={certificate['successful_relation_prime']}")
    if record["status"].startswith("BLOCKED"):
        print(f"blocker={record['blocker']}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()

