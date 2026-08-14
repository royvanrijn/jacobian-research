#!/usr/bin/env python3
"""Conditional Delta=11/5 explicit-formula closure for exact rank-16 T=62/35.

The input certificate proves algebraic rank at least 16 unconditionally.
This diagnostic applies the independently audited Bober sinc-squared formula
to every prime power in the Delta=11/5 support and adds the pinned 0.001
numerical allowance.  Under GRH the value is an analytic-rank upper bound;
BSD is additionally required to identify analytic and algebraic rank.  No
unconditional algebraic-rank upper bound is claimed.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import mpmath
import os
from pathlib import Path
import platform
import shlex
import subprocess
import sys
from typing import Any

from explicit_formula_rank20_t5081_delta22 import (
    archimedean_closed_form,
    explicit_formula_upper,
    gp_program,
    parse_prime_sum,
    prime_limit,
)
from pari_bridge import pari_version


Q = Fraction
ROOT = Path(__file__).resolve().parents[2]
CERTIFICATE = (
    ROOT
    / "artifacts/generated-results/"
    "elliptic_mestre_02557104116148_t62_35_rank16_certificate.json"
)
EXPECTED_CERTIFICATE_SHA256 = (
    "2c6d918546548227ac8f83287b3242e8d4261a98facd2665d506a8308f4c9fc7"
)
SOURCE_ENGINE = (
    ROOT / "elliptic-curves/cas/explicit_formula_rank20_t5081_delta22.py"
)
EXPECTED_SOURCE_ENGINE_SHA256 = (
    "c33da285159adb34bd14dcd2caf5ebdd77bddf5454d0eb226ca44da8543a80b7"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/"
    "elliptic_mestre_02557104116148_t62_35_explicit_formula_delta22.json"
)
DELTA = Q(11, 5)
EXPECTED_PRIME_LIMIT = 1_007_525
NUMERICAL_ALLOWANCE = "0.001"
EXPECTED_CONSERVATIVE_UPPER_PREFIX = "17.386729017908180130769670272467018284184848304"


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_prime_sum(
    model: list[int], *, timeout: float, stack_bytes: int
) -> dict[str, Any]:
    support = prime_limit(DELTA)
    if support != EXPECTED_PRIME_LIMIT:
        raise AssertionError("the Delta=11/5 support changed")
    process = subprocess.run(
        ["gp", "-q", "-s", str(stack_bytes)],
        input=gp_program(model, delta=DELTA, support_limit=support),
        text=True,
        capture_output=True,
        timeout=timeout,
    )
    if process.returncode != 0 or "***" in process.stderr:
        raise RuntimeError(
            f"PARI explicit-formula call failed: {process.stderr.strip()}"
        )
    return parse_prime_sum(process.stdout)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--timeout", type=float, default=120.0)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    return parser


def exclusive_write(path: Path, artifact: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    with os.fdopen(descriptor, "w") as stream:
        json.dump(artifact, stream, indent=2, sort_keys=True)
        stream.write("\n")


def main() -> None:
    args = build_parser().parse_args()
    if not 0 < args.timeout <= 180 or args.stack_bytes < 64_000_000:
        raise SystemExit("invalid explicit-formula resource bounds")
    if args.output.exists():
        raise SystemExit("refusing to overwrite the Delta=11/5 diagnostic")
    if file_sha256(CERTIFICATE) != EXPECTED_CERTIFICATE_SHA256:
        raise AssertionError("the exact rank-16 certificate changed")
    if file_sha256(SOURCE_ENGINE) != EXPECTED_SOURCE_ENGINE_SHA256:
        raise AssertionError("the audited explicit-formula engine changed")
    certificate = json.loads(CERTIFICATE.read_text())
    curve = certificate["curve"]
    if (
        certificate["theorem"]["certified_algebraic_rank_lower_bound"] != 16
        or curve["parameter_T"] != "62/35"
        or curve["root_number"] != 1
    ):
        raise AssertionError("the exact rank/root checkpoint changed")
    prime_sum = run_prime_sum(
        curve["minimal_model"],
        timeout=args.timeout,
        stack_bytes=args.stack_bytes,
    )
    if prime_sum["conductor"] != curve["conductor"]:
        raise AssertionError("the exact conductor replay changed")
    with mpmath.workdps(60):
        _, archimedean = archimedean_closed_form(DELTA)
        upper = explicit_formula_upper(
            log_conductor=str(prime_sum["log_conductor"]),
            prime_sum=str(prime_sum["prime_sum"]),
            delta=DELTA,
            archimedean_contribution=archimedean,
            numerical_allowance=mpmath.mpf(NUMERICAL_ALLOWANCE),
        )
        raw_value = upper - mpmath.mpf(NUMERICAL_ALLOWANCE)
        upper_text = mpmath.nstr(upper, 60)
        if not upper_text.startswith(EXPECTED_CONSERVATIVE_UPPER_PREFIX):
            raise AssertionError("the pinned conservative formula value changed")
        if upper >= 18:
            raise AssertionError("the conservative formula no longer closes below 18")
        numerical = {
            "archimedean_closed_form_contribution": mpmath.nstr(
                archimedean, 60
            ),
            "explicit_formula_value_before_allowance": mpmath.nstr(
                raw_value, 60
            ),
            "numerical_allowance": "0.001",
            "conservative_explicit_formula_upper": upper_text,
            "strict_margin_below_18": mpmath.nstr(18 - upper, 60),
        }
    script = Path(__file__).resolve()
    artifact = {
        "schema_version": 1,
        "status": "conditional fixed-fiber rank-16 closure",
        "target_hit": False,
        "curve": {
            "roots": curve["family_roots"],
            "parameter_T": curve["parameter_T"],
            "minimal_model": curve["minimal_model"],
            "conductor": curve["conductor"],
            "log_conductor": curve["log_conductor"],
            "root_number": curve["root_number"],
            "unconditional_algebraic_rank_lower_bound": 16,
        },
        "explicit_formula": {
            "formula": "Bober equation (3), sinc-squared test function",
            "delta": "11/5",
            "support_prime_limit": EXPECTED_PRIME_LIMIT,
            "prime_limit_definition": "floor(exp(2*pi*Delta))",
            **prime_sum,
            **numerical,
        },
        "conclusion": {
            "unconditional": (
                "The exact finite-reduction certificate proves algebraic rank "
                "at least 16; this diagnostic proves no rank upper bound."
            ),
            "under_grh": (
                "The conservative explicit-formula value is below 18. Root "
                "number +1 forces even analytic rank, hence analytic rank is "
                "at most 16."
            ),
            "under_bsd_and_grh": "The algebraic rank is exactly 16.",
            "search_priority": (
                "Conditionally close hidden-point searches on this fixed fiber "
                "and search disjoint nearby specializations instead."
            ),
        },
        "provenance": {
            "rank16_certificate": str(CERTIFICATE.relative_to(ROOT)),
            "rank16_certificate_sha256": EXPECTED_CERTIFICATE_SHA256,
            "audited_formula_engine": str(SOURCE_ENGINE.relative_to(ROOT)),
            "audited_formula_engine_sha256": EXPECTED_SOURCE_ENGINE_SHA256,
            "script": str(script.relative_to(ROOT)),
            "script_sha256": file_sha256(script),
            "command": " ".join(
                shlex.quote(part) for part in [sys.executable, *sys.argv]
            ),
            "python": platform.python_version(),
            "pari_gp": pari_version(),
            "single_foreground_pari_call": True,
            "no_retry": True,
        },
    }
    exclusive_write(args.output, artifact)
    print(
        "conservative upper="
        + artifact["explicit_formula"]["conservative_explicit_formula_upper"]
        + " < 18",
        flush=True,
    )
    print(f"wrote {args.output}", flush=True)


if __name__ == "__main__":
    main()
