#!/usr/bin/env python3
"""Delta=11/5 explicit-formula diagnostic for the exact rank-15 T=2 fiber.

The exact finite-reduction certificate supplies algebraic rank at least 15,
conductor, minimal model, and root number -1.  This script evaluates Bober's
sinc-squared explicit formula using every prime power in its Delta=11/5
support and the audited closed-form archimedean term.  It adds the same 0.001
numerical allowance as the pinned rank-20 diagnostic.

Under GRH, a conservative value below 17 and root number -1 force analytic
rank at most 15.  BSD plus GRH would then identify exact algebraic rank 15.
No unconditional rank upper bound is claimed.  The artifact is still written
if the value fails to close below 17, so that the decision to search the fixed
fiber or its parameter neighborhood is explicit and reproducible.
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
    "elliptic_mestre_02136217261290_t2_rank15_certificate.json"
)
EXPECTED_CERTIFICATE_SHA256 = (
    "35abefefab42b19f49fad074f0c2cd65b039e8f36c398fbe7b46f68a0c2f09ea"
)
EXPECTED_CERTIFICATE_RESULT_SHA256 = (
    "c1de5071cf9ac8bb993345804bb0ab6f96656c72912c294f8e5fe097d002a77b"
)
SOURCE_ENGINE = ROOT / "elliptic-curves/cas/explicit_formula_rank20_t5081_delta22.py"
EXPECTED_SOURCE_ENGINE_SHA256 = (
    "c33da285159adb34bd14dcd2caf5ebdd77bddf5454d0eb226ca44da8543a80b7"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/"
    "elliptic_mestre_02136217261290_t2_explicit_formula_delta22.json"
)
DELTA = Q(11, 5)
EXPECTED_PRIME_LIMIT = 1_007_525
NUMERICAL_ALLOWANCE = mpmath.mpf("0.001")


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_prime_sum(model: list[int], *, timeout: float, stack_bytes: int) -> dict[str, Any]:
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
        raise RuntimeError(f"PARI explicit-formula call failed: {process.stderr.strip()}")
    return parse_prime_sum(process.stdout)


def exclusive_write(path: Path, artifact: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    with os.fdopen(descriptor, "w") as stream:
        json.dump(artifact, stream, indent=2, sort_keys=True)
        stream.write("\n")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--timeout", type=float, default=120.0)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if not 0 < args.timeout <= 120 or args.stack_bytes != 512_000_000:
        raise SystemExit("the explicit-formula resource bounds are pinned")
    if args.output.exists():
        raise SystemExit("refusing to overwrite the Delta=11/5 diagnostic")
    if file_sha256(CERTIFICATE) != EXPECTED_CERTIFICATE_SHA256:
        raise AssertionError("the rank-15 certificate changed")
    if file_sha256(SOURCE_ENGINE) != EXPECTED_SOURCE_ENGINE_SHA256:
        raise AssertionError("the audited explicit-formula engine changed")
    certificate = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    if certificate["result_sha256"] != EXPECTED_CERTIFICATE_RESULT_SHA256:
        raise AssertionError("the rank-15 certificate result changed")
    curve = certificate["curve"]
    if (
        certificate["theorem"]["certified_algebraic_rank_lower_bound"] != 15
        or curve["root_number"] != -1
    ):
        raise AssertionError("the certified lower bound or root number changed")

    prime_sum = run_prime_sum(
        curve["minimal_model"], timeout=args.timeout, stack_bytes=args.stack_bytes
    )
    if prime_sum["conductor"] != curve["conductor"]:
        raise AssertionError("the conductor replay changed")
    with mpmath.workdps(60):
        _, archimedean = archimedean_closed_form(DELTA)
        upper = explicit_formula_upper(
            log_conductor=str(prime_sum["log_conductor"]),
            prime_sum=str(prime_sum["prime_sum"]),
            delta=DELTA,
            archimedean_contribution=archimedean,
            numerical_allowance=NUMERICAL_ALLOWANCE,
        )
        raw_value = upper - NUMERICAL_ALLOWANCE
        closes_below_17 = bool(upper < 17)
        numerical = {
            "archimedean_closed_form_contribution": mpmath.nstr(archimedean, 50),
            "explicit_formula_value_before_allowance": mpmath.nstr(raw_value, 50),
            "numerical_allowance": "0.001",
            "conservative_explicit_formula_upper": mpmath.nstr(upper, 50),
            "strict_margin_below_17": mpmath.nstr(17 - upper, 50),
            "conservative_upper_strictly_below_17": closes_below_17,
        }

    if closes_below_17:
        under_grh = (
            "The conservative explicit-formula value is below 17. Root number -1 "
            "forces odd analytic rank, hence analytic rank is at most 15."
        )
        under_bsd_grh = "The algebraic rank is exactly 15."
        priority = (
            "Conditionally close the fixed fiber and prioritize disjoint nearby "
            "rational specializations."
        )
    else:
        under_grh = (
            "This Delta=11/5 value does not force analytic rank below 17; no "
            "fixed-fiber rank upper bound follows from this diagnostic."
        )
        under_bsd_grh = (
            "BSD plus GRH does not identify the exact algebraic rank from this "
            "diagnostic."
        )
        priority = "Prioritize a bounded hidden-point search on the fixed fiber."

    script = Path(__file__).resolve()
    artifact = {
        "schema_version": 1,
        "status": "conditional explicit-formula diagnostic complete",
        "target_hit": False,
        "curve": {
            "roots": curve["family_roots"],
            "parameter": curve["parameter_T"],
            "minimal_model": curve["minimal_model"],
            "conductor": curve["conductor"],
            "log_conductor": curve["log_conductor"],
            "root_number": curve["root_number"],
            "unconditional_algebraic_rank_lower_bound": 15,
        },
        "explicit_formula": {
            "delta": "11/5",
            "support_prime_limit": EXPECTED_PRIME_LIMIT,
            **prime_sum,
            **numerical,
        },
        "conclusion": {
            "unconditional": (
                "The exact finite-reduction certificate proves algebraic rank at "
                "least 15; this diagnostic proves no unconditional rank upper bound."
            ),
            "under_grh": under_grh,
            "under_bsd_and_grh": under_bsd_grh,
            "search_priority": priority,
        },
        "provenance": {
            "rank15_certificate": str(CERTIFICATE.relative_to(ROOT)),
            "rank15_certificate_sha256": EXPECTED_CERTIFICATE_SHA256,
            "rank15_certificate_result_sha256": EXPECTED_CERTIFICATE_RESULT_SHA256,
            "audited_formula_engine": str(SOURCE_ENGINE.relative_to(ROOT)),
            "audited_formula_engine_sha256": EXPECTED_SOURCE_ENGINE_SHA256,
            "script": str(script.relative_to(ROOT)),
            "script_sha256": file_sha256(script),
            "command": " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv]),
            "python": platform.python_version(),
            "pari_gp": pari_version(),
            "single_foreground_pari_call": True,
            "same_stage_retries": 0,
            "owned_processes_remaining": 0,
        },
    }
    exclusive_write(args.output, artifact)
    comparison = "< 17" if closes_below_17 else ">= 17"
    print(
        "upper="
        + artifact["explicit_formula"]["conservative_explicit_formula_upper"]
        + f" {comparison}",
        flush=True,
    )
    print(f"wrote {args.output}", flush=True)


if __name__ == "__main__":
    main()
