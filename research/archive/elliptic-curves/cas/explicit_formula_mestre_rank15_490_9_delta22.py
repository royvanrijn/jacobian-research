#!/usr/bin/env python3
"""Conditional explicit-formula closure for the exact rank-15 fiber T=490/9.

The exact finite-reduction certificate supplies algebraic rank at least 15,
conductor, minimal model, and root number -1.  Here Bober's sinc-squared
explicit formula is evaluated at Delta=11/5 using every prime power in its
support.  The closed-form archimedean term and the same conservative 0.001
numerical allowance used by the independently audited rank-20 diagnostic are
applied verbatim.

Under GRH the resulting value bounds analytic rank.  Since it is below 17 and
the root number is -1, analytic rank is then at most 15.  BSD plus GRH would
therefore identify the exact algebraic rank as 15.  No unconditional rank
upper bound is claimed.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import mpmath
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
    ROOT / "artifacts" / "generated-results" / "elliptic_mestre_rank15_490_9.json"
)
EXPECTED_CERTIFICATE_SHA256 = "50b2b9c8bd24bcb5533534446af6404f3a9a761b5f33e0e28e04dc572227f950"
SOURCE_ENGINE = ROOT / "elliptic-curves" / "cas" / "explicit_formula_rank20_t5081_delta22.py"
EXPECTED_SOURCE_ENGINE_SHA256 = "c33da285159adb34bd14dcd2caf5ebdd77bddf5454d0eb226ca44da8543a80b7"
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "elliptic_mestre_rank15_490_9_explicit_formula_delta22.json"
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


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--timeout", type=float, default=120.0)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if not 0 < args.timeout <= 180 or args.stack_bytes < 64_000_000:
        raise SystemExit("invalid explicit-formula resource bounds")
    if file_sha256(CERTIFICATE) != EXPECTED_CERTIFICATE_SHA256:
        raise AssertionError("the rank-15 certificate changed")
    if file_sha256(SOURCE_ENGINE) != EXPECTED_SOURCE_ENGINE_SHA256:
        raise AssertionError("the audited explicit-formula engine changed")
    certificate = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    curve = certificate["curve"]
    if certificate["claim"]["certified_algebraic_rank_lower_bound"] != 15:
        raise AssertionError("the algebraic lower bound changed")
    if curve["root_number"] != -1:
        raise AssertionError("the root number changed")
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
        if upper >= 17:
            raise AssertionError("the conservative formula no longer closes below 17")
        margin = 17 - upper
        numerical = {
            "archimedean_closed_form_contribution": mpmath.nstr(archimedean, 50),
            "explicit_formula_value_before_allowance": mpmath.nstr(raw_value, 50),
            "numerical_allowance": "0.001",
            "conservative_explicit_formula_upper": mpmath.nstr(upper, 50),
            "strict_margin_below_17": mpmath.nstr(margin, 50),
        }
    script = Path(__file__).resolve()
    artifact = {
        "schema_version": 1,
        "status": "conditional_fixed_fiber_rank15_closure",
        "target_hit": False,
        "curve": {
            "roots": curve["roots"],
            "parameter": curve["parameter"],
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
                "The exact finite-reduction certificate proves algebraic rank at least 15; "
                "this diagnostic proves no rank upper bound."
            ),
            "under_grh": (
                "The conservative explicit-formula value is below 17. Root number -1 "
                "forces odd analytic rank, hence analytic rank is at most 15."
            ),
            "under_bsd_and_grh": "The algebraic rank is exactly 15.",
            "search_priority": (
                "Deprioritize hidden-point searches on this fixed fiber and search nearby "
                "specializations instead."
            ),
        },
        "provenance": {
            "rank15_certificate": str(CERTIFICATE.relative_to(ROOT)),
            "rank15_certificate_sha256": EXPECTED_CERTIFICATE_SHA256,
            "audited_formula_engine": str(SOURCE_ENGINE.relative_to(ROOT)),
            "audited_formula_engine_sha256": EXPECTED_SOURCE_ENGINE_SHA256,
            "script": str(script.relative_to(ROOT)),
            "script_sha256": file_sha256(script),
            "command": " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv]),
            "python": platform.python_version(),
            "pari_gp": pari_version(),
            "single_foreground_pari_call": True,
            "no_retry": True,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    print(
        "upper="
        + artifact["explicit_formula"]["conservative_explicit_formula_upper"]
        + " < 17",
        flush=True,
    )
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
