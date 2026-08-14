#!/usr/bin/env python3
"""Sharpen the imported Fermigier rank-20 near miss at Delta=11/5.

The exact rank-lower-bound/conductor artifact imported under
``artifacts/generated-results/elliptic-curves`` is treated as immutable input.
We replay the conductor and Bober prime-power sum through every prime at most
``floor(exp(22*pi/5))``.  The conclusion is conditional exactly as in the
section-7 rank-20 diagnostic: GRH is needed for the analytic upper bound, and
BSD is additionally needed to turn it into an algebraic-rank equality.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
import platform
import shlex
import subprocess
import sys
from typing import Any

import mpmath

from explicit_formula_rank20_t5081_delta22 import (
    DELTA,
    DELTA2_DIAGNOSTIC,
    DELTA2_DIAGNOSTIC_SHA256,
    EXPECTED_PRIME_LIMIT,
    explicit_formula_upper,
    gp_program,
    parse_prime_sum,
    prime_limit,
    quadrature_record,
    sha256_file,
)
from explicit_formula_rank_diagnostic import BOBER_SOURCE
from pari_bridge import pari_version


Q = Fraction
REPOSITORY = Path(__file__).resolve().parents[2]
IMPORTED_CERTIFICATE = (
    REPOSITORY
    / "artifacts/generated-results/elliptic-curves/fermigier_rank20_near_miss_v1.json"
)
IMPORTED_CERTIFICATE_SHA256 = (
    "8416e835887236e9e4eafcb01384a710ce4f1be0628701a97f4a7d7a07fe63b1"
)
EXPECTED_PARAMETER = "28917/20"
EXPECTED_MINIMAL_MODEL = [
    1,
    1,
    1,
    -4437412060110743641525245114305,
    3586842216822165612930264910099076801587288127,
]
EXPECTED_CONDUCTOR = (
    2876153493562761211278364526603564191699143885403233935132057708367930
)
NUMERICAL_ALLOWANCE_TEXT = "0.001"
PUBLISHED_E20_UPPER_TEXT = "21.70"
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/explicit_formula_fermigier_rank20_28917_20_delta22.py"
)


def run_prime_sum(*, timeout: float, stack_bytes: int) -> dict[str, str | int]:
    if not 0 < timeout <= 120 or stack_bytes < 64_000_000:
        raise ValueError("invalid PARI resource bounds")
    support = prime_limit(DELTA)
    if support != EXPECTED_PRIME_LIMIT:
        raise AssertionError("the Delta=11/5 support cutoff changed")
    result = subprocess.run(
        ["gp", "-q", "-s", str(stack_bytes)],
        input=gp_program(
            EXPECTED_MINIMAL_MODEL, delta=DELTA, support_limit=support
        ),
        text=True,
        capture_output=True,
        timeout=timeout,
    )
    if result.returncode != 0 or "***" in result.stderr:
        raise RuntimeError(f"PARI explicit-formula sum failed: {result.stderr.strip()}")
    record = parse_prime_sum(result.stdout)
    if int(record["conductor"]) != EXPECTED_CONDUCTOR:
        raise AssertionError("the exact conductor replay changed")
    return record


def load_inputs() -> tuple[dict[str, Any], dict[str, Any]]:
    if sha256_file(IMPORTED_CERTIFICATE) != IMPORTED_CERTIFICATE_SHA256:
        raise AssertionError("the imported rank-20 certificate changed")
    if sha256_file(DELTA2_DIAGNOSTIC) != DELTA2_DIAGNOSTIC_SHA256:
        raise AssertionError("the pinned Delta=2 diagnostic changed")
    certificate = json.loads(IMPORTED_CERTIFICATE.read_text(encoding="utf-8"))
    curve = certificate["global_curve"]
    if certificate["family"]["adapter_parameter"] != EXPECTED_PARAMETER:
        raise AssertionError("the imported parameter changed")
    if curve["minimal_model"] != EXPECTED_MINIMAL_MODEL:
        raise AssertionError("the imported minimal model changed")
    if int(curve["conductor"]) != EXPECTED_CONDUCTOR or curve["root_number"] != 1:
        raise AssertionError("the imported conductor/root number changed")
    if certificate["point_cloud"]["selected_count"] != 20:
        raise AssertionError("the imported finite-reduction certificate changed")
    delta2 = json.loads(DELTA2_DIAGNOSTIC.read_text(encoding="utf-8"))
    return certificate, delta2


def build_artifact(args: argparse.Namespace) -> dict[str, Any]:
    certificate, delta2 = load_inputs()
    candidate_prime_sum = run_prime_sum(
        timeout=args.timeout, stack_bytes=args.stack_bytes
    )
    delta2_quadrature = quadrature_record(Q(2))
    delta22_quadrature = quadrature_record(DELTA)
    with mpmath.workdps(70):
        numerical_allowance = mpmath.mpf(NUMERICAL_ALLOWANCE_TEXT)
        published_e20_upper = mpmath.mpf(PUBLISHED_E20_UPPER_TEXT)
        calibration = delta2["curves"]["E20"]
        calibration_value = explicit_formula_upper(
            log_conductor=calibration["log_conductor"],
            prime_sum=calibration["prime_sum"],
            delta=Q(2),
            archimedean_contribution=mpmath.mpf(
                delta2_quadrature["closed_form_contribution_to_explicit_formula"]
            ),
            numerical_allowance=mpmath.mpf(0),
        )
        calibration_upper = calibration_value + numerical_allowance
        candidate_value = explicit_formula_upper(
            log_conductor=str(candidate_prime_sum["log_conductor"]),
            prime_sum=str(candidate_prime_sum["prime_sum"]),
            delta=DELTA,
            archimedean_contribution=mpmath.mpf(
                delta22_quadrature["closed_form_contribution_to_explicit_formula"]
            ),
            numerical_allowance=mpmath.mpf(0),
        )
        candidate_upper = candidate_value + numerical_allowance
        if not calibration_upper < published_e20_upper:
            raise AssertionError("the direct Delta=2 calibration missed Bober E20")
        if not candidate_upper < 22:
            raise AssertionError("the conditional upper bound is no longer below 22")

    return {
        "schema_version": 1,
        "status": "imported Fermigier rank-20 conditional closure complete",
        "method": {
            "formula": "Bober equation (3), sinc-squared test function",
            "delta": "11/5",
            "prime_limit": EXPECTED_PRIME_LIMIT,
            "prime_limit_definition": "floor(exp(2*pi*Delta))",
            "bober_source": BOBER_SOURCE,
            "declared_final_numerical_allowance": NUMERICAL_ALLOWANCE_TEXT,
        },
        "input": {
            "imported_certificate_path": str(IMPORTED_CERTIFICATE),
            "imported_certificate_sha256": IMPORTED_CERTIFICATE_SHA256,
            "adapter_parameter": EXPECTED_PARAMETER,
            "certified_algebraic_rank_lower_bound": 20,
            "minimal_model": EXPECTED_MINIMAL_MODEL,
            "conductor": str(EXPECTED_CONDUCTOR),
            "root_number": 1,
            "bounded_search_limitations": certificate["limitations"],
        },
        "direct_delta2_calibration": {
            "curve": "Bober E20",
            "published_upper_value_rounded_up": PUBLISHED_E20_UPPER_TEXT,
            "quadrature": delta2_quadrature,
            "direct_value_before_allowance": mpmath.nstr(calibration_value, 55),
            "direct_conservative_upper": mpmath.nstr(calibration_upper, 55),
            "passes_published_value": True,
        },
        "delta_11_over_5": {
            "prime_sum": candidate_prime_sum,
            "quadrature": delta22_quadrature,
            "explicit_formula_value_before_allowance": mpmath.nstr(
                candidate_value, 60
            ),
            "conservative_explicit_formula_upper": mpmath.nstr(
                candidate_upper, 60
            ),
            "strictly_less_than_22": True,
        },
        "interpretation": {
            "under_grh": (
                "the analytic rank is at most 21; root number +1 forces even "
                "analytic order, hence analytic rank at most 20"
            ),
            "under_bsd_and_grh": (
                "the exact algebraic lower bound 20 and conditional analytic "
                "upper bound 20 force algebraic and analytic rank exactly 20"
            ),
            "unconditional": (
                "this supplies no algebraic-rank upper bound; the imported "
                "unconditional statement remains rank at least 20"
            ),
            "search_priority": (
                "deprioritize additional fixed-fiber point searches and use the "
                "fiber as an anchor for parameter-neighborhood searches"
            ),
        },
        "declared_budget": {
            "pari_timeout_seconds": args.timeout,
            "pari_stack_bytes": args.stack_bytes,
        },
        "software": {
            "python": platform.python_version(),
            "pari_gp": pari_version(),
            "mpmath": mpmath.__version__,
        },
        "reproducing_command": REPRODUCING_COMMAND,
        "actual_command": " ".join(
            shlex.quote(part) for part in [sys.executable, *sys.argv]
        ),
        "script_sha256": sha256_file(Path(__file__).resolve()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--timeout", type=float, default=60.0)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            REPOSITORY
            / "artifacts/generated-results/elliptic_fermigier_rank20_28917_20_explicit_formula_delta22.json"
        ),
    )
    args = parser.parse_args()
    artifact = build_artifact(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(
        f"wrote {args.output}: conservative upper="
        f"{artifact['delta_11_over_5']['conservative_explicit_formula_upper']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
