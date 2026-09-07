#!/usr/bin/env python3
"""Calibrated Delta=2 explicit-formula diagnostic for Nagao's rank-20 fiber.

The common archimedean term is eliminated by comparison with Bober's E20
calibration curve.  The prime-power sum is exact through
``floor(exp(4*pi))=286751``.  The result is conditional on GRH as an analytic
rank bound, and comparison with algebraic rank is additionally conditional on
BSD; it is never promoted to an unconditional rank upper bound.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import shlex
import subprocess
import sys
from typing import Any

import mpmath

from certify_nagao_rank20_t5081 import CONSTRUCTION, PARAMETER_T
from explicit_formula_rank_diagnostic import (
    BOBER_SOURCE,
    DELTA,
    E20_MODEL,
    PRIME_LIMIT,
    REFERENCE_BOUND,
    gp_program,
)
from nagao_1994 import PRIMARY_SOURCE, short_jacobian_coefficients
from pari_bridge import pari_version


Q = Fraction
REPOSITORY = Path(__file__).resolve().parents[2]
RANK20_CERTIFICATE = (
    REPOSITORY
    / "artifacts/generated-results/elliptic_nagao_rank20_t5081_rank20_certificate.json"
)
RANK20_CERTIFICATE_SHA256 = (
    "466946076dc0c3fa02d0c5edd90b947d5ee3d10a4fb8cb16567049ab4380f88d"
)
T5081_MODEL = short_jacobian_coefficients(CONSTRUCTION, PARAMETER_T)
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/explicit_formula_rank20_t5081.py"
)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_prime_sums(*, timeout: float, stack_bytes: int) -> dict[str, dict[str, Any]]:
    if timeout <= 0 or timeout > 120 or stack_bytes < 64_000_000:
        raise ValueError("invalid PARI resource bounds")
    result = subprocess.run(
        ["gp", "-q", "-s", str(stack_bytes)],
        input=gp_program((("E20", E20_MODEL), ("T5081", T5081_MODEL))),
        text=True,
        capture_output=True,
        timeout=timeout,
    )
    if result.returncode != 0 or "***" in result.stderr:
        raise RuntimeError(f"PARI explicit-formula sum failed: {result.stderr.strip()}")
    records: dict[str, dict[str, Any]] = {}
    for line in result.stdout.splitlines():
        if not line.startswith("ROW|"):
            continue
        _, label, conductor, log_conductor, prime_sum, limit = line.split("|")
        records[label] = {
            "conductor": conductor,
            "log_conductor": log_conductor,
            "prime_sum": prime_sum,
            "prime_limit": int(limit),
        }
    if set(records) != {"E20", "T5081"}:
        raise RuntimeError("PARI omitted an explicit-formula record")
    if any(record["prime_limit"] != PRIME_LIMIT for record in records.values()):
        raise AssertionError("the explicit-formula support cutoff changed")
    return records


def calibrated_bound(records: dict[str, dict[str, Any]]) -> dict[str, str]:
    mpmath.mp.dps = 80
    reference = records["E20"]
    candidate = records["T5081"]
    log_difference = (
        mpmath.mpf(candidate["log_conductor"])
        - mpmath.mpf(reference["log_conductor"])
    ) / (2 * mpmath.pi * DELTA)
    prime_difference = -(
        mpmath.mpf(candidate["prime_sum"])
        - mpmath.mpf(reference["prime_sum"])
    ) / (mpmath.pi * DELTA)
    total_difference = log_difference + prime_difference
    upper = REFERENCE_BOUND + total_difference
    return {
        "log_conductor_term_difference": mpmath.nstr(log_difference, 70),
        "prime_term_difference": mpmath.nstr(prime_difference, 70),
        "total_difference": mpmath.nstr(total_difference, 70),
        "calibrated_upper_value": mpmath.nstr(upper, 70),
    }


def build_artifact(args: argparse.Namespace) -> dict[str, Any]:
    if sha256_file(RANK20_CERTIFICATE) != RANK20_CERTIFICATE_SHA256:
        raise AssertionError("the pinned rank-20 certificate changed")
    certificate = json.loads(RANK20_CERTIFICATE.read_text(encoding="utf-8"))
    if certificate["candidate"]["root_number"] != 1:
        raise AssertionError("the certified curve's root number changed")
    if (
        certificate["exact_rank_certificate"][
            "certified_algebraic_rank_lower_bound"
        ]
        != 20
    ):
        raise AssertionError("the certified algebraic lower bound changed")

    records = run_prime_sums(timeout=args.timeout, stack_bytes=args.stack_bytes)
    if records["T5081"]["conductor"] != certificate["candidate"]["conductor"]:
        raise AssertionError("the explicit-formula conductor changed")
    bound = calibrated_bound(records)
    upper = mpmath.mpf(bound["calibrated_upper_value"])
    if not mpmath.mpf("22.20") < upper < mpmath.mpf("22.22"):
        raise AssertionError("the calibrated value left its pinned interval")

    return {
        "schema_version": 1,
        "status": "conditional explicit-formula diagnostic complete",
        "method": {
            "delta": DELTA,
            "prime_limit": PRIME_LIMIT,
            "prime_limit_definition": "floor(exp(2*pi*Delta))",
            "reference_curve": "Bober E20",
            "reference_published_upper_value": str(REFERENCE_BOUND),
            "common_archimedean_term_cancels_in_difference": True,
            "bober_source": BOBER_SOURCE,
            "nagao_source": PRIMARY_SOURCE,
        },
        "curves": records,
        "input": {
            "rank20_certificate_path": str(RANK20_CERTIFICATE),
            "rank20_certificate_sha256": RANK20_CERTIFICATE_SHA256,
            "certified_algebraic_rank_lower_bound": 20,
            "constructor_parameter_T": rational_to_string(PARAMETER_T),
        },
        "comparison": {
            **bound,
            "calibrated_upper_value_less_than_23": True,
            "candidate_root_number": 1,
            "even_functional_equation_forces_even_analytic_order": True,
        },
        "interpretation": {
            "under_grh": (
                "the calibrated value below 23 bounds analytic rank; the even "
                "functional-equation sign restricts the analytic rank to at most 22"
            ),
            "under_bsd_and_grh": (
                "the exact algebraic lower bound 20 and conditional analytic "
                "upper bound 22 leave ranks 20 and 22 as the parity-compatible cases"
            ),
            "unconditional": (
                "this computation supplies no algebraic-rank upper bound and does "
                "not weaken or strengthen the exact rank-at-least-20 certificate"
            ),
        },
        "software": {
            "python": platform.python_version(),
            "pari_gp": pari_version(),
            "mpmath": mpmath.__version__,
        },
        "declared_budget": {
            "timeout_seconds": args.timeout,
            "pari_stack_bytes": args.stack_bytes,
        },
        "reproducing_command": REPRODUCING_COMMAND,
        "actual_command": " ".join(
            shlex.quote(part) for part in [sys.executable, *sys.argv]
        ),
        "script_sha256": sha256_file(Path(__file__).resolve()),
    }


def rational_to_string(value: Fraction) -> str:
    value = Q(value)
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--timeout", type=float, default=60.0)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            REPOSITORY
            / "artifacts/generated-results/elliptic_nagao_rank20_t5081_explicit_formula.json"
        ),
    )
    args = parser.parse_args()
    artifact = build_artifact(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(
        f"wrote {args.output}: calibrated upper="
        f"{artifact['comparison']['calibrated_upper_value']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
