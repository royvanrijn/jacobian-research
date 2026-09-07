#!/usr/bin/env python3
"""Conditional explicit-formula diagnostics for two section-7 rank-17 fibers.

The exact global-search artifact certifies rank at least 17 for ``T=599/2``
and ``T=426``.  This script compares both curves with Bober's E20 calibration
at ``Delta=2``.  Under GRH the resulting values are analytic-rank upper
bounds.  The values are below 19 and both curves have root number -1, so GRH
forces analytic rank at most 17.  Identifying algebraic and analytic rank also
uses BSD; no unconditional rank upper bound is claimed here.
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

from explicit_formula_rank_diagnostic import (
    BOBER_SOURCE,
    DELTA,
    E20_MODEL,
    PRIME_LIMIT,
    gp_program,
)
from nagao_1994 import PRIMARY_SOURCE, short_jacobian_coefficients
from nagao_1994_section7 import SECTION7_CONSTRUCTION
from pari_bridge import pari_version


Q = Fraction
REPOSITORY = Path(__file__).resolve().parents[2]
GLOBAL_ARTIFACT = (
    REPOSITORY
    / "artifacts/generated-results/elliptic_nagao_section7_global.json"
)
EXPECTED_GLOBAL_ARTIFACT_SHA256 = (
    "c86c2b39acfe278802d3b654e134d3031772013e984d81e2b78073eca1f53568"
)
REFERENCE_BOUND_TEXT = "21.70"
CANDIDATES = (("T599", Q(599, 2)), ("T426", Q(426)))
EXPECTED_BOUND_INTERVALS = {
    "T599": (mpmath.mpf("17.95"), mpmath.mpf("17.96")),
    "T426": (mpmath.mpf("18.73"), mpmath.mpf("18.74")),
}
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/explicit_formula_section7_rank17.py"
)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rational_to_string(value: Fraction) -> str:
    value = Q(value)
    return str(value.numerator) if value.denominator == 1 else str(value)


def load_exact_inputs() -> dict[str, dict[str, Any]]:
    if sha256_file(GLOBAL_ARTIFACT) != EXPECTED_GLOBAL_ARTIFACT_SHA256:
        raise AssertionError("the pinned section-7 global artifact changed")
    data = json.loads(GLOBAL_ARTIFACT.read_text(encoding="utf-8"))
    checkpoints = data["exact_checkpoints_stable_numerical_rank_at_least_17"]
    final_by_parameter = {
        row["constructor_parameter_T"]: row for row in data["final_frontier"]
    }
    answer = {}
    for label, parameter in CANDIDATES:
        identifier = (
            f"section7-global-{parameter.numerator}-{parameter.denominator}"
        )
        checkpoint = checkpoints[identifier]
        final = final_by_parameter[rational_to_string(parameter)]
        conductor = final["conductor"]
        if (
            checkpoint["status"] != "certified"
            or checkpoint["certified_algebraic_rank_lower_bound"] != 17
            or checkpoint["combined_exact_rank_over_F2"] != 17
            or int(conductor["root_number"]) != -1
        ):
            raise AssertionError(f"the exact rank-17 input changed for {label}")
        answer[label] = {
            "identifier": identifier,
            "constructor_parameter_T": rational_to_string(parameter),
            "conductor": str(conductor["conductor"]),
            "log_conductor": conductor["log_conductor"],
            "root_number": int(conductor["root_number"]),
            "certified_algebraic_rank_lower_bound": 17,
            "saturated_basis_sha256": checkpoint["saturated_point_sha256"],
        }
    return answer


def run_prime_sums(*, timeout: float, stack_bytes: int) -> dict[str, dict[str, Any]]:
    if timeout <= 0 or timeout > 120 or stack_bytes < 64_000_000:
        raise ValueError("invalid PARI resource bounds")
    curves = (
        ("E20", E20_MODEL),
        *(
            (
                label,
                short_jacobian_coefficients(SECTION7_CONSTRUCTION, parameter),
            )
            for label, parameter in CANDIDATES
        ),
    )
    result = subprocess.run(
        ["gp", "-q", "-s", str(stack_bytes)],
        input=gp_program(curves),
        text=True,
        capture_output=True,
        timeout=timeout,
    )
    if result.returncode != 0 or "***" in result.stderr:
        raise RuntimeError(
            f"PARI explicit-formula sum failed: {result.stderr.strip()}"
        )
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
    if set(records) != {"E20", "T599", "T426"}:
        raise RuntimeError("PARI omitted an explicit-formula record")
    if any(record["prime_limit"] != PRIME_LIMIT for record in records.values()):
        raise AssertionError("the explicit-formula support cutoff changed")
    return records


def calibrated_bounds(
    records: dict[str, dict[str, Any]]
) -> dict[str, dict[str, str]]:
    mpmath.mp.dps = 80
    reference_bound = mpmath.mpf(REFERENCE_BOUND_TEXT)
    reference = records["E20"]
    answer = {}
    for label, _ in CANDIDATES:
        candidate = records[label]
        log_difference = (
            mpmath.mpf(candidate["log_conductor"])
            - mpmath.mpf(reference["log_conductor"])
        ) / (2 * mpmath.pi * DELTA)
        prime_difference = -(
            mpmath.mpf(candidate["prime_sum"])
            - mpmath.mpf(reference["prime_sum"])
        ) / (mpmath.pi * DELTA)
        upper = reference_bound + log_difference + prime_difference
        lower_expected, upper_expected = EXPECTED_BOUND_INTERVALS[label]
        if not lower_expected < upper < upper_expected:
            raise AssertionError(f"the calibrated {label} value left its interval")
        answer[label] = {
            "log_conductor_term_difference": mpmath.nstr(log_difference, 70),
            "prime_term_difference": mpmath.nstr(prime_difference, 70),
            "calibrated_upper_value": mpmath.nstr(upper, 70),
        }
    return answer


def build_artifact(args: argparse.Namespace) -> dict[str, Any]:
    exact_inputs = load_exact_inputs()
    records = run_prime_sums(
        timeout=args.timeout, stack_bytes=args.stack_bytes
    )
    comparisons = calibrated_bounds(records)
    for label, _ in CANDIDATES:
        if records[label]["conductor"] != exact_inputs[label]["conductor"]:
            raise AssertionError(f"the explicit-formula conductor changed for {label}")
    return {
        "schema_version": 1,
        "status": "conditional explicit-formula diagnostics complete",
        "method": {
            "delta": DELTA,
            "prime_limit": PRIME_LIMIT,
            "prime_limit_definition": "floor(exp(2*pi*Delta))",
            "reference_curve": "Bober E20",
            "reference_published_upper_value": REFERENCE_BOUND_TEXT,
            "common_archimedean_term_cancels_in_difference": True,
            "bober_source": BOBER_SOURCE,
            "nagao_source": PRIMARY_SOURCE,
        },
        "curves": records,
        "exact_inputs": exact_inputs,
        "comparisons": comparisons,
        "interpretation": {
            "under_grh": (
                "both calibrated values are below 19; root number -1 forces "
                "odd analytic order, hence analytic rank at most 17"
            ),
            "under_bsd_and_grh": (
                "the exact algebraic lower bounds and analytic upper bounds "
                "force algebraic rank exactly 17 for both curves"
            ),
            "unconditional": (
                "the global artifact certifies rank at least 17 for both "
                "curves; this diagnostic supplies no algebraic-rank upper bound"
            ),
        },
        "input": {
            "global_artifact": str(GLOBAL_ARTIFACT),
            "global_artifact_sha256": EXPECTED_GLOBAL_ARTIFACT_SHA256,
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


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--timeout", type=float, default=60.0)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            REPOSITORY
            / "artifacts/generated-results/elliptic_nagao_section7_rank17_explicit_formula.json"
        ),
    )
    args = parser.parse_args()
    artifact = build_artifact(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(
        "wrote "
        f"{args.output}: T599={artifact['comparisons']['T599']['calibrated_upper_value']} "
        f"T426={artifact['comparisons']['T426']['calibrated_upper_value']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
