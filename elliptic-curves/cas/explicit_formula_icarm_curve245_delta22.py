#!/usr/bin/env python3
"""Conditional Delta=11/5 explicit-formula closure for ICARM curve 245.

This is not an unconditional rank upper bound.  Under GRH the sinc-squared
explicit formula bounds the analytic rank.  The exact root number then supplies
parity; BSD is additionally required to identify analytic and algebraic rank.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import subprocess

import mpmath

import icarm_curve245 as curve245
from explicit_formula_rank20_t5081_delta22 import (
    EXPECTED_PRIME_LIMIT,
    archimedean_closed_form,
    explicit_formula_upper,
    gp_program,
    parse_prime_sum,
    prime_limit,
    quadrature_record,
)
from pari_bridge import pari_version


Q = Fraction
DELTA = Q(11, 5)
ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "icarm_curve245_explicit_formula_delta22_v1.json"
)
HELPER = Path(__file__).with_name("explicit_formula_rank20_t5081_delta22.py")


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_gp(program: str, *, timeout: float, stack_bytes: int) -> str:
    process = subprocess.run(
        ["gp", "-q", "-s", str(stack_bytes)],
        input=program,
        text=True,
        capture_output=True,
        timeout=timeout,
    )
    if process.returncode != 0 or "***" in process.stderr:
        raise RuntimeError(f"PARI/GP failed: {process.stderr.strip()}")
    return process.stdout


def build_payload(*, timeout: float, stack_bytes: int) -> dict:
    if not 0 < timeout <= 120 or stack_bytes < 64_000_000:
        raise ValueError("invalid PARI resource bounds")
    support = prime_limit(DELTA)
    if support != EXPECTED_PRIME_LIMIT:
        raise AssertionError("the Delta=11/5 support cutoff changed")

    model = tuple(int(value) for value in curve245.GENERAL_WEIERSTRASS_COEFFICIENTS)
    prime_output = run_gp(
        gp_program(model, delta=DELTA, support_limit=support),
        timeout=timeout,
        stack_bytes=stack_bytes,
    )
    prime_record = parse_prime_sum(prime_output)
    if int(prime_record["conductor"]) != curve245.CONDUCTOR:
        raise AssertionError("the exact conductor replay changed")
    # Wall-clock timing is useful in the terminal but not deterministic data.
    prime_record.pop("pari_milliseconds_after_curve_setup", None)

    vector = ",".join(str(value) for value in model)
    root_output = run_gp(
        f'E=ellinit([{vector}]);print("ROOT|",ellrootno(E));quit\n',
        timeout=min(timeout, 30),
        stack_bytes=stack_bytes,
    )
    root_rows = [line for line in root_output.splitlines() if line.startswith("ROOT|")]
    if root_rows != ["ROOT|1"]:
        raise AssertionError("the exact root number is no longer +1")

    with mpmath.workdps(60):
        _integral, archimedean = archimedean_closed_form(DELTA)
        upper = explicit_formula_upper(
            log_conductor=str(prime_record["log_conductor"]),
            prime_sum=str(prime_record["prime_sum"]),
            delta=DELTA,
            archimedean_contribution=archimedean,
        )
        if not mpmath.mpf(20) < upper < mpmath.mpf(22):
            raise AssertionError("the conditional rank-20 enclosure changed")
        upper_text = mpmath.nstr(upper, 60)

    return {
        "schema_version": 1,
        "artifact_kind": "conditional_explicit_formula_rank_diagnostic",
        "status": "conditional_fixed_fiber_closure",
        "curve": {
            "source": "ICARM curve 245",
            "ainvs": [str(value) for value in model],
            "conductor": str(curve245.CONDUCTOR),
            "rank_lower_bound": 20,
            "exact_root_number": 1,
        },
        "explicit_formula": {
            "test_function": "Bober sinc-squared",
            "delta": "11/5",
            "conservative_upper": upper_text,
            "strictly_below_22": True,
            "prime_sum": prime_record,
            "archimedean_cross_check": quadrature_record(DELTA),
        },
        "conclusion": {
            "unconditional": "rank(E(Q)) >= 20 only; no rank upper bound",
            "under_grh": (
                "analytic rank <= 20: the explicit-formula upper is below 22 "
                "and root number +1 forces even analytic rank"
            ),
            "under_grh_and_bsd": "algebraic rank = 20",
            "search_consequence": (
                "redirect next-point work from this fixed fiber to nearby "
                "specializations or a different construction"
            ),
        },
        "reproduction": {
            "command": (
                "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
                "elliptic-curves/cas/explicit_formula_icarm_curve245_delta22.py --check"
            ),
            "pari_gp": list(pari_version()),
            "script_sha256": sha256_file(Path(__file__)),
            "helper_sha256": sha256_file(HELPER),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--timeout", type=float, default=120.0)
    parser.add_argument("--stack-bytes", type=int, default=1_000_000_000)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()

    payload = build_payload(timeout=args.timeout, stack_bytes=args.stack_bytes)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.write:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        if not args.output.exists() or args.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit("stale or missing explicit-formula artifact")
    print(
        "IC245EF|status=PASS|upper="
        f'{payload["explicit_formula"]["conservative_upper"]}'
        f"|mode={'write' if args.write else 'check'}"
    )


if __name__ == "__main__":
    main()
