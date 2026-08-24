#!/usr/bin/env python3
"""Calibrate a Delta=2 explicit-formula rank bound for the T=6793/64 lead.

The archimedean term in Bober's explicit formula is the same for every
elliptic curve over Q.  We therefore avoid a second numerical quadrature by
computing the *difference* between the new curve and Bober's published E20
calibration curve.  At fixed Delta this difference is

    (log(N_new)-log(N_ref))/(2*pi*Delta)
      - (prime_sum_new-prime_sum_ref)/(pi*Delta).

Bober reports an upper value 21.70 for E20 at Delta=2.  The prime sums here
run through floor(exp(2*pi*Delta)) exactly as in his formula.  The resulting
statement is conditional on GRH; identifying analytic and algebraic rank is
additionally conditional on BSD.  It is a diagnostic, never an unconditional
rank upper bound.
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
from typing import Any, Sequence

import mpmath

from nagao_1994 import RANK21_CONSTRUCTION, short_jacobian_coefficients
from pari_bridge import pari_version


Q = Fraction
DELTA = 2
PRIME_LIMIT = 286_751
REFERENCE_BOUND = mpmath.mpf("21.70")
PARAMETER_T = Q(6793, 64)
REPOSITORY = Path(__file__).resolve().parents[2]
RANK19_CERTIFICATE = (
    REPOSITORY
    / "artifacts/generated-results/elliptic_nagao_rank21_t6793_rank19_certificate.json"
)
RANK19_CERTIFICATE_SHA256 = (
    "0df1ea176d66c0446de16774ac5129253ade2379d6b6af9a591c444c8822bb6e"
)
E20_MODEL = (
    Q(1),
    Q(0),
    Q(0),
    Q(-431092980766333677958362095891166),
    Q(5156283555366643659035652799871176909391533088196),
)
T6793_MODEL = short_jacobian_coefficients(RANK21_CONSTRUCTION, PARAMETER_T)
BOBER_SOURCE = "https://antsmath.org/ANTSX/bober/paper.pdf"
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/explicit_formula_rank_diagnostic.py"
)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def gp_rational(value: Fraction) -> str:
    value = Q(value)
    if value.denominator == 1:
        return str(value.numerator)
    return f"({value.numerator}/{value.denominator})"


def good_prime_power_sums(trace: int, prime: int, count: int) -> tuple[int, ...]:
    """Return alpha^k+beta^k for alpha+beta=trace, alpha*beta=prime."""

    if prime < 2 or count < 0:
        raise ValueError("invalid prime-power request")
    if count == 0:
        return ()
    previous, current = 2, int(trace)
    answer = [current]
    for _ in range(2, count + 1):
        previous, current = current, trace * current - prime * previous
        answer.append(current)
    return tuple(answer)


def gp_program(models: Sequence[tuple[str, Sequence[Fraction]]]) -> str:
    lines = [
        "default(realprecision,80);",
        f"D={DELTA};",
        f"LIM={PRIME_LIMIT};",
    ]
    for label, coefficients in models:
        vector = ",".join(gp_rational(Q(value)) for value in coefficients)
        lines.extend(
            (
                f"E=ellminimalmodel(ellinit([{vector}]));",
                "N=ellglobalred(E)[1];",
                "PS=0;",
                (
                    "forprime(p=2,LIM,"
                    "a=ellap(E,p);K=floor(2*Pi*D/log(p));"
                    "if(valuation(E.disc,p)>0,"
                    "sk=1;for(k=1,K,sk*=a;"
                    "PS+=log(p)*sk/p^k*(1-k*log(p)/(2*Pi*D))),"
                    "s0=2;s1=a;for(k=1,K,"
                    "if(k==1,sk=s1,sk=a*s1-p*s0;s0=s1;s1=sk);"
                    "PS+=log(p)*sk/p^k*(1-k*log(p)/(2*Pi*D))))"
                    ");"
                ),
                f'print("ROW|{label}|",N,"|",log(N),"|",PS,"|",LIM);',
            )
        )
    lines.append("quit")
    return "\n".join(lines) + "\n"


def run_prime_sums(*, timeout: float, stack_bytes: int) -> dict[str, dict[str, Any]]:
    if timeout <= 0 or stack_bytes < 64_000_000:
        raise ValueError("invalid PARI resource bounds")
    models = (("E20", E20_MODEL), ("T6793", T6793_MODEL))
    result = subprocess.run(
        ["gp", "-q", "-s", str(stack_bytes)],
        input=gp_program(models),
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
    if set(records) != {"E20", "T6793"}:
        raise RuntimeError("PARI omitted an explicit-formula record")
    if any(record["prime_limit"] != PRIME_LIMIT for record in records.values()):
        raise AssertionError("the explicit-formula support cutoff changed")
    return records


def calibrated_bound(records: dict[str, dict[str, Any]]) -> dict[str, str]:
    mpmath.mp.dps = 80
    reference = records["E20"]
    candidate = records["T6793"]
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
    if sha256_file(RANK19_CERTIFICATE) != RANK19_CERTIFICATE_SHA256:
        raise AssertionError("the pinned rank-19 certificate changed")
    rank_certificate = json.loads(RANK19_CERTIFICATE.read_text(encoding="utf-8"))
    if rank_certificate["candidate"]["root_number"] != -1:
        raise AssertionError("the certified curve's root number changed")
    if (
        rank_certificate["exact_rank_certificate"][
            "certified_algebraic_rank_lower_bound"
        ]
        != 19
    ):
        raise AssertionError("the certified algebraic lower bound changed")
    records = run_prime_sums(timeout=args.timeout, stack_bytes=args.stack_bytes)
    bound = calibrated_bound(records)
    upper = mpmath.mpf(bound["calibrated_upper_value"])
    if not upper < 21:
        raise AssertionError("the calibrated value no longer lies below 21")
    script_path = Path(__file__).resolve()
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
            "source": BOBER_SOURCE,
        },
        "curves": records,
        "input": {
            "rank19_certificate_path": str(RANK19_CERTIFICATE),
            "rank19_certificate_sha256": RANK19_CERTIFICATE_SHA256,
            "certified_algebraic_rank_lower_bound": 19,
        },
        "comparison": {
            **bound,
            "calibrated_upper_value_less_than_21": True,
            "candidate_root_number": -1,
            "odd_functional_equation_forces_odd_analytic_order": True,
        },
        "interpretation": {
            "under_grh": (
                "the explicit-formula value bounds analytic rank; being below 21 "
                "and odd forces analytic rank at most 19"
            ),
            "under_bsd_and_grh": (
                "the independent algebraic lower bound 19 and the analytic upper "
                "bound 19 predict algebraic rank exactly 19"
            ),
            "unconditional": (
                "this computation supplies no algebraic-rank upper bound and does "
                "not weaken the exact rank-at-least-19 certificate"
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
        "actual_command": " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv]),
        "script_sha256": sha256_file(script_path),
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
            / "artifacts/generated-results/elliptic_nagao_rank21_t6793_explicit_formula.json"
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
