#!/usr/bin/env python3
"""Sharpen the conditional rank diagnostic for the section-7 rank-20 curve.

Bober's sinc-squared explicit formula is evaluated at ``Delta=11/5``.  The
prime-power sum uses every prime through ``floor(exp(22*pi/5))``.  The common
archimedean contribution is evaluated from a closed dilogarithm identity and
cross-checked by two finite-interval quadratures whose omitted tail is bounded
from the convergent series for ``Re(psi(1+i*t))``.

The same direct implementation is first checked at ``Delta=2`` against
Bober's published E20 value.  As in the earlier diagnostic, the conclusion is
conditional: GRH makes the explicit-formula value an analytic-rank bound, and
BSD is additionally needed to identify analytic and algebraic rank.
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

from certify_nagao_rank20_t5081 import (
    EXPECTED_CONDUCTOR,
    EXPECTED_MINIMAL_MODEL,
    PARAMETER_T,
)
from ek_k3 import rational_to_string
from explicit_formula_rank_diagnostic import BOBER_SOURCE
from pari_bridge import pari_version


Q = Fraction
REPOSITORY = Path(__file__).resolve().parents[2]
DELTA = Q(11, 5)
EXPECTED_PRIME_LIMIT = 1_007_525
ARCHIMEDEAN_TRUNCATION = 32
QUADRATURE_DIGITS = 40
NUMERICAL_ALLOWANCE = mpmath.mpf("0.001")
PUBLISHED_E20_UPPER = mpmath.mpf("21.70")
RANK20_CERTIFICATE = (
    REPOSITORY
    / "artifacts/generated-results/elliptic-curves/elliptic_nagao_rank20_t5081_rank20_certificate.json"
)
RANK20_CERTIFICATE_SHA256 = (
    "466946076dc0c3fa02d0c5edd90b947d5ee3d10a4fb8cb16567049ab4380f88d"
)
DELTA2_DIAGNOSTIC = (
    REPOSITORY
    / "archive/elliptic-curves/artifacts/generated-results/elliptic_nagao_rank20_t5081_explicit_formula.json"
)
DELTA2_DIAGNOSTIC_SHA256 = (
    "2f799421101235c2092956a045c3b0e2cca0afccef26112297406d26e6432485"
)
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/explicit_formula_rank20_t5081_delta22.py"
)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prime_limit(delta: Fraction) -> int:
    with mpmath.workdps(80):
        return int(mpmath.floor(mpmath.exp(2 * mpmath.pi * Q(delta).numerator / Q(delta).denominator)))


def gp_program(
    model: Sequence[int], *, delta: Fraction, support_limit: int
) -> str:
    delta = Q(delta)
    vector = ",".join(str(int(value)) for value in model)
    return f"""default(realprecision,80);
D={delta.numerator}/{delta.denominator};
LIM={support_limit};
E=ellminimalmodel(ellinit([{vector}]));
N=ellglobalred(E)[1];
PS=0;
PSGOOD=0;
PSBAD=0;
NPC=0;
NTERMS=0;
MAXP=0;
gettime();
forprime(p=2,LIM,NPC++;MAXP=p;a=ellap(E,p);K=floor(2*Pi*D/log(p));if(valuation(E.disc,p)>0,sk=1;for(k=1,K,sk*=a;z=log(p)*sk/p^k*(1-k*log(p)/(2*Pi*D));PS+=z;PSBAD+=z;NTERMS++),s0=2;s1=a;for(k=1,K,if(k==1,sk=s1,sk=a*s1-p*s0;s0=s1;s1=sk);z=log(p)*sk/p^k*(1-k*log(p)/(2*Pi*D));PS+=z;PSGOOD+=z;NTERMS++)));
print("ROW|",N,"|",log(N),"|",PS,"|",PSGOOD,"|",PSBAD,"|",LIM,"|",NPC,"|",NTERMS,"|",MAXP,"|",gettime());
quit
"""


def parse_prime_sum(output: str) -> dict[str, str | int]:
    rows = [line for line in output.splitlines() if line.startswith("ROW|")]
    if len(rows) != 1:
        raise RuntimeError("PARI omitted or duplicated the prime-sum row")
    (
        _,
        conductor,
        log_conductor,
        prime_sum,
        good_prime_sum,
        bad_prime_sum,
        support,
        prime_count,
        term_count,
        maximum_prime,
        milliseconds,
    ) = rows[0].split("|")
    return {
        "conductor": conductor,
        "log_conductor": log_conductor,
        "prime_sum": prime_sum,
        "good_reduction_prime_sum": good_prime_sum,
        "bad_reduction_prime_sum": bad_prime_sum,
        "prime_limit": int(support),
        "prime_count": int(prime_count),
        "prime_power_term_count": int(term_count),
        "maximum_prime": int(maximum_prime),
        "pari_milliseconds_after_curve_setup": int(milliseconds),
    }


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
    if record["prime_limit"] != EXPECTED_PRIME_LIMIT:
        raise AssertionError("PARI used the wrong prime support")
    return record


def archimedean_integrand(t_value: mpmath.mpf, delta: mpmath.mpf) -> mpmath.mpf:
    if t_value == 0:
        test_value = mpmath.mpf(1)
    else:
        test_value = (
            mpmath.sin(mpmath.pi * delta * t_value)
            / (mpmath.pi * delta * t_value)
        ) ** 2
    return 2 * mpmath.re(mpmath.digamma(1 + 1j * t_value)) * test_value


def finite_archimedean_integral(
    delta: Fraction, *, method: str, truncation: int = ARCHIMEDEAN_TRUNCATION
) -> mpmath.mpf:
    if method not in {"tanh-sinh", "gauss-legendre"}:
        raise ValueError("unsupported quadrature method")
    delta_real = mpmath.mpf(Q(delta).numerator) / Q(delta).denominator
    return mpmath.fsum(
        mpmath.quad(
            lambda value: archimedean_integrand(value, delta_real),
            [interval, interval + 1],
            method=method,
        )
        for interval in range(truncation)
    )


def archimedean_tail_absolute_bound(
    delta: Fraction, *, truncation: int = ARCHIMEDEAN_TRUNCATION
) -> mpmath.mpf:
    """Bound the two-sided omitted tail beyond ``abs(t)=truncation``.

    From

      Re psi(1+it) = -gamma + sum_{n>=1} t^2/[n(n^2+t^2)],

    splitting at ``floor(t)`` gives ``abs(Re psi(1+it)) <= log(t)+3``
    for ``t>=2``.  Combining this with ``sinc^2 <= 1/(pi*Delta*t)^2``
    and integrating explicitly gives the returned bound.
    """

    if truncation < 2:
        raise ValueError("the tail estimate requires truncation at least two")
    delta_real = mpmath.mpf(Q(delta).numerator) / Q(delta).denominator
    truncation_real = mpmath.mpf(truncation)
    return (
        2
        / (mpmath.pi**2 * delta_real**2)
        * (mpmath.log(truncation_real) + 4)
        / truncation_real
    )


def archimedean_closed_form(delta: Fraction) -> tuple[mpmath.mpf, mpmath.mpf]:
    """Return the full digamma integral and its ``1/pi`` formula contribution.

    Combining the integral representation of ``Re psi(1+it)`` with the
    triangular Fourier transform of the sinc-squared test function gives

      A_D = 1/(pi*D) * [-gamma
            + (pi^2/6 - Li_2(exp(-2*pi*D)))/(2*pi*D)].

    Here ``A_D`` is the archimedean contribution already including the
    ``1/pi`` factor in Bober's equation (3).
    """

    delta_real = mpmath.mpf(Q(delta).numerator) / Q(delta).denominator
    contribution = 1 / (mpmath.pi * delta_real) * (
        -mpmath.euler
        + (
            mpmath.pi**2 / 6
            - mpmath.polylog(2, mpmath.exp(-2 * mpmath.pi * delta_real))
        )
        / (2 * mpmath.pi * delta_real)
    )
    return mpmath.pi * contribution, contribution


def explicit_formula_upper(
    *,
    log_conductor: str,
    prime_sum: str,
    delta: Fraction,
    archimedean_contribution: mpmath.mpf,
    numerical_allowance: mpmath.mpf = NUMERICAL_ALLOWANCE,
) -> mpmath.mpf:
    delta_real = mpmath.mpf(Q(delta).numerator) / Q(delta).denominator
    return (
        mpmath.mpf(log_conductor) / (2 * mpmath.pi * delta_real)
        - mpmath.log(2 * mpmath.pi) / (mpmath.pi * delta_real)
        + archimedean_contribution
        - mpmath.mpf(prime_sum) / (mpmath.pi * delta_real)
        + numerical_allowance
    )


def quadrature_record(delta: Fraction) -> dict[str, str]:
    with mpmath.workdps(QUADRATURE_DIGITS):
        tanh_sinh = finite_archimedean_integral(delta, method="tanh-sinh")
        gauss_legendre = finite_archimedean_integral(delta, method="gauss-legendre")
        disagreement = abs(tanh_sinh - gauss_legendre)
        if disagreement >= mpmath.mpf("1e-30"):
            raise AssertionError("the independent finite quadratures disagreed")
        tail = archimedean_tail_absolute_bound(delta)
        closed_integral, closed_contribution = archimedean_closed_form(delta)
        distance = abs(closed_integral - tanh_sinh)
        if distance > tail:
            raise AssertionError("the closed form escaped the analytic tail enclosure")
        return {
            "delta": rational_to_string(Q(delta)),
            "finite_interval": f"[-{ARCHIMEDEAN_TRUNCATION},{ARCHIMEDEAN_TRUNCATION}]",
            "tanh_sinh_value": mpmath.nstr(tanh_sinh, 38),
            "gauss_legendre_value": mpmath.nstr(gauss_legendre, 38),
            "quadrature_disagreement": mpmath.nstr(disagreement, 10),
            "analytic_absolute_tail_bound": mpmath.nstr(tail, 38),
            "closed_form_distance_from_truncated_integral": mpmath.nstr(distance, 38),
            "closed_form_full_integral": mpmath.nstr(closed_integral, 38),
            "closed_form_contribution_to_explicit_formula": mpmath.nstr(
                closed_contribution, 38
            ),
            "closed_form_inside_tail_enclosure": True,
        }


def build_artifact(args: argparse.Namespace) -> dict[str, Any]:
    if sha256_file(RANK20_CERTIFICATE) != RANK20_CERTIFICATE_SHA256:
        raise AssertionError("the pinned rank-20 certificate changed")
    if sha256_file(DELTA2_DIAGNOSTIC) != DELTA2_DIAGNOSTIC_SHA256:
        raise AssertionError("the pinned Delta=2 diagnostic changed")
    certificate = json.loads(RANK20_CERTIFICATE.read_text(encoding="utf-8"))
    if (
        certificate["candidate"]["root_number"] != 1
        or certificate["exact_rank_certificate"][
            "certified_algebraic_rank_lower_bound"
        ]
        != 20
    ):
        raise AssertionError("the pinned rank/root checkpoint changed")
    delta2 = json.loads(DELTA2_DIAGNOSTIC.read_text(encoding="utf-8"))
    candidate_prime_sum = run_prime_sum(
        timeout=args.timeout, stack_bytes=args.stack_bytes
    )

    delta2_quadrature = quadrature_record(Q(2))
    delta22_quadrature = quadrature_record(DELTA)
    with mpmath.workdps(60):
        calibration = delta2["curves"]["E20"]
        calibration_value = explicit_formula_upper(
            log_conductor=calibration["log_conductor"],
            prime_sum=calibration["prime_sum"],
            delta=Q(2),
            archimedean_contribution=mpmath.mpf(
                delta2_quadrature[
                    "closed_form_contribution_to_explicit_formula"
                ]
            ),
            numerical_allowance=mpmath.mpf(0),
        )
        calibration_upper = calibration_value + NUMERICAL_ALLOWANCE
        candidate_value = explicit_formula_upper(
            log_conductor=str(candidate_prime_sum["log_conductor"]),
            prime_sum=str(candidate_prime_sum["prime_sum"]),
            delta=DELTA,
            archimedean_contribution=mpmath.mpf(
                delta22_quadrature[
                    "closed_form_contribution_to_explicit_formula"
                ]
            ),
            numerical_allowance=mpmath.mpf(0),
        )
        candidate_upper = candidate_value + NUMERICAL_ALLOWANCE
        if not calibration_upper < PUBLISHED_E20_UPPER:
            raise AssertionError("the direct Delta=2 calibration missed Bober E20")
        if not candidate_upper < 22:
            raise AssertionError("the sharpened conditional bound is no longer below 22")

    return {
        "schema_version": 1,
        "status": "sharpened conditional explicit-formula diagnostic complete",
        "method": {
            "formula": "Bober equation (3), sinc-squared test function",
            "delta": rational_to_string(DELTA),
            "prime_limit": EXPECTED_PRIME_LIMIT,
            "prime_limit_definition": "floor(exp(2*pi*Delta))",
            "bober_source": BOBER_SOURCE,
            "quadrature_digits": QUADRATURE_DIGITS,
            "declared_final_numerical_allowance": str(NUMERICAL_ALLOWANCE),
            "archimedean_closed_form": (
                "1/(pi*Delta)*[-EulerGamma+(pi^2/6-"
                "Li_2(exp(-2*pi*Delta)))/(2*pi*Delta)]"
            ),
            "tail_bound_derivation": (
                "the digamma series gives |Re psi(1+it)|<=log(t)+3 for t>=2; "
                "combine with sinc^2<=1/(pi*Delta*t)^2 and integrate"
            ),
        },
        "input": {
            "rank20_certificate_path": str(RANK20_CERTIFICATE),
            "rank20_certificate_sha256": RANK20_CERTIFICATE_SHA256,
            "delta2_diagnostic_path": str(DELTA2_DIAGNOSTIC),
            "delta2_diagnostic_sha256": DELTA2_DIAGNOSTIC_SHA256,
            "constructor_parameter_T": rational_to_string(PARAMETER_T),
            "certified_algebraic_rank_lower_bound": 20,
            "root_number": 1,
        },
        "direct_delta2_calibration": {
            "curve": "Bober E20",
            "published_upper_value_rounded_up": str(PUBLISHED_E20_UPPER),
            "quadrature": delta2_quadrature,
            "direct_value_before_allowance": mpmath.nstr(calibration_value, 50),
            "direct_conservative_upper": mpmath.nstr(calibration_upper, 50),
            "passes_published_value": True,
        },
        "delta_11_over_5": {
            "prime_sum": candidate_prime_sum,
            "quadrature": delta22_quadrature,
            "explicit_formula_value_before_allowance": mpmath.nstr(
                candidate_value, 50
            ),
            "conservative_explicit_formula_upper": mpmath.nstr(candidate_upper, 50),
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
                "this diagnostic supplies no algebraic-rank upper bound; the "
                "unconditional statement remains rank at least 20"
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
            / "archive/elliptic-curves/artifacts/generated-results/elliptic_nagao_rank20_t5081_explicit_formula_delta22.json"
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
