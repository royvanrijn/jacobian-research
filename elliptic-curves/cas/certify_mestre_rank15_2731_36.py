#!/usr/bin/env python3
"""Certify rank at least 15 for the Mestre fiber T=2731/36.

The sole discovery input is the pinned thirteen-family rational-frontier
artifact.  This script performs no parameter or point search.  It replays the
stored H=10^6 fifteen-point subset on the exact short Weierstrass curve, asks
PARI once for saturation at primes strictly below 20, checks every returned
point exactly, and proves independence of the saturated basis by exact
finite reductions modulo 3.  The independence theorem does not rely on
PARI's finite-index hypothesis for ``ellsaturation``.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import json
from math import factorial
import os
from pathlib import Path
import platform
import shlex
import sys
from typing import Any, Sequence

from extend_nagao_u42_frontier import saturate_exact_basis
from search_mestre_root_tuple_scale import (
    capped_minimal_curve_data,
    point_digest,
    point_on_short_curve,
    sha256_file,
)
from search_mestre_root_tuple_scale_max100 import stable_json_digest
from search_mestre_root_tuple_scale_max200 import mod3_independence_certificate


Q = Fraction
if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

ROOTS = (0, 7, 93, 154, 161, 191)
PARAMETER = Q(2731, 36)
STAGE = "H1000000"
FRONTIER_FILENAME = "elliptic_mestre_rank13_multifamily_rational.json"
EXPECTED_FRONTIER_SHA256 = (
    "0f664e937b9983bd7fa1cfb80269b5c734faddbf6d02dbe4dfca0e3b573ac41f"
)
EXPECTED_FRONTIER_RESULT_SHA256 = (
    "cd191a414ed35235169349e8b3a38929da2f494aa775a54ba7a6b5f52580576d"
)
EXPECTED_FRONTIER_SCRIPT_SHA256 = (
    "cc16c8aeaa2319eb009be5ea6e6b58f041ba646e147b6786ac7ace5e11895335"
)
EXPECTED_INPUT_POINT_SHA256 = (
    "423d2e5c126501d0b0ed12f1bf8fadff6584ffffb20bc4e147e07d132ce8810e"
)
EXPECTED_COEFFICIENTS = (
    Q(0), Q(0), Q(0),
    Q(-7243875957312210121029861327889, 80621568),
    Q(15439913501945598734077174101842443821189557113, 1880739938304),
)
SATURATION_BOUND = 20
SATURATION_TIMEOUT = 60.0
CONDUCTOR_TIMEOUT = 20.0
CERTIFICATE_PRIME_BOUND = 499
STACK_BYTES = 512_000_000
LOG_CONDUCTOR_TARGET = Q(18272, 100)
LOG_TEN_UPPER_BOUND = Q(231, 100)


def rational_string(value: Fraction) -> str:
    value = Q(value)
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def evaluate_polynomial(coefficients: Sequence[int], value: Fraction) -> Fraction:
    answer = Q(0)
    for coefficient in reversed(coefficients):
        answer = answer * value + coefficient
    return answer


def exact_log_conductor_certificate(conductor: int) -> dict[str, Any]:
    """Prove log(N)<182.72 by a rational decimal-digit estimate."""

    digits = len(str(conductor))
    if conductor <= 0 or conductor >= 10**digits:
        raise AssertionError("the exact decimal conductor bound failed")
    exponential_partial_sum = sum(
        LOG_TEN_UPPER_BOUND**degree / factorial(degree)
        for degree in range(8)
    )
    if exponential_partial_sum <= 10:
        raise AssertionError("the rational proof of log(10)<2.31 failed")
    upper = digits * LOG_TEN_UPPER_BOUND
    if upper >= LOG_CONDUCTOR_TARGET:
        raise AssertionError("the digit estimate does not prove the strict target")
    return {
        "conductor_less_than_power_of_ten": f"10^{digits}",
        "decimal_digit_count": digits,
        "exp_231_over_100_degree_7_partial_sum": str(exponential_partial_sum),
        "partial_sum_greater_than_10": True,
        "deduced_log_10_upper_bound": str(LOG_TEN_UPPER_BOUND),
        "deduced_log_conductor_upper_bound": str(upper),
        "strict_target_as_rational": str(LOG_CONDUCTOR_TARGET),
        "strict_target_proved_exactly": True,
    }


def exclusive_write(path: Path, artifact: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    with os.fdopen(descriptor, "w") as stream:
        json.dump(artifact, stream, indent=2, sort_keys=True)
        stream.write("\n")


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    generated = root / "artifacts/generated-results"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frontier", type=Path, default=generated / FRONTIER_FILENAME)
    parser.add_argument(
        "--output", type=Path,
        default=generated / "elliptic_mestre_rank15_2731_36.json",
    )
    parser.add_argument("--saturation-bound", type=int, default=SATURATION_BOUND)
    parser.add_argument("--saturation-timeout", type=float, default=SATURATION_TIMEOUT)
    parser.add_argument("--conductor-timeout", type=float, default=CONDUCTOR_TIMEOUT)
    parser.add_argument("--certificate-prime-bound", type=int, default=CERTIFICATE_PRIME_BOUND)
    parser.add_argument("--stack-bytes", type=int, default=STACK_BYTES)
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if (
        args.saturation_bound != SATURATION_BOUND
        or args.saturation_timeout != SATURATION_TIMEOUT
        or args.conductor_timeout != CONDUCTOR_TIMEOUT
        or args.certificate_prime_bound != CERTIFICATE_PRIME_BOUND
        or args.stack_bytes != STACK_BYTES
    ):
        raise SystemExit("the certificate resource bounds are pinned")


def main() -> None:
    args = build_parser().parse_args()
    validate_args(args)
    if args.output.exists():
        raise SystemExit("refusing to overwrite the rank-15 certificate artifact")
    script_path = Path(__file__).resolve()
    root = script_path.parents[2]
    frontier_script = script_path.with_name(
        "search_mestre_rank13_multifamily_rational.py"
    )
    saturation_engine = script_path.with_name("extend_nagao_u42_frontier.py")
    finite_engine = script_path.with_name("search_mestre_root_tuple_scale_max200.py")
    if sha256_file(args.frontier) != EXPECTED_FRONTIER_SHA256:
        raise AssertionError("the pinned multifamily frontier artifact changed")
    if sha256_file(frontier_script) != EXPECTED_FRONTIER_SCRIPT_SHA256:
        raise AssertionError("the pinned multifamily frontier script changed")
    frontier = json.loads(args.frontier.read_text())
    if frontier["result_sha256"] != EXPECTED_FRONTIER_RESULT_SHA256:
        raise AssertionError("the multifamily frontier result digest changed")

    family = next(
        row for row in frontier["families"] if tuple(row["roots"]) == ROOTS
    )
    a_coefficients = tuple(int(value) for value in family["A_coefficients_ascending"])
    b_coefficients = tuple(int(value) for value in family["B_coefficients_ascending"])
    coefficients = (
        Q(0), Q(0), Q(0),
        evaluate_polynomial(a_coefficients, PARAMETER),
        evaluate_polynomial(b_coefficients, PARAMETER),
    )
    if coefficients != EXPECTED_COEFFICIENTS:
        raise AssertionError("the exact short Weierstrass model changed")

    record = next(
        row for row in frontier["selected_records"]
        if tuple(row["roots"]) == ROOTS
        and row["numerator"] == PARAMETER.numerator
        and row["denominator"] == PARAMETER.denominator
    )
    stage = record["point_stages"][STAGE]
    if stage["status"] != "completed" or stage["stable_numerical_rank"] != 15:
        raise AssertionError("the frozen H1m leader changed")
    input_points = tuple(
        (Q(row["x"]), Q(row["y"])) for row in stage["numerical_subset"]
    )
    if len(input_points) != 15 or any(
        not point_on_short_curve(coefficients, point) for point in input_points
    ):
        raise AssertionError("the proposed exact input point set changed")
    if point_digest(input_points) != EXPECTED_INPUT_POINT_SHA256:
        raise AssertionError("the H1m subset digest changed")

    saturated, saturation = saturate_exact_basis(
        coefficients,
        input_points,
        prime_bound=args.saturation_bound,
        timeout=args.saturation_timeout,
        stack_bytes=args.stack_bytes,
    )
    if len(saturated) != 15 or any(
        not point_on_short_curve(coefficients, point) for point in saturated
    ):
        raise AssertionError("small-prime saturation did not return 15 exact points")
    certificate = mod3_independence_certificate(
        coefficients, saturated, prime_bound=args.certificate_prime_bound
    )
    if certificate["certified_algebraic_rank_lower_bound"] != 15:
        raise AssertionError("finite reductions did not certify all 15 points")

    source_conductor = record["conductor_phase"]
    replay = capped_minimal_curve_data(
        coefficients,
        timeout=args.conductor_timeout,
        stack_bytes=args.stack_bytes,
    )
    if (
        replay["conductor"] != source_conductor["conductor"]
        or replay["minimal_model"] != source_conductor["minimal_model"]
        or replay["minimal_discriminant"] != source_conductor["minimal_discriminant"]
        or replay["root_number"] != source_conductor["root_number"]
    ):
        raise AssertionError("the direct conductor/minimal-model replay changed")
    exact_log_bound = exact_log_conductor_certificate(int(replay["conductor"]))

    artifact: dict[str, Any] = {
        "schema_version": 1,
        "status": "certified exact algebraic rank lower bound 15 after small-prime saturation",
        "curve": {
            "roots": list(ROOTS),
            "parameter": rational_string(PARAMETER),
            "sign_equivalent_parameter": rational_string(-PARAMETER),
            "weierstrass_coefficients": [rational_string(value) for value in coefficients],
            "minimal_model": replay["minimal_model"],
            "minimal_discriminant": replay["minimal_discriminant"],
            "conductor": replay["conductor"],
            "log_conductor": replay["log_conductor"],
            "root_number": replay["root_number"],
            "strict_log_conductor_target": "182.72",
            "exact_log_conductor_bound": exact_log_bound,
        },
        "point_source": {
            "frontier_stage": STAGE,
            "height_bound": stage["height_bound"],
            "mapping_truncated": stage["mapping_truncated"],
            "input_point_count": len(input_points),
            "input_point_sha256": point_digest(input_points),
            "exact_curve_membership_replayed": True,
            "selection_statement": (
                "the certificate replays the fixed H1m subset and performs no "
                "new parameter or point search"
            ),
        },
        "input_points": [
            {"x": rational_string(x_value), "y": rational_string(y_value)}
            for x_value, y_value in input_points
        ],
        "small_prime_saturation": saturation,
        "saturated_basis": [
            {"x": rational_string(x_value), "y": rational_string(y_value)}
            for x_value, y_value in saturated
        ],
        "finite_reduction_certificate": certificate,
        "claim": {
            "certified_algebraic_rank_lower_bound": 15,
            "independence_uses_numerical_heights": False,
            "independence_depends_on_ellsaturation_finite_index_hypothesis": False,
            "does_not_claim_exact_mordell_weil_rank": True,
            "does_not_hit_rank21_or_rank30_target": True,
        },
        "provenance": {
            "script_path": str(script_path.relative_to(root)),
            "script_sha256": sha256_file(script_path),
            "frontier_path": str(args.frontier.relative_to(root)),
            "frontier_sha256": EXPECTED_FRONTIER_SHA256,
            "frontier_result_sha256": EXPECTED_FRONTIER_RESULT_SHA256,
            "frontier_script_sha256": EXPECTED_FRONTIER_SCRIPT_SHA256,
            "saturation_engine_path": str(saturation_engine.relative_to(root)),
            "saturation_engine_sha256": sha256_file(saturation_engine),
            "finite_reduction_engine_path": str(finite_engine.relative_to(root)),
            "finite_reduction_engine_sha256": sha256_file(finite_engine),
            "saturation_prime_bound_strict_upper_limit": SATURATION_BOUND,
            "certificate_prime_bound": CERTIFICATE_PRIME_BOUND,
            "direct_conductor_replay": True,
            "external_calls": 2,
            "reproducing_command": " ".join(
                shlex.quote(part) for part in [sys.executable, *sys.argv]
            ),
        },
        "software": {"python": platform.python_version()},
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    artifact["result_sha256"] = stable_json_digest(
        {
            "curve": artifact["curve"],
            "point_source": artifact["point_source"],
            "input_points": artifact["input_points"],
            "saturation": artifact["small_prime_saturation"],
            "saturated_basis": artifact["saturated_basis"],
            "certificate": artifact["finite_reduction_certificate"],
            "claim": artifact["claim"],
        }
    )
    exclusive_write(args.output, artifact)
    print(
        f"certified rank>={certificate['certified_algebraic_rank_lower_bound']} "
        f"at T={rational_string(PARAMETER)} saturated_sha={point_digest(saturated)} "
        f"output={args.output}",
        flush=True,
    )


if __name__ == "__main__":
    main()
