#!/usr/bin/env python3
"""Certify four Nagao specializations of rank at least 17 below the target.

The input candidates were selected by numerical height matrices.  Selection
is not certification.  For each candidate this script asks PARI only for a
small-prime saturation candidate, checks every returned point exactly, and
then proves independence with finite reduction signatures.  Consequently the
rank lower bounds do not depend on numerical heights or on the finite-index
hypothesis in the documentation of ``ellsaturation``.
"""

from __future__ import annotations

import argparse
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
from math import factorial
from pathlib import Path
import platform
import shlex
import sys
from typing import Any

from extend_nagao_u42_frontier import saturate_exact_basis
from mod2_reduction_independence import (
    combined_mod2_rank,
    find_mod2_reduction_certificate,
    find_two_torsion_certificate_prime,
)
from nagao_1994 import rank13_base_changed_short_jacobian_coefficients
from pari_bridge import minimal_curve_data, pari_version
from triage_nagao_rank13_finalists import point_digest, point_on_short_curve


Q = Fraction
CANDIDATE_PARAMETERS = (Q(135, 2), Q(471, 11), Q(42), Q(74))
LOG_CONDUCTOR_TARGET = Decimal("182.72")
LOG_CONDUCTOR_TARGET_RATIONAL = Q(18272, 100)
LOG_TEN_UPPER_BOUND = Q(231, 100)


def exact_log_conductor_certificate(conductor: int) -> dict[str, Any]:
    """Prove ``log(conductor) < 182.72`` using rational inequalities only."""

    if conductor <= 0:
        raise ValueError("the conductor must be positive")
    decimal_digits = len(str(conductor))
    if not conductor < 10**decimal_digits:
        raise AssertionError("the decimal digit bound failed")

    # The positive Taylor series gives exp(231/100) greater than this partial
    # sum.  Once the partial sum exceeds 10, monotonicity gives
    # log(10)<231/100 without any floating-point logarithm.
    exponential_partial_sum = sum(
        (LOG_TEN_UPPER_BOUND**degree) / factorial(degree)
        for degree in range(8)
    )
    if exponential_partial_sum <= 10:
        raise AssertionError("the rational exponential lower bound failed")
    rational_log_upper_bound = decimal_digits * LOG_TEN_UPPER_BOUND
    if rational_log_upper_bound >= LOG_CONDUCTOR_TARGET_RATIONAL:
        raise AssertionError("the exact logarithmic conductor bound is too weak")
    return {
        "conductor_less_than_power_of_ten": f"10^{decimal_digits}",
        "decimal_digit_count": decimal_digits,
        "exp_231_over_100_degree_7_partial_sum": str(exponential_partial_sum),
        "partial_sum_greater_than_10": True,
        "deduced_log_10_upper_bound": str(LOG_TEN_UPPER_BOUND),
        "deduced_log_conductor_upper_bound": str(rational_log_upper_bound),
        "strict_target_as_rational": str(LOG_CONDUCTOR_TARGET_RATIONAL),
        "strict_target_proved_exactly": True,
    }


def load_candidates(paths: tuple[Path, ...]) -> dict[Fraction, dict[str, Any]]:
    records: dict[Fraction, dict[str, Any]] = {}
    for path in paths:
        data = json.loads(path.read_text())
        for record in data["escalation_box"]["records"]:
            parameter = Q(record["parameter_u"])
            if parameter in records:
                continue
            copied = dict(record)
            copied["certificate_search_source"] = path.name
            records[parameter] = copied
    missing = set(CANDIDATE_PARAMETERS) - set(records)
    if missing:
        raise AssertionError(f"rank-gain artifact omitted candidates {sorted(missing)}")
    return {parameter: records[parameter] for parameter in CANDIDATE_PARAMETERS}


def input_subset(record: dict[str, Any]) -> tuple[tuple[Fraction, Fraction], ...]:
    points = tuple(
        (Q(item["jacobian_x"]), Q(item["jacobian_y"]))
        for item in record["explicit_numerically_independent_subset"]
    )
    if len(points) != 17 or int(record["stable_pool_numerical_rank"]) != 17:
        raise AssertionError("a frontier input does not contain its declared 17-point subset")
    return points


def certify_candidate(
    parameter_u: Fraction,
    record: dict[str, Any],
    *,
    saturation_bound: int,
    saturation_timeout: float,
    stack_bytes: int,
    certificate_prime_bound: int,
    conductor_timeout: float,
) -> dict[str, Any]:
    coefficients = rank13_base_changed_short_jacobian_coefficients(parameter_u)
    points = input_subset(record)
    if any(not point_on_short_curve(coefficients, point) for point in points):
        raise AssertionError("an input point failed exact membership")
    saturated, saturation = saturate_exact_basis(
        coefficients,
        points,
        prime_bound=saturation_bound,
        timeout=saturation_timeout,
        stack_bytes=stack_bytes,
    )
    if len(saturated) != 17:
        raise AssertionError("small-prime saturation changed the point count")
    if any(not point_on_short_curve(coefficients, point) for point in saturated):
        raise AssertionError("a returned saturation point failed exact membership")

    signatures = find_mod2_reduction_certificate(
        coefficients, saturated, prime_bound=certificate_prime_bound
    )
    binary_rank = combined_mod2_rank(signatures, len(saturated))
    if binary_rank != 17:
        raise AssertionError("the bounded finite-reduction search did not certify rank 17")
    two_torsion_prime = find_two_torsion_certificate_prime(
        coefficients, prime_bound=200
    )

    conductor = record["conductor_probe"]
    if conductor["status"] != "completed":
        raise AssertionError("the candidate lacks a completed conductor computation")
    log_conductor = Decimal(str(conductor["log_conductor"]))
    if log_conductor >= LOG_CONDUCTOR_TARGET:
        raise AssertionError("the candidate does not meet the strict conductor target")

    # Replay the exact conductor from the rational short model.  The discovery
    # artifact remains provenance for candidate selection, but is not trusted
    # as the sole source of the conductor claim in this certificate.
    conductor_replay = minimal_curve_data(
        coefficients,
        timeout=conductor_timeout,
        stack_bytes=stack_bytes,
    )
    if conductor_replay["conductor"] != int(conductor["conductor"]):
        raise AssertionError("the directly replayed conductor changed")
    if tuple(conductor_replay["minimal_model"]) != tuple(conductor["minimal_model"]):
        raise AssertionError("the directly replayed minimal model changed")
    if conductor_replay["root_number"] != int(conductor["root_number"]):
        raise AssertionError("the directly replayed root number changed")
    exact_log_bound = exact_log_conductor_certificate(conductor_replay["conductor"])

    return {
        "parameter_u": str(parameter_u),
        "parameter_t": str(Q(record["parameter_t"])),
        "short_weierstrass_coefficients": [str(value) for value in coefficients],
        "minimal_model": conductor["minimal_model"],
        "conductor": str(conductor["conductor"]),
        "log_conductor": str(conductor["log_conductor"]),
        "root_number": int(conductor["root_number"]),
        "below_strict_log_conductor_target": True,
        "exact_log_conductor_bound": exact_log_bound,
        "direct_conductor_replay": {
            "minimal_model_matches_discovery_artifact": True,
            "conductor_matches_discovery_artifact": True,
            "root_number_matches_discovery_artifact": True,
            "log_conductor": conductor_replay["log_conductor"],
        },
        "input_subset": {
            "point_count": len(points),
            "sha256": point_digest(points),
            "selection_basis": "stable numerical height rank at two precisions",
            "source": record["certificate_search_source"],
        },
        "small_prime_saturation": saturation,
        "saturated_basis": [
            {"jacobian_x": str(point[0]), "jacobian_y": str(point[1])}
            for point in saturated
        ],
        "finite_reduction_certificate": {
            "two_torsion_certificate_prime": two_torsion_prime,
            "reduced_2_division_cubic_has_no_root": True,
            "certificate_primes": [signature.prime for signature in signatures],
            "signatures": [
                {
                    "prime": signature.prime,
                    "group_order": signature.group_order,
                    "doubled_subgroup_order": signature.doubled_subgroup_order,
                    "quotient_dimension": signature.quotient_dimension,
                    "rows": [list(row) for row in signature.rows],
                }
                for signature in signatures
            ],
            "combined_exact_rank_over_F2": binary_rank,
            "certified_algebraic_rank_lower_bound": 17,
        },
    }


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    generated = root / "artifacts/generated-results/elliptic-curves"
    archived = root / "archive/elliptic-curves/artifacts/generated-results"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=archived / "elliptic_nagao_rank13_rank_gain_search.json",
    )
    parser.add_argument(
        "--mutation-input",
        type=Path,
        default=archived / "elliptic_nagao_rank13_rank_gain_mutations.json",
    )
    parser.add_argument("--saturation-bound", type=int, default=20)
    parser.add_argument("--saturation-timeout", type=float, default=20.0)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument("--certificate-prime-bound", type=int, default=500)
    parser.add_argument("--conductor-timeout", type=float, default=10.0)
    parser.add_argument(
        "--output",
        type=Path,
        default=generated / "elliptic_nagao_rank17_frontier_certificate.json",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if (
        args.saturation_bound < 3
        or args.saturation_timeout <= 0
        or args.conductor_timeout <= 0
    ):
        raise SystemExit("invalid saturation bounds")
    if args.certificate_prime_bound < 3:
        raise SystemExit("certificate prime bound must be at least 3")
    input_paths = (args.input, args.mutation_input)
    records = load_candidates(input_paths)
    certificates = []
    for parameter in CANDIDATE_PARAMETERS:
        certificate = certify_candidate(
            parameter,
            records[parameter],
            saturation_bound=args.saturation_bound,
            saturation_timeout=args.saturation_timeout,
            stack_bytes=args.stack_bytes,
            certificate_prime_bound=args.certificate_prime_bound,
            conductor_timeout=args.conductor_timeout,
        )
        certificates.append(certificate)
        print(
            f"u={parameter}: exact rank >=17, "
            f"logN={certificate['log_conductor']}",
            flush=True,
        )

    script_path = Path(__file__).resolve()
    engine_path = script_path.with_name("mod2_reduction_independence.py")
    saturation_path = script_path.with_name("extend_nagao_u42_frontier.py")
    artifact = {
        "schema_version": 1,
        "status": "four exact unconditional rank-at-least-17 certificates",
        "strict_log_conductor_target": str(LOG_CONDUCTOR_TARGET),
        "certificates": certificates,
        "theorem": {
            "argument": (
                "full column rank in a product of E(F_p)/2E(F_p) forces every "
                "integral relation coefficient to be even; trivial E(Q)[2] permits "
                "iteration, hence all coefficients vanish"
            ),
            "uses_numerical_heights_for_certification": False,
            "uses_BSD": False,
            "uses_parity_conjecture": False,
            "uses_full_2_descent": False,
            "depends_on_ellsaturation_finite_index_hypothesis": False,
            "uses_floating_logarithms_for_strict_target": False,
        },
        "frontier": {
            "certified_curve_count": len(certificates),
            "maximum_certified_rank_lower_bound": 17,
            "smallest_certified_log_conductor": min(
                (item["log_conductor"] for item in certificates), key=Decimal
            ),
            "rank21_log_conductor_target_certified": False,
            "rank30_target_certified": False,
            "remaining_independent_points_needed": 4,
        },
        "inputs": [
            {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
            for path in input_paths
        ],
        "parameters": {
            "saturation_prime_bound": args.saturation_bound,
            "saturation_timeout_seconds": args.saturation_timeout,
            "pari_stack_bytes": args.stack_bytes,
            "certificate_prime_bound": args.certificate_prime_bound,
            "conductor_timeout_seconds": args.conductor_timeout,
        },
        "software": {"python": platform.python_version(), "pari_gp": pari_version()},
        "reproducing_command": " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv]),
        "script_sha256": hashlib.sha256(script_path.read_bytes()).hexdigest(),
        "certificate_engine_sha256": hashlib.sha256(engine_path.read_bytes()).hexdigest(),
        "saturation_engine_sha256": hashlib.sha256(saturation_path.read_bytes()).hexdigest(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
