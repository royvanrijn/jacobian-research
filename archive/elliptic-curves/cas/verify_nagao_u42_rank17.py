#!/usr/bin/env python3
"""Exact rank-at-least-17 certificate for Nagao's ``u=42`` specialization.

This verifier does not use numerical heights, BSD, root-number parity, or a
full 2-descent.  It proves that 17 exact rational points have independent
images in a product of groups ``E(F_p)/2E(F_p)`` and separately proves that
the rational curve has no rational 2-torsion.  Infinite descent on the
coefficients of an integral point relation then proves rank at least 17.
"""

from __future__ import annotations

import argparse
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import shlex
import sys
from typing import Any

from mod2_reduction_independence import (
    combined_mod2_rank,
    mod2_reduction_signature,
    short_curve_has_no_rational_2_torsion_modular_certificate,
)
from nagao_1994 import rank13_base_changed_short_jacobian_coefficients
from triage_nagao_rank13_finalists import point_digest, point_on_short_curve


Q = Fraction
PARAMETER_U = 42
CERTIFICATE_PRIMES = (11, 23, 29, 41, 43, 53, 59, 61, 67, 79, 83, 97, 107, 131, 137)
TWO_TORSION_CERTIFICATE_PRIME = 31
LOG_CONDUCTOR_TARGET = "182.72"


def load_basis(path: Path) -> tuple[tuple[Fraction, Fraction], ...]:
    data = json.loads(path.read_text())
    if int(data["candidate"]["parameter_u"]) != PARAMETER_U:
        raise AssertionError("the point artifact is not the u=42 specialization")
    points = tuple(
        (Q(record["jacobian_x"]), Q(record["jacobian_y"]))
        for record in data["small_prime_saturation"]["saturated_basis"]
    )
    if len(points) != 17:
        raise AssertionError("the pinned basis does not contain 17 points")
    if point_digest(points) != data["small_prime_saturation"]["saturated_basis_sha256"]:
        raise AssertionError("the exact point digest changed")
    return points


def load_conductor_record(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text())
    record = next(
        item
        for item in data["final_conductor_candidates"]
        if int(item["parameter_u"]) == PARAMETER_U
    )
    if Q(record["parameter_t"]) != Q(3631, 14):
        raise AssertionError("the u=42 base parameter changed")
    return record


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    generated = root / "artifacts/generated-results"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--point-input",
        type=Path,
        default=generated / "elliptic_nagao_u42_height_10000000.json",
    )
    parser.add_argument(
        "--conductor-input",
        type=Path,
        default=generated / "elliptic_nagao_rank13_integer_u.json",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=generated / "elliptic_nagao_u42_rank17_certificate.json",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    points = load_basis(args.point_input)
    conductor_record = load_conductor_record(args.conductor_input)
    coefficients = rank13_base_changed_short_jacobian_coefficients(Q(PARAMETER_U))
    if any(not point_on_short_curve(coefficients, point) for point in points):
        raise AssertionError("an exact certificate point is not on the rational curve")

    signatures = tuple(
        mod2_reduction_signature(coefficients, points, prime)
        for prime in CERTIFICATE_PRIMES
    )
    binary_rank = combined_mod2_rank(signatures, len(points))
    if binary_rank != len(points):
        raise AssertionError("the combined reduction signatures lost full rank")
    no_two_torsion = short_curve_has_no_rational_2_torsion_modular_certificate(
        coefficients, TWO_TORSION_CERTIFICATE_PRIME
    )
    if not no_two_torsion:
        raise AssertionError("the modular 2-torsion certificate failed")

    log_conductor = str(conductor_record["pari"]["log_conductor"])
    if not Decimal(log_conductor) < Decimal(LOG_CONDUCTOR_TARGET):
        raise AssertionError("the certified curve no longer meets the conductor bound")

    script_path = Path(__file__).resolve()
    engine_path = script_path.with_name("mod2_reduction_independence.py")
    artifact = {
        "schema_version": 1,
        "status": "exact unconditional rank lower-bound certificate",
        "candidate": {
            "parameter_u": PARAMETER_U,
            "parameter_t": str(Q(conductor_record["parameter_t"])),
            "short_weierstrass_coefficients": [str(value) for value in coefficients],
            "conductor": str(conductor_record["pari"]["conductor"]),
            "log_conductor": str(log_conductor),
            "strict_log_conductor_target": LOG_CONDUCTOR_TARGET,
            "below_strict_log_conductor_target": True,
        },
        "exact_points": {
            "count": len(points),
            "sha256": point_digest(points),
            "all_exact_curve_membership_checks": True,
        },
        "rational_two_torsion": {
            "certificate_prime": TWO_TORSION_CERTIFICATE_PRIME,
            "reduced_2_division_cubic_has_no_root": True,
            "conclusion": "E(Q)[2] is trivial by irreducibility modulo 31",
        },
        "mod2_reduction_signatures": [
            {
                "prime": signature.prime,
                "group_order": signature.group_order,
                "doubled_subgroup_order": signature.doubled_subgroup_order,
                "quotient_dimension": signature.quotient_dimension,
                "rows": [list(row) for row in signature.rows],
            }
            for signature in signatures
        ],
        "combined_binary_matrix": {
            "row_count": sum(signature.quotient_dimension for signature in signatures),
            "column_count": len(points),
            "exact_rank_over_F2": binary_rank,
        },
        "theorem": {
            "certified_algebraic_rank_lower_bound": 17,
            "argument": (
                "a rational integral relation reduces to every E(F_p)/2E(F_p); "
                "full column rank makes all coefficients even. Dividing by two "
                "uses E(Q)[2]=0 and iterating forces every integer coefficient to vanish"
            ),
            "uses_numerical_heights": False,
            "uses_BSD": False,
            "uses_parity_conjecture": False,
            "uses_full_2_descent": False,
            "depends_on_ellsaturation_finite_index_hypothesis": False,
        },
        "target_status": {
            "rank21_log_conductor_target_certified": False,
            "rank30_target_certified": False,
            "remaining_independent_points_needed_for_rank21": 4,
        },
        "inputs": [
            {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
            for path in (args.point_input, args.conductor_input)
        ],
        "software": {"python": platform.python_version()},
        "reproducing_command": " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv]),
        "script_sha256": hashlib.sha256(script_path.read_bytes()).hexdigest(),
        "certificate_engine_sha256": hashlib.sha256(engine_path.read_bytes()).hexdigest(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(f"exact combined F2 rank: {binary_rank}/{len(points)}")
    print("certified algebraic rank >= 17")
    print(f"log conductor: {log_conductor} < {LOG_CONDUCTOR_TARGET}")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
