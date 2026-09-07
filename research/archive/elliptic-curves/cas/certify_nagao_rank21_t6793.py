#!/usr/bin/env python3
"""Replay the exact rank-19 lower-bound certificate at ``T=6793/64``.

The 19 input points were selected by the unbiased height-1,000,000 search.
This verifier pins that source artifact by SHA-256, checks every point on the
exact Jacobian, repeats small-prime saturation, and proves independence by
finite reductions modulo good primes.  Numerical heights select the input;
they are not used in the rank proof.  The conductor, minimal model, root
number, and a rational proof of the strict logarithmic bound are replayed.
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

from certify_nagao_rank17_frontier import exact_log_conductor_certificate
from ek_k3 import rational_to_string
from extend_nagao_u42_frontier import saturate_exact_basis
from mod2_reduction_independence import (
    combined_mod2_rank,
    find_mod2_reduction_certificate,
    find_two_torsion_certificate_prime,
)
from nagao_1994 import PRIMARY_SOURCE, RANK21_CONSTRUCTION, short_jacobian_coefficients
from pari_bridge import minimal_curve_data, pari_version
from triage_nagao_rank13_finalists import point_digest, point_on_short_curve


Q = Fraction
PARAMETER_T = Q(6793, 64)
SOURCE_CANDIDATE_ID = "unbiased-6793-64"
SOURCE_STAGE_HEIGHT = 1_000_000
EXPECTED_SOURCE_SHA256 = (
    "5bf7406855af5ec39b269fa4105c9225adb4a10d13fab5480b15264cc3e8fe1d"
)
EXPECTED_SOURCE_SCRIPT_SHA256 = (
    "f58e093a059f386a0ee80f7e73d44e647f3b90662109cbd9b864bad508ac1f11"
)
EXPECTED_INPUT_POINT_SHA256 = (
    "22ddb4c3f21237f0045dec6d6ad325ac74e8e319ec7d5b0b90f40b4cc89223aa"
)
EXPECTED_SATURATED_BASIS_SHA256 = (
    "cdb7328683b523f49aac5efe3588631e82e9f18dd603b1329bc4e9f7c89e44dd"
)
EXPECTED_CERTIFICATE_PRIMES = (
    17,
    19,
    29,
    43,
    47,
    53,
    59,
    61,
    73,
    89,
    97,
    101,
    103,
    113,
    127,
    131,
)
EXPECTED_CONDUCTOR = (
    736590680363874295586257187445083659346347463739287340897386953172230
)
EXPECTED_MINIMAL_MODEL = (
    1,
    0,
    0,
    -1772157768312348985364211445680,
    954109168232910161122740568758658458094442752,
)
EXPECTED_MINIMAL_DISCRIMINANT = (
    -37065670903236142858397191464229556726205264807933586145063935448714902343475493731303424000
)
TARGET_LOG_CONDUCTOR = Decimal("182.72")
CERTIFIED_RANK_LOWER_BOUND = 19
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/certify_nagao_rank21_t6793.py"
)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_pinned_input(path: Path) -> tuple[tuple[tuple[Fraction, Fraction], ...], dict[str, Any]]:
    """Load only the pinned H=1m subset and return its search provenance."""

    actual_sha256 = sha256_file(path)
    if actual_sha256 != EXPECTED_SOURCE_SHA256:
        raise AssertionError(
            "the unbiased source artifact changed: "
            f"expected {EXPECTED_SOURCE_SHA256}, got {actual_sha256}"
        )
    data = json.loads(path.read_text())
    if data.get("script_sha256") != EXPECTED_SOURCE_SCRIPT_SHA256:
        raise AssertionError("the source artifact's generating-script hash changed")
    stages = tuple(data["point_stages"])
    stage = next(
        item
        for item in stages
        if int(item["quartic_naive_height_bound"]) == SOURCE_STAGE_HEIGHT
    )
    record = next(
        item
        for item in stage["ranked_population"]
        if item["candidate_id"] == SOURCE_CANDIDATE_ID
    )
    if record["constructor_parameter"] != rational_to_string(PARAMETER_T):
        raise AssertionError("the source parameter changed")
    rank_record = record["height_rank"]
    if (
        rank_record.get("status") != "completed"
        or int(rank_record["stable_numerical_rank"]) != CERTIFIED_RANK_LOWER_BOUND
    ):
        raise AssertionError("the source no longer contains its numerical rank-19 subset")
    points = tuple(
        (Q(item["jacobian_x"]), Q(item["jacobian_y"]))
        for item in record["exact_selected_points"]
    )
    if len(points) != CERTIFIED_RANK_LOWER_BOUND:
        raise AssertionError("the source subset does not contain 19 points")
    if point_digest(points) != EXPECTED_INPUT_POINT_SHA256:
        raise AssertionError("the pinned source point digest changed")

    stage_ranks = {}
    for source_stage in stages:
        source_record = next(
            (
                item
                for item in source_stage["ranked_population"]
                if item["candidate_id"] == SOURCE_CANDIDATE_ID
            ),
            None,
        )
        if source_record is not None:
            stage_ranks[str(source_stage["quartic_naive_height_bound"])] = {
                "stable_numerical_rank": source_record["height_rank"].get(
                    "stable_numerical_rank"
                ),
                "signed_points": source_record["point_search"]["signed_points"],
                "new_distinct_jacobian_sign_pairs": source_record["point_search"][
                    "new_distinct_jacobian_sign_pairs"
                ],
            }
    provenance = {
        "source_artifact": str(path),
        "source_artifact_sha256": actual_sha256,
        "source_script_sha256": data["script_sha256"],
        "source_candidate_id": SOURCE_CANDIDATE_ID,
        "uniform_height_bound": SOURCE_STAGE_HEIGHT,
        "stage_point_and_rank_history": stage_ranks,
        "selected_subset_indices_one_based": rank_record[
            "selected_subset_indices_one_based"
        ],
        "input_point_count": len(points),
        "input_point_sha256": point_digest(points),
        "stable_numerical_rank": rank_record["stable_numerical_rank"],
        "numerical_rank_is_selection_only": True,
    }
    return points, provenance


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=(
            root
            / "artifacts/generated-results/elliptic_nagao_rank21_unbiased.json"
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            root
            / "artifacts/generated-results/elliptic_nagao_rank21_t6793_rank19_certificate.json"
        ),
    )
    parser.add_argument("--saturation-timeout", type=float, default=30.0)
    parser.add_argument("--conductor-timeout", type=float, default=30.0)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument("--certificate-prime-bound", type=int, default=1000)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if not 0 < args.saturation_timeout <= 120:
        raise SystemExit("--saturation-timeout must be in (0,120]")
    if not 0 < args.conductor_timeout <= 120:
        raise SystemExit("--conductor-timeout must be in (0,120]")
    if args.stack_bytes < 64_000_000:
        raise SystemExit("--stack-bytes is too small")
    if not 131 <= args.certificate_prime_bound <= 10_000:
        raise SystemExit("--certificate-prime-bound must be in [131,10000]")

    input_points, provenance = load_pinned_input(args.input.resolve())
    coefficients = short_jacobian_coefficients(RANK21_CONSTRUCTION, PARAMETER_T)
    if any(not point_on_short_curve(coefficients, point) for point in input_points):
        raise AssertionError("a pinned input point is not on the exact Jacobian")

    saturated_basis, saturation = saturate_exact_basis(
        coefficients,
        input_points,
        prime_bound=20,
        timeout=args.saturation_timeout,
        stack_bytes=args.stack_bytes,
    )
    if len(saturated_basis) != CERTIFIED_RANK_LOWER_BOUND:
        raise AssertionError("small-prime saturation changed the basis length")
    if any(not point_on_short_curve(coefficients, point) for point in saturated_basis):
        raise AssertionError("a saturated point is not on the exact Jacobian")
    saturated_digest = point_digest(saturated_basis)
    if saturated_digest != EXPECTED_SATURATED_BASIS_SHA256:
        raise AssertionError("the pinned saturated basis changed")

    signatures = find_mod2_reduction_certificate(
        coefficients,
        saturated_basis,
        prime_bound=args.certificate_prime_bound,
    )
    exact_binary_rank = combined_mod2_rank(signatures, len(saturated_basis))
    if exact_binary_rank != CERTIFIED_RANK_LOWER_BOUND:
        raise AssertionError("finite reductions did not certify all 19 points")
    certificate_primes = tuple(signature.prime for signature in signatures)
    if certificate_primes != EXPECTED_CERTIFICATE_PRIMES:
        raise AssertionError("the deterministic finite-reduction certificate changed")
    two_torsion_prime = find_two_torsion_certificate_prime(coefficients)

    conductor = minimal_curve_data(
        coefficients,
        timeout=args.conductor_timeout,
        stack_bytes=args.stack_bytes,
    )
    if int(conductor["conductor"]) != EXPECTED_CONDUCTOR:
        raise AssertionError("the exact conductor changed")
    if tuple(conductor["minimal_model"]) != EXPECTED_MINIMAL_MODEL:
        raise AssertionError("the exact minimal model changed")
    if int(conductor["minimal_discriminant"]) != EXPECTED_MINIMAL_DISCRIMINANT:
        raise AssertionError("the exact minimal discriminant changed")
    if int(conductor["root_number"]) != -1:
        raise AssertionError("the exact root number changed")
    if Decimal(conductor["log_conductor"]) >= TARGET_LOG_CONDUCTOR:
        raise AssertionError("the candidate crossed the strict conductor target")
    exact_log_bound = exact_log_conductor_certificate(EXPECTED_CONDUCTOR)

    script_path = Path(__file__).resolve()
    artifact: dict[str, Any] = {
        "schema_version": 1,
        "status": "exact_rank_at_least_19_certificate_complete",
        "theorem": (
            "The Nagao rank-21-family specialization T=6793/64 has "
            "Mordell-Weil rank at least 19 over Q and log conductor below 182.72."
        ),
        "candidate": {
            "parameter_t": rational_to_string(PARAMETER_T),
            "short_weierstrass_coefficients": [
                rational_to_string(value) for value in coefficients
            ],
            "minimal_model": list(conductor["minimal_model"]),
            "minimal_discriminant": str(conductor["minimal_discriminant"]),
            "conductor": str(conductor["conductor"]),
            "log_conductor": conductor["log_conductor"],
            "root_number": conductor["root_number"],
            "strict_log_conductor_target": str(TARGET_LOG_CONDUCTOR),
            "below_strict_log_conductor_target": True,
            "exact_log_conductor_bound": exact_log_bound,
        },
        "search_provenance": provenance,
        "exact_rank_certificate": {
            "argument": (
                "Exact reduction maps at the listed good primes give a combined "
                "19-dimensional image in E(F_p)/2E(F_p).  The separate "
                "2-division-cubic certificate proves E(Q)[2]=0, so the 19 "
                "displayed rational points are independent."
            ),
            "small_prime_saturation": saturation,
            "saturated_basis_sha256": saturated_digest,
            "saturated_basis": [
                {
                    "jacobian_x": rational_to_string(point[0]),
                    "jacobian_y": rational_to_string(point[1]),
                    "exact_jacobian_membership_checked": True,
                }
                for point in saturated_basis
            ],
            "two_torsion_certificate_prime": two_torsion_prime,
            "finite_reduction_signatures": [
                {
                    "prime": signature.prime,
                    "group_order": signature.group_order,
                    "doubled_subgroup_order": signature.doubled_subgroup_order,
                    "quotient_dimension": signature.quotient_dimension,
                    "rows": [list(row) for row in signature.rows],
                }
                for signature in signatures
            ],
            "combined_exact_rank_over_F2": exact_binary_rank,
            "certified_algebraic_rank_lower_bound": exact_binary_rank,
            "height_matrices_not_used_in_certificate": True,
        },
        "interpretation": {
            "target_rank": 21,
            "target_reached": False,
            "rank_upper_bound_not_claimed": True,
            "root_number_parity_is_not_used_in_the_rank_certificate": True,
        },
        "primary_source": PRIMARY_SOURCE,
        "software": {
            "python": platform.python_version(),
            "pari_gp": pari_version(),
            "platform": platform.platform(),
        },
        "reproducing_command": REPRODUCING_COMMAND,
        "actual_command": " ".join(
            shlex.quote(part) for part in [sys.executable, *sys.argv]
        ),
        "script_sha256": sha256_file(script_path),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(
        f"wrote {args.output}: exact_rank>={exact_binary_rank} "
        f"logN={conductor['log_conductor']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
