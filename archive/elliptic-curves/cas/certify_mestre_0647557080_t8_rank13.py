#!/usr/bin/env python3
"""Certify rank at least 13 for the max-root-100 Mestre lead.

The pinned max-root-100 artifact found stable numerical rank 13 at roots
``(0,6,47,55,70,80)`` and ``T=8``.  This script exactly reconstructs its
H=5000 quartic pool, separates the twelve displayed points from the first
three accidental directions, uses PARI saturation only to propose a better
basis, and then proves independence of that basis by exact finite reductions.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from fractions import Fraction
import json
import os
from pathlib import Path
import platform
import shlex
import sys
from typing import Any

from ek_k3 import rational_to_string
from extend_nagao_u42_frontier import saturate_exact_basis
from mestre_root_tuples import SixRootMestreConstruction
from search_mestre_root_tuple_scale import (
    bounded_quartic_points,
    canonical_signless_points,
    finite_reduction_attempt,
    height_matrix_replay,
    numerical_subset,
    point_digest,
    point_record,
    point_on_short_curve,
    primitive_visible_points,
    quartic_point_to_jacobian,
    quartic_value,
    sha256_file,
)


Q = Fraction
ROOTS = (0, 6, 47, 55, 70, 80)
PARAMETER = Q(8)
HEIGHT_BOUND = 5_000
STACK_BYTES = 512_000_000
INPUT_ARTIFACT = "elliptic_mestre_root_tuple_scale_max100.json"
EXPECTED_INPUT_ARTIFACT_SHA256 = (
    "63dcd39555ad8b39c7b584a16663164bf73e6c6c59906b6a230bfa9b9f65a3bb"
)
EXPECTED_POOL_SHA256 = (
    "5cbb40eac95c97f6da9fcd6e6c7be57da78fd9b5c2f9b7daed29924f05596681"
)
EXPECTED_SATURATED_BASIS_SHA256 = (
    "1a038326d2caff9bac0310cc484f21921270666c20f9bfdbb2d48cf8abd7975f"
)
EXPECTED_CERTIFICATE_PRIMES = (5, 13, 19, 23, 29, 31, 41, 47, 67, 71, 73)


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--point-timeout", type=float, default=15.0)
    parser.add_argument("--height-timeout", type=float, default=15.0)
    parser.add_argument("--saturation-timeout", type=float, default=30.0)
    parser.add_argument(
        "--output",
        type=Path,
        default=root
        / "artifacts"
        / "generated-results"
        / "elliptic_mestre_0647557080_t8_rank13_certificate.json",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if min(args.point_timeout, args.height_timeout, args.saturation_timeout) <= 0:
        raise SystemExit("subprocess caps must be positive")
    if max(args.point_timeout, args.height_timeout, args.saturation_timeout) > 30:
        raise SystemExit("subprocess caps must not exceed 30 seconds")
    if args.output.exists():
        raise SystemExit("refusing to overwrite the rank-13 certificate")
    root = Path(__file__).resolve().parents[2]
    input_path = root / "artifacts" / "generated-results" / INPUT_ARTIFACT
    if sha256_file(input_path) != EXPECTED_INPUT_ARTIFACT_SHA256:
        raise AssertionError("the pinned max-root-100 artifact changed")
    input_artifact = json.loads(input_path.read_text())
    lead = next(
        record
        for record in input_artifact["specialization_screen"]["h5000_records"]
        if record["identifier"] == "r0_6_47_55_70_80_t8"
    )

    construction = SixRootMestreConstruction(tuple(Q(root) for root in ROOTS))
    coefficients = construction.primitive_jacobian_coefficients(PARAMETER)
    quartic_coefficients = construction.primitive_quartic_coefficients(PARAMETER)
    visible_quartic = primitive_visible_points(construction, PARAMETER)
    visible_jacobian = tuple(
        quartic_point_to_jacobian(construction, PARAMETER, point)
        for point in visible_quartic
    )
    if len(visible_quartic) != 12 or len({point[0] for point in visible_quartic}) != 12:
        raise AssertionError("the twelve displayed points changed")
    if any(not point_on_short_curve(coefficients, point) for point in visible_jacobian):
        raise AssertionError("a displayed point missed the exact Jacobian")

    raw = bounded_quartic_points(
        quartic_coefficients,
        height_bound=HEIGHT_BOUND,
        timeout=args.point_timeout,
        stack_bytes=STACK_BYTES,
    )
    signless = canonical_signless_points(raw)
    if any(
        point[1] ** 2 != quartic_value(quartic_coefficients, point[0])
        for point in signless
    ):
        raise AssertionError("an H=5000 point missed the exact quartic")
    searched_jacobian = tuple(
        quartic_point_to_jacobian(construction, PARAMETER, point)
        for point in signless
    )
    pool = list(visible_jacobian)
    pool_sources: list[dict[str, Any]] = [
        {
            "source": "displayed generic point",
            "quartic_point": point_record(point),
        }
        for point in visible_quartic
    ]
    seen_x = {point[0] for point in pool}
    for quartic_point, jacobian_point in zip(signless, searched_jacobian, strict=True):
        if jacobian_point[0] in seen_x:
            continue
        seen_x.add(jacobian_point[0])
        pool.append(jacobian_point)
        pool_sources.append(
            {
                "source": "H=5000 accidental point",
                "quartic_point": point_record(quartic_point),
            }
        )
    pool_tuple = tuple(pool)
    if len(signless) != 61 or len(pool_tuple) != 61:
        raise AssertionError("the exact H=5000 pool cardinality changed")
    if point_digest(pool_tuple) != EXPECTED_POOL_SHA256:
        raise AssertionError("the exact H=5000 pool digest changed")

    prefix_records = []
    for prefix_count in (12, 13, 14, 15, 61):
        height = height_matrix_replay(
            coefficients,
            pool_tuple[:prefix_count],
            precisions=(72, 120),
            timeout=args.height_timeout,
            stack_bytes=STACK_BYTES,
        )
        prefix_records.append(
            {
                "pool_prefix_count": prefix_count,
                "height_matrix_runs": list(height),
                "stable_numerical_rank": int(height[-1]["numerical_rank"]),
                "subset_indices_one_based": height[-1]["subset_indices_one_based"],
            }
        )
    if [record["stable_numerical_rank"] for record in prefix_records] != [
        10,
        11,
        12,
        13,
        13,
    ]:
        raise AssertionError("the displayed/accidental rank split changed")
    full_height = prefix_records[-1]["height_matrix_runs"]
    proposed_basis = numerical_subset(pool_tuple, full_height)
    if len(proposed_basis) != 13:
        raise AssertionError("the numerical basis proposal changed rank")

    saturated_basis, saturation = saturate_exact_basis(
        coefficients,
        proposed_basis,
        prime_bound=50,
        timeout=args.saturation_timeout,
        stack_bytes=STACK_BYTES,
    )
    if len(saturated_basis) != 13:
        raise AssertionError("small-prime saturation changed basis length")
    if point_digest(saturated_basis) != EXPECTED_SATURATED_BASIS_SHA256:
        raise AssertionError("the saturated basis digest changed")
    certificate = finite_reduction_attempt(
        coefficients, saturated_basis, prime_bound=100
    )
    if (
        certificate["certified_algebraic_rank_lower_bound"] != 13
        or certificate["combined_exact_rank_over_F2"] != 13
        or tuple(certificate["certificate_primes"]) != EXPECTED_CERTIFICATE_PRIMES
        or certificate["two_torsion_certificate_prime"] != 11
    ):
        raise AssertionError("the exact finite-reduction certificate changed")

    symmetry_panel = tuple(
        Q(value) for value in (1, 2, 3, 4, 5, 6, 7, 8)
    ) + (Q(1, 2), Q(3, 2), Q(7, 3))
    if any(
        construction.primitive_quartic_coefficients(parameter)
        != construction.primitive_quartic_coefficients(-parameter)
        or construction.primitive_jacobian_coefficients(parameter)
        != construction.primitive_jacobian_coefficients(-parameter)
        for parameter in symmetry_panel
    ):
        raise AssertionError("the exact T -> -T symmetry failed its replay panel")

    artifact = {
        "schema_version": 1,
        "status": "exact algebraic rank lower bound 13 certified",
        "input": {
            "max100_artifact": INPUT_ARTIFACT,
            "max100_artifact_sha256": EXPECTED_INPUT_ARTIFACT_SHA256,
            "roots": list(ROOTS),
            "parameter": 8,
            "conductor_phase": lead["conductor_phase"],
        },
        "symmetry": {
            "exact_quotient": "T is identified with -T",
            "structural_reason": (
                "the root multiset of q(X-T)q(X+T) is unchanged by T -> -T; "
                "the unique square approximant and remainder are therefore even, "
                "and division by T^2 preserves equality"
            ),
            "quartic_and_jacobian_coefficients_replayed_equal_on_panel": [
                rational_to_string(parameter) for parameter in symmetry_panel
            ],
            "visible_point_sets_are_exchanged_up_to_ordinate_sign": True,
        },
        "H5000_reconstruction": {
            "height_bound": HEIGHT_BOUND,
            "signed_points_returned": len(raw),
            "signless_quartic_point_count": len(signless),
            "displayed_point_count": len(visible_quartic),
            "pool_point_count_modulo_inverse": len(pool_tuple),
            "pool_point_sha256": point_digest(pool_tuple),
            "prefix_rank_replay": prefix_records,
            "displayed_points_span_stable_numerical_rank": 10,
            "first_three_accidental_points_raise_stable_rank_successively": [
                11,
                12,
                13,
            ],
            "first_three_accidental_directions": [
                {
                    **pool_sources[index],
                    "pool_index_one_based": index + 1,
                    "jacobian_point": point_record(pool_tuple[index]),
                }
                for index in (12, 13, 14)
            ],
            "all_points_checked_exactly": True,
        },
        "basis_proposal": {
            "source": "precision-stable numerical height subset; not itself a proof",
            "point_count": len(proposed_basis),
            "point_sha256": point_digest(proposed_basis),
            "points": [point_record(point) for point in proposed_basis],
        },
        "small_prime_saturation": saturation,
        "exact_finite_reduction_certificate": certificate,
        "conclusion": {
            "certified_algebraic_rank_lower_bound": 13,
            "depends_on_PARIs_finite_index_saturation_hypothesis": False,
            "reason": (
                "PARI saturation only proposed exact rational points; exact curve "
                "membership and the full finite-reduction rank prove independence"
            ),
        },
        "parameters": {
            "point_timeout_seconds": args.point_timeout,
            "height_timeout_seconds": args.height_timeout,
            "saturation_timeout_seconds": args.saturation_timeout,
            "saturation_prime_bound_strict_upper_limit": 50,
            "finite_reduction_prime_bound": 100,
            "stack_bytes": STACK_BYTES,
            "no_retries": True,
        },
        "provenance": {
            "script": str(Path(__file__).resolve().relative_to(root)),
            "script_sha256": sha256_file(Path(__file__).resolve()),
            "reproducing_command": " ".join(
                shlex.quote(part) for part in [sys.executable, *sys.argv]
            ),
            "owned_processes_remaining": 0,
        },
        "software": {
            "python": platform.python_version(),
            "platform": platform.platform(),
        },
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(args.output, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    with os.fdopen(descriptor, "w") as stream:
        json.dump(artifact, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(
        json.dumps(
            {
                "certified_algebraic_rank_lower_bound": 13,
                "pool_sha256": point_digest(pool_tuple),
                "saturated_basis_sha256": point_digest(saturated_basis),
                "certificate_primes": certificate["certificate_primes"],
            },
            sort_keys=True,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
