#!/usr/bin/env python3
"""Exact rank-17 certificate for the rank-21-family lead ``T=5777/32``.

The script replays one uniform height-1,000,000 quartic search, checks all
returned points exactly, and reconstructs the pinned 52-point Jacobian pool.
A two-precision height matrix only selects seventeen points.  PARI then
proposes their small-prime saturation, after which exact finite reductions
prove independence without relying on heights or on a full descent.

The conductor, root number, and the engineered local profile at 5, 7, 11,
and 13 are recomputed from the rational short model.  In particular, this
also records the corrected positive-parameter residue ``T mod 13 = 3``;
``10`` belongs to the equivalent negative-parameter convention.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import shlex
import sys
from typing import Any

from certify_nagao_rank17_frontier import exact_log_conductor_certificate
from ek_k3 import rational_to_string, valuation
from extend_nagao_u42_frontier import saturate_exact_basis
from mod2_reduction_independence import (
    combined_mod2_rank,
    find_mod2_reduction_certificate,
    find_two_torsion_certificate_prime,
)
from nagao_1994 import (
    PRIMARY_SOURCE,
    RANK21_CONSTRUCTION,
    primitive_quartic_coefficients,
    primitive_visible_points,
    quartic_point_to_short_jacobian,
    quartic_value,
    short_jacobian_coefficients,
)
from pari_bridge import minimal_curve_data, pari_version
from search_extra_points import signless_quartic_points
from search_nagao_u42_skew_height import run_mobius_charts
from triage_nagao_rank13_finalists import (
    height_matrix_replay,
    point_digest,
    point_on_short_curve,
    stable_height_rank,
)


Q = Fraction
PARAMETER_T = Q(5777, 32)
UNIFORM_HEIGHT = 1_000_000
EXPECTED_SIGNED_POINTS = 104
EXPECTED_NEW_IMAGES = 40
EXPECTED_POOL_SIZE = 52
EXPECTED_HEIGHT_SUBSET = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 16, 17, 19, 21)
LOCAL_PRIMES = (5, 7, 11, 13)
EXPECTED_RESIDUES = {5: 1, 7: 4, 11: 9, 13: 3}
EXPECTED_PRIMITIVE_DISCRIMINANT_VALUATIONS = {5: 7, 7: 6, 11: 5, 13: 4}
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/verify_nagao_rank21_t5777.py"
)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--search-timeout", type=float, default=90.0)
    parser.add_argument("--height-timeout", type=float, default=30.0)
    parser.add_argument("--saturation-timeout", type=float, default=20.0)
    parser.add_argument("--conductor-timeout", type=float, default=30.0)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            root
            / "artifacts/generated-results/elliptic_nagao_rank21_t5777_rank17_certificate.json"
        ),
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    for name in (
        "search_timeout",
        "height_timeout",
        "saturation_timeout",
        "conductor_timeout",
    ):
        if not 0 < getattr(args, name) <= 120:
            raise SystemExit(f"--{name.replace('_', '-')} must be in (0,120]")
    if args.stack_bytes < 64_000_000:
        raise SystemExit("--stack-bytes is too small")

    quartic_coefficients = primitive_quartic_coefficients(
        RANK21_CONSTRUCTION, PARAMETER_T
    )
    short_coefficients = short_jacobian_coefficients(
        RANK21_CONSTRUCTION, PARAMETER_T
    )
    raw_by_chart, pari_milliseconds, search_wall = run_mobius_charts(
        quartic_coefficients,
        (("identity", (1, 0, 0, 1)),),
        height_bound=UNIFORM_HEIGHT,
        timeout=args.search_timeout,
        stack_bytes=args.stack_bytes,
    )
    raw_points = raw_by_chart["identity"]
    if len(raw_points) != EXPECTED_SIGNED_POINTS:
        raise AssertionError("the pinned uniform signed-point count changed")
    if any(
        ordinate**2 != quartic_value(quartic_coefficients, parameter)
        for parameter, ordinate in raw_points
    ):
        raise AssertionError("PARI returned a point off the exact quartic")
    signless = signless_quartic_points(raw_points)
    visible_quartic = primitive_visible_points(RANK21_CONSTRUCTION, PARAMETER_T)
    visible_x = {point[0] for point in visible_quartic}
    visible_images = tuple(
        quartic_point_to_short_jacobian(
            RANK21_CONSTRUCTION, PARAMETER_T, quartic_point
        )
        for quartic_point in visible_quartic
    )
    seen_image_x = {point[0] for point in visible_images}
    new_records = []
    new_images = []
    for quartic_point in signless:
        if quartic_point[0] in visible_x or quartic_point[1] == 0:
            continue
        image = quartic_point_to_short_jacobian(
            RANK21_CONSTRUCTION, PARAMETER_T, quartic_point
        )
        if not point_on_short_curve(short_coefficients, image):
            raise AssertionError("a uniform-search image missed the Jacobian")
        if image[0] in seen_image_x:
            continue
        seen_image_x.add(image[0])
        new_images.append(image)
        new_records.append(
            {
                "quartic_x": rational_to_string(quartic_point[0]),
                "quartic_z": rational_to_string(quartic_point[1]),
                "jacobian_x": rational_to_string(image[0]),
                "jacobian_y": rational_to_string(image[1]),
                "exact_quartic_and_jacobian_membership_checked": True,
            }
        )
    if len(new_images) != EXPECTED_NEW_IMAGES:
        raise AssertionError("the pinned uniform new-image count changed")
    pool = visible_images + tuple(new_images)
    if len(pool) != EXPECTED_POOL_SIZE:
        raise AssertionError("the exact Jacobian pool changed size")

    height_runs = height_matrix_replay(
        short_coefficients,
        pool,
        precisions=(72, 120),
        timeout=args.height_timeout,
        stack_bytes=args.stack_bytes,
    )
    if stable_height_rank(height_runs) != 17:
        raise AssertionError("the pool's stable numerical rank changed")
    selected_indices = tuple(height_runs[-1]["subset_indices_one_based"])
    if selected_indices != EXPECTED_HEIGHT_SUBSET:
        raise AssertionError("the stable numerical subset changed")
    selected = tuple(pool[index - 1] for index in selected_indices)
    saturated_basis, saturation = saturate_exact_basis(
        short_coefficients,
        selected,
        prime_bound=20,
        timeout=args.saturation_timeout,
        stack_bytes=args.stack_bytes,
    )
    if len(saturated_basis) != 17:
        raise AssertionError("small-prime saturation changed the basis length")
    signatures = find_mod2_reduction_certificate(
        short_coefficients, saturated_basis, prime_bound=500
    )
    exact_rank = combined_mod2_rank(signatures, len(saturated_basis))
    if exact_rank != 17:
        raise AssertionError("finite reductions did not certify rank 17")
    two_torsion_prime = find_two_torsion_certificate_prime(short_coefficients)

    primitive_discriminant = RANK21_CONSTRUCTION.primitive_discriminant_value(
        PARAMETER_T
    )
    residues = {
        prime: (
            PARAMETER_T.numerator
            * pow(PARAMETER_T.denominator, -1, prime)
            % prime
        )
        for prime in LOCAL_PRIMES
    }
    discriminant_valuations = {
        prime: valuation(primitive_discriminant, prime)
        for prime in LOCAL_PRIMES
    }
    if residues != EXPECTED_RESIDUES:
        raise AssertionError("the positive-T local residue profile changed")
    if discriminant_valuations != EXPECTED_PRIMITIVE_DISCRIMINANT_VALUATIONS:
        raise AssertionError("the primitive discriminant valuations changed")
    conductor = minimal_curve_data(
        short_coefficients,
        timeout=args.conductor_timeout,
        local_primes=LOCAL_PRIMES,
        stack_bytes=args.stack_bytes,
    )
    for prime in LOCAL_PRIMES:
        local = conductor["local_reduction"][str(prime)]
        if (
            local["conductor_exponent"] != 1
            or local["minimal_c4_valuation"] != 0
            or local["minimal_discriminant_valuation"]
            != EXPECTED_PRIMITIVE_DISCRIMINANT_VALUATIONS[prime]
            or local["ellap"] != 1
        ):
            raise AssertionError(f"the exact local profile changed at {prime}")
    exact_log_bound = exact_log_conductor_certificate(conductor["conductor"])

    script_path = Path(__file__).resolve()
    artifact: dict[str, Any] = {
        "schema_version": 1,
        "status": "exact_rank17_certificate_complete",
        "theorem": (
            "the Nagao rank-21-family specialization T=5777/32 has "
            "Mordell-Weil rank at least 17 and log conductor below 182.72"
        ),
        "candidate": {
            "parameter_t": rational_to_string(PARAMETER_T),
            "short_weierstrass_coefficients": [
                rational_to_string(value) for value in short_coefficients
            ],
            "minimal_model": list(conductor["minimal_model"]),
            "conductor": str(conductor["conductor"]),
            "log_conductor": conductor["log_conductor"],
            "root_number": conductor["root_number"],
            "below_strict_log_conductor_target": True,
            "exact_log_conductor_bound": exact_log_bound,
        },
        "primary_source": PRIMARY_SOURCE,
        "uniform_search": {
            "height_bound": UNIFORM_HEIGHT,
            "signed_point_count": len(raw_points),
            "signless_point_count": len(signless),
            "visible_point_count": len(visible_images),
            "new_distinct_jacobian_images": len(new_images),
            "raw_point_sha256": point_digest(raw_points),
            "new_image_sha256": point_digest(new_images),
            "exact_pool_point_count": len(pool),
            "exact_pool_sha256": point_digest(pool),
            "all_memberships_checked_exactly": True,
            "pari_reported_milliseconds": pari_milliseconds["identity"],
            "wall_seconds": search_wall,
            "new_points": new_records,
        },
        "height_selection": {
            "runs": list(height_runs),
            "stable_numerical_rank": 17,
            "selected_pool_indices_one_based": list(selected_indices),
            "selection_is_not_certification": True,
        },
        "exact_rank_certificate": {
            "small_prime_saturation": saturation,
            "saturated_basis_sha256": point_digest(saturated_basis),
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
            "combined_exact_rank_over_F2": exact_rank,
            "certified_algebraic_rank_lower_bound": exact_rank,
        },
        "engineered_local_profile": {
            "positive_parameter_residues": {
                str(prime): residues[prime] for prime in LOCAL_PRIMES
            },
            "sign_convention_correction": (
                "for positive T, T mod 13 is 3; residue 10 is -T mod 13"
            ),
            "primitive_discriminant_valuations": {
                str(prime): discriminant_valuations[prime]
                for prime in LOCAL_PRIMES
            },
            "minimal_local_reduction": conductor["local_reduction"],
            "all_four_split_multiplicative_with_conductor_exponent_one": True,
        },
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
        f"wrote {args.output}: exact_rank>={exact_rank} "
        f"logN={conductor['log_conductor']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
