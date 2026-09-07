#!/usr/bin/env python3
"""Certify rank at least 20 for Nagao's section-7 specialization.

Nagao prints the parameter ``t=5081/94`` for the six roots
``(346,260,255,146,55,0)``.  The repository's symmetric constructor uses
``T=2t=5081/47``.  This verifier replays one uniform quartic search, constructs
48 deterministic cross-ratio charts from its exact abscissas, maps all points
exactly to the Jacobian, and uses numerical heights only to select 20 points.
Small-prime saturation and finite-reduction signatures then prove their
independence.  The minimal model and conductor are recomputed by PARI.
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
from mestre_root_tuples import SixRootMestreConstruction
from mod2_reduction_independence import (
    combined_mod2_rank,
    find_mod2_reduction_certificate,
    find_two_torsion_certificate_prime,
)
from nagao_1994 import (
    PRIMARY_SOURCE,
    primitive_quartic_coefficients,
    primitive_visible_points,
    quartic_point_to_short_jacobian,
    quartic_value,
    short_jacobian_coefficients,
)
from pari_bridge import minimal_curve_data, pari_version
from search_extra_points import signless_quartic_points
from search_nagao_rank21_t6793_skew import (
    map_run_points,
    optimized_cross_ratio_charts,
)
from search_nagao_u42_skew_height import run_mobius_charts
from triage_nagao_rank13_finalists import (
    height_matrix_replay,
    point_digest,
    point_on_short_curve,
    stable_height_rank,
)


Q = Fraction
ROOTS = (346, 260, 255, 146, 55, 0)
CONSTRUCTION = SixRootMestreConstruction(tuple(Q(root) for root in ROOTS))
PAPER_PARAMETER = Q(5081, 94)
PARAMETER_T = 2 * PAPER_PARAMETER
UNIFORM_HEIGHT = 1_000_000
CHART_HEIGHT = 50_000
CHART_COUNT = 48
EXPECTED_UNIFORM_SIGNED_POINTS = 50
EXPECTED_UNIFORM_ABSCISSAS = 25
EXPECTED_TOTAL_ABSCISSAS = 34
EXPECTED_POOL_SIZE = 34
EXPECTED_POOL_SHA256 = (
    "ebc619800df687d7efaa3c5540492fdf22d874c78c1f86d835058bbfc281758d"
)
EXPECTED_SELECTED_INDICES = (
    1,
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    9,
    10,
    11,
    13,
    14,
    15,
    17,
    18,
    19,
    20,
    23,
    25,
)
EXPECTED_SELECTED_SHA256 = (
    "e951f0eee8e59502cfb949e4b91c92fea28061bcddcb37c662e69e9ffe48074c"
)
EXPECTED_SATURATED_BASIS_SHA256 = (
    "cae79a6eb2ab158e601c7536d24b8bded78e572c5e925d9ee9e0f53dc531eeac"
)
EXPECTED_CERTIFICATE_PRIMES = (
    11,
    19,
    41,
    53,
    59,
    67,
    71,
    79,
    97,
    101,
    103,
    109,
    113,
    131,
    173,
)
EXPECTED_CONDUCTOR = int(
    "4739512365768104141634183882739432010081578062727282562384233196211945306960"
)
EXPECTED_MINIMAL_MODEL = (
    0,
    -1,
    0,
    -47433564031723813622493745045480,
    124574716166660957649866283198474133374724238272,
)
EXPECTED_MINIMAL_DISCRIMINANT = int(
    "126112692807661787461041931534709893261155274569961052591643175490028008470229912366239440000000"
)
TARGET_LOG_CONDUCTOR = Decimal("182.72")
CERTIFIED_RANK_LOWER_BOUND = 20
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/certify_nagao_rank20_t5081.py"
)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exact_curve_data() -> tuple[
    tuple[Fraction, ...],
    tuple[tuple[Fraction, Fraction], ...],
    tuple[Fraction, ...],
]:
    quartic = primitive_quartic_coefficients(CONSTRUCTION, PARAMETER_T)
    visible_quartic = primitive_visible_points(CONSTRUCTION, PARAMETER_T)
    if len(visible_quartic) != 12 or len({point[0] for point in visible_quartic}) != 12:
        raise AssertionError("the twelve visible quartic points collided")
    if any(point[1] ** 2 != quartic_value(quartic, point[0]) for point in visible_quartic):
        raise AssertionError("a visible point missed the primitive quartic")
    visible_jacobian = tuple(
        quartic_point_to_short_jacobian(CONSTRUCTION, PARAMETER_T, point)
        for point in visible_quartic
    )
    coefficients = short_jacobian_coefficients(CONSTRUCTION, PARAMETER_T)
    if any(not point_on_short_curve(coefficients, point) for point in visible_jacobian):
        raise AssertionError("a visible point missed the exact Jacobian")
    return quartic, visible_jacobian, coefficients


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--uniform-timeout", type=float, default=40.0)
    parser.add_argument("--chart-timeout", type=float, default=60.0)
    parser.add_argument("--height-timeout", type=float, default=30.0)
    parser.add_argument("--saturation-timeout", type=float, default=30.0)
    parser.add_argument("--conductor-timeout", type=float, default=30.0)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument("--certificate-prime-bound", type=int, default=1000)
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            root
            / "artifacts/generated-results/elliptic-curves/elliptic_nagao_rank20_t5081_rank20_certificate.json"
        ),
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    for name in (
        "uniform_timeout",
        "chart_timeout",
        "height_timeout",
        "saturation_timeout",
        "conductor_timeout",
    ):
        if not 0 < getattr(args, name) <= 120:
            raise SystemExit(f"--{name.replace('_', '-')} must be in (0,120]")
    if args.stack_bytes < 64_000_000:
        raise SystemExit("--stack-bytes is too small")
    if not 173 <= args.certificate_prime_bound <= 10_000:
        raise SystemExit("--certificate-prime-bound must be in [173,10000]")

    quartic, visible_jacobian, coefficients = exact_curve_data()
    uniform_by_chart, uniform_ms, uniform_wall = run_mobius_charts(
        quartic,
        (("identity", (1, 0, 0, 1)),),
        height_bound=UNIFORM_HEIGHT,
        timeout=args.uniform_timeout,
        stack_bytes=args.stack_bytes,
    )
    uniform_raw = uniform_by_chart["identity"]
    uniform_points = signless_quartic_points(uniform_raw)
    if len(uniform_raw) != EXPECTED_UNIFORM_SIGNED_POINTS:
        raise AssertionError("the pinned uniform signed-point count changed")
    if len(uniform_points) != EXPECTED_UNIFORM_ABSCISSAS:
        raise AssertionError("the pinned uniform abscissa count changed")
    if any(point[1] ** 2 != quartic_value(quartic, point[0]) for point in uniform_points):
        raise AssertionError("a uniform point missed the exact quartic")

    charts = optimized_cross_ratio_charts(
        tuple(point[0] for point in uniform_points), count=CHART_COUNT
    )
    raw_by_chart, chart_ms, chart_wall = run_mobius_charts(
        quartic,
        tuple((chart.identifier, chart.matrix) for chart in charts),
        height_bound=CHART_HEIGHT,
        timeout=args.chart_timeout,
        stack_bytes=args.stack_bytes,
    )
    quartic_by_x = {point[0]: point for point in uniform_points}
    chart_records = []
    for chart in charts:
        mapped = map_run_points(
            quartic, raw_by_chart[chart.identifier], chart.matrix
        )
        before = len(quartic_by_x)
        for point in mapped:
            quartic_by_x.setdefault(point[0], point)
        chart_records.append(
            {
                "chart_id": chart.identifier,
                "matrix": list(chart.matrix),
                "raw_signed_points": len(raw_by_chart[chart.identifier]),
                "mapped_distinct_abscissas": len(mapped),
                "new_global_abscissas": len(quartic_by_x) - before,
                "pari_milliseconds": chart_ms[chart.identifier],
            }
        )
    if len(quartic_by_x) != EXPECTED_TOTAL_ABSCISSAS:
        raise AssertionError("the cross-ratio abscissa union changed")

    pool = list(visible_jacobian)
    seen_jacobian_x = {point[0] for point in visible_jacobian}
    new_images = []
    for quartic_point in sorted(quartic_by_x.values()):
        if quartic_point[1] == 0:
            continue
        image = quartic_point_to_short_jacobian(
            CONSTRUCTION, PARAMETER_T, quartic_point
        )
        if not point_on_short_curve(coefficients, image):
            raise AssertionError("a searched point missed the exact Jacobian")
        if image[0] in seen_jacobian_x:
            continue
        seen_jacobian_x.add(image[0])
        new_images.append(image)
        pool.append(image)
    pool_tuple = tuple(pool)
    if len(pool_tuple) != EXPECTED_POOL_SIZE or point_digest(pool_tuple) != EXPECTED_POOL_SHA256:
        raise AssertionError("the pinned exact Jacobian pool changed")

    height_runs = height_matrix_replay(
        coefficients,
        pool_tuple,
        precisions=(72, 120),
        timeout=args.height_timeout,
        stack_bytes=args.stack_bytes,
    )
    if stable_height_rank(height_runs) != CERTIFIED_RANK_LOWER_BOUND:
        raise AssertionError("the stable numerical selection rank changed")
    selected_indices = tuple(height_runs[-1]["subset_indices_one_based"])
    if selected_indices != EXPECTED_SELECTED_INDICES:
        raise AssertionError("the pinned numerical subset changed")
    selected = tuple(pool_tuple[index - 1] for index in selected_indices)
    if point_digest(selected) != EXPECTED_SELECTED_SHA256:
        raise AssertionError("the selected-point digest changed")

    saturated_basis, saturation = saturate_exact_basis(
        coefficients,
        selected,
        prime_bound=20,
        timeout=args.saturation_timeout,
        stack_bytes=args.stack_bytes,
    )
    saturated_digest = point_digest(saturated_basis)
    if (
        len(saturated_basis) != CERTIFIED_RANK_LOWER_BOUND
        or saturated_digest != EXPECTED_SATURATED_BASIS_SHA256
    ):
        raise AssertionError("the pinned saturated basis changed")
    if any(not point_on_short_curve(coefficients, point) for point in saturated_basis):
        raise AssertionError("a saturated basis point missed the exact curve")

    signatures = find_mod2_reduction_certificate(
        coefficients,
        saturated_basis,
        prime_bound=args.certificate_prime_bound,
    )
    exact_binary_rank = combined_mod2_rank(signatures, len(saturated_basis))
    certificate_primes = tuple(signature.prime for signature in signatures)
    if exact_binary_rank != CERTIFIED_RANK_LOWER_BOUND:
        raise AssertionError("finite reductions did not certify all 20 points")
    if certificate_primes != EXPECTED_CERTIFICATE_PRIMES:
        raise AssertionError("the deterministic reduction certificate changed")
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
    if int(conductor["root_number"]) != 1:
        raise AssertionError("the exact root number changed")
    if Decimal(conductor["log_conductor"]) >= TARGET_LOG_CONDUCTOR:
        raise AssertionError("the specialization crossed the conductor target")
    exact_log_bound = exact_log_conductor_certificate(EXPECTED_CONDUCTOR)

    artifact: dict[str, Any] = {
        "schema_version": 1,
        "status": "exact_rank_at_least_20_certificate_complete",
        "theorem": (
            "Nagao's section-7 specialization with roots "
            "(346,260,255,146,55,0) and constructor T=5081/47 has "
            "Mordell-Weil rank at least 20 over Q and log conductor below 182.72."
        ),
        "candidate": {
            "paper_parameter_t": rational_to_string(PAPER_PARAMETER),
            "constructor_parameter_T": rational_to_string(PARAMETER_T),
            "factor_two_convention_checked": PARAMETER_T == 2 * PAPER_PARAMETER,
            "roots": list(ROOTS),
            "primitive_quartic_square_scale": rational_to_string(
                CONSTRUCTION.quartic_square_scale
            ),
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
        "point_search": {
            "uniform_height": UNIFORM_HEIGHT,
            "uniform_signed_points": len(uniform_raw),
            "uniform_distinct_abscissas": len(uniform_points),
            "uniform_pari_milliseconds": uniform_ms["identity"],
            "uniform_wall_seconds": uniform_wall,
            "cross_ratio_chart_count": len(charts),
            "cross_ratio_height": CHART_HEIGHT,
            "cross_ratio_wall_seconds": chart_wall,
            "chart_records": chart_records,
            "total_distinct_quartic_abscissas": len(quartic_by_x),
            "visible_jacobian_points": len(visible_jacobian),
            "new_jacobian_sign_pairs": len(new_images),
            "exact_jacobian_pool_count": len(pool_tuple),
            "exact_jacobian_pool_sha256": point_digest(pool_tuple),
            "all_memberships_checked_exactly": True,
        },
        "height_selection": {
            "runs": list(height_runs),
            "stable_numerical_rank": CERTIFIED_RANK_LOWER_BOUND,
            "selected_pool_indices_one_based": list(selected_indices),
            "selected_point_sha256": point_digest(selected),
            "selection_is_not_certification": True,
        },
        "exact_rank_certificate": {
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
        "script_sha256": sha256_file(Path(__file__).resolve()),
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
