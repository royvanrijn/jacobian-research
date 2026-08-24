#!/usr/bin/env python3
"""Numerically audit the displayed rank-eleven subgroup at diameter 235.

This is deliberately a seed-fibre audit, not a Shioda calculation.  At the
regular fibre ``p=-294, T=2`` the first ten visible points and the first
selected affine point are exactly independent by the finite-reduction
certificate in the component verifier.  Here we replay the real height matrix
at two precisions and make one bounded PARI ``ellsaturation`` call.  The
result is diagnostic: PARI documents ``ellsaturation`` under a finite-index
hypothesis, and neither calculation determines the generic Mordell--Weil
group, intersections, or reducible-fibre corrections.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
from typing import Any

from extend_nagao_u42_frontier import saturate_exact_basis
from mestre_root_tuples import SixRootMestreConstruction
from screen_mestre_fermigier_two_section_escape import square_root
from search_mestre_root_tuple_scale import (
    height_matrix_replay,
    primitive_visible_points,
    quartic_point_to_jacobian,
    quartic_value,
)
from verify_mestre_diameter235_eight_companion_component import (
    ROOTS_AT_SEED,
    SEED_PARAMETER,
    SEED_SECTIONS,
)


Q = Fraction
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/audit_mestre_diameter235_displayed_lattice.py"
)


def displayed_basis() -> tuple[tuple[Q, ...], tuple[tuple[Q, Q], ...]]:
    """Return the exact independent eleven-point seed basis.

    Independence is not inferred here: the component verifier supplies its
    exact finite-reduction certificate.  This helper only reproduces the same
    points before calling numerical lattice routines.
    """

    construction = SixRootMestreConstruction(ROOTS_AT_SEED)
    parameter = Q(2)
    visible = tuple(
        quartic_point_to_jacobian(construction, parameter, point)
        for point in primitive_visible_points(construction, parameter)
    )
    intercept, slope = SEED_SECTIONS[0]
    abscissa = intercept + slope * parameter
    ordinate = square_root(
        quartic_value(construction.primitive_quartic_coefficients(parameter), abscissa)
    )
    if ordinate is None:
        raise AssertionError("the selected affine point lost its exact ordinate")
    affine = quartic_point_to_jacobian(construction, parameter, (abscissa, ordinate))
    return construction.primitive_jacobian_coefficients(parameter), (*visible[:10], affine)


def audit(
    *,
    precisions: tuple[int, ...],
    height_timeout: float,
    saturation_bound: int,
    saturation_timeout: float,
    stack_bytes: int,
) -> dict[str, Any]:
    coefficients, points = displayed_basis()
    heights = height_matrix_replay(
        coefficients,
        points,
        precisions=precisions,
        timeout=height_timeout,
        stack_bytes=stack_bytes,
    )
    _, saturation = saturate_exact_basis(
        coefficients,
        points,
        prime_bound=saturation_bound,
        timeout=saturation_timeout,
        stack_bytes=stack_bytes,
    )
    return {
        "status": "completed seed-fibre numerical lattice audit",
        "reproducing_command": REPRODUCING_COMMAND,
        "specialization": {"p": str(SEED_PARAMETER), "T": "2"},
        "basis": {
            "point_count": len(points),
            "description": "first ten visible points followed by the first selected affine point",
            "exact_independence_source": "EC-MD235 finite-reduction certificate",
        },
        "height_matrix_replay": list(heights),
        "bounded_pari_saturation": {
            key: saturation[key]
            for key in (
                "prime_bound_strict_upper_limit",
                "input_point_count",
                "returned_point_count",
                "exact_returned_points_on_curve",
                "original_height_determinant",
                "saturated_height_determinant",
                "height_determinant_ratio",
                "saturated_basis_sha256",
                "scope_warning",
            )
        },
        "not_established": [
            "a saturated Mordell-Weil basis at this fibre or generically",
            "a full Shioda Gram matrix or the affine-pair intersection number",
            "generic rank at least 12 or 14, or an upper bound for the full Mordell-Weil group",
        ],
    }


def parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    generated = root / "artifacts" / "generated-results"
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--precisions", default="72,120")
    result.add_argument("--height-timeout", type=float, default=60.0)
    result.add_argument("--saturation-bound", type=int, default=20)
    result.add_argument("--saturation-timeout", type=float, default=120.0)
    result.add_argument("--stack-bytes", type=int, default=1_024_000_000)
    result.add_argument(
        "--output",
        type=Path,
        default=generated / "elliptic_mestre_diameter235_displayed_lattice_seed_audit.json",
    )
    return result


def main() -> None:
    args = parser().parse_args()
    precisions = tuple(int(value) for value in args.precisions.split(",") if value)
    if len(precisions) < 2 or min(precisions) < 16:
        raise SystemExit("--precisions needs at least two values of at least 16 digits")
    if min(args.height_timeout, args.saturation_timeout) <= 0:
        raise SystemExit("timeouts must be positive")
    if args.saturation_bound < 3 or args.stack_bytes < 8_000_000:
        raise SystemExit("invalid saturation bound or PARI stack size")
    rendered = json.dumps(
        audit(
            precisions=precisions,
            height_timeout=args.height_timeout,
            saturation_bound=args.saturation_bound,
            saturation_timeout=args.saturation_timeout,
            stack_bytes=args.stack_bytes,
        ),
        indent=2,
        sort_keys=True,
    ) + "\n"
    args.output.write_text(rendered)


if __name__ == "__main__":
    main()
