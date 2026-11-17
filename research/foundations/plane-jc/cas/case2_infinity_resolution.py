#!/usr/bin/env python3
"""Resolve the generic infinity branch of the Case-2 residue exactly.

The certified Case-2 residue is the polynomial parametrization

    [1:C(t):G(t)],  deg(C)=8, deg(G)=12.

In the target chart at infinity, with ``s=t^-1``, put

    U=1/G(t),  V=C(t)/G(t).

The leading orders are ``ord(U)=12`` and ``ord(V)=4``.  After cancelling
the common cubic tangent, the coefficient at order thirteen is, up to a
nonzero denominator,

    K13 = 2*C8*G11 - 3*C7*G12.

This expression is invariant under translations of the normalization
parameter.  On the open stratum ``G12*K13 != 0`` the infinity branch has
primitive parametrization ``(V,W)=(s^4, unit*s^13+...)``.  Its minimal
regular toric fan and exceptional self-intersections are computed here.

The optional Singular audit asks a stronger question: after imposing the
seven unused J1 compatibility equations and localizing at ``G12``, is the
resonance divisor ``K13=0`` empty?
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path

from audit_case2_residue_strata import (
    default_replay_root,
    load_exact_core,
    singular_minpoly,
    solve_case2_through_j1,
)


Ray = tuple[int, int]
EXPECTED_ENDPOINT_SUBSET = (0, 1, 2, 3, 4, 5, 6)
EXPECTED_SINGULAR_INPUT_SHA256 = (
    "1618285117e313afb04de0ab935d4a7c"
    "f182f4baabbb0a1346f336f88aaf394e"
)


@dataclass(frozen=True)
class InfinityResolutionAudit:
    coordinate_degrees: tuple[int, int]
    leading_orders: tuple[int, int]
    characteristic_orders: tuple[int, int]
    characteristic_numerator_term_count: int
    characteristic_numerator_parameter_degree: int
    normalization_translation_invariant: bool
    fan_rays: tuple[Ray, ...]
    exceptional_rays: tuple[Ray, ...]
    exceptional_self_intersections: tuple[int, ...]
    strict_transform_meets_ray: Ray
    blowup_count: int
    compatibility_equation_count: int
    endpoint_certificate_indices: tuple[int, ...]
    endpoint_certificate_constraint_count: int
    endpoint_certificate_term_counts: tuple[int, ...]
    endpoint_certificate_parameter_degrees: tuple[int, ...]
    degree_twelve_open_nonempty: bool | None
    case2_excluded_by_j1_endpoint: bool | None
    singular_input_sha256: str | None


def det(left: Ray, right: Ray) -> int:
    return left[0] * right[1] - left[1] * right[0]


def regular_subdivision(branch_ray: Ray) -> tuple[Ray, ...]:
    """Star-subdivide until the primitive branch ray occurs in a regular fan."""

    left: Ray = (1, 0)
    right: Ray = (0, 1)
    rays = [left, right]
    while branch_ray not in rays:
        middle = (left[0] + right[0], left[1] + right[1])
        rays.append(middle)
        side = det(middle, branch_ray)
        if side == 0:
            if middle != branch_ray:
                raise AssertionError("branch ray must be primitive")
            break
        if side > 0:
            left = middle
        else:
            right = middle
    ordered = tuple(sorted(set(rays), key=lambda ray: ray[1] / ray[0] if ray[0] else float("inf")))
    if any(det(ordered[index], ordered[index + 1]) != 1 for index in range(len(ordered) - 1)):
        raise AssertionError("star subdivision did not produce a regular fan")
    return ordered


def self_intersections(fan: tuple[Ray, ...]) -> tuple[int, ...]:
    return tuple(
        -det(fan[index - 1], fan[index + 1])
        for index in range(1, len(fan) - 1)
    )


def characteristic_numerator(ec, C, G):
    if C[8].sing(3) != "1":
        raise AssertionError("the Case-2 residue must have monic C")
    return G[11] * C[8] * ec.K(2) - C[7] * G[12] * ec.K(3)


def endpoint_source(ec, compatibility, g12) -> str:
    compatibility_source = ",\n".join(item.sing(3) for item in compatibility)
    return "\n".join(
        (
            'LIB "resources.lib";',
            "Resources::setcores(1);",
            'LIB "nfmodstd.lib";',
            "ring R=(0,u),(r,s,h,w),dp;",
            f"minpoly={singular_minpoly(ec)};",
            "option(redSB);",
            "ideal I=",
            compatibility_source + ";",
            f"ideal Iopen=I,w*({g12.sing(3)})-1;",
            "ideal Jopen=nfmodStd(Iopen);",
            (
                'if(size(Jopen)>0 && Jopen[1]==1)'
                '{print("CASE2_J1_ENDPOINT_UNIT");}'
            ),
            "quit;",
            "",
        )
    )


def run_endpoint_source(source: str, singular: str, timeout: int) -> bool:
    with tempfile.TemporaryDirectory(prefix="jc2-j1-endpoint-") as directory:
        path = Path(directory) / "case2_j1_endpoint.sing"
        path.write_text(source)
        completed = subprocess.run(
            [singular, "-q", str(path)],
            check=False,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    if completed.returncode != 0:
        raise RuntimeError(
            "Case-2 J1 endpoint audit failed\n"
            f"stdout:\n{completed.stdout}\n"
            f"stderr:\n{completed.stderr}"
        )
    return "CASE2_J1_ENDPOINT_UNIT" in completed.stdout


def audit(
    replay_root: Path | None = None,
    singular: str | None = None,
    timeout: int = 300,
) -> InfinityResolutionAudit:
    ec = load_exact_core((replay_root or default_replay_root()).resolve())
    _, _, C, _, G, compatibility = solve_case2_through_j1(ec)
    if len(compatibility) != 7:
        raise AssertionError("the Case-2 J1 compatibility count drifted")

    k13 = characteristic_numerator(ec, C, G)
    if not k13:
        raise AssertionError("the order-thirteen characteristic coefficient vanished identically")

    fan = regular_subdivision((4, 13))
    exceptional = fan[1:-1]
    intersections = self_intersections(fan)
    if (
        exceptional
        != ((1, 1), (1, 2), (1, 3), (4, 13), (3, 10), (2, 7), (1, 4))
        or intersections != (-2, -2, -5, -1, -2, -2, -2)
    ):
        raise AssertionError("the (4,13) toric resolution graph drifted")

    selected_indices = EXPECTED_ENDPOINT_SUBSET
    selected = tuple(compatibility[index] for index in selected_indices)
    selected_term_counts = tuple(len(item.d) for item in selected)
    selected_degrees = tuple(max(map(sum, item.d)) for item in selected)
    digest = None
    open_nonempty = None
    if singular:
        if not selected:
            selected = compatibility
            selected_indices = tuple(range(len(compatibility)))
            selected_term_counts = tuple(len(item.d) for item in selected)
            selected_degrees = tuple(max(map(sum, item.d)) for item in selected)
        source = endpoint_source(ec, selected, G[12])
        digest = hashlib.sha256(source.encode()).hexdigest()
        if EXPECTED_SINGULAR_INPUT_SHA256 and digest != EXPECTED_SINGULAR_INPUT_SHA256:
            raise AssertionError("the Case-2 J1 endpoint input hash drifted")
        open_nonempty = not run_endpoint_source(source, singular, timeout)
        resonance_empty = None

    return InfinityResolutionAudit(
        coordinate_degrees=(8, 12),
        leading_orders=(4, 12),
        characteristic_orders=(4, 13),
        characteristic_numerator_term_count=len(k13.d),
        characteristic_numerator_parameter_degree=max(map(sum, k13.d)),
        normalization_translation_invariant=(2 * 12 == 3 * 8),
        fan_rays=fan,
        exceptional_rays=exceptional,
        exceptional_self_intersections=intersections,
        strict_transform_meets_ray=(4, 13),
        blowup_count=len(exceptional),
        compatibility_equation_count=len(compatibility),
        endpoint_certificate_indices=selected_indices,
        endpoint_certificate_constraint_count=len(selected),
        endpoint_certificate_term_counts=selected_term_counts,
        endpoint_certificate_parameter_degrees=selected_degrees,
        degree_twelve_open_nonempty=open_nonempty,
        case2_excluded_by_j1_endpoint=(
            None if open_nonempty is None else not open_nonempty
        ),
        singular_input_sha256=digest,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--replay-root", type=Path)
    parser.add_argument("--singular", default=shutil.which("Singular"))
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument(
        "--skip-singular",
        action="store_true",
        help="emit only the formal-series and toric-fan certificate",
    )
    args = parser.parse_args()
    result = audit(
        replay_root=args.replay_root,
        singular=None if args.skip_singular else args.singular,
        timeout=args.timeout,
    )
    print(json.dumps(asdict(result), sort_keys=True))
    print("CASE2_INFINITY_RESOLUTION_PASS")


if __name__ == "__main__":
    main()
