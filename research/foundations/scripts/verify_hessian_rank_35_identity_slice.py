#!/usr/bin/env python3
"""Certify the rank-35 nonhomogeneous HN identity-slice witness.

The homogeneous 22-variable ``qb+x2s`` cubic collision has an identity
output.  Restricting that coordinate to one gives a 21-variable
nonhomogeneous map ``X+K``.  Its cotangent potential ``y.K(x)`` has 42
variables and degree at most four.

This script certifies the two generic ranks over characteristic zero with
Singular syzygies and exact good-prime lower specializations.  The all-order
HN and Vanishing arguments are written in
``extended-geometry/HESSIAN_RANK_35_IDENTITY_SLICE.md``.
"""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path

from audit_fixed_rank_hessian_witness import (
    cotangent_hessian,
    decode_h,
    deterministic_point,
    exact_singular_certificate,
    specialization_profile,
)


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "hessian_rank_reduced_bcw_22_counterexample.json"
)
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "hessian_rank_35_identity_slice_counterexample.json"
)
PRIME = 1_000_003
SEEDS = (20_260_728, 20_260_729, 20_260_730)


def identity_slice(
    homogeneous_h: list[dict[tuple[int, ...], Fraction]],
) -> list[dict[tuple[int, ...], Fraction]]:
    """Set the final source coordinate to one and remove its zero output."""

    assert homogeneous_h[-1] == {}
    sliced: list[dict[tuple[int, ...], Fraction]] = []
    for component in homogeneous_h[:-1]:
        result: dict[tuple[int, ...], Fraction] = {}
        for exponent, coefficient in component.items():
            reduced = exponent[:-1]
            result[reduced] = result.get(reduced, Fraction(0)) + coefficient
        sliced.append(
            {
                exponent: coefficient
                for exponent, coefficient in result.items()
                if coefficient
            }
        )
    return sliced


def evaluate(
    component: dict[tuple[int, ...], Fraction],
    point: list[Fraction],
) -> Fraction:
    answer = Fraction(0)
    for exponent, coefficient in component.items():
        monomial = coefficient
        for value, power in zip(point, exponent):
            monomial *= value**power
        answer += monomial
    return answer


def encode(
    components: list[dict[tuple[int, ...], Fraction]],
) -> list[list[dict[str, object]]]:
    return [
        [
            {
                "coefficient": str(coefficient),
                "monomial": [
                    [index, power]
                    for index, power in enumerate(exponent)
                    if power
                ],
            }
            for exponent, coefficient in sorted(component.items(), reverse=True)
        ]
        for component in components
    ]


def main() -> None:
    stored = json.loads(SOURCE.read_text())
    assert stored["dimension"] == 22
    assert stored["H"][-1] == []
    assert stored["statistics"]["generic_cotangent_hessian_rank_over_QQ_xy"] == 37
    assert stored["statistics"]["nilpotency_index_JH"] == 18

    points = [
        [Fraction(value) for value in point]
        for point in stored["collision_points"]
    ]
    assert [point[-1] for point in points] == [1, 1, 1]
    sliced_points = [point[:-1] for point in points]
    assert len(set(map(tuple, sliced_points))) == 3

    homogeneous_h = decode_h(stored)
    sliced_h = identity_slice(homogeneous_h)
    assert len(sliced_h) == 21
    degrees = sorted(
        {
            sum(exponent)
            for component in sliced_h
            for exponent in component
        }
    )
    assert degrees == [1, 2, 3]

    images = [
        [
            coordinate + evaluate(component, point)
            for coordinate, component in zip(point, sliced_h)
        ]
        for point in sliced_points
    ]
    assert images[0] == images[1] == images[2]
    assert [point[0] for point in sliced_points] == [0, 1, -1]

    jacobian, upper_left, hessian = cotangent_hessian(sliced_h)
    profiles = [
        specialization_profile(jacobian, upper_left, hessian, seed)
        for seed in SEEDS
    ]
    assert profiles == [(17, 14, 35)] * len(SEEDS)

    jacobian_point = deterministic_point(21, PRIME, SEEDS[0])
    jacobian_syzygies, jacobian_rank, jacobian_kernel_check = (
        exact_singular_certificate(
            jacobian,
            jacobian_point,
            expected_kernel_rank=4,
        )
    )
    assert (jacobian_syzygies, jacobian_rank, jacobian_kernel_check) == (
        4,
        17,
        True,
    )

    hessian_point = deterministic_point(42, PRIME, SEEDS[0])
    hessian_syzygies, hessian_rank, hessian_kernel_check = (
        exact_singular_certificate(
            hessian,
            hessian_point,
            expected_kernel_rank=7,
        )
    )
    assert (hessian_syzygies, hessian_rank, hessian_kernel_check) == (
        13,
        35,
        True,
    )

    artifact = {
        "format": "hessian-rank-35-identity-slice-v1",
        "field": "QQ for the slice and QQ(I) for the HN change",
        "source_artifact": str(SOURCE.relative_to(ROOT)),
        "source_dimension": 22,
        "identity_coordinate": 21,
        "identity_slice_value": "1",
        "slice_dimension": 21,
        "slice_map": "X+K(X), where K(X)=H_0,...,H_20 evaluated at x_21=1",
        "K": encode(sliced_h),
        "collision_points": [
            [str(value) for value in point] for point in sliced_points
        ],
        "common_image": [str(value) for value in images[0]],
        "distinguished_collision_coordinate": 0,
        "distinguished_values": ["0", "1", "-1"],
        "cotangent_potential": {
            "variables": "x_0,...,x_20,y_0,...,y_20",
            "definition": "p(x,y)=sum_i y_i K_i(x)",
            "dimension": 42,
            "degrees": [2, 3, 4],
            "generic_hessian_rank": 35,
        },
        "HN_potential": {
            "variables": "u_0,...,u_20,v_0,...,v_20",
            "definition": (
                "P(u,v)=(1/2)*(u-I*v).K(u+I*v), with the slice value "
                "x_21=1 already applied"
            ),
            "dimension": 42,
            "degrees": [2, 3, 4],
            "generic_hessian_rank": 35,
            "certificate": (
                "Hess(P) is nilpotent; I+grad(P) is noninjective; hence "
                "Delta^m(P^(m+1)) is nonzero for infinitely many m"
            ),
            "homogeneity_warning": (
                "P is not homogeneous, so this does not lower the "
                "homogeneous-quartic rank or dimension frontiers"
            ),
        },
        "exact_rank_certificates": {
            "generic_rank_JK": {
                "rank": 17,
                "syzygy_generators": jacobian_syzygies,
                "independent_kernel_columns": 4,
                "exact_specialization_rank": jacobian_rank,
            },
            "generic_rank_Hess_p": {
                "rank": 35,
                "syzygy_generators": hessian_syzygies,
                "independent_kernel_columns": 7,
                "exact_specialization_rank": hessian_rank,
            },
        },
        "good_prime_profiles": {
            "prime": PRIME,
            "seeds": list(SEEDS),
            "block_ranks": [list(profile) for profile in profiles],
            "block_order": ["rank_JK", "rank_sum_y_Hess_K", "rank_Hess_p"],
        },
        "status": (
            "exact nonhomogeneous HN rank witness; no change to the "
            "homogeneous quartic HN rank-37 or dimension-42 endpoints"
        ),
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")
    print("PASS rank-35 slice: exact 21D collision with identity slice x_21=1")
    print("PASS rank-35 slice: generic rank(JK)=17 over QQ(x)")
    print("PASS rank-35 slice: generic rank Hess(y.K)=35 over QQ(x,y)")
    print("PASS rank-35 slice: 42-variable degree-2/3/4 nonhomogeneous HN witness")
    print(f"PASS rank-35 slice: wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
