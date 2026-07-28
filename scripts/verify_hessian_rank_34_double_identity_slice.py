#!/usr/bin/env python3
"""Certify the 40-variable rank-34 nonhomogeneous HN witness.

The rank-35 identity slice has the exact output relation

    K_9 - 3*K_1 - K_6 = 0.

Hence ``-3*x_1-x_6+x_9`` is a second identity output.  Its three collision
points lie in the zero slice.  Eliminating ``x_9=3*x_1+x_6`` gives a
20-variable nilpotent-Jacobian collision ``X+L``.  The cotangent potential
``y.L(x)`` has 40 variables and exact generic Hessian rank 34.
"""

from __future__ import annotations

from fractions import Fraction
from math import comb
import json
from pathlib import Path

from audit_fixed_rank_hessian_witness import (
    constant_kernel_coefficient_rank,
    cotangent_hessian,
    deterministic_point,
    exact_singular_certificate,
    specialization_profile,
)


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "hessian_rank_35_identity_slice_counterexample.json"
)
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "hessian_rank_34_double_identity_slice_counterexample.json"
)
PRIME = 1_000_003
SEEDS = (20_260_728, 20_260_729, 20_260_730)
REMOVED_COORDINATE = 9
KEPT_COORDINATES = tuple(
    index for index in range(21) if index != REMOVED_COORDINATE
)

Exponent = tuple[int, ...]
Polynomial = dict[Exponent, Fraction]


def add_term(
    polynomial: Polynomial,
    exponent: Exponent,
    coefficient: Fraction,
) -> None:
    value = polynomial.get(exponent, Fraction(0)) + coefficient
    if value:
        polynomial[exponent] = value
    else:
        polynomial.pop(exponent, None)


def decode(components: list[list[dict[str, object]]]) -> list[Polynomial]:
    result: list[Polynomial] = []
    for component in components:
        polynomial: Polynomial = {}
        for term in component:
            exponent = [0] * 21
            for variable, power in term["monomial"]:
                exponent[variable] = power
            add_term(
                polynomial,
                tuple(exponent),
                Fraction(term["coefficient"]),
            )
        result.append(polynomial)
    return result


def linear_combination(
    terms: tuple[tuple[Fraction, Polynomial], ...],
) -> Polynomial:
    result: Polynomial = {}
    for scalar, polynomial in terms:
        for exponent, coefficient in polynomial.items():
            add_term(result, exponent, scalar * coefficient)
    return result


def eliminate_x9(polynomial: Polynomial) -> Polynomial:
    """Substitute x_9=3*x_1+x_6 and drop coordinate nine."""

    result: Polynomial = {}
    for exponent, coefficient in polynomial.items():
        power = exponent[REMOVED_COORDINATE]
        reduced = [exponent[index] for index in KEPT_COORDINATES]
        for x6_power in range(power + 1):
            expanded = reduced[:]
            expanded[1] += power - x6_power
            expanded[6] += x6_power
            scalar = (
                coefficient
                * comb(power, x6_power)
                * 3 ** (power - x6_power)
            )
            add_term(result, tuple(expanded), scalar)
    return result


def evaluate(polynomial: Polynomial, point: list[Fraction]) -> Fraction:
    result = Fraction(0)
    for exponent, coefficient in polynomial.items():
        term = coefficient
        for value, power in zip(point, exponent):
            term *= value**power
        result += term
    return result


def encode(components: list[Polynomial]) -> list[list[dict[str, object]]]:
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
    assert stored["format"] == "hessian-rank-35-identity-slice-v1"
    assert stored["slice_dimension"] == 21
    assert stored["exact_rank_certificates"]["generic_rank_JK"]["rank"] == 17
    assert (
        stored["exact_rank_certificates"]["generic_rank_Hess_p"]["rank"]
        == 35
    )
    k = decode(stored["K"])
    assert len(k) == 21

    relation = linear_combination(
        (
            (Fraction(1), k[9]),
            (Fraction(-3), k[1]),
            (Fraction(-1), k[6]),
        )
    )
    assert relation == {}

    points = [
        [Fraction(value) for value in point]
        for point in stored["collision_points"]
    ]
    identity_values = [
        -3 * point[1] - point[6] + point[9] for point in points
    ]
    assert identity_values == [0, 0, 0]
    sliced_points = [
        [point[index] for index in KEPT_COORDINATES] for point in points
    ]
    assert len(set(map(tuple, sliced_points))) == 3

    ell = [eliminate_x9(k[index]) for index in KEPT_COORDINATES]
    assert len(ell) == 20
    degrees = sorted(
        {
            sum(exponent)
            for component in ell
            for exponent in component
        }
    )
    assert degrees == [1, 2, 3]

    images = [
        [
            coordinate + evaluate(component, point)
            for coordinate, component in zip(point, ell)
        ]
        for point in sliced_points
    ]
    assert images[0] == images[1] == images[2]
    assert [point[0] for point in sliced_points] == [0, 1, -1]

    jacobian, upper_left, hessian = cotangent_hessian(ell)
    profiles = [
        specialization_profile(jacobian, upper_left, hessian, seed)
        for seed in SEEDS
    ]
    assert profiles == [(17, 12, 34)] * len(SEEDS)

    jacobian_syzygies, jacobian_rank, jacobian_kernel_check = (
        exact_singular_certificate(
            jacobian,
            deterministic_point(20, PRIME, SEEDS[0]),
            expected_kernel_rank=3,
        )
    )
    assert (jacobian_syzygies, jacobian_rank, jacobian_kernel_check) == (
        3,
        17,
        True,
    )

    hessian_syzygies, hessian_rank, hessian_kernel_check = (
        exact_singular_certificate(
            hessian,
            deterministic_point(40, PRIME, SEEDS[0]),
            expected_kernel_rank=6,
        )
    )
    assert (hessian_syzygies, hessian_rank, hessian_kernel_check) == (
        12,
        34,
        True,
    )
    assert constant_kernel_coefficient_rank(hessian) == 40

    artifact = {
        "format": "hessian-rank-34-double-identity-slice-v1",
        "field": "QQ for the slice and QQ(I) for the HN change",
        "source_artifact": str(SOURCE.relative_to(ROOT)),
        "source_dimension": 21,
        "identity_linear_form": "-3*x_1-x_6+x_9",
        "identity_output_relation": "K_9-3*K_1-K_6=0",
        "collision_slice_value": "0",
        "elimination": "x_9=3*x_1+x_6",
        "removed_coordinate": REMOVED_COORDINATE,
        "kept_source_coordinates": list(KEPT_COORDINATES),
        "slice_dimension": 20,
        "slice_map": "X+L(X), obtained by the displayed elimination",
        "L": encode(ell),
        "collision_points": [
            [str(value) for value in point] for point in sliced_points
        ],
        "common_image": [str(value) for value in images[0]],
        "distinguished_collision_coordinate": 0,
        "distinguished_values": ["0", "1", "-1"],
        "nilpotent_jacobian_certificate": (
            "after the rational coordinate change, JL is the upper-left "
            "block of the exact nilpotent Jacobian JK on the identity slice"
        ),
        "cotangent_potential": {
            "variables": "x_0,...,x_19,y_0,...,y_19",
            "definition": "p(x,y)=sum_i y_i L_i(x)",
            "dimension": 40,
            "degrees": [2, 3, 4],
            "generic_hessian_rank": 34,
            "constant_hessian_kernel_dimension": 0,
        },
        "HN_potential": {
            "variables": "u_0,...,u_19,v_0,...,v_19",
            "definition": "P(u,v)=(1/2)*(u-I*v).L(u+I*v)",
            "dimension": 40,
            "degrees": [2, 3, 4],
            "generic_hessian_rank": 34,
            "certificate": (
                "Hess(P) is nilpotent; I+grad(P) is noninjective; hence "
                "Delta^m(P^(m+1)) is nonzero for infinitely many m"
            ),
            "homogeneity_warning": (
                "P is not homogeneous, so this does not lower the "
                "homogeneous-quartic rank or dimension frontiers"
            ),
        },
        "ordinary_laplacian_GVC": {
            "dimension": 40,
            "multiplier": "the surviving coordinate x_0",
            "status": (
                "alternative rank-optimized witness at the existing "
                "ordinary-Laplacian dimension bound 40"
            ),
        },
        "exact_rank_certificates": {
            "generic_rank_JL": {
                "rank": 17,
                "syzygy_generators": jacobian_syzygies,
                "independent_kernel_columns": 3,
                "exact_specialization_rank": jacobian_rank,
            },
            "generic_rank_Hess_p": {
                "rank": 34,
                "syzygy_generators": hessian_syzygies,
                "independent_kernel_columns": 6,
                "exact_specialization_rank": hessian_rank,
                "cotangent_kernel_excess": 0,
            },
        },
        "good_prime_profiles": {
            "prime": PRIME,
            "seeds": list(SEEDS),
            "block_ranks": [list(profile) for profile in profiles],
            "block_order": ["rank_JL", "rank_sum_y_Hess_L", "rank_Hess_p"],
        },
        "status": (
            "exact 40-variable nonhomogeneous HN rank-34 witness; no change "
            "to the homogeneous quartic HN rank-37 or dimension-42 endpoints"
        ),
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")
    print("PASS rank-34 slice: exact output relation K_9-3*K_1-K_6=0")
    print("PASS rank-34 slice: exact 20D collision on -3*x_1-x_6+x_9=0")
    print("PASS rank-34 slice: generic rank(JL)=17 over QQ(x)")
    print("PASS rank-34 slice: generic rank Hess(y.L)=34 and excess zero")
    print("PASS rank-34 slice: no further constant Hessian-kernel direction")
    print(f"PASS rank-34 slice: wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
