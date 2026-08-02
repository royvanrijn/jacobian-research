#!/usr/bin/env python3
"""Exact bounded replays for consequences of the homogeneous GVC(3) witness.

The written companion note proves the all-order statements. This checker
verifies, with integer/rational arithmetic only:

* radial padding from Delta^6 to Delta^k for k >= 6;
* the complete positive-phase multiplier ladder x^(2 ell);
* the corresponding rank-one SIC contractions;
* homogeneous Gaussian and higher-sphere transfer formulas.
"""

from __future__ import annotations

import json
from fractions import Fraction
from math import comb, factorial
from pathlib import Path

from verify_gvc3_homogeneous_counterexample import (
    ZERO,
    X,
    Y,
    T,
    add,
    apolar_scalar,
    double_factorial,
    factorial_product,
    gaussian_expectation,
    homogeneous_degree,
    monomial,
    multiply,
    power,
    scale,
)

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "gvc3_homogeneous_spillovers.json"
)


def radial_even_moment(dimension: int, half_degree: int) -> int:
    """E[||G||^(2*half_degree)] for G standard Gaussian in R^dimension."""
    result = 1
    for j in range(half_degree):
        result *= dimension + 2 * j
    return result


def endpoint_constant(order: int) -> Fraction:
    return Fraction(
        2 ** (2 * order) * factorial(2 * order),
        double_factorial(4 * order + 1),
    )


def expected_gaussian_mixed(k: int, order: int, ell: int) -> int:
    value = (
        Fraction(double_factorial(2 * k * order + 2 * ell + 1), 1)
        * comb(order - 1, ell - 1)
        * endpoint_constant(order)
    )
    assert value.denominator == 1
    return value.numerator


def expected_detector(k: int, order: int, ell: int) -> int:
    return (
        2 ** (k * order + ell)
        * factorial(k * order + ell)
        * expected_gaussian_mixed(k, order, ell)
    )


def main() -> None:
    rho = add(monomial((0, 0, 2)), monomial((1, 1, 0)))
    a = add(rho, monomial((2, 0, 0)))
    c = add(
        multiply(Y, power(rho, 2)),
        scale(multiply(multiply(X, power(T, 2)), rho), -2),
        scale(monomial((3, 0, 2)), -1),
    )
    base_p = multiply(a, power(c, 2))
    delta = add(monomial((1, 1, 0), 4), monomial((0, 0, 2), 1))

    checks = []
    for k in range(6, 10):
        p_k = multiply(power(rho, k - 6), base_p)
        assert homogeneous_degree(p_k) == 2 * k
        assert all(sum(exponent) % 2 == 0 for exponent in p_k)

        p_power = {ZERO: 1}
        for order in range(1, 4):
            p_power = multiply(p_power, p_k)
            pure = apolar_scalar(power(delta, k * order), p_power)
            assert pure == 0
            assert gaussian_expectation(p_power) == 0

            mixed = []
            for ell in range(1, min(3, order) + 1):
                multiplier = monomial((2 * ell, 0, 0))
                mixed_input = multiply(multiplier, p_power)
                gaussian_value = gaussian_expectation(mixed_input)
                expected_gaussian = expected_gaussian_mixed(k, order, ell)
                assert gaussian_value == expected_gaussian
                assert gaussian_value != 0

                detector = apolar_scalar(
                    power(delta, k * order + ell),
                    mixed_input,
                )
                expected = expected_detector(k, order, ell)
                assert detector == expected
                assert detector != 0

                sphere_values = {}
                for dimension in range(3, 9):
                    value = Fraction(
                        gaussian_value,
                        radial_even_moment(
                            dimension,
                            k * order + ell,
                        ),
                    )
                    assert value != 0
                    sphere_values[str(dimension)] = str(value)

                mixed.append(
                    {
                        "ell": ell,
                        "gaussian_moment": gaussian_value,
                        "detector": detector,
                        "sphere_moments_dimensions_3_through_8": sphere_values,
                    }
                )

            checks.append(
                {
                    "k": k,
                    "m": order,
                    "degree_P_k": 2 * k,
                    "term_count_P_k": len(p_k),
                    "pure": pure,
                    "mixed_ladder": mixed,
                }
            )

    # The binomial functions form a basis: a nonzero positive-phase
    # multiplier A(u)=sum_{ell>=1} a_ell u^ell cannot have identically zero
    # mixed sequence. Replay one nontrivial cancellation-prone example.
    coefficients = {1: 2, 2: -3, 4: 1}
    values = []
    for order in range(4, 12):
        value = sum(
            coefficient * comb(order - 1, ell - 1)
            for ell, coefficient in coefficients.items()
        )
        values.append(value)
    assert any(values)

    artifact = {
        "format": "gvc3-homogeneous-spillovers-v1",
        "field": "characteristic zero",
        "base_polynomial": "P=(rho+x^2)*(y*rho^2-2*x*t^2*rho-x^3*t^2)^2",
        "rho": "t^2+x*y",
        "family": {
            "P_k": "rho^(k-6)*P",
            "operator": "Lambda_k=(4*d_x*d_y+d_t^2)^k",
            "range": "every integer k>=6",
            "multiplier_ladder": "Q_ell=x^(2*ell), m>=ell",
            "sphere_coefficient": (
                "binom(m-1,ell-1)*2^(2m)*(2m)!/(4m+1)!!"
            ),
            "detector": (
                "Delta^(k*m+ell)(x^(2*ell)*P_k^m)="
                "2^((k+2)m+ell)*(k*m+ell)!*(2m)!*"
                "(2*k*m+2*ell+1)!!/(4m+1)!!*binom(m-1,ell-1)"
            ),
        },
        "bounded_replay": checks,
        "positive_phase_polynomial_regression": {
            "coefficients_by_ell": coefficients,
            "orders": list(range(4, 12)),
            "binomial_sequence_values": values,
        },
        "consequences_checked": [
            "pure GVC identities for k=6..9 and m=1..3",
            "positive-phase mixed ladder through ell=3",
            "rank-one SIC contractions (same apolar identities)",
            "homogeneous Gaussian moments",
            "nonzero sphere moments in ambient dimensions 3..8",
        ],
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")

    print("PASS radial padding Delta^k for k=6..9")
    print("PASS positive-phase multiplier ladder through ell=3")
    print("PASS homogeneous Gaussian and sphere transfer dimensions 3..8")
    print("PASS rank-one SIC contractions")
    print(f"PASS wrote {OUTPUT}")


if __name__ == "__main__":
    main()
