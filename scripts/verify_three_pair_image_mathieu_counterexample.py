#!/usr/bin/env python3
"""Exact replay of the three-pair Image-Mathieu counterexample.

With contraction pairs (tau,t), (w,z), (v,y), the witness is

    f = tau^3 (t+z) (w*t^3 - v*y*(t+y)^2),
    g = z.

The written proof gives, for every m >= 1,

    E(f^m) = 0,
    [t]E(g*f^m) = (3m+1)! m!.

This dependency-free checker expands the eight-term polynomial and verifies
the identities by exact sparse contraction through the declared cutoff.
"""

from __future__ import annotations

import json
from math import comb, factorial
from pathlib import Path


Index = tuple[int, int, int]
Monomial = tuple[Index, Index]
Polynomial = dict[Monomial, int]
CUTOFF = 10
ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "three_pair_image_mathieu_counterexample.json"
)


def multiply(left: Polynomial, right: Polynomial) -> Polynomial:
    result: Polynomial = {}
    for (left_zeta, left_z), left_coefficient in left.items():
        for (right_zeta, right_z), right_coefficient in right.items():
            zeta_exponent = tuple(
                left_zeta[index] + right_zeta[index] for index in range(3)
            )
            z_exponent = tuple(
                left_z[index] + right_z[index] for index in range(3)
            )
            monomial = (zeta_exponent, z_exponent)
            result[monomial] = (
                result.get(monomial, 0)
                + left_coefficient * right_coefficient
            )
    return {
        monomial: coefficient
        for monomial, coefficient in result.items()
        if coefficient
    }


def powers(polynomial: Polynomial, cutoff: int) -> list[Polynomial]:
    result: list[Polynomial] = [{((0, 0, 0), (0, 0, 0)): 1}]
    for _ in range(cutoff):
        result.append(multiply(result[-1], polynomial))
    return result


def contraction(polynomial: Polynomial) -> dict[Index, int]:
    """Return E(polynomial) in the variable order (t,z,y)."""
    result: dict[Index, int] = {}
    for (zeta_exponent, z_exponent), coefficient in polynomial.items():
        if any(
            zeta_exponent[index] > z_exponent[index]
            for index in range(3)
        ):
            continue
        residual = tuple(
            z_exponent[index] - zeta_exponent[index]
            for index in range(3)
        )
        value = coefficient
        for derivative, degree in zip(zeta_exponent, z_exponent):
            value *= factorial(degree) // factorial(degree - derivative)
        result[residual] = result.get(residual, 0) + value
    return {
        exponent: coefficient
        for exponent, coefficient in result.items()
        if coefficient
    }


def witness() -> tuple[Polynomial, Polynomial]:
    # First summand: tau^3*w*t^3*(t+z).
    f: Polynomial = {
        ((3, 1, 0), (4, 0, 0)): 1,
        ((3, 1, 0), (3, 1, 0)): 1,
    }

    # Second summand: -tau^3*v*y*(t+z)*(t+y)^2.
    for outer in ((1, 0, 0), (0, 1, 0)):
        for inner, coefficient in (
            ((2, 0, 0), -1),
            ((1, 0, 1), -2),
            ((0, 0, 2), -1),
        ):
            z_exponent = (
                outer[0] + inner[0],
                outer[1] + inner[1],
                1 + outer[2] + inner[2],
            )
            monomial = ((3, 0, 1), z_exponent)
            f[monomial] = f.get(monomial, 0) + coefficient

    g: Polynomial = {((0, 0, 0), (0, 1, 0)): 1}
    return f, g


def main() -> None:
    f, g = witness()
    assert len(f) == 8
    assert all(sum(zeta) == sum(z) == 4 for zeta, z in f)
    assert all(sum(zeta) == 0 and sum(z) == 1 for zeta, z in g)

    f_powers = powers(f, CUTOFF)
    for order in range(1, CUTOFF + 1):
        assert contraction(f_powers[order]) == {}
        mixed = contraction(multiply(g, f_powers[order]))
        assert mixed.get((1, 0, 0)) == (
            factorial(3 * order + 1) * factorial(order)
        )

    # Independent finite replay of the two binomial identities used in the
    # all-order proof of the underlying two-pair circular moments.
    for order in range(1, 100):
        pure_sum = sum(
            (-1) ** index * comb(order, index)
            for index in range(order + 1)
        )
        mixed_sum = sum(
            (-1) ** index * comb(order, index + 1)
            for index in range(order)
        )
        assert pure_sum == 0
        assert mixed_sum == 1

    artifact = {
        "format": "three-pair-image-mathieu-counterexample-v1",
        "field": "characteristic zero",
        "contraction_pairs": [["tau", "t"], ["w", "z"], ["v", "y"]],
        "f": "tau^3*(t+z)*(w*t^3-v*y*(t+y)^2)",
        "g": "z",
        "expanded_f_term_count": len(f),
        "bidegrees": {"f": [4, 4], "g": [0, 1]},
        "all_order_identities": {
            "E(f^m)": "0",
            "[t]E(g*f^m)": "(3m+1)!*m!",
        },
        "exact_sparse_replay_cutoff": CUTOFF,
        "written_proof": (
            "extended-geometry/"
            "THREE_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md"
        ),
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")

    print(
        "PASS SIC(3): E(f^m)=0 and "
        "[t]E(g*f^m)=(3m+1)!m! through m=10"
    )
    print("PASS SIC(3): all-order binomial identities through m=99")
    print("PASS SIC(3): f has eight terms and bidegree (4,4)")
    print(f"PASS SIC(3): wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
