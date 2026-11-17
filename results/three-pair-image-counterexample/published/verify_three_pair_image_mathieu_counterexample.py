#!/usr/bin/env python3
"""Exact replay of the four-term three-pair Image-Mathieu counterexample.

With contraction pairs (tau,t), (w,z), (v,y), the witness is

    f = tau (t-y) (w*z + v*t),
    g = y.

The written proof gives, for every m >= 1,

    E(f^m) = 0,
    [t]E(g*f^m) = (-1)^(m-1) (m+1)! m!.

This dependency-free checker expands the four-term polynomial and verifies
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
SCRIPT_PATH = Path(__file__).resolve()
IN_REPOSITORY = SCRIPT_PATH.parent.name == "scripts"
ROOT = SCRIPT_PATH.parents[1] if IN_REPOSITORY else SCRIPT_PATH.parent
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "three_pair_image_mathieu_counterexample.json"
    if IN_REPOSITORY
    else ROOT / "three_pair_image_mathieu_counterexample.json"
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
    f: Polynomial = {
        # tau*v*t^2 - tau*v*t*y
        ((1, 0, 1), (2, 0, 0)): 1,
        ((1, 0, 1), (1, 0, 1)): -1,
        # tau*w*t*z - tau*w*y*z
        ((1, 1, 0), (1, 1, 0)): 1,
        ((1, 1, 0), (0, 1, 1)): -1,
    }
    g: Polynomial = {((0, 0, 0), (0, 0, 1)): 1}
    return f, g


def main() -> None:
    f, g = witness()
    assert len(f) == 4
    assert all(sum(zeta) == sum(z) == 2 for zeta, z in f)
    assert all(sum(zeta) == 0 and sum(z) == 1 for zeta, z in g)

    f_powers = powers(f, CUTOFF)
    for order in range(1, CUTOFF + 1):
        assert contraction(f_powers[order]) == {}
        mixed = contraction(multiply(g, f_powers[order]))
        assert mixed.get((1, 0, 0)) == (
            (-1) ** (order - 1)
            * factorial(order + 1)
            * factorial(order)
        )

    # Independent finite replay of the two binomial identities used in the
    # all-order proof of the underlying two-pair circular moments.
    for order in range(1, 100):
        pure_sum = sum(
            (-1) ** index * comb(order, index)
            for index in range(order + 1)
        )
        mixed_sum = sum(
            (-1) ** index * comb(order, index)
            for index in range(order)
        )
        assert pure_sum == 0
        assert mixed_sum == (-1) ** (order - 1)

    artifact = {
        "format": "three-pair-image-mathieu-counterexample-v2",
        "field": "characteristic zero",
        "contraction_pairs": [["tau", "t"], ["w", "z"], ["v", "y"]],
        "f": "tau*(t-y)*(w*z+v*t)",
        "g": "y",
        "expanded_f_term_count": len(f),
        "bidegrees": {"f": [2, 2], "g": [0, 1]},
        "all_order_identities": {
            "E(f^m)": "0",
            "[t]E(g*f^m)": "(-1)^(m-1)*(m+1)!*m!",
        },
        "exact_sparse_replay_cutoff": CUTOFF,
        "written_proof": (
            "extended-geometry/"
            "THREE_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md"
            if IN_REPOSITORY
            else "main.tex"
        ),
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")

    print(
        "PASS SIC(3): E(f^m)=0 and "
        "[t]E(g*f^m)=(-1)^(m-1)(m+1)!m! through m=10"
    )
    print("PASS SIC(3): all-order binomial identities through m=99")
    print("PASS SIC(3): f has four terms and bidegree (2,2)")
    print(f"PASS SIC(3): wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
