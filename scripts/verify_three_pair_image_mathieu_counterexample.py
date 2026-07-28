#!/usr/bin/env python3
"""Exact replay of the four-term SIC(3) and four-real GMC counterexamples.

With contraction pairs (tau,t), (w,z), (v,y), the witness is

    f = tau (t-y) (w*z + v*t),
    g = y.

The written proof gives, for every m >= 1,

    E(f^m) = 0,
    [t]E(g*f^m) = (-1)^(m-1) (m+1)! m!.

This dependency-free checker expands the four-term polynomial and verifies
the identities by exact sparse contraction through the declared cutoff.  It
also reads the two-pair seed as a polynomial in two independent circular
complex Gaussians and verifies its GMC moments by a separate Wick
calculation.
"""

from __future__ import annotations

import json
from math import comb, factorial
from pathlib import Path


Index = tuple[int, int, int]
Monomial = tuple[Index, Index]
Polynomial = dict[Monomial, int]
GaussianMonomial = tuple[int, int, int, int]  # W1, Z1, W2, Z2
GaussianPolynomial = dict[GaussianMonomial, int]
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


def gaussian_multiply(
    left: GaussianPolynomial, right: GaussianPolynomial
) -> GaussianPolynomial:
    result: GaussianPolynomial = {}
    for left_monomial, left_coefficient in left.items():
        for right_monomial, right_coefficient in right.items():
            monomial = tuple(
                left_monomial[index] + right_monomial[index]
                for index in range(4)
            )
            result[monomial] = (
                result.get(monomial, 0)
                + left_coefficient * right_coefficient
            )
    return {
        monomial: coefficient
        for monomial, coefficient in result.items()
        if coefficient
    }


def gaussian_expectation(polynomial: GaussianPolynomial) -> int:
    """Apply the Wick rule for two independent circular complex Gaussians."""
    result = 0
    for (w1, z1, w2, z2), coefficient in polynomial.items():
        if w1 == z1 and w2 == z2:
            result += coefficient * factorial(w1) * factorial(w2)
    return result


def gaussian_witness() -> tuple[GaussianPolynomial, GaussianPolynomial]:
    # P=(1-Z2)(W1*Z1+W2), Q=Z2.
    p: GaussianPolynomial = {
        (1, 1, 0, 0): 1,
        (0, 0, 1, 0): 1,
        (1, 1, 0, 1): -1,
        (0, 0, 1, 1): -1,
    }
    q: GaussianPolynomial = {(0, 0, 0, 1): 1}
    return p, q


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

    gaussian_p, gaussian_q = gaussian_witness()
    assert len(gaussian_p) == 4
    gaussian_power: GaussianPolynomial = {(0, 0, 0, 0): 1}
    for order in range(1, CUTOFF + 1):
        gaussian_power = gaussian_multiply(gaussian_power, gaussian_p)
        assert gaussian_expectation(gaussian_power) == 0
        assert gaussian_expectation(
            gaussian_multiply(gaussian_q, gaussian_power)
        ) == (-1) ** (order - 1) * factorial(order)

    artifact = {
        "format": "three-pair-image-mathieu-counterexample-v3",
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
        "four_real_gaussian_corollary": {
            "real_gaussians": ["X1", "Y1", "X2", "Y2"],
            "circular_coordinates": {
                "Zj": "(Xj+i*Yj)/sqrt(2)",
                "Wj": "(Xj-i*Yj)/sqrt(2)",
            },
            "P": "(1-Z2)*(W1*Z1+W2)",
            "expanded_P": "W1*Z1+W2-W1*Z1*Z2-W2*Z2",
            "Q": "Z2",
            "expanded_P_term_count": len(gaussian_p),
            "total_degree": 3,
            "all_order_identities": {
                "E(P^m)": "0",
                "E(Q*P^m)": "(-1)^(m-1)*m!",
            },
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
    print(
        "PASS GMC(4): four-term cubic has pure moments zero and "
        "mixed moments (-1)^(m-1)m! through m=10"
    )
    print(f"PASS SIC(3): wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
