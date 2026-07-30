#!/usr/bin/env python3
"""Verify the compact projective-resultant proof for (72,108) Case 2.

The pinned characteristic-zero residuals are reduced at the good place
``101, U=55``.  Cramer elimination on the two projective (r,s)-charts leaves
degree-eight univariate pairs.  Explicit extended-gcd identities over F_101
show that all relevant resultants are nonzero.
"""

from __future__ import annotations

import hashlib
import pickle
from pathlib import Path

from flint import nmod_poly


PRIME = 101
FIELD_ROOT = 55
CERTIFICATE_SHA256 = (
    "cfbc3c39d7a28013671144f43ef76f0498542eaf6d562dd624bba3311194e4aa"
)

ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = (
    ROOT
    / "external"
    / "zenodo-21479814"
    / "bilLkarkariy-jc2-72-108-exact-certificates-d9ea4fd"
    / "release_bundle"
    / "exact_replay"
    / "case2_exact_certificate.pkl"
)


def reduce_rational(value: tuple[int, int]) -> int:
    numerator, denominator = map(int, value)
    if denominator % PRIME == 0:
        raise RuntimeError("the chosen place meets a certificate denominator")
    return numerator % PRIME * pow(denominator, -1, PRIME) % PRIME


def reduce_field(serialized: dict[int, tuple[int, int]]) -> int:
    return sum(
        reduce_rational(value) * pow(FIELD_ROOT, int(degree), PRIME)
        for degree, value in serialized.items()
    ) % PRIME


def reduce_polynomial(serialized):
    result = {}
    for monomial, coefficient in serialized.items():
        value = reduce_field(coefficient)
        if value:
            result[tuple(map(int, monomial))] = value
    return result


def polynomial(coefficients) -> nmod_poly:
    return nmod_poly(list(coefficients), PRIME)


def dense(coefficient_map: dict[int, int]) -> nmod_poly:
    return polynomial(
        coefficient_map.get(degree, 0)
        for degree in range(max(coefficient_map, default=-1) + 1)
    )


def linear_parts(generator, chart: str):
    """Write a cubic residual as A(x)h+B(x)q+C(x)."""
    variable_index = 0 if chart == "s" else 1
    parts = ({}, {}, {})
    for monomial, coefficient in generator.items():
        r_degree, s_degree, h_degree = monomial[:3]
        radial_degree = r_degree + s_degree
        x_degree = monomial[variable_index]
        if radial_degree == 1 and h_degree == 1:
            part = parts[0]
        elif radial_degree == 3 and h_degree == 0:
            part = parts[1]
        elif radial_degree == 1 and h_degree == 0:
            part = parts[2]
        else:
            raise RuntimeError(f"unexpected cubic monomial {monomial}")
        part[x_degree] = (part.get(x_degree, 0) + coefficient) % PRIME
    return tuple(dense(part) for part in parts)


def cramer_numerator(generator, chart, determinant, h_numerator, q_numerator):
    """Substitute h=H/D and q=Q/D, then multiply by D^2."""
    variable_index = 0 if chart == "s" else 1
    result = polynomial([])
    for monomial, coefficient in generator.items():
        r_degree, s_degree, h_degree = monomial[:3]
        radial_degree = r_degree + s_degree
        if radial_degree not in (0, 2, 4):
            raise RuntimeError(f"unexpected quartic monomial {monomial}")
        q_degree = radial_degree // 2
        denominator_degree = h_degree + q_degree
        if denominator_degree > 2:
            raise RuntimeError(f"unexpected Cramer degree in {monomial}")
        x_degree = monomial[variable_index]
        x_power = polynomial([0] * x_degree + [1])
        result += (
            coefficient
            * x_power
            * h_numerator**h_degree
            * q_numerator**q_degree
            * determinant ** (2 - denominator_degree)
        )
    return result


def verify_bezout(left: nmod_poly, right: nmod_poly, label: str) -> None:
    gcd, left_multiplier, right_multiplier = left.xgcd(right)
    if gcd != 1:
        raise RuntimeError(f"{label}: resultant vanishes modulo {PRIME}")
    if left_multiplier * left + right_multiplier * right != 1:
        raise RuntimeError(f"{label}: extended-gcd identity failed")
    print(
        f"{label}: degrees ({left.degree()},{right.degree()}), "
        f"Bezout degrees ({left_multiplier.degree()},{right_multiplier.degree()})"
    )


def origin_polynomial(generator) -> nmod_poly:
    coefficients = {}
    for monomial, coefficient in generator.items():
        r_degree, s_degree, h_degree = monomial[:3]
        if r_degree == s_degree == 0:
            coefficients[h_degree] = coefficient
    return dense(coefficients)


def main() -> None:
    raw = CERTIFICATE.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    if digest != CERTIFICATE_SHA256:
        raise RuntimeError(f"unexpected certificate SHA-256: {digest}")
    payload = pickle.loads(raw)
    if payload.get("format") != "jc2-case2-nullstellensatz-v1":
        raise RuntimeError("unsupported certificate format")

    minimal_polynomial = polynomial(
        reduce_rational(value)
        for value in payload["field_minpoly_coefficients"]
    )
    if minimal_polynomial(FIELD_ROOT) != 0:
        raise RuntimeError("U=55 is not a root of H modulo 101")

    generators = [reduce_polynomial(item) for item in payload["generators"]]
    if tuple(map(len, generators)) != (8, 8, 14, 12):
        raise RuntimeError("unexpected reduced generator supports")

    for chart in ("s", "r"):
        a0, b0, c0 = linear_parts(generators[0], chart)
        a1, b1, c1 = linear_parts(generators[1], chart)
        determinant = a0 * b1 - a1 * b0
        h_numerator = b0 * c1 - b1 * c0
        q_numerator = a1 * c0 - a0 * c1
        if (
            determinant.degree(),
            h_numerator.degree(),
            q_numerator.degree(),
        ) != (4, 4, 2):
            raise RuntimeError(f"{chart}-chart Cramer degrees dropped")

        # If D=0, consistency of the first two affine-linear equations
        # requires Q=0.  Their nonzero resultant removes this branch.
        verify_bezout(
            determinant, q_numerator, f"{chart}-chart singular Cramer branch"
        )

        numerator_7 = cramer_numerator(
            generators[2],
            chart,
            determinant,
            h_numerator,
            q_numerator,
        )
        numerator_9 = cramer_numerator(
            generators[3],
            chart,
            determinant,
            h_numerator,
            q_numerator,
        )
        if (numerator_7.degree(), numerator_9.degree()) != (8, 8):
            raise RuntimeError(f"{chart}-chart eliminant degrees dropped")
        verify_bezout(
            numerator_7, numerator_9, f"{chart}-chart Cramer eliminants"
        )

    verify_bezout(
        origin_polynomial(generators[2]),
        origin_polynomial(generators[3]),
        "origin",
    )
    print("good place: p=101, U=55")
    print("CASE2_PROJECTIVE_RESULTANT_PASS")


if __name__ == "__main__":
    main()
