#!/usr/bin/env python3
"""Independently verify the explicit (72,108) Case-2 Bezout identity.

The checker deliberately imports neither the equation generator nor
``exact_core`` and does not invoke Singular.  Its only mathematical input is
the pinned serialized certificate.  Arithmetic is performed directly in
Q[U]/(H) with python-flint.
"""

from __future__ import annotations

import hashlib
import pickle
from pathlib import Path

from flint import fmpq, fmpq_poly


CERTIFICATE_SHA256 = (
    "cfbc3c39d7a28013671144f43ef76f0498542eaf6d562dd624bba3311194e4aa"
)
EXPECTED_GENERATOR_DEGREES = (3, 3, 4, 4)
EXPECTED_GENERATOR_TERMS = (8, 8, 14, 14)
EXPECTED_MULTIPLIER_DEGREES = (5, 5, 4, 4)
EXPECTED_MULTIPLIER_TERMS = (23, 20, 14, 11)
EXPECTED_GENERATOR_INDICES = (0, 1, 7, 9)

ROOT = Path(__file__).resolve().parents[1]
REPLAY = (
    ROOT
    / "external"
    / "zenodo-21479814"
    / "bilLkarkariy-jc2-72-108-exact-certificates-d9ea4fd"
    / "release_bundle"
    / "exact_replay"
)
CERTIFICATE = REPLAY / "case2_exact_certificate.pkl"


def total_degree(polynomial: dict[tuple[int, ...], object]) -> int:
    return max((sum(monomial) for monomial in polynomial), default=-1)


def decode_rational(pair: tuple[int, int]) -> fmpq:
    return fmpq(int(pair[0]), int(pair[1]))


def decode_field(
    serialized: dict[int, tuple[int, int]], modulus: fmpq_poly
) -> fmpq_poly:
    coefficients = [fmpq(0) for _ in range(modulus.degree())]
    for degree, value in serialized.items():
        degree = int(degree)
        if not 0 <= degree < modulus.degree():
            raise RuntimeError("certificate coefficient is not field-reduced")
        coefficients[degree] = decode_rational(value)
    return fmpq_poly(coefficients)


def decode_polynomial(
    serialized: dict[tuple[int, ...], dict[int, tuple[int, int]]],
    modulus: fmpq_poly,
) -> dict[tuple[int, ...], fmpq_poly]:
    return {
        tuple(int(exponent) for exponent in monomial): decode_field(
            coefficient, modulus
        )
        for monomial, coefficient in serialized.items()
    }


def add(
    left: dict[tuple[int, ...], fmpq_poly],
    right: dict[tuple[int, ...], fmpq_poly],
) -> dict[tuple[int, ...], fmpq_poly]:
    result = dict(left)
    for monomial, coefficient in right.items():
        value = result.get(monomial, fmpq_poly()) + coefficient
        if value:
            result[monomial] = value
        else:
            result.pop(monomial, None)
    return result


def multiply(
    left: dict[tuple[int, ...], fmpq_poly],
    right: dict[tuple[int, ...], fmpq_poly],
    modulus: fmpq_poly,
) -> dict[tuple[int, ...], fmpq_poly]:
    result: dict[tuple[int, ...], fmpq_poly] = {}
    for left_monomial, left_coefficient in left.items():
        for right_monomial, right_coefficient in right.items():
            if len(left_monomial) != len(right_monomial):
                raise RuntimeError("inconsistent monomial arity")
            monomial = tuple(
                a + b for a, b in zip(left_monomial, right_monomial)
            )
            product = (left_coefficient * right_coefficient) % modulus
            result[monomial] = (
                result.get(monomial, fmpq_poly()) + product
            ) % modulus
    return {monomial: value for monomial, value in result.items() if value}


def main() -> None:
    raw = CERTIFICATE.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    if digest != CERTIFICATE_SHA256:
        raise RuntimeError(f"unexpected certificate SHA-256: {digest}")

    payload = pickle.loads(raw)
    if payload.get("format") != "jc2-case2-nullstellensatz-v1":
        raise RuntimeError("unsupported certificate format")
    if tuple(payload.get("generator_indices", ())) != EXPECTED_GENERATOR_INDICES:
        raise RuntimeError("unexpected residual-generator indices")
    if payload.get("macaulay_degree") != 8:
        raise RuntimeError("unexpected Macaulay degree")
    if len(payload.get("row_monomials", ())) != 165:
        raise RuntimeError("unexpected Macaulay row count")
    if len(payload.get("columns", ())) != 182:
        raise RuntimeError("unexpected Macaulay column count")
    if len(payload.get("pivot_columns", ())) != 151:
        raise RuntimeError("unexpected exact Macaulay rank")

    modulus = fmpq_poly(
        [decode_rational(value) for value in payload["field_minpoly_coefficients"]]
    )
    if modulus.degree() != 35 or modulus.leading_coefficient() != 1:
        raise RuntimeError("the stored field polynomial is not monic of degree 35")
    unit, factors = modulus.factor()
    if (
        unit == 0
        or len(factors) != 1
        or factors[0][0].degree() != 35
        or factors[0][1] != 1
    ):
        raise RuntimeError("the stored degree-35 field polynomial is reducible")

    generators = [
        decode_polynomial(polynomial, modulus)
        for polynomial in payload["generators"]
    ]
    multipliers = [
        decode_polynomial(polynomial, modulus)
        for polynomial in payload["multipliers"]
    ]
    if tuple(total_degree(q) for q in generators) != EXPECTED_GENERATOR_DEGREES:
        raise RuntimeError("unexpected generator degrees")
    if tuple(len(q) for q in generators) != EXPECTED_GENERATOR_TERMS:
        raise RuntimeError("unexpected generator term counts")
    if tuple(total_degree(q) for q in multipliers) != EXPECTED_MULTIPLIER_DEGREES:
        raise RuntimeError("unexpected multiplier degrees")
    if tuple(len(q) for q in multipliers) != EXPECTED_MULTIPLIER_TERMS:
        raise RuntimeError("unexpected multiplier term counts")

    monomial_arity = len(next(iter(generators[0])))
    zero_monomial = (0,) * monomial_arity
    residual = {zero_monomial: fmpq_poly([-1])}
    for generator, multiplier in zip(generators, multipliers):
        residual = add(residual, multiply(generator, multiplier, modulus))
    if residual:
        raise RuntimeError(
            f"Bezout identity has {len(residual)} nonzero residual coefficients"
        )

    print("field degree=35 irreducible")
    print("Macaulay rows=165 columns=182 rank=151 degree=8")
    print("multiplier terms=23,20,14,11")
    print("1 = T0*R0 + T1*R1 + T7*R7 + T9*R9")
    print("CASE2_EXPLICIT_SYZYGY_PASS")


if __name__ == "__main__":
    main()
