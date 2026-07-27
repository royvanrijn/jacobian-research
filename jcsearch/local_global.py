"""Constructive coefficient weak approximation for local polynomial data.

The module can derive a universal discriminant-based stability radius for a
monic p-integral local model, or consume a sharper caller-supplied radius.
It then performs the global prime-power CRT step exactly.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math
from typing import Mapping, Sequence

import sympy as sp


@dataclass(frozen=True)
class CoefficientCRTLift:
    """One rational vector produced by the denominator-compatible CRT grid."""

    coefficients: tuple[Fraction, ...]
    common_denominator: int
    base_denominator: int
    crt_modulus: int
    multiplier: int


@dataclass(frozen=True)
class MonicPolynomialCRTLift:
    """A synthesized monic polynomial and its coefficient-grid certificate."""

    polynomial: sp.Poly
    certificate: CoefficientCRTLift
    automatic_local_certificates: tuple["LocalAlgebraStabilityCertificate", ...] = ()


@dataclass(frozen=True)
class LocalAlgebraStabilityCertificate:
    """A universal discriminant-radius certificate for one local model."""

    prime: int
    polynomial: sp.Poly
    discriminant: sp.Rational
    discriminant_valuation: int
    coefficient_precision: int


def rational_valuation(value: Fraction, prime: int) -> int | None:
    """Return ``v_p(value)``, with ``None`` representing positive infinity."""
    prime = int(prime)
    if not sp.isprime(prime):
        raise ValueError("prime must be a rational prime")
    value = Fraction(value)
    if value == 0:
        return None

    def integer_valuation(integer: int) -> int:
        integer = abs(integer)
        exponent = 0
        while integer % prime == 0:
            integer //= prime
            exponent += 1
        return exponent

    return integer_valuation(value.numerator) - integer_valuation(value.denominator)


def automatic_local_stability_certificate(
    polynomial: sp.Poly | sp.Expr,
    variable: sp.Symbol,
    prime: int,
) -> LocalAlgebraStabilityCertificate:
    """Certify a coefficient radius preserving a finite étale Q_p-algebra.

    If ``f`` is monic, p-integral, and squarefree, put
    ``D=v_p(Disc(f))``.  Every monic same-degree ``g`` satisfying
    ``g == f (mod p**(2*D+1))`` defines an isomorphic finite étale
    ``Q_p``-algebra.  The bound is uniform and deliberately conservative.
    """
    prime = int(prime)
    if not sp.isprime(prime):
        raise ValueError("prime must be a rational prime")
    local_model = sp.Poly(polynomial, variable, domain=sp.QQ)
    if local_model.degree() < 1 or local_model.LC() != 1:
        raise ValueError("the local model must be monic of positive degree")
    for coefficient in local_model.all_coeffs():
        value = rational_valuation(Fraction(coefficient), prime)
        if value is not None and value < 0:
            raise ValueError("the local model must be p-integral")

    discriminant = sp.Rational(sp.discriminant(local_model.as_expr(), variable))
    if not discriminant:
        raise ValueError("the local model must be squarefree")
    discriminant_valuation = rational_valuation(Fraction(discriminant), prime)
    assert discriminant_valuation is not None and discriminant_valuation >= 0
    precision = 2 * discriminant_valuation + 1
    return LocalAlgebraStabilityCertificate(
        prime=prime,
        polynomial=local_model,
        discriminant=discriminant,
        discriminant_valuation=discriminant_valuation,
        coefficient_precision=precision,
    )


def certify_polynomial_in_local_stability_ball(
    candidate: sp.Poly | sp.Expr,
    certificate: LocalAlgebraStabilityCertificate,
) -> None:
    """Verify the exact congruences required by a stability certificate."""
    variable = certificate.polynomial.gens[0]
    polynomial = sp.Poly(candidate, variable, domain=sp.QQ)
    if (
        polynomial.degree() != certificate.polynomial.degree()
        or polynomial.LC() != 1
    ):
        raise ValueError("candidate must be monic of the certified degree")
    for index in range(polynomial.degree()):
        difference = (
            polynomial.coeff_monomial(variable**index)
            - certificate.polynomial.coeff_monomial(variable**index)
        )
        value = rational_valuation(Fraction(difference), certificate.prime)
        if value is not None and value < certificate.coefficient_precision:
            raise ValueError("candidate lies outside the certified coefficient ball")


def _prime_power_crt(residues: Sequence[tuple[int, int]]) -> tuple[int, int]:
    """Combine residues modulo pairwise-coprime prime powers."""
    modulus = math.prod(local_modulus for local_modulus, _ in residues)
    value = sum(
        (residue % local_modulus)
        * (modulus // local_modulus)
        * pow(modulus // local_modulus, -1, local_modulus)
        for local_modulus, residue in residues
    )
    return value % modulus, modulus


def constructive_coefficient_crt(
    local_balls: Mapping[int, tuple[int, Sequence[Fraction | int]]],
    real_intervals: Sequence[tuple[Fraction | int, Fraction | int]],
) -> CoefficientCRTLift:
    """Meet certified p-adic coefficient balls and one rational real box.

    ``local_balls[p]`` is ``(m, center)``, representing the coordinatewise
    conditions ``v_p(c_i-center_i) >= m``.  Centers must be rational and all
    vectors must have the same dimension as ``real_intervals``.  The
    precision ``m`` is required to be positive; any open p-adic ball may be
    shrunk to this form.

    The result uses a common denominator ``D=D_0(1+kM)``.  Here ``D_0``
    clears the local centers and ``M`` includes the extra denominator
    valuations needed to preserve every ball.  The translated lattices have
    mesh ``M/D``, which eventually fits inside every real interval.
    """
    intervals = tuple(
        (Fraction(left), Fraction(right)) for left, right in real_intervals
    )
    if not intervals or any(left >= right for left, right in intervals):
        raise ValueError("real_intervals must be a nonempty rational open box")
    dimension = len(intervals)

    normalized: dict[int, tuple[int, tuple[Fraction, ...]]] = {}
    for raw_prime, (raw_precision, raw_center) in local_balls.items():
        prime = int(raw_prime)
        precision = int(raw_precision)
        center = tuple(Fraction(value) for value in raw_center)
        if not sp.isprime(prime):
            raise ValueError("local_balls keys must be rational primes")
        if prime in normalized:
            raise ValueError("local_balls keys must be distinct")
        if precision <= 0:
            raise ValueError("local coefficient precisions must be positive")
        if len(center) != dimension:
            raise ValueError("every local center must match the real-box dimension")
        normalized[prime] = (precision, center)

    if not normalized:
        midpoint = tuple((left + right) / 2 for left, right in intervals)
        denominator = math.lcm(*(value.denominator for value in midpoint))
        return CoefficientCRTLift(midpoint, denominator, denominator, 1, 0)

    base_denominator = math.lcm(
        *(
            value.denominator
            for _, center in normalized.values()
            for value in center
        )
    )

    local_moduli: dict[int, int] = {}
    for prime, (precision, _) in normalized.items():
        denominator_valuation = rational_valuation(Fraction(base_denominator), prime)
        assert denominator_valuation is not None
        exponent = precision + denominator_valuation
        local_moduli[prime] = prime**exponent

    coefficient_residues: list[int] = []
    crt_modulus = math.prod(local_moduli.values())
    for coordinate in range(dimension):
        congruences = []
        for prime, (_, center) in normalized.items():
            scaled_center = base_denominator * center[coordinate]
            assert scaled_center.denominator == 1
            congruences.append((local_moduli[prime], scaled_center.numerator))
        residue, modulus = _prime_power_crt(congruences)
        assert modulus == crt_modulus
        coefficient_residues.append(residue)

    minimum_width = min(right - left for left, right in intervals)
    guaranteed_multiplier = 0
    while Fraction(
        crt_modulus,
        base_denominator * (1 + guaranteed_multiplier * crt_modulus),
    ) >= minimum_width:
        guaranteed_multiplier += 1

    for multiplier in range(guaranteed_multiplier + 1):
        denominator = base_denominator * (1 + multiplier * crt_modulus)
        result = []
        for residue, (lower, upper) in zip(coefficient_residues, intervals):
            shift = math.floor((lower * denominator - residue) / crt_modulus) + 1
            candidate = Fraction(residue + shift * crt_modulus, denominator)
            if not lower < candidate < upper:
                break
            result.append(candidate)
        else:
            coefficients = tuple(result)
            for prime, (precision, center) in normalized.items():
                for coefficient, local_center in zip(coefficients, center):
                    value = rational_valuation(coefficient - local_center, prime)
                    assert value is None or value >= precision
            return CoefficientCRTLift(
                coefficients=coefficients,
                common_denominator=denominator,
                base_denominator=base_denominator,
                crt_modulus=crt_modulus,
                multiplier=multiplier,
            )

    raise AssertionError("guaranteed coefficient CRT mesh bound failed")


def synthesize_monic_polynomial(
    local_models: Mapping[int, tuple[int | None, sp.Poly | sp.Expr]],
    real_coefficient_intervals: Sequence[
        tuple[Fraction | int, Fraction | int]
    ],
    variable: sp.Symbol,
) -> MonicPolynomialCRTLift:
    """Synthesize a monic polynomial from local models and a real box.

    Coefficient vectors are ordered from the constant term through degree
    ``N-1``.  Every local model must be monic of the common degree ``N``.
    Set a model's precision to ``None`` to use the automatic
    ``2*v_p(Disc)+1`` certificate.
    """
    degree = len(real_coefficient_intervals)
    coefficient_balls = {}
    normalized_models = {}
    normalized_precisions = {}
    automatic_certificates = []
    for prime, (precision, model) in local_models.items():
        polynomial = sp.Poly(model, variable, domain=sp.QQ)
        if polynomial.degree() != degree or polynomial.LC() != 1:
            raise ValueError("every local model must be monic of the real-box degree")
        if precision is None:
            automatic = automatic_local_stability_certificate(
                polynomial, variable, int(prime)
            )
            precision = automatic.coefficient_precision
            automatic_certificates.append(automatic)
        center = tuple(
            Fraction(polynomial.coeff_monomial(variable**index))
            for index in range(degree)
        )
        coefficient_balls[int(prime)] = (int(precision), center)
        normalized_models[int(prime)] = polynomial
        normalized_precisions[int(prime)] = int(precision)

    certificate = constructive_coefficient_crt(
        coefficient_balls, real_coefficient_intervals
    )
    expression = variable**degree + sum(
        sp.Rational(value.numerator, value.denominator) * variable**index
        for index, value in enumerate(certificate.coefficients)
    )
    polynomial = sp.Poly(expression, variable, domain=sp.QQ)

    for prime in local_models:
        normalized_prime = int(prime)
        precision = normalized_precisions[normalized_prime]
        model = normalized_models[normalized_prime]
        for index in range(degree):
            difference = (
                polynomial.coeff_monomial(variable**index)
                - model.coeff_monomial(variable**index)
            )
            value = rational_valuation(Fraction(difference), normalized_prime)
            assert value is None or value >= precision

    for automatic in automatic_certificates:
        certify_polynomial_in_local_stability_ball(polynomial, automatic)

    return MonicPolynomialCRTLift(
        polynomial=polynomial,
        certificate=certificate,
        automatic_local_certificates=tuple(automatic_certificates),
    )
