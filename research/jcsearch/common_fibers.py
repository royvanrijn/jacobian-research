"""Local synthesis inside the fixed common-fiber Keller pencil."""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Mapping

import sympy as sp

from jcsearch.local_global import (
    CoefficientCRTLift,
    LocalAlgebraStabilityCertificate,
    automatic_local_stability_certificate,
    certify_polynomial_in_local_stability_ball,
    constructive_coefficient_crt,
)


@dataclass(frozen=True)
class FixedPencilLocalCertificate:
    """One local parameter and its full Q_p-algebra stability certificate."""

    degree: int
    prime: int
    local_parameter: Fraction
    local_polynomial: sp.Poly
    stability: LocalAlgebraStabilityCertificate


@dataclass(frozen=True)
class FixedCommonFiberSynthesis:
    """A rational parameter meeting local and real fixed-pencil conditions."""

    degree: int
    parameter: Fraction
    polynomial: sp.Poly
    local_certificates: tuple[FixedPencilLocalCertificate, ...]
    crt_certificate: CoefficientCRTLift
    real_interval: tuple[Fraction, Fraction]


def fixed_common_fiber_polynomial(
    degree: int,
    parameter: int | Fraction | sp.Rational,
    variable: sp.Symbol,
) -> sp.Poly:
    """Return ``T^N+T^3-2T^2+T+u`` over Q."""
    degree = int(degree)
    if degree < 4:
        raise ValueError("the fixed common-fiber pair starts in degree four")
    if isinstance(parameter, Fraction):
        value = sp.Rational(parameter.numerator, parameter.denominator)
    else:
        value = sp.Rational(parameter)
    return sp.Poly(
        variable**degree + variable**3 - 2 * variable**2 + variable + value,
        variable,
        domain=sp.QQ,
    )


def fixed_pencil_local_certificate(
    degree: int,
    prime: int,
    local_parameter: int | Fraction | sp.Rational,
    variable: sp.Symbol,
) -> FixedPencilLocalCertificate:
    """Derive the automatic parameter radius at one local witness."""
    parameter = Fraction(local_parameter)
    polynomial = fixed_common_fiber_polynomial(degree, parameter, variable)
    stability = automatic_local_stability_certificate(
        polynomial, variable, int(prime)
    )
    return FixedPencilLocalCertificate(
        degree=int(degree),
        prime=int(prime),
        local_parameter=parameter,
        local_polynomial=polynomial,
        stability=stability,
    )


def certify_fixed_pencil_parameter(
    parameter: int | Fraction | sp.Rational,
    certificate: FixedPencilLocalCertificate,
    variable: sp.Symbol,
) -> None:
    """Verify that a global parameter preserves one certified local algebra."""
    candidate = fixed_common_fiber_polynomial(
        certificate.degree, parameter, variable
    )
    certify_polynomial_in_local_stability_ball(candidate, certificate.stability)


def synthesize_fixed_common_fiber_parameter(
    degree: int,
    local_parameters: Mapping[int, int | Fraction | sp.Rational],
    real_interval: tuple[Fraction | int, Fraction | int],
    variable: sp.Symbol,
) -> FixedCommonFiberSynthesis:
    """Construct one rational parameter in all selected local balls."""
    interval = (Fraction(real_interval[0]), Fraction(real_interval[1]))
    certificates = tuple(
        fixed_pencil_local_certificate(
            degree, int(prime), local_parameter, variable
        )
        for prime, local_parameter in local_parameters.items()
    )
    local_balls = {
        certificate.prime: (
            certificate.stability.coefficient_precision,
            (certificate.local_parameter,),
        )
        for certificate in certificates
    }
    crt = constructive_coefficient_crt(local_balls, (interval,))
    parameter = crt.coefficients[0]
    polynomial = fixed_common_fiber_polynomial(degree, parameter, variable)
    for certificate in certificates:
        certify_fixed_pencil_parameter(parameter, certificate, variable)
    return FixedCommonFiberSynthesis(
        degree=int(degree),
        parameter=parameter,
        polynomial=polynomial,
        local_certificates=certificates,
        crt_certificate=crt,
        real_interval=interval,
    )
