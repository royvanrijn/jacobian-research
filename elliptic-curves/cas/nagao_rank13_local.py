#!/usr/bin/env python3
"""Exact p-adic root balls for Nagao's base-changed rank-13 family.

For

``T=(23550-u^2)/(2u)``,

let ``G(u)=(2u)^20 disc(R_T)`` for the primitive Mestre quartic.  At an odd
prime where ``u`` is a unit, the denominator is a unit, so roots of ``G``
are precisely bad-reduction conditions for the associated short Jacobian.

This module discovers those conditions rather than relying on a handwritten
residue table.  A multiplicative classification is uniform on a returned
ball: ``c4`` is a unit and the split symbol depends only on ``u mod p``.
Additive balls are retained but explicitly left unresolved and are not used
by the default CRT search.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from fractions import Fraction
from math import comb
from typing import Any, Sequence

from ek_k3 import legendre_symbol
from multiple_root_lifting import (
    RootBall,
    affine_variable_coefficients,
    all_roots_mod_prime_power,
    fixed_divisor_valuation,
)
from nagao_1994 import (
    RANK13_BASE_CHANGE_CONSTANT,
    rank13_base_changed_discriminant_numerator,
)


Q = Fraction
Polynomial = tuple[int, ...]

DEFAULT_TARGET_EXPONENTS = {
    7: 5,
    11: 3,
    13: 3,
    17: 3,
    19: 3,
    23: 3,
    31: 3,
}
DEFAULT_CRT_PRIMES = (7, 11, 13, 19, 31)


def _trim(coefficients: list[int]) -> Polynomial:
    while len(coefficients) > 1 and coefficients[-1] == 0:
        coefficients.pop()
    return tuple(coefficients)


def polynomial_add(*polynomials: Sequence[int]) -> Polynomial:
    degree = max((len(polynomial) for polynomial in polynomials), default=1)
    return _trim(
        [
            sum(
                polynomial[index] if index < len(polynomial) else 0
                for polynomial in polynomials
            )
            for index in range(degree)
        ]
    )


def polynomial_scale(polynomial: Sequence[int], scalar: int) -> Polynomial:
    return _trim([scalar * coefficient for coefficient in polynomial])


def polynomial_multiply(*polynomials: Sequence[int]) -> Polynomial:
    answer = (1,)
    for polynomial in polynomials:
        product = [0] * (len(answer) + len(polynomial) - 1)
        for left_degree, left in enumerate(answer):
            for right_degree, right in enumerate(polynomial):
                product[left_degree + right_degree] += left * right
        answer = _trim(product)
    return answer


def polynomial_value_mod(polynomial: Sequence[int], value: int, modulus: int) -> int:
    answer = 0
    for coefficient in reversed(polynomial):
        answer = (answer * value + coefficient) % modulus
    return answer


def _base_change_clear(polynomial: Sequence[int], denominator_degree: int) -> Polynomial:
    """Return ``(2u)^D f((23550-u^2)/(2u))`` for ``D=denominator_degree``."""

    answer = [0] * (2 * denominator_degree + 1)
    constant = RANK13_BASE_CHANGE_CONSTANT
    for source_degree, coefficient in enumerate(polynomial):
        if coefficient == 0:
            continue
        if source_degree > denominator_degree:
            raise ValueError("the clearing degree is too small")
        scalar = coefficient * 2 ** (denominator_degree - source_degree)
        base_degree = denominator_degree - source_degree
        for chosen in range(source_degree + 1):
            target_degree = base_degree + 2 * chosen
            answer[target_degree] += (
                scalar
                * comb(source_degree, chosen)
                * constant ** (source_degree - chosen)
                * (-1) ** chosen
            )
    return _trim(answer)


def rank13_base_changed_invariant_numerators() -> dict[str, Polynomial]:
    """Return cleared numerators of ``I,J,a4,a6,c4`` as polynomials in u."""

    # Ascending coefficients in T, transcribed from Nagao's printed quartic.
    quartic_e = (4_156_297_690_000, 0, 891_699_592, 0, -159_200, 0, 9)
    quartic_d = (-284_435_346_600, 0, -29_575_350, 0, 2_700)
    quartic_c = (6_706_476_489, 0, 396_150, 0, -18)
    quartic_b = (-63_901_710, 0, -2_700)
    quartic_a = (211_950, 0, 9)
    invariant_i = polynomial_add(
        polynomial_scale(polynomial_multiply(quartic_a, quartic_e), 12),
        polynomial_scale(polynomial_multiply(quartic_b, quartic_d), -3),
        polynomial_multiply(quartic_c, quartic_c),
    )
    invariant_j = polynomial_add(
        polynomial_scale(
            polynomial_multiply(quartic_a, quartic_c, quartic_e), 72
        ),
        polynomial_scale(
            polynomial_multiply(quartic_b, quartic_c, quartic_d), 9
        ),
        polynomial_scale(
            polynomial_multiply(quartic_a, quartic_d, quartic_d), -27
        ),
        polynomial_scale(
            polynomial_multiply(quartic_b, quartic_b, quartic_e), -27
        ),
        polynomial_scale(
            polynomial_multiply(quartic_c, quartic_c, quartic_c), -2
        ),
    )
    cleared_i = _base_change_clear(invariant_i, 8)
    cleared_j = _base_change_clear(invariant_j, 12)
    return {
        "I": cleared_i,
        "J": cleared_j,
        "a4": polynomial_scale(cleared_i, -27),
        "a6": polynomial_scale(cleared_j, -27),
        "c4": polynomial_scale(cleared_i, 1296),
    }


BASE_CHANGED_INVARIANTS = rank13_base_changed_invariant_numerators()
BASE_CHANGED_DISCRIMINANT = tuple(
    int(coefficient) for coefficient in rank13_base_changed_discriminant_numerator()
)


@dataclass(frozen=True)
class NagaoLocalBall:
    prime: int
    requested_discriminant_exponent: int
    exponent: int
    residue: int
    forced_discriminant_valuation: int
    reduction: str
    conductor_exponent: int | None
    split_multiplicative: bool | None
    proof: str

    @property
    def modulus(self) -> int:
        return self.prime**self.exponent

    def as_json(self) -> dict[str, Any]:
        answer = asdict(self)
        answer["modulus"] = self.modulus
        return answer


@dataclass(frozen=True)
class NagaoLocalDiscovery:
    prime: int
    requested_discriminant_exponent: int
    level_counts: tuple[int, ...]
    fixed_divisor_valuation: int
    balls: tuple[NagaoLocalBall, ...]

    def as_json(self) -> dict[str, Any]:
        return {
            "prime": self.prime,
            "requested_discriminant_exponent": self.requested_discriminant_exponent,
            "level_counts": list(self.level_counts),
            "fixed_divisor_valuation": self.fixed_divisor_valuation,
            "unit_balls": [ball.as_json() for ball in self.balls],
        }


def _short_coefficients_mod(residue: int, prime: int) -> tuple[int, int]:
    residue %= prime
    if residue == 0:
        raise ValueError("the local base-change parameter must be a p-adic unit")
    denominator_a4 = pow(2 * residue, 8, prime)
    denominator_a6 = pow(2 * residue, 12, prime)
    coefficient_a = (
        polynomial_value_mod(BASE_CHANGED_INVARIANTS["a4"], residue, prime)
        * pow(denominator_a4, -1, prime)
        % prime
    )
    coefficient_b = (
        polynomial_value_mod(BASE_CHANGED_INVARIANTS["a6"], residue, prime)
        * pow(denominator_a6, -1, prime)
        % prime
    )
    return coefficient_a, coefficient_b


def classify_root_ball(
    ball: RootBall, requested_discriminant_exponent: int
) -> NagaoLocalBall:
    """Classify one unit root ball, proving multiplicative cases uniformly."""

    if ball.residue % ball.prime == 0:
        raise ValueError("classification is restricted to p-adic unit balls")
    transformed_discriminant = affine_variable_coefficients(
        BASE_CHANGED_DISCRIMINANT, ball.residue, ball.modulus
    )
    forced_valuation = fixed_divisor_valuation(
        transformed_discriminant, ball.prime
    )
    if forced_valuation < requested_discriminant_exponent:
        raise AssertionError("the root ball does not force the requested valuation")

    coefficient_a, coefficient_b = _short_coefficients_mod(
        ball.residue, ball.prime
    )
    c4 = (-48 * coefficient_a) % ball.prime
    if c4:
        double_root = (
            -3
            * coefficient_b
            * pow(2 * coefficient_a, -1, ball.prime)
            % ball.prime
        )
        tangent_symbol = legendre_symbol(3 * double_root, ball.prime)
        if tangent_symbol == 0:
            raise AssertionError("a multiplicative fiber became cuspidal")
        split = tangent_symbol == 1
        return NagaoLocalBall(
            prime=ball.prime,
            requested_discriminant_exponent=requested_discriminant_exponent,
            exponent=ball.exponent,
            residue=ball.residue,
            forced_discriminant_valuation=forced_valuation,
            reduction="split multiplicative" if split else "nonsplit multiplicative",
            conductor_exponent=1,
            split_multiplicative=split,
            proof=(
                "u is a unit, the short-model c4 is a unit modulo p throughout "
                "the ball, and the tangent symbol is fixed by u modulo p"
            ),
        )
    return NagaoLocalBall(
        prime=ball.prime,
        requested_discriminant_exponent=requested_discriminant_exponent,
        exponent=ball.exponent,
        residue=ball.residue,
        forced_discriminant_valuation=forced_valuation,
        reduction="additive or unresolved after minimalization",
        conductor_exponent=None,
        split_multiplicative=None,
        proof=(
            "the presented c4 is divisible by p; this coarse ball is retained "
            "as an exact discriminant condition but is not used by default CRT"
        ),
    )


def discover_local_conditions(
    prime: int,
    requested_discriminant_exponent: int,
    *,
    max_roots: int = 100_000,
) -> NagaoLocalDiscovery:
    """Discover every unit ball forcing ``p^requested_exponent | G(u)``."""

    result = all_roots_mod_prime_power(
        BASE_CHANGED_DISCRIMINANT,
        prime,
        requested_discriminant_exponent,
        max_roots=max_roots,
    )
    balls = tuple(
        classify_root_ball(ball, requested_discriminant_exponent)
        for ball in result.maximal_balls()
        if ball.residue % prime
    )
    return NagaoLocalDiscovery(
        prime=prime,
        requested_discriminant_exponent=requested_discriminant_exponent,
        level_counts=result.level_counts,
        fixed_divisor_valuation=fixed_divisor_valuation(
            BASE_CHANGED_DISCRIMINANT, prime
        ),
        balls=balls,
    )


def default_local_discoveries() -> tuple[NagaoLocalDiscovery, ...]:
    return tuple(
        discover_local_conditions(prime, exponent)
        for prime, exponent in DEFAULT_TARGET_EXPONENTS.items()
    )


def default_crt_balls() -> dict[int, tuple[NagaoLocalBall, ...]]:
    """Return the exponent-one split balls used by the bounded CRT search."""

    discoveries = {
        discovery.prime: discovery for discovery in default_local_discoveries()
    }
    answer = {
        prime: tuple(
            ball
            for ball in discoveries[prime].balls
            if ball.exponent == 1 and ball.reduction == "split multiplicative"
        )
        for prime in DEFAULT_CRT_PRIMES
    }
    if any(not balls for balls in answer.values()):
        raise AssertionError("a default CRT prime has no certified split ball")
    return answer


def rational_discriminant_valuation(
    numerator: int, denominator: int, prime: int
) -> int:
    """Return ``v_p(G(numerator/denominator))`` for a p-unit denominator."""

    if denominator % prime == 0:
        raise ValueError("the rational denominator must be a p-adic unit")
    degree = len(BASE_CHANGED_DISCRIMINANT) - 1
    homogenized = sum(
        coefficient * numerator**index * denominator ** (degree - index)
        for index, coefficient in enumerate(BASE_CHANGED_DISCRIMINANT)
    )
    if homogenized == 0:
        raise ValueError("the specialization is singular")
    value = abs(homogenized)
    valuation = 0
    while value % prime == 0:
        value //= prime
        valuation += 1
    return valuation

