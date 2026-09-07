"""Exact small-prime local tables for integral Weierstrass models."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class LocalData:
    prime: int
    reduction: str
    discriminant_mod_prime: int
    trace: int | None
    point_count: int | None
    local_euler_coefficient: int | None


def legendre_symbol(value: int, prime: int) -> int:
    """Return the quadratic character modulo an odd caller-supplied prime."""

    if prime < 3 or prime % 2 == 0:
        raise ValueError("an odd prime is required")
    value %= prime
    if value == 0:
        return 0
    result = pow(value, (prime - 1) // 2, prime)
    return -1 if result == prime - 1 else result


def weierstrass_local_data(
    coefficients: tuple[int, int, int, int, int], prime: int
) -> LocalData:
    """Classify an integral Weierstrass equation modulo ``p >= 5``.

    For good reduction, the trace is computed by exhaustive character sum.
    When ``c4`` is a unit, the bad special fiber is clean multiplicative and
    this records its local Euler coefficient ``+1`` or ``-1``.  When both
    ``c4`` and the discriminant vanish modulo ``p``, residue data alone cannot
    distinguish additive reduction of a minimal elliptic model from a
    nonminimal model or a globally singular input; that branch is deliberately
    left unresolved.  Final candidates must be globally nonsingular and passed
    through minimalization.
    """

    if prime < 5:
        raise ValueError("this clean short-Weierstrass classifier requires p >= 5")
    a1, a2, a3, a4, a6 = (
        coefficient % prime for coefficient in coefficients
    )
    b2 = (a1 * a1 + 4 * a2) % prime
    b4 = (2 * a4 + a1 * a3) % prime
    b6 = (a3 * a3 + 4 * a6) % prime
    b8 = (
        a1 * a1 * a6
        + 4 * a2 * a6
        - a1 * a3 * a4
        + a2 * a3 * a3
        - a4 * a4
    ) % prime
    discriminant = (
        -b2 * b2 * b8
        - 8 * b4**3
        - 27 * b6**2
        + 9 * b2 * b4 * b6
    ) % prime
    if discriminant:
        character_sum = sum(
            legendre_symbol(
                (a1 * x + a3) ** 2
                + 4 * (x**3 + a2 * x**2 + a4 * x + a6),
                prime,
            )
            for x in range(prime)
        )
        trace = -character_sum
        return LocalData(
            prime=prime,
            reduction="good",
            discriminant_mod_prime=discriminant,
            trace=trace,
            point_count=prime + 1 - trace,
            local_euler_coefficient=trace,
        )
    c4 = (b2 * b2 - 24 * b4) % prime
    if c4:
        c6 = (-b2**3 + 36 * b2 * b4 - 216 * b6) % prime
        minus_c6 = -c6
        split = legendre_symbol(minus_c6, prime) == 1
        return LocalData(
            prime=prime,
            reduction="split_multiplicative" if split else "nonsplit_multiplicative",
            discriminant_mod_prime=0,
            trace=None,
            point_count=None,
            local_euler_coefficient=1 if split else -1,
        )
    return LocalData(
        prime=prime,
        reduction="unresolved_bad",
        discriminant_mod_prime=0,
        trace=None,
        point_count=None,
        local_euler_coefficient=None,
    )


def short_weierstrass_local_data(a4: int, a6: int, prime: int) -> LocalData:
    """Classify ``y^2=x^3+a4*x+a6`` at a small prime ``p >= 5``."""

    return weierstrass_local_data((0, 0, 0, a4, a6), prime)


def calibration_family_local_data(parameter: int, prime: int) -> LocalData:
    """Local data for ``y^2=x^3-t^2*x+t^2`` at ``t mod p``."""

    parameter %= prime
    parameter_square = parameter * parameter
    return short_weierstrass_local_data(
        -parameter_square,
        parameter_square,
        prime,
    )
