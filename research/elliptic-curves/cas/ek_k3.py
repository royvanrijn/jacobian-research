#!/usr/bin/env python3
"""Exact arithmetic for the Elkies--Klagsbrun rank-nine K3 family.

The family is [Elkies--Klagsbrun, equations (3)--(4)]

    y^2 = x^3 + 2 A(u,t) x^2 + B(u,t) x,
    B(u,t) = product(B_i(u,t), i=1,...,8).

All global calculations use :class:`fractions.Fraction`.  Modular routines
reject primes at which a rational coefficient is not integral.  Nothing in
this module estimates or certifies the rank of a specialization.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import gcd, isqrt
from typing import Iterable, Sequence


Q = Fraction


def parse_rational(value: str) -> Fraction:
    """Parse an integer or ``numerator/denominator`` string."""

    if "/" in value:
        numerator, denominator = value.split("/", 1)
        return Q(int(numerator), int(denominator))
    return Q(int(value), 1)


def rational_to_string(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def fraction_mod(value: Fraction, modulus: int) -> int:
    """Return a rational number modulo ``modulus`` when its denominator is a unit."""

    denominator = value.denominator % modulus
    if gcd(denominator, modulus) != 1:
        raise ValueError(
            f"denominator {value.denominator} is not invertible modulo {modulus}"
        )
    return (value.numerator % modulus) * pow(denominator, -1, modulus) % modulus


def valuation(value: Fraction, prime: int) -> int:
    """Return the normalized p-adic valuation of a nonzero rational number."""

    if value == 0:
        raise ValueError("the valuation of zero is infinite")

    def integer_valuation(integer: int) -> int:
        integer = abs(integer)
        answer = 0
        while integer % prime == 0:
            integer //= prime
            answer += 1
        return answer

    return integer_valuation(value.numerator) - integer_valuation(value.denominator)


def legendre_symbol(value: int, prime: int) -> int:
    """Return the Legendre symbol for an odd prime."""

    value %= prime
    if value == 0:
        return 0
    symbol = pow(value, (prime - 1) // 2, prime)
    return -1 if symbol == prime - 1 else symbol


def rational_square_root(value: Fraction) -> Fraction | None:
    """Return the nonnegative rational square root, or ``None`` if none exists."""

    if value < 0:
        return None
    numerator = isqrt(value.numerator)
    denominator = isqrt(value.denominator)
    if numerator * numerator != value.numerator:
        return None
    if denominator * denominator != value.denominator:
        return None
    return Q(numerator, denominator)


@dataclass(frozen=True)
class LocalData:
    prime: int
    residue: int
    good_reduction: bool
    point_count: int | None
    trace: int | None
    vanishing_factors: tuple[int, ...]
    split_multiplicative: bool | None


@dataclass(frozen=True)
class PowerRoot:
    prime: int
    exponent: int
    modulus: int
    residue: int
    factor_index: int
    split_multiplicative: bool


@dataclass(frozen=True)
class EKK3Family:
    """A fixed rational member ``u`` of the two-parameter K3 family."""

    u: Fraction
    m: Fraction | None = None

    def __post_init__(self) -> None:
        if self.m is not None:
            expected = 2 * (self.m**2 - self.m - 1) / (self.m**2 + 1)
            if expected != self.u:
                raise ValueError("u and m do not satisfy the Shimura-curve parametrization")

    @property
    def excluded_primes(self) -> tuple[int, ...]:
        denominators = [self.u.denominator]
        if self.m is not None:
            denominators.append(self.m.denominator)
        product = 2
        for denominator in denominators:
            product *= denominator
        return tuple(p for p in prime_factors(product))

    def a(self, t: Fraction) -> Fraction:
        u = self.u
        return (
            (u**8 - 18 * u**6 + 163 * u**4 - 1152 * u**2 + 4096) * t**4
            + (3 * u**7 - 35 * u**5 - 120 * u**3 + 1536 * u) * t**3
            + (u**8 - 13 * u**6 + 32 * u**4 - 152 * u**2 + 1536) * t**2
            + (u**7 + 3 * u**5 - 156 * u**3 + 672 * u) * t
            + (3 * u**6 - 33 * u**4 + 112 * u**2 - 80)
        )

    def factor_coefficients(self) -> tuple[tuple[Fraction, Fraction], ...]:
        """Return ``(slope, intercept)`` for the eight linear factors of B."""

        u = self.u
        return (
            (u**2 + u - 8, -u + 2),
            (u**2 - u - 8, -u - 2),
            (u**2 - u - 8, u**2 + u - 10),
            (u**2 + u - 8, -u**2 + u + 10),
            (u**2 - 7 * u + 8, -u**2 + u + 2),
            (u**2 + 7 * u + 8, u**2 + u - 2),
            (u**2 + 5 * u + 8, u**2 + 3 * u + 2),
            (u**2 - 5 * u + 8, -u**2 + 3 * u - 2),
        )

    def b_factors(self, t: Fraction) -> tuple[Fraction, ...]:
        return tuple(slope * t + intercept for slope, intercept in self.factor_coefficients())

    def b(self, t: Fraction) -> Fraction:
        answer = Q(1)
        for factor in self.b_factors(t):
            answer *= factor
        return answer

    def coefficients(self, t: Fraction) -> tuple[Fraction, ...]:
        """Return the extended Weierstrass vector ``(a1,a2,a3,a4,a6)``."""

        return (Q(0), 2 * self.a(t), Q(0), self.b(t), Q(0))

    def invariants(self, t: Fraction) -> dict[str, Fraction]:
        a_value = self.a(t)
        b_value = self.b(t)
        return {
            "c4": 16 * (4 * a_value**2 - 3 * b_value),
            "c6": 64 * a_value * (-8 * a_value**2 + 9 * b_value),
            "discriminant": 64 * b_value**2 * (a_value**2 - b_value),
        }

    def is_nonsingular(self, t: Fraction) -> bool:
        return self.invariants(t)["discriminant"] != 0

    def _check_prime(self, prime: int) -> None:
        if prime < 3:
            raise ValueError("local engineering is restricted to odd primes")
        if any(prime % excluded == 0 for excluded in self.excluded_primes):
            raise ValueError(f"prime {prime} divides a fixed family denominator")

    def a_mod(self, residue: int, prime: int) -> int:
        self._check_prime(prime)
        return fraction_mod(self.a(Q(residue)), prime)

    def factors_mod(self, residue: int, prime: int) -> tuple[int, ...]:
        self._check_prime(prime)
        return tuple(fraction_mod(value, prime) for value in self.b_factors(Q(residue)))

    def b_mod(self, residue: int, prime: int) -> int:
        answer = 1
        for factor in self.factors_mod(residue, prime):
            answer = answer * factor % prime
        return answer

    def discriminant_mod(self, residue: int, prime: int) -> int:
        a_value = self.a_mod(residue, prime)
        b_value = self.b_mod(residue, prime)
        return 64 * b_value * b_value * (a_value * a_value - b_value) % prime

    def local_data(self, residue: int, prime: int) -> LocalData:
        """Compute exact finite-field data by direct point counting.

        At bad residues only the clean case where exactly one B-factor vanishes
        and A is nonzero is classified as split/non-split multiplicative.
        """

        self._check_prime(prime)
        residue %= prime
        factors = self.factors_mod(residue, prime)
        vanishing = tuple(index + 1 for index, value in enumerate(factors) if value == 0)
        a_value = self.a_mod(residue, prime)
        b_value = 1
        for factor in factors:
            b_value = b_value * factor % prime
        discriminant = 64 * b_value * b_value * (a_value * a_value - b_value) % prime
        if discriminant == 0:
            split = None
            if len(vanishing) == 1 and a_value != 0:
                split = legendre_symbol(2 * a_value, prime) == 1
            return LocalData(prime, residue, False, None, None, vanishing, split)

        character_sum = 0
        for x_value in range(prime):
            rhs = (
                x_value**3 + 2 * a_value * x_value**2 + b_value * x_value
            ) % prime
            character_sum += legendre_symbol(rhs, prime)
        point_count = prime + 1 + character_sum
        trace = -character_sum
        return LocalData(prime, residue, True, point_count, trace, (), None)

    def power_roots(
        self, prime: int, exponent: int, *, split_only: bool = False
    ) -> tuple[PowerRoot, ...]:
        """Lift all clean roots of the linear B-factors modulo ``prime**exponent``.

        A clean root is unique among the eight factors and has A nonzero mod p.
        It therefore forces multiplicative reduction, provided the rational
        denominator of a later specialization remains prime to p.
        """

        self._check_prime(prime)
        if exponent < 1:
            raise ValueError("the lifting exponent must be positive")
        modulus = prime**exponent
        roots: dict[tuple[int, int], PowerRoot] = {}
        for index, (slope, intercept) in enumerate(self.factor_coefficients(), start=1):
            slope_mod_p = fraction_mod(slope, prime)
            if slope_mod_p == 0:
                continue
            root_mod_p = -fraction_mod(intercept, prime) * pow(slope_mod_p, -1, prime) % prime
            local = self.local_data(root_mod_p, prime)
            if local.vanishing_factors != (index,) or local.split_multiplicative is None:
                continue

            slope_modulus = fraction_mod(slope, modulus)
            intercept_modulus = fraction_mod(intercept, modulus)
            root = -intercept_modulus * pow(slope_modulus, -1, modulus) % modulus
            if (slope_modulus * root + intercept_modulus) % modulus != 0:
                raise AssertionError("linear Hensel lift failed")
            if split_only and not local.split_multiplicative:
                continue
            record = PowerRoot(
                prime=prime,
                exponent=exponent,
                modulus=modulus,
                residue=root,
                factor_index=index,
                split_multiplicative=bool(local.split_multiplicative),
            )
            roots[(root, index)] = record
        return tuple(sorted(roots.values(), key=lambda item: (item.residue, item.factor_index)))

    def known_points(self, t: Fraction) -> tuple[tuple[Fraction, Fraction], ...]:
        """Return the nine published section specializations when ``m`` is known.

        This verifies only that the specialized points lie on the curve.  It
        does not assert that they remain independent at this specialization.
        """

        if self.m is None:
            raise ValueError("the ninth section requires a Shimura parameter m")
        factors = self.b_factors(t)
        b1, b2, b3, b4, b5, b6, b7, b8 = factors
        x_values = (
            -b1 * b2 * b3 * b6,
            -b1 * b2 * b4 * b5,
            4 * b1 * b2 * b5 * b6,
            b1 * b3 * b4 * b6,
            -b1 * b3 * b4 * b7,
            b1 * b3 * b4 * b8,
            b1 * b3 * b5 * b6,
            -b1 * b5 * b6 * b7,
            -(self.m - 1) ** 2 * b1 * b2 * b3 * b8,
        )
        a_value = self.a(t)
        b_value = self.b(t)
        points: list[tuple[Fraction, Fraction]] = []
        for x_value in x_values:
            rhs = x_value**3 + 2 * a_value * x_value**2 + b_value * x_value
            y_value = rational_square_root(rhs)
            if y_value is None:
                raise AssertionError("a published section failed the curve equation")
            points.append((x_value, y_value))
        return tuple(points)


def prime_factors(integer: int) -> Iterable[int]:
    integer = abs(integer)
    divisor = 2
    while divisor * divisor <= integer:
        if integer % divisor == 0:
            yield divisor
            while integer % divisor == 0:
                integer //= divisor
        divisor = 3 if divisor == 2 else divisor + 2
    if integer > 1:
        yield integer


def primes_up_to(bound: int) -> tuple[int, ...]:
    """Return all primes at most ``bound`` by an elementary sieve."""

    if bound < 2:
        return ()
    sieve = bytearray(b"\x01") * (bound + 1)
    sieve[0:2] = b"\x00\x00"
    for prime in range(2, isqrt(bound) + 1):
        if sieve[prime]:
            sieve[prime * prime : bound + 1 : prime] = b"\x00" * (
                (bound - prime * prime) // prime + 1
            )
    return tuple(index for index, is_prime in enumerate(sieve) if is_prime)


def polynomial_eval(coefficients: Sequence[int], value: int, modulus: int) -> int:
    """Evaluate ascending-order integral coefficients modulo ``modulus``."""

    answer = 0
    for coefficient in reversed(coefficients):
        answer = (answer * value + coefficient) % modulus
    return answer


def polynomial_derivative(coefficients: Sequence[int]) -> tuple[int, ...]:
    return tuple(index * coefficient for index, coefficient in enumerate(coefficients[1:], 1))


def hensel_lift_simple_root(
    coefficients: Sequence[int], root: int, prime: int, exponent: int
) -> int:
    """Lift a simple root modulo p to a root modulo p**exponent.

    The Newton step is performed at growing moduli.  The derivative must be a
    unit modulo p; multiple-root lifting is deliberately outside this pilot.
    """

    if exponent < 1:
        raise ValueError("the lifting exponent must be positive")
    root %= prime
    if polynomial_eval(coefficients, root, prime) != 0:
        raise ValueError("the supplied value is not a root modulo p")
    derivative = polynomial_derivative(coefficients)
    derivative_mod_p = polynomial_eval(derivative, root, prime)
    if derivative_mod_p == 0:
        raise ValueError("Hensel lifting requires a simple root")

    modulus = prime
    lifted = root
    for _ in range(1, exponent):
        next_modulus = modulus * prime
        value = polynomial_eval(coefficients, lifted, next_modulus)
        quotient = (value // modulus) % prime
        correction = -quotient * pow(derivative_mod_p, -1, prime) % prime
        lifted += correction * modulus
        modulus = next_modulus
    if polynomial_eval(coefficients, lifted, prime**exponent) != 0:
        raise AssertionError("Hensel lift verification failed")
    return lifted % (prime**exponent)

