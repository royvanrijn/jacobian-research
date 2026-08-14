#!/usr/bin/env python3
"""Exact digit lifting for roots of integral polynomials modulo ``p^k``.

The usual Newton--Hensel formula only applies when the derivative is a unit
modulo ``p``.  Discriminant polynomials frequently have multiple roots, for
which a root at one level can have either no children or all ``p`` children
at the next level.  This module enumerates those branches exactly.

The root cap is a safety limit, not a truncation: exceeding it raises an
exception, so a returned result always contains *all* roots modulo the
requested prime power.  ``maximal_balls`` then compresses complete sibling
sets.  A ball ``r (mod p^j)`` means that every extension modulo ``p^k`` is a
root, exposing the true congruence cost of repeated-root conditions.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import comb, isqrt
from typing import Sequence


def _is_prime(integer: int) -> bool:
    if integer < 2:
        return False
    if integer % 2 == 0:
        return integer == 2
    divisor = 3
    while divisor <= isqrt(integer):
        if integer % divisor == 0:
            return False
        divisor += 2
    return True


def _polynomial_value(coefficients: Sequence[int], value: int) -> int:
    answer = 0
    for coefficient in reversed(coefficients):
        answer = answer * value + coefficient
    return answer


def _polynomial_value_mod(
    coefficients: Sequence[int], value: int, modulus: int
) -> int:
    answer = 0
    for coefficient in reversed(coefficients):
        answer = (answer * value + coefficient) % modulus
    return answer


def _derivative(coefficients: Sequence[int]) -> tuple[int, ...]:
    return tuple(
        index * coefficient
        for index, coefficient in enumerate(coefficients[1:], start=1)
    )


@dataclass(frozen=True, order=True)
class RootBall:
    """A residue class on which every refinement is a requested root."""

    prime: int
    exponent: int
    residue: int

    @property
    def modulus(self) -> int:
        return self.prime**self.exponent


@dataclass(frozen=True)
class PrimePowerRootResult:
    """The complete roots modulo a prime power and their lift profile."""

    prime: int
    exponent: int
    modulus: int
    roots: tuple[int, ...]
    level_counts: tuple[int, ...]
    candidate_digits_checked: int

    @property
    def density(self) -> Fraction:
        """Fraction of all residues modulo ``p^k`` that are roots."""

        return Fraction(len(self.roots), self.modulus)

    @property
    def reciprocal_density(self) -> Fraction | None:
        """Average search cost of accepting any root; ``None`` if there are none."""

        if not self.roots:
            return None
        return Fraction(self.modulus, len(self.roots))

    def maximal_balls(self) -> tuple[RootBall, ...]:
        """Compress the root set into disjoint, maximal ``p``-adic balls.

        The output covers exactly ``roots`` modulo ``p^exponent``.  Exponent
        zero is allowed: it denotes the unique class modulo one, i.e. the
        polynomial value is always divisible by the requested prime power.
        """

        balls: set[tuple[int, int]] = {
            (self.exponent, root) for root in self.roots
        }
        for child_exponent in range(self.exponent, 0, -1):
            child_modulus = self.prime**child_exponent
            parent_modulus = child_modulus // self.prime
            parents = {
                residue % parent_modulus
                for exponent, residue in balls
                if exponent == child_exponent
            }
            for parent in parents:
                children = {
                    (
                        child_exponent,
                        parent + digit * parent_modulus,
                    )
                    for digit in range(self.prime)
                }
                if children <= balls:
                    balls.difference_update(children)
                    balls.add((child_exponent - 1, parent))

        answer = tuple(
            sorted(
                (
                    RootBall(self.prime, exponent, residue)
                    for exponent, residue in balls
                ),
                key=lambda ball: (ball.exponent, ball.residue),
            )
        )
        covered = sum(
            self.prime ** (self.exponent - ball.exponent) for ball in answer
        )
        if covered != len(self.roots):
            raise AssertionError("compressed root balls do not preserve the root set")
        return answer


class RootLiftCapExceeded(RuntimeError):
    """Raised instead of returning an incomplete prime-power root set."""

    def __init__(
        self,
        *,
        prime: int,
        requested_exponent: int,
        reached_exponent: int,
        cap: int,
    ) -> None:
        self.prime = prime
        self.requested_exponent = requested_exponent
        self.reached_exponent = reached_exponent
        self.cap = cap
        super().__init__(
            f"more than {cap} roots occur modulo {prime}^{reached_exponent}; "
            f"cannot exactly enumerate roots modulo {prime}^{requested_exponent}"
        )


def verify_prime_power_roots(
    coefficients: Sequence[int], result: PrimePowerRootResult
) -> None:
    """Verify normalization, distinctness, and every modular root identity."""

    roots = result.roots
    if tuple(sorted(set(roots))) != roots:
        raise AssertionError("roots must be sorted and distinct")
    if any(root < 0 or root >= result.modulus for root in roots):
        raise AssertionError("a root is outside the canonical residue range")
    for root in roots:
        if _polynomial_value_mod(coefficients, root, result.modulus) != 0:
            raise AssertionError(
                f"{root} is not a root modulo {result.prime}^{result.exponent}"
            )


def all_roots_mod_prime_power(
    coefficients: Sequence[int],
    prime: int,
    exponent: int,
    *,
    max_roots: int | None = 100_000,
    verify: bool = True,
) -> PrimePowerRootResult:
    """Return every root of an integral polynomial modulo ``prime**exponent``.

    Starting from all roots modulo ``p``, each root ``r (mod p^j)`` is lifted
    through the digits ``r + d*p^j``.  The exact first-order congruence

    ``f(r+d*p^j)/p^j = f(r)/p^j + d*f'(r) (mod p)``

    chooses the valid digits.  It remains valid when ``f'(r)=0 (mod p)``:
    then either every digit works or none does.

    ``max_roots`` bounds the number of surviving roots at any level.  The
    function raises :class:`RootLiftCapExceeded` rather than returning a
    partial set.  Pass ``None`` only when the possible ``p``-fold explosion
    is known to be manageable.
    """

    coefficients = tuple(coefficients)
    if not coefficients:
        raise ValueError("the coefficient sequence must not be empty")
    if not all(isinstance(coefficient, int) for coefficient in coefficients):
        raise TypeError("all polynomial coefficients must be integers")
    if not any(coefficients):
        raise ValueError("the zero polynomial has too many roots to enumerate")
    if not _is_prime(prime):
        raise ValueError("the modulus base must be prime")
    if exponent < 1:
        raise ValueError("the lifting exponent must be positive")
    if max_roots is not None and max_roots < 1:
        raise ValueError("max_roots must be positive or None")

    derivative = _derivative(coefficients)
    roots = [
        residue
        for residue in range(prime)
        if _polynomial_value_mod(coefficients, residue, prime) == 0
    ]
    if max_roots is not None and len(roots) > max_roots:
        raise RootLiftCapExceeded(
            prime=prime,
            requested_exponent=exponent,
            reached_exponent=1,
            cap=max_roots,
        )
    counts = [len(roots)]
    candidate_digits_checked = prime
    modulus = prime

    for reached_exponent in range(2, exponent + 1):
        next_modulus = modulus * prime
        next_roots: list[int] = []
        for root in roots:
            value_mod_next = _polynomial_value_mod(
                coefficients, root, next_modulus
            )
            if value_mod_next % modulus != 0:
                raise AssertionError("the previous-level root invariant failed")
            quotient = (value_mod_next // modulus) % prime
            derivative_mod_prime = _polynomial_value_mod(
                derivative, root, prime
            )

            if derivative_mod_prime:
                digits = (
                    -quotient * pow(derivative_mod_prime, -1, prime) % prime,
                )
                candidate_digits_checked += 1
            elif quotient:
                digits = ()
                candidate_digits_checked += prime
            else:
                digits = range(prime)
                candidate_digits_checked += prime

            for digit in digits:
                next_roots.append(root + digit * modulus)
                if max_roots is not None and len(next_roots) > max_roots:
                    raise RootLiftCapExceeded(
                        prime=prime,
                        requested_exponent=exponent,
                        reached_exponent=reached_exponent,
                        cap=max_roots,
                    )
        roots = sorted(next_roots)
        modulus = next_modulus
        counts.append(len(roots))

    result = PrimePowerRootResult(
        prime=prime,
        exponent=exponent,
        modulus=modulus,
        roots=tuple(roots),
        level_counts=tuple(counts),
        candidate_digits_checked=candidate_digits_checked,
    )
    if verify:
        verify_prime_power_roots(coefficients, result)
    return result


def scaled_variable_coefficients(
    coefficients: Sequence[int], scale: int
) -> tuple[int, ...]:
    """Return the coefficients of ``f(scale*x)`` in ascending order."""

    return affine_variable_coefficients(coefficients, 0, scale)


def affine_variable_coefficients(
    coefficients: Sequence[int], offset: int, scale: int
) -> tuple[int, ...]:
    """Return the ascending coefficients of ``f(offset + scale*x)``."""

    coefficients = tuple(coefficients)
    if not coefficients:
        raise ValueError("the coefficient sequence must not be empty")
    if not all(isinstance(coefficient, int) for coefficient in coefficients):
        raise TypeError("all polynomial coefficients must be integers")
    if not isinstance(offset, int) or not isinstance(scale, int):
        raise TypeError("the affine parameters must be integers")
    answer = [0] * len(coefficients)
    for source_degree, coefficient in enumerate(coefficients):
        for target_degree in range(source_degree + 1):
            answer[target_degree] += (
                coefficient
                * comb(source_degree, target_degree)
                * offset ** (source_degree - target_degree)
                * scale**target_degree
            )
    return tuple(answer)


def fixed_divisor_valuation(coefficients: Sequence[int], prime: int) -> int:
    """Return ``v_p(gcd(f(n): n in Z))`` for a nonzero integral polynomial.

    For degree ``d``, the values at ``0,...,d`` suffice.  This follows from
    the integer-valued Newton expansion in binomial polynomials.  The result
    can prove automatic prime-power divisibility that coefficient content
    alone does not reveal.
    """

    coefficients = tuple(coefficients)
    if not coefficients or not any(coefficients):
        raise ValueError("the polynomial must be nonzero")
    if not all(isinstance(coefficient, int) for coefficient in coefficients):
        raise TypeError("all polynomial coefficients must be integers")
    if not _is_prime(prime):
        raise ValueError("the valuation base must be prime")
    degree = max(index for index, coefficient in enumerate(coefficients) if coefficient)

    def integer_valuation(integer: int) -> int:
        if integer == 0:
            return 10**9
        value = abs(integer)
        answer = 0
        while value % prime == 0:
            value //= prime
            answer += 1
        return answer

    return min(
        integer_valuation(_polynomial_value(coefficients, integer))
        for integer in range(degree + 1)
    )
