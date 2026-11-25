#!/usr/bin/env python3
"""Replay the exact finite prefix and contour endpoint for 1 plus XY."""

from bisect import bisect_right
from fractions import Fraction
from math import factorial, isqrt


LAST_FINITE_TERMINAL = 99_999


def primes_through(limit):
    """List the primes through limit by the sieve of Eratosthenes."""
    sieve = bytearray(b"\1") * (limit + 1)
    sieve[0:2] = b"\0\0"
    for prime in range(2, isqrt(limit) + 1):
        if sieve[prime]:
            start = prime * prime
            sieve[start::prime] = b"\0" * ((limit - start) // prime + 1)
    return [value for value in range(2, limit + 1) if sieve[value]]


def reciprocal_mod(order, prime):
    """Compute R_order over the given prime field, constant term first."""
    coefficient = 1
    polynomial = [coefficient]
    for k in range(1, order + 1):
        coefficient = coefficient * (order - k + 1) * k % prime
        polynomial.append(coefficient)
    polynomial.reverse()
    return polynomial


def coprime_mod(left, right, prime):
    """Test coprimality by Euclidean division over a prime field."""
    while right:
        inverse = pow(right[-1], -1, prime)
        while len(left) >= len(right):
            scale = left[-1] * inverse % prime
            offset = len(left) - len(right)
            for index, coefficient in enumerate(right):
                left[offset + index] = (
                    left[offset + index] - scale * coefficient
                ) % prime
            while left and left[-1] == 0:
                left.pop()
        left, right = right, left
    return len(left) == 1


def check_pair(order, prime):
    return coprime_mod(
        reciprocal_mod(order, prime),
        reciprocal_mod(order - 1, prime),
        prime,
    )


def verify_finite_prefix():
    primes = primes_through(LAST_FINITE_TERMINAL)
    prime_terminals = near_prime_terminals = 0
    direct_moduli = {4: 5, 10: 11}
    alternate_primes = {50_533: 50_513}
    direct_cases = []
    alternate_cases = []
    largest_primary_gap = 0

    for terminal in range(2, LAST_FINITE_TERMINAL + 1):
        index = bisect_right(primes, terminal) - 1
        prime = primes[index]
        if prime == terminal:
            prime_terminals += 1
            continue

        remainder = terminal - prime
        largest_primary_gap = max(largest_primary_gap, remainder)
        if 2 * remainder < prime - 2 and check_pair(remainder, prime):
            near_prime_terminals += 1
            continue

        if terminal in direct_moduli:
            modulus = direct_moduli[terminal]
            if modulus not in primes or not check_pair(terminal, modulus):
                raise ArithmeticError((terminal, modulus))
            direct_cases.append((terminal, modulus))
            continue

        alternative = alternate_primes.get(terminal)
        if alternative is None or alternative not in primes:
            raise ArithmeticError((terminal, prime, remainder))
        remainder = terminal - alternative
        if not (0 < 2 * remainder < alternative - 2):
            raise ArithmeticError((terminal, alternative, remainder))
        if not check_pair(remainder, alternative):
            raise ArithmeticError((terminal, alternative, remainder))
        alternate_cases.append((terminal, alternative, remainder))

    checked = (
        prime_terminals + near_prime_terminals + len(direct_cases)
        + len(alternate_cases)
    )
    if checked != LAST_FINITE_TERMINAL - 1:
        raise ArithmeticError(checked)

    print(f"Terminals 2 through {LAST_FINITE_TERMINAL}: {checked}")
    print(f"Prime terminals: {prime_terminals}")
    print(f"Primary near-prime terminals: {near_prime_terminals}")
    print(f"Largest gap from the preceding prime: {largest_primary_gap}")
    print(f"Alternate near-prime certificates: {alternate_cases}")
    print(f"Direct modular certificates: {direct_cases}")
    print("Unresolved terminals: 0")


def verify_contour_endpoint():
    """Prove the squared real-endpoint-to-saddle ratio is below 24/25."""
    c = Fraction(71, 100)
    checked = 0
    for k in range(530, 1310):
        low, high = Fraction(k, 1000), Fraction(k + 1, 1000)
        distance = max(Fraction(0), 1 - high**2, low**2 - 1)
        if distance > c:
            continue

        t = (c**2 - distance**2) / (2 * (1 + low**2))
        x = 2 * high
        exp_upper = sum(x**j / factorial(j) for j in range(17))
        exp_upper += x**17 / factorial(17) / (1 - x / 18)
        numerator = exp_upper * (
            (2 - low**2)**2 + (4 + 2 * high**2) * t + t**2
        )
        if numerator >= Fraction(24, 25) * 4 * (1 + low)**2:
            raise ArithmeticError(k)
        checked += 1

    print(f"Exact contour endpoint intervals checked: {checked}")


if __name__ == "__main__":
    verify_contour_endpoint()
    verify_finite_prefix()
