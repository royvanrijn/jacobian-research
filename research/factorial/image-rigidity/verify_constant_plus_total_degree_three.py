#!/usr/bin/env python3
"""Verify the degree-three contour margins and finite terminal prefix."""

from bisect import bisect_right
from fractions import Fraction as Q

from verify_constant_plus_product_of_two_variables import coprime_mod, primes_through


LAST_FINITE_TERMINAL = 99_999
DIRECT_MODULI = {
    (1, 1, 1): {4: 5, 9: 11, 10: 11},
    (2, 1): {
        2: 3, 3: 7, 4: 7, 5: 11, 6: 13, 7: 13,
        9: 17, 10: 19, 11: 23, 13: 29, 17: 37,
    },
}


def reciprocal_mod(order, prime, profile):
    """Compute the monic reciprocal over the given prime field."""
    coefficient = 1
    polynomial = [coefficient]
    for k in range(1, order + 1):
        if profile == (1, 1, 1):
            factorial_step = k * k
        else:
            factorial_step = (2 * k) * (2 * k - 1)
        coefficient = coefficient * (order - k + 1) * factorial_step % prime
        polynomial.append(coefficient)
    polynomial.reverse()
    return polynomial


def check_pair(order, prime, profile):
    return coprime_mod(
        reciprocal_mod(order, prime, profile),
        reciprocal_mod(order - 1, prime, profile),
        prime,
    )


def verify_finite_prefix(profile, primes):
    prime_count = near_count = alternate_count = 0
    direct = []
    largest_checked_remainder = 0
    for terminal in range(2, LAST_FINITE_TERMINAL + 1):
        index = bisect_right(primes, terminal) - 1
        largest = primes[index]
        if profile == (1, 1, 1) and largest == terminal:
            prime_count += 1
            continue

        for shift in range(min(4, index + 1)):
            prime = primes[index - shift]
            remainder = terminal - prime
            if profile == (1, 1, 1):
                admissible = 0 < 2 * remainder < prime - 2
            else:
                admissible = 0 < 3 * remainder < prime - 3
            if not admissible:
                continue
            if check_pair(remainder, prime, profile):
                near_count += 1
                alternate_count += int(shift != 0)
                largest_checked_remainder = max(largest_checked_remainder, remainder)
                break
        else:
            modulus = DIRECT_MODULI[profile].get(terminal)
            if modulus is None or modulus not in primes:
                raise ArithmeticError((profile, terminal, largest))
            if not check_pair(terminal, modulus, profile):
                raise ArithmeticError((profile, terminal, modulus))
            direct.append((terminal, modulus))

    if prime_count + near_count + len(direct) != LAST_FINITE_TERMINAL - 1:
        raise ArithmeticError(profile)

    print(f"Profile {profile}: {LAST_FINITE_TERMINAL - 1} terminals")
    print(f"Prime terminals: {prime_count}")
    print(f"Near-prime certificates: {near_count}, with {alternate_count} alternatives")
    print(f"Largest degree in a successful near-prime gcd: {largest_checked_remainder}")
    print(f"Direct modular certificates: {direct}")
    print("Unresolved terminals: 0")


def verify_rational_margins():
    """Check the scalar comparisons used in the quantitative saddle."""
    minimum_radius = Q(101, 40)
    g = Q(51, 500)
    third_derivative = (
        6 / minimum_radius**3
        + Q(15, 4) * g / (1 - g)
        + 3 * Q(19, 10) * Q(6, 5) * g**2 / (1 - g)**2
        + 2 * (Q(6, 5) * g)**3 / (1 - g)**3
    )
    comparisons = (
        Q(163, 60)**3 > 20,
        Q(3, 2) * Q(37, 40)**2 > Q(6, 5),
        Q(3, 20) * (1 + Q(43, 20)) < Q(1, 2),
        Q(81, 50) / Q(3, 2)**3 < Q(1, 2),
        Q(3, 40) + Q(3, 40)**2 / Q(37, 40) < Q(1, 10),
        Q(3, 11) + Q(5, 6)**3 / (1 - Q(27, 800)) < Q(22, 25),
        Q(81, 50) / minimum_radius**3 < g,
        3 / minimum_radius < Q(6, 5),
        12 / minimum_radius**2 < Q(19, 10),
        60 / minimum_radius**3 < Q(15, 4),
        third_derivative < Q(9, 10),
        (8 * Q(25, 24)**2)**3 < Q(11, 2)**4,
        Q(459, 384) * (Q(3, 2) + Q(11, 10)) < 4,
        Q(4, 223) + Q(10, 2**50) < Q(1, 40),
        Q(10, 316) + Q(2, 25) + Q(10 * 317, 2**20) < Q(1, 8),
        Q(10, 9) * Q(101, 100) / Q(7, 6) < Q(97, 100),
        Q(11, 4)**4 / 400 < Q(3, 20),
        Q(44, 45) / Q(99, 100) < Q(99, 100),
        Q(111, 40) / (3 * Q(7, 6)) > Q(3, 4),
    )
    if not all(comparisons):
        raise ArithmeticError(("rational margins", comparisons))
    print(f"Exact contour scalar comparisons: {len(comparisons)}")


if __name__ == "__main__":
    verify_rational_margins()
    prime_table = primes_through(LAST_FINITE_TERMINAL)
    for exponents in DIRECT_MODULI:
        verify_finite_prefix(exponents, prime_table)
