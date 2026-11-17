#!/usr/bin/env python3
"""Exact regressions for the cancellation-parameter irreducibility results."""
from __future__ import annotations

import bisect
import math

from flint import nmod_poly
import sympy as sp

from master_cancellation import (
    DISPLAYED_INSTANCES,
    parameter_polynomial,
    raw_parameter_polynomial,
)


q = sp.symbols("q")
x = sp.symbols("x")


def valuation(value: int, prime: int) -> int:
    order = 0
    while value % prime == 0:
        value //= prime
        order += 1
    return order


def prime_power(value: int) -> tuple[int, int] | None:
    factors = sp.factorint(value)
    if len(factors) != 1:
        return None
    prime, exponent = next(iter(factors.items()))
    return int(prime), int(exponent)


def assert_eisenstein(m: int, r: int, prime: int) -> None:
    polynomial = sp.Poly(parameter_polynomial(m, r, q), q, domain=sp.ZZ)
    coefficients = polynomial.all_coeffs()
    assert coefficients[0] == 1
    assert all(int(coefficient) % prime == 0 for coefficient in coefficients[1:])
    assert int(coefficients[-1]) % (prime * prime) != 0


def subset_sum_mask(degrees: list[int], maximum: int) -> int:
    """Encode attainable subset sums through ``maximum`` as a bit mask."""
    sums = 1
    mask = (1 << (maximum + 1)) - 1
    for degree in degrees:
        if degree <= maximum:
            sums |= (sums << degree) & mask
    return sums


def modular_factor_degrees(polynomial: sp.Poly, prime: int) -> list[int]:
    """Factor a monic integral polynomial over ``F_prime`` using FLINT."""
    coefficients = [
        int(polynomial.nth(index)) % prime
        for index in range(polynomial.degree() + 1)
    ]
    factors = nmod_poly(coefficients, prime).factor()[1]
    degrees: list[int] = []
    for factor, multiplicity in factors:
        degrees.extend([len(factor) - 1] * multiplicity)
    return degrees


def binomial_mod_prime(total: int, index: int, prime: int) -> int:
    """Return binomial(total,index) modulo a prime via Lucas's theorem."""
    result = 1
    while total or index:
        total_digit = total % prime
        index_digit = index % prime
        if index_digit > total_digit:
            return 0
        result = result * math.comb(total_digit, index_digit) % prime
        total //= prime
        index //= prime
    return result


def modular_parameter_polynomial(
    degree: int, total: int, prime: int
) -> nmod_poly:
    """Construct M_(m,r) modulo ``prime`` without huge integer coefficients."""
    coefficients = [
        (-1) ** (degree - index)
        * binomial_mod_prime(total, degree - index, prime)
        % prime
        for index in range(degree + 1)
    ]
    return nmod_poly(coefficients, prime)


def flint_factor_degrees(polynomial: nmod_poly) -> list[int]:
    """Return irreducible factor degrees with multiplicity."""
    factors = polynomial.factor()[1]
    degrees: list[int] = []
    for factor, multiplicity in factors:
        degrees.extend([len(factor) - 1] * multiplicity)
    return degrees


def modular_parameter_factor_degrees(
    degree: int, total: int, prime: int
) -> list[int]:
    """Factor M_(m,r) modulo ``prime``."""
    return flint_factor_degrees(
        modular_parameter_polynomial(degree, total, prime)
    )


def modular_small_factor_degrees(
    degree: int, total: int, prime: int, maximum: int
) -> list[int] | None:
    """Return all factor degrees through ``maximum`` at a squarefree prime.

    The degrees are recovered by distinct-degree factorization: the degree of
    gcd(f,x^(p^d)-x) is the sum of e*n_e over e|d, where n_e is the number of
    irreducible degree-e factors.  ``None`` marks a nonsquarefree reduction.
    """
    polynomial = modular_parameter_polynomial(degree, total, prime)
    if len(polynomial.gcd(polynomial.derivative())) > 1:
        return None

    variable = nmod_poly([0, 1], prime)
    frobenius = variable
    factor_counts = [0] * (maximum + 1)
    degrees: list[int] = []
    for current_degree in range(1, maximum + 1):
        frobenius = frobenius.pow_mod(prime, polynomial)
        divisor_degree = len((frobenius - variable).gcd(polynomial)) - 1
        known_degree = sum(
            divisor * factor_counts[divisor]
            for divisor in range(1, current_degree)
            if current_degree % divisor == 0
        )
        count = (divisor_degree - known_degree) // current_degree
        assert count >= 0
        assert known_degree + current_degree * count == divisor_degree
        factor_counts[current_degree] = count
        degrees.extend([current_degree] * count)
    return degrees


certified: list[tuple[int, int, int, int]] = []
for m in range(1, 21):
    for r in range(1, 13):
        n = m * r
        total = n + r + 1
        power = prime_power(total)
        if power is None:
            continue
        prime, exponent = power
        if valuation(n, prime) != exponent - 1:
            continue
        assert_eisenstein(m, r, prime)
        certified.append((m, r, prime, exponent))

displayed_pairs = {(instance.m, instance.r) for instance in DISPLAYED_INSTANCES}
assert displayed_pairs <= {(m, r) for m, r, _, _ in certified}

cyclotomic: list[tuple[int, int, int, int]] = []
for m in range(1, 13):
    for r in range(1, 13):
        n = m * r
        total = n + r + 1
        prime = total + 1
        ell = n + 1
        if not (sp.isprime(prime) and sp.isprime(ell)):
            continue
        if sp.n_order(prime % ell, ell) != ell - 1:
            continue

        polynomial = sp.Poly(parameter_polynomial(m, r, q), q, modulus=prime)
        geometric = sp.Poly(sum(q**j for j in range(n + 1)), q, modulus=prime)
        assert polynomial == geometric
        assert polynomial.is_irreducible
        cyclotomic.append((m, r, prime, ell))

leading_prime: list[tuple[int, int, int]] = []
for m in range(1, 31):
    for r in range(1, 13):
        n = m * r
        leading = int(sp.binomial(n + r, r))
        if not sp.isprime(leading):
            continue

        transformed = sp.cancel(
            x**n * raw_parameter_polynomial(m, r, q).subs(q, 1 - 1 / x)
        )
        binomial_sum = sum(sp.binomial(j + r, r) * x**j for j in range(n + 1))
        scale = sp.factorial(n) * sp.factorial(r) / sp.factorial(n + r + 1)
        assert sp.cancel(transformed - scale * binomial_sum) == 0
        assert sp.Poly(binomial_sum, x, domain=sp.QQ).is_irreducible
        leading_prime.append((m, r, leading))

# Translation to the classical truncated binomial P_{N,k}.  The cited
# Khanduja--Khassa--Laishram theorem applies uniformly when m=1; the bounded
# loop checks the exact reciprocal identification and representative cases.
truncated_binomial_column: list[int] = []
for r in range(1, 31):
    n = r
    total = 2 * r + 1
    truncated = sum(sp.binomial(total, j) * x**j for j in range(n + 1))
    reciprocal = sp.cancel(q**n * truncated.subs(x, -1 / q))
    assert sp.expand(reciprocal - parameter_polynomial(1, r, q)) == 0
    assert 2 <= 2 * n <= total < (n + 1) ** 3
    assert sp.Poly(truncated, x, domain=sp.QQ).is_irreducible
    truncated_binomial_column.append(r)

# A prime p=N-u in (k,N) gives the V-shaped Newton polygon with vertices
# (0,1), (u,0), (k,1) after translating P_{N,k}(x) by x -> x-1.  Two such
# primes give incompatible factor-degree pairs on the diagonal when m>=2.
two_prime_newton_pairs = 0
for m in range(2, 9):
    for r in range(1, 11):
        n = m * r
        total = n + r + 1
        primes = list(sp.primerange(n + 1, total))
        if len(primes) < 2:
            continue

        truncated = sum(sp.binomial(total, j) * x**j for j in range(n + 1))
        translated = sp.Poly(sp.expand(truncated.subs(x, x - 1)), x, domain=sp.ZZ)
        degree_pairs: list[tuple[int, int]] = []
        for prime in primes[:2]:
            u = total - int(prime)
            assert 1 <= u <= r <= n // 2
            for j in range(n + 1):
                expected = 0 if j == u else 1
                assert valuation(int(translated.nth(j)), int(prime)) == expected
            degree_pairs.append((u, n - u))
        assert degree_pairs[0] != degree_pairs[1]
        two_prime_newton_pairs += 1

# Dusart's explicit interval theorem gives a prime in
# (x,x(1+1/log(x)^3)] for x>=89693.  Two iterations lie below
# x(1001/1000)^2, since log(x)>10.  Thus for every m<=300 the
# two-prime theorem handles every k=mr at or above this cutoff.
DUSART_THRESHOLD = 89_693
COMPLETE_COLUMN_MAX = 300
assert DUSART_THRESHOLD > 3**10
assert 1001**2 * COMPLETE_COLUMN_MAX < (COMPLETE_COLUMN_MAX + 1) * 1000**2
finite_prime_limit = DUSART_THRESHOLD + (DUSART_THRESHOLD - 1) // 2 + 3
finite_primes = [
    int(prime) for prime in sp.primerange(2, finite_prime_limit)
]
finite_prime_flags = bytearray(finite_prime_limit)
for prime in finite_primes:
    finite_prime_flags[prime] = 1
finite_prime_prefix = [0] * finite_prime_limit
for value in range(1, finite_prime_limit):
    finite_prime_prefix[value] = (
        finite_prime_prefix[value - 1] + finite_prime_flags[value]
    )

# Exhaust the finite part of the 299 columns m=2,...,300.  Prime N is the
# prime-total Eisenstein case; two interval primes invoke the theorem above.
# The remainder is small enough for exact modular degree certificates.
column_residual_targets: list[tuple[int, int, int | None]] = []
column_zero_prime_targets = 0
column_one_prime_targets = 0
column_residual_max_degree = 0
column_max_forced_degree = 0
for m in range(2, COMPLETE_COLUMN_MAX + 1):
    for r in range(1, (DUSART_THRESHOLD - 1) // m + 1):
        degree = m * r
        total = degree + r + 1
        interval_prime_count = (
            finite_prime_prefix[total - 1] - finite_prime_prefix[degree]
        )
        if finite_prime_flags[total] or interval_prime_count >= 2:
            continue

        forced_degree: int | None = None
        if interval_prime_count == 1:
            prime_index = bisect.bisect_right(finite_primes, degree)
            interval_prime = finite_primes[prime_index]
            assert degree < interval_prime < total
            forced_degree = total - interval_prime
            column_one_prime_targets += 1
            column_max_forced_degree = max(
                column_max_forced_degree, forced_degree
            )
        else:
            column_zero_prime_targets += 1

        column_residual_targets.append((m, r, forced_degree))
        column_residual_max_degree = max(column_residual_max_degree, degree)

assert len(column_residual_targets) == 2_598
assert column_zero_prime_targets == 843
assert column_one_prime_targets == 1_755
assert column_residual_max_degree == 9_555
assert column_max_forced_degree == 28

# The same elementary Dusart estimate controls the tails through m=499.
# Record exactly what remains to enlarge the finite replay from 300 to that
# natural cutoff.
DUSART_COARSE_MAX = 499
assert 1001**2 * DUSART_COARSE_MAX < (DUSART_COARSE_MAX + 1) * 1000**2
frontier_residual_count = 0
frontier_zero_prime_count = 0
frontier_one_prime_count = 0
frontier_max_degree = 0
for m in range(COMPLETE_COLUMN_MAX + 1, DUSART_COARSE_MAX + 1):
    for r in range(1, (DUSART_THRESHOLD - 1) // m + 1):
        degree = m * r
        total = degree + r + 1
        interval_prime_count = (
            finite_prime_prefix[total - 1] - finite_prime_prefix[degree]
        )
        if finite_prime_flags[total] or interval_prime_count >= 2:
            continue
        frontier_residual_count += 1
        frontier_zero_prime_count += interval_prime_count == 0
        frontier_one_prime_count += interval_prime_count == 1
        frontier_max_degree = max(frontier_max_degree, degree)

assert frontier_residual_count == 2_192
assert frontier_zero_prime_count == 816
assert frontier_one_prime_count == 1_376
assert frontier_max_degree == 31_395

auxiliary_primes = [int(prime) for prime in sp.primerange(2, 300)]

# Cross-check the Lucas/FLINT fast path against direct construction before it
# is used for the large finite exhaustion.
for m in range(1, 5):
    for r in range(1, 5):
        polynomial = sp.Poly(parameter_polynomial(m, r, q), q, domain=sp.ZZ)
        degree = m * r
        total = degree + r + 1
        for prime in (2, 3, 5, 7):
            direct_degrees = sorted(modular_factor_degrees(polynomial, prime))
            assert sorted(
                modular_parameter_factor_degrees(degree, total, prime)
            ) == direct_degrees
            small_degrees = modular_small_factor_degrees(
                degree, total, prime, min(degree, 6)
            )
            if small_degrees is not None:
                assert sorted(small_degrees) == [
                    value for value in direct_degrees if value <= 6
                ]

degree_sieve_max_prime = 0
degree_sieve_max_steps = 0
degree_sieve_pairs = 0
for m, r, forced_degree in column_residual_targets:
    degree = m * r
    total = degree + r + 1
    steps = 0
    certified_irreducible = False

    if forced_degree is not None:
        # One interval prime forces a reducible polynomial to have an
        # irreducible factor of this degree.  Distinct-degree factorization
        # only through that small degree is enough to exclude it.
        for prime in auxiliary_primes:
            degrees = modular_small_factor_degrees(
                degree, total, prime, forced_degree
            )
            if degrees is None:
                continue
            steps += 1
            modular_degrees = subset_sum_mask(degrees, forced_degree)
            if not (modular_degrees >> forced_degree) & 1:
                certified_irreducible = True
                degree_sieve_max_prime = max(degree_sieve_max_prime, prime)
                degree_sieve_max_steps = max(degree_sieve_max_steps, steps)
                break
    else:
        possible_degrees = (1 << (degree + 1)) - 1
        for prime in auxiliary_primes:
            polynomial = modular_parameter_polynomial(degree, total, prime)
            if len(polynomial.gcd(polynomial.derivative())) > 1:
                continue
            degrees = flint_factor_degrees(polynomial)
            possible_degrees &= subset_sum_mask(degrees, degree)
            steps += 1
            if possible_degrees == 1 | (1 << degree):
                certified_irreducible = True
                degree_sieve_max_prime = max(degree_sieve_max_prime, prime)
                degree_sieve_max_steps = max(degree_sieve_max_steps, steps)
                break

    if certified_irreducible:
        if forced_degree is None:
            assert possible_degrees == 1 | (1 << degree)
    else:
        raise AssertionError(f"missing degree certificate for {(m, r)}")
    degree_sieve_pairs += 1

# Retain the earlier bounded complete-degree checks.
legacy_degree_sieve_pairs = 0
degree_sieve_targets = {
    (m, r)
    for m in range(1, 31)
    for r in range(1, 31)
    if m * r <= 30
}
degree_sieve_targets |= {
    (m, r)
    for m in range(2, 7)
    for r in range(1, 118)
    if m * r < 118
}
for m, r in sorted(degree_sieve_targets):
    polynomial = sp.Poly(parameter_polynomial(m, r, q), q, domain=sp.ZZ)
    degree = polynomial.degree()
    possible_degrees = (1 << (degree + 1)) - 1
    steps = 0
    for prime in auxiliary_primes:
        # Monicity preserves every rational factor degree under reduction,
        # even at a discriminant prime, so no good-reduction filter is needed.
        degrees = modular_factor_degrees(polynomial, prime)
        possible_degrees &= subset_sum_mask(degrees, degree)
        steps += 1
        if possible_degrees == 1 | (1 << degree):
            degree_sieve_max_prime = max(degree_sieve_max_prime, prime)
            degree_sieve_max_steps = max(degree_sieve_max_steps, steps)
            break
    assert possible_degrees == 1 | (1 << degree)
    legacy_degree_sieve_pairs += 1

prime_cases = sum(1 for _, _, _, exponent in certified if exponent == 1)
print(
    "PASS: prime-power Eisenstein certificates for "
    f"{len(certified)} pairs ({prime_cases} prime-total cases)"
)
print("PASS: every displayed parameter polynomial has an Eisenstein certificate")
print(f"PASS: irreducible cyclotomic reductions for {len(cyclotomic)} pairs")
print(f"PASS: unit-disk leading-prime certificates for {len(leading_prime)} pairs")
print("PASS: truncated-binomial theorem translation for the full m=1 column")
print(
    "PASS: translated two-prime Newton polygons for "
    f"{two_prime_newton_pairs} representative diagonal pairs"
)
print(
    "PASS: Dusart tail plus exact finite certificates prove all "
    "1<=m<=300 columns"
)
print(
    f"PASS: {degree_sieve_pairs} residual pairs below {DUSART_THRESHOLD} "
    f"({column_zero_prime_targets} zero-prime, "
    f"{column_one_prime_targets} one-prime; maximum degree "
    f"{column_residual_max_degree})"
)
print(
    "PASS: separated Dusart frontier 301<=m<=499 has exactly "
    f"{frontier_residual_count} residual pairs "
    f"({frontier_zero_prime_count} zero-prime, "
    f"{frontier_one_prime_count} one-prime; maximum degree "
    f"{frontier_max_degree})"
)
print(
    f"PASS: retained modular degree-sieve certificates for "
    f"{legacy_degree_sieve_pairs} earlier bounded pairs "
    f"(largest prime {degree_sieve_max_prime}, at most {degree_sieve_max_steps} steps)"
)
