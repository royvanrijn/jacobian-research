#!/usr/bin/env python3
"""Exact regressions for the C24 prime-power Eisenstein theorem."""
from __future__ import annotations

import warnings

import sympy as sp
from sympy.utilities.exceptions import SymPyDeprecationWarning

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


def subset_sums(degrees: list[int]) -> set[int]:
    sums = {0}
    for degree in degrees:
        sums |= {value + degree for value in tuple(sums)}
    return sums


def modular_factor_degrees(polynomial: sp.Poly, prime: int) -> list[int]:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", SymPyDeprecationWarning)
        factors = sp.factor_list(polynomial.as_expr(), modulus=prime)[1]
    degrees: list[int] = []
    for factor, multiplicity in factors:
        degrees.extend([int(sp.degree(factor, q))] * multiplicity)
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

degree_sieve_max_prime = 0
degree_sieve_max_steps = 0
degree_sieve_pairs = 0
for m in range(1, 31):
    for r in range(1, 31):
        if m * r > 30:
            continue
        polynomial = sp.Poly(parameter_polynomial(m, r, q), q, domain=sp.ZZ)
        possible_degrees = set(range(polynomial.degree() + 1))
        steps = 0
        for prime in sp.primerange(2, 100):
            degrees = modular_factor_degrees(polynomial, int(prime))
            possible_degrees &= subset_sums(degrees)
            steps += 1
            if possible_degrees == {0, polynomial.degree()}:
                degree_sieve_max_prime = max(degree_sieve_max_prime, int(prime))
                degree_sieve_max_steps = max(degree_sieve_max_steps, steps)
                break
        assert possible_degrees == {0, polynomial.degree()}
        degree_sieve_pairs += 1

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
    f"PASS: modular degree-sieve certificates for all {degree_sieve_pairs} pairs "
    "with mr<=30 "
    f"(largest prime {degree_sieve_max_prime}, at most {degree_sieve_max_steps} steps)"
)
