#!/usr/bin/env python3
"""Exact regressions for finite-character digit separation.

For a Laurent polynomial ``f`` and fixed base-p digits ``n_j``, put

    N = n_0 + n_1*p + ... + n_s*p**s.

If ``p`` is larger than every exponent that can occur in the small
blocks ``f**n_j``, the freshman's-dream identity and uniqueness of
balanced signed base-p digits give

    CT(f**N) = product_j CT(f**n_j)**(p**j)  (mod p).

At a completely split prime the Frobenius powers disappear.  Applying
this to every character component of a torsion--torus trace and using
repeated equal digits produces all power sums of the component moments.
Newton identities then force every component moment to vanish.

The script checks the digit factorization directly in one and two free
Laurent variables, checks its compatibility with finite trace sums, and
replays the Newton-identity endpoint.  The characteristic-zero proof
also uses specialization to a number field and infinitely many completely
split good primes; those standard arithmetic steps are documented in the
canonical note and are not replaced by this finite regression.
"""

from __future__ import annotations

import itertools
from collections.abc import Iterable

import sympy as sp


Exponent = tuple[int, ...]
Laurent = dict[Exponent, int]


def add_term(
    polynomial: Laurent,
    exponent: Exponent,
    coefficient: int,
    modulus: int | None = None,
) -> None:
    value = polynomial.get(exponent, 0) + coefficient
    if modulus is not None:
        value %= modulus
    if value:
        polynomial[exponent] = value
    elif exponent in polynomial:
        del polynomial[exponent]


def multiply(
    left: Laurent,
    right: Laurent,
    modulus: int | None = None,
) -> Laurent:
    answer: Laurent = {}
    for left_exponent, left_value in left.items():
        for right_exponent, right_value in right.items():
            exponent = tuple(
                left_coordinate + right_coordinate
                for left_coordinate, right_coordinate in zip(
                    left_exponent,
                    right_exponent,
                    strict=True,
                )
            )
            add_term(
                answer,
                exponent,
                left_value * right_value,
                modulus,
            )
    return answer


def power(
    polynomial: Laurent,
    exponent: int,
    modulus: int | None = None,
) -> Laurent:
    dimension = len(next(iter(polynomial)))
    answer: Laurent = {(0,) * dimension: 1}
    base = polynomial
    remaining = exponent
    while remaining:
        if remaining & 1:
            answer = multiply(answer, base, modulus)
        remaining //= 2
        if remaining:
            base = multiply(base, base, modulus)
    return answer


def constant_term(polynomial: Laurent) -> int:
    dimension = len(next(iter(polynomial)))
    return polynomial.get((0,) * dimension, 0)


def digit_index(digits: Iterable[int], prime: int) -> int:
    return sum(
        digit * prime**position
        for position, digit in enumerate(digits)
    )


def block_bound(
    polynomial: Laurent,
    digits: Iterable[int],
) -> int:
    largest_coordinate = max(
        abs(coordinate)
        for exponent in polynomial
        for coordinate in exponent
    )
    return largest_coordinate * max(digits)


def verify_digit_factorization(
    polynomial: Laurent,
    digits: tuple[int, ...],
    prime: int,
) -> tuple[int, int]:
    assert prime > block_bound(polynomial, digits)
    moment = digit_index(digits, prime)
    direct = constant_term(power(polynomial, moment, prime)) % prime
    factors = tuple(
        constant_term(power(polynomial, digit, prime)) % prime
        for digit in digits
    )
    predicted = 1
    for position, factor in enumerate(factors):
        predicted *= pow(factor, prime**position, prime)
        predicted %= prime
    assert direct == predicted
    return direct, predicted


def verify_trace_factorization(
    polynomials: tuple[Laurent, ...],
    digits: tuple[int, ...],
    prime: int,
) -> tuple[int, int]:
    direct = 0
    predicted = 0
    moment = digit_index(digits, prime)
    for polynomial in polynomials:
        assert prime > block_bound(polynomial, digits)
        direct += constant_term(power(polynomial, moment, prime))
        factors = (
            constant_term(power(polynomial, digit, prime)) % prime
            for digit in digits
        )
        product = 1
        for position, factor in enumerate(factors):
            product *= pow(factor, prime**position, prime)
            product %= prime
        predicted += product
    direct %= prime
    predicted %= prime
    assert direct == predicted
    return direct, predicted


def verify_newton_endpoint(component_count: int) -> None:
    values = sp.symbols(f"x0:{component_count}")
    power_sums = {
        degree: sum(value**degree for value in values)
        for degree in range(1, component_count + 1)
    }
    elementary = {0: sp.Integer(1)}
    for degree in range(1, component_count + 1):
        numerator = sum(
            (-1) ** (index - 1)
            * elementary[degree - index]
            * power_sums[index]
            for index in range(1, degree + 1)
        )
        elementary[degree] = sp.expand(numerator / degree)
        expected = sum(
            sp.prod(chosen)
            for chosen in itertools.combinations(
                values,
                degree,
            )
        )
        assert sp.expand(elementary[degree] - expected) == 0

    # After imposing the first q power sums, Newton's recursion makes
    # every elementary symmetric function zero, so the monic polynomial
    # with the component moments as roots is T**q.
    zero_substitution = {
        power_sums[degree]: 0
        for degree in range(1, component_count + 1)
    }
    reduced_elementary = {0: sp.Integer(1)}
    for degree in range(1, component_count + 1):
        numerator = sum(
            (-1) ** (index - 1)
            * reduced_elementary[degree - index]
            * zero_substitution[power_sums[index]]
            for index in range(1, degree + 1)
        )
        reduced_elementary[degree] = sp.expand(numerator / degree)
        assert reduced_elementary[degree] == 0


def verify() -> None:
    one_variable: Laurent = {
        (-1,): 2,
        (0,): 3,
        (1,): -1,
    }
    two_variables: Laurent = {
        (-1, 1): 1,
        (-1, -1): -2,
        (0, 0): 3,
        (1, -1): 1,
        (1, 1): 2,
    }
    companion: Laurent = {
        (-1, 0): -1,
        (0, -1): 2,
        (0, 0): -2,
        (1, 1): 1,
        (1, 0): 1,
    }

    prime = 5
    digit_rows = (
        verify_digit_factorization(
            one_variable,
            (2, 1, 1),
            prime,
        ),
        verify_digit_factorization(
            two_variables,
            (1, 2, 1),
            prime,
        ),
    )
    trace_row = verify_trace_factorization(
        (two_variables, companion),
        (2, 1, 2),
        prime,
    )
    verify_newton_endpoint(2)
    verify_newton_endpoint(3)
    verify_newton_endpoint(5)

    print(f"digit rows mod {prime}: {digit_rows}")
    print(f"two-component trace row mod {prime}: {trace_row}")
    print(
        "PASS signed base-p digit uniqueness and finite-trace "
        "factorization"
    )
    print(
        "PASS Newton endpoint for 2, 3, and 5 character components"
    )
    print(
        "STATUS: exact regression for the characteristic-zero "
        "finite-character separation theorem"
    )


if __name__ == "__main__":
    verify()
