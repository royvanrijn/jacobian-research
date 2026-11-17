#!/usr/bin/env python3
"""Independent pure-Python audit of prime-endpoint rigidity.

This implementation shares no symbolic-algebra code with
verify_two_real_gmc_three_level_rigidity.py.  It checks:

1. the two prime-endpoint coefficient congruences;
2. the polynomial identities H_p=C^p and H_(2p)=C^(2p)+2D^p mod p;
3. the normalized factorial endpoint in both U-adic order cases.
"""

from __future__ import annotations

from math import factorial


Polynomial = dict[int, int]


def add(left: Polynomial, right: Polynomial) -> Polynomial:
    result = dict(left)
    for degree, coefficient in right.items():
        result[degree] = result.get(degree, 0) + coefficient
        if result[degree] == 0:
            del result[degree]
    return result


def scale(polynomial: Polynomial, scalar: int) -> Polynomial:
    return {
        degree: scalar * coefficient
        for degree, coefficient in polynomial.items()
        if scalar * coefficient
    }


def multiply(left: Polynomial, right: Polynomial) -> Polynomial:
    result: Polynomial = {}
    for left_degree, left_coefficient in left.items():
        for right_degree, right_coefficient in right.items():
            degree = left_degree + right_degree
            result[degree] = (
                result.get(degree, 0)
                + left_coefficient * right_coefficient
            )
    return {degree: value for degree, value in result.items() if value}


def power(polynomial: Polynomial, exponent: int) -> Polynomial:
    result: Polynomial = {0: 1}
    base = polynomial
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = multiply(result, base)
        base = multiply(base, base)
        remaining //= 2
    return result


def coefficient(moment: int, circuit_power: int) -> int:
    return (
        factorial(moment)
        // factorial(circuit_power) ** 2
        // factorial(moment - 2 * circuit_power)
    )


def transform(moment: int, centered: Polynomial, circuit: Polynomial) -> Polynomial:
    result: Polynomial = {}
    for circuit_power in range(moment // 2 + 1):
        term = multiply(
            power(centered, moment - 2 * circuit_power),
            power(circuit, circuit_power),
        )
        result = add(
            result,
            scale(term, coefficient(moment, circuit_power)),
        )
    return result


def reduce_mod(polynomial: Polynomial, prime: int) -> Polynomial:
    return {
        degree: coefficient % prime
        for degree, coefficient in polynomial.items()
        if coefficient % prime
    }


def factorial_functional(polynomial: Polynomial) -> int:
    return sum(
        coefficient * factorial(degree)
        for degree, coefficient in polynomial.items()
    )


def main() -> None:
    primes = (3, 5, 7, 11)

    for prime in primes:
        assert all(
            coefficient(prime, circuit_power) % prime == 0
            for circuit_power in range(1, prime // 2 + 1)
        )
        assert all(
            coefficient(2 * prime, circuit_power) % prime == 0
            for circuit_power in range(1, prime)
        )
        assert coefficient(2 * prime, prime) % prime == 2

        centered_mod = {0: 2, 1: 3, 3: 1}
        circuit_mod = {1: 4, 2: 1, 4: 2}
        assert reduce_mod(
            transform(prime, centered_mod, circuit_mod), prime
        ) == reduce_mod(power(centered_mod, prime), prime)
        assert reduce_mod(
            transform(2 * prime, centered_mod, circuit_mod), prime
        ) == reduce_mod(
            add(
                power(centered_mod, 2 * prime),
                scale(power(circuit_mod, prime), 2),
            ),
            prime,
        )

        for centered_order in range(4):
            centered = {
                centered_order: 1,
                centered_order + 1: 2,
                centered_order + 3: 3,
            }
            for circuit_order in range(5):
                circuit = {
                    circuit_order: 1,
                    circuit_order + 2: 5,
                    circuit_order + 3: 7,
                }
                if circuit_order >= 2 * centered_order:
                    moment = prime
                    baseline = centered_order * prime
                    expected = 1
                else:
                    moment = 2 * prime
                    baseline = circuit_order * prime
                    expected = 2

                value = factorial_functional(
                    transform(moment, centered, circuit)
                )
                divisor = factorial(baseline)
                assert value % divisor == 0
                assert value // divisor % prime == expected

    print("PASS independent audit: prime-endpoint lemma")
    print("PASS independent audit: normalized factorial endpoints")


if __name__ == "__main__":
    main()
