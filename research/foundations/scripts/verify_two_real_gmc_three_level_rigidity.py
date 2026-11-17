#!/usr/bin/env python3
"""Exact arithmetic regression for all-degree three-level rigidity.

The proof is mathematical: at an odd prime p, the Bessel coefficients of
H_p have only the C^p endpoint modulo p, while those of H_(2p) have only the
C^(2p) and binomial(2p,p)D^p endpoints.  Comparing the U-adic orders of C
and D then isolates one endpoint after applying L(U^j)=j!.

This script checks the coefficient divisibility and the resulting normalized
factorial congruences over a representative range of primes and orders.
"""

from __future__ import annotations

import math

import sympy as sp


def bessel_coefficient(moment: int, circuit_power: int) -> int:
    return (
        math.factorial(moment)
        // math.factorial(circuit_power) ** 2
        // math.factorial(moment - 2 * circuit_power)
    )


def factorial_functional(expression: sp.Expr, radial: sp.Symbol) -> int:
    polynomial = sp.Poly(sp.expand(expression), radial)
    return sum(
        int(coefficient) * math.factorial(exponent[0])
        for exponent, coefficient in polynomial.terms()
    )


def bessel_polynomial(
    moment: int, centered: sp.Expr, circuit: sp.Expr
) -> sp.Expr:
    return sp.expand(
        sum(
            bessel_coefficient(moment, power)
            * centered ** (moment - 2 * power)
            * circuit**power
            for power in range(moment // 2 + 1)
        )
    )


def main() -> None:
    radial = sp.symbols("U")
    primes = (3, 5, 7, 11, 13)

    for prime in primes:
        assert all(
            bessel_coefficient(prime, power) % prime == 0
            for power in range(1, prime // 2 + 1)
        )
        assert all(
            bessel_coefficient(2 * prime, power) % prime == 0
            for power in range(1, prime)
        )
        assert bessel_coefficient(2 * prime, prime) % prime == 2

        for centered_order in range(4):
            centered = radial**centered_order * (
                1 + 2 * radial + 3 * radial**2
            )
            for circuit_order in range(5):
                circuit = radial**circuit_order * (
                    1 + 5 * radial + 7 * radial**2
                )

                if circuit_order >= 2 * centered_order:
                    moment = prime
                    baseline = centered_order * prime
                    expected = 1
                else:
                    moment = 2 * prime
                    baseline = circuit_order * prime
                    expected = 2

                value = factorial_functional(
                    bessel_polynomial(moment, centered, circuit),
                    radial,
                )
                baseline_factorial = math.factorial(baseline)
                assert value % baseline_factorial == 0
                assert value // baseline_factorial % prime == expected

    print("PASS three-level rigidity: prime Bessel endpoint coefficients")
    print("PASS three-level rigidity: both U-adic order cases")


if __name__ == "__main__":
    main()
