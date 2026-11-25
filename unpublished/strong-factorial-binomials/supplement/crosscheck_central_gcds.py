#!/usr/bin/env python3
"""Recompute the exported central gcd pairs over a second prime field."""

from math import isqrt, prod
import sys

from flint import nmod_poly


PRIME = 2**31 - 1


def reciprocal(order, lower_exponent, gap, inverses):
    coefficients = [1]
    binomial = factorial_ratio = 1
    for k in range(1, order + 1):
        binomial = binomial * (order - k + 1) * inverses[k] % PRIME
        start = lower_exponent * order + gap * (k - 1)
        factorial_ratio *= prod(range(start + 1, start + gap + 1))
        factorial_ratio %= PRIME
        coefficients.append(binomial * factorial_ratio % PRIME)
    return nmod_poly(coefficients[::-1], PRIME)


def main():
    if any(PRIME % divisor == 0 for divisor in range(2, isqrt(PRIME) + 1)):
        raise ValueError("Composite modulus")
    pairs = [
        tuple(map(int, line.split()[2:]))
        for line in sys.stdin
        if line.startswith("Replay pair ")
    ]
    if len(pairs) != 64 or len(set(pairs)) != 64:
        raise ValueError("Expected 64 distinct exported pairs")
    inverses = [0] + [pow(k, -1, PRIME) for k in range(1, 3150)]
    for n, a, d in pairs:
        if not (25 <= n < 3150 and a >= 5 and 2 <= d < 3 * a and (a + d) * n < PRIME):
            raise ValueError(f"Invalid support {(n, a, d)}")
        lower = reciprocal(n - 1, a, d, inverses)
        upper = reciprocal(n, a, d, inverses)
        if lower.degree() != n - 1 or upper.degree() != n or lower.gcd(upper).degree() != 0:
            raise ValueError(f"Coprimality check failed for {(n, a, d)}")
    print(f"FLINT verified all {len(pairs)} exported pairs modulo {PRIME}.")


if __name__ == "__main__":
    main()
