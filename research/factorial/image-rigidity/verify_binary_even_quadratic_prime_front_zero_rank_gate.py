"""Replay the finite prime-front gcd table. The valuation reductions are written proofs."""

from math import comb, factorial

from sympy import Poly, Symbol, primerange


X = Symbol("X")


def factorial_moment(order, prime):
    coefficients = [comb(order, i) * factorial(2 * i) % prime for i in range(order + 1)]
    return Poly.from_list(coefficients[::-1], X, modulus=prime)


def gate_gcd(prime, position):
    h = (prime - 1) // 2
    upper = factorial_moment(h + 1, prime)
    lower = factorial_moment(h, prime)
    x = Poly(X, X, modulus=prime)
    if position == "second":
        left = lower
        right = 3 * upper * (1 + 2 * x - 4 * upper**2) - 2
    else:
        w = x - 2 * upper**2
        left = 15 * upper * (1 + 4 * w + 24 * w**2) - 8
        right = -12 * upper**3 + 3 * upper + 6 * x * upper - 2 + 6 * x * lower
    return left.gcd(right).monic()


if __name__ == "__main__":
    for position, cutoff, expected_count, expected in (
        ("second", 550, 98, {19: [1, 2]}),
        ("first", 200, 43, {11: [1, 6], 37: [1, 9]}),
    ):
        primes = list(primerange(7, cutoff))
        survivors = {}
        for prime in primes:
            gcd = gate_gcd(prime, position)
            if gcd.degree() > 0:
                survivors[prime] = [int(value) % prime for value in gcd.all_coeffs()]
        assert len(primes) == expected_count
        assert survivors == expected
        print(f"Prime in {position} position: {len(primes)} exact gcds, surviving factors {survivors}.")
