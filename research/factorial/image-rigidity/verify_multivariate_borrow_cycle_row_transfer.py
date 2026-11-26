"""Exact rational controls for borrowed exponent rows and quotient cycles."""

from fractions import Fraction
from math import factorial, gcd, prod


def unit_congruent_to_one(value, modulus):
    return (value.numerator - value.denominator) % modulus == 0 and gcd(value.denominator, modulus) == 1


def unit_factorial(n, prime):
    return prod(k for k in range(1, n + 1) if k % prime)


def inflation_scale(exponent, prime):
    return prod(prime ** (prime ** (h - 1)) * unit_factorial(prime**h, prime)
                for h in range(1, exponent + 1))


def dot(row, counts):
    return sum(entry * count for entry, count in zip(row, counts, strict=True))


def normalize_borrow(prime, A, B, counts, e, s):
    x = dot(A, counts)
    y = dot(B, counts)
    diagonal = Fraction(inflation_scale(e, prime), inflation_scale(s, prime)) ** x
    scalar = Fraction((-1) ** (y - 1), prime ** (e - s))
    actual = Fraction(factorial(prime**e * x - y), 1) / (diagonal * scalar)
    virtual = Fraction(factorial(prime**s * x - 1), factorial(y - 1))
    return actual, virtual


def main():
    local_controls = 0
    for prime, depth, e, s in ((2, 3, 7, 2), (3, 3, 5, 2), (5, 2, 4, 1)):
        for counts in ((1, 0, 0), (0, 1, 0), (0, 0, 1), (2, 1, 0), (1, 1, 1)):
            A, B = (1, 2, 3), (5, 2, 4)
            actual, virtual = normalize_borrow(prime, A, B, counts, e, s)
            assert unit_congruent_to_one(actual / virtual, prime**depth)
            local_controls += 1

    prime, depth = 2, 3
    vertices = ((1, 3), (4, 2), (2, 5))
    scales = (2, 3, 2)
    cycle = ((0, 1, 10), (1, 2, 9), (2, 0, 8))
    graph_controls = 0
    for counts in ((1, 0), (0, 1), (1, 1), (2, 1), (1, 2)):
        actuals = []
        virtuals = []
        for source, target, exponent in cycle:
            B = tuple(prime ** scales[target] * value for value in vertices[target])
            actual, virtual = normalize_borrow(
                prime, vertices[source], B, counts, exponent, scales[source]
            )
            actuals.append(actual)
            virtuals.append(virtual)
        assert prod(virtuals) == 1
        assert unit_congruent_to_one(prod(actuals), prime**depth)
        graph_controls += 1

        start = factorial(prime ** scales[0] * dot(vertices[0], counts) - 1)
        end = factorial(prime ** scales[2] * dot(vertices[2], counts) - 1)
        path = Fraction(start, end)
        assert prod(virtuals[:2]) == path
        assert unit_congruent_to_one(prod(actuals[:2]) / path, prime**depth)
        graph_controls += 1

    actual, virtual = normalize_borrow(2, (1, 2), (1, 2), (0, 1), 2, 0)
    assert not unit_congruent_to_one(actual / virtual, 8)
    print(f"{local_controls} local units, {graph_controls} graph controls, one rejected depth")


if __name__ == "__main__":
    main()
