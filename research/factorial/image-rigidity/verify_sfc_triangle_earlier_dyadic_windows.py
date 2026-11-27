"""Exact arithmetic checks for the quadratic triangle's earlier dyadic windows."""

from sage.all import GF, PolynomialRing, matrix


def check_factorial_layers():
    ring = PolynomialRing(GF(2), names=("a", "b", "c"))
    a, b, c = ring.gens()
    s, r = a + b + c, a * b * c
    complements = (b * c, c * a, a * b)
    for k in range(3, 9):
        Q = 2**k
        digits = [value.bit_count() for value in range(Q)]
        fifth = s**(Q - 5) + sum(
            s**(Q - 2**(j + 1) + 3) * sum(x**(2**j - 4) for x in complements)
            for j in range(3, k)
        )
        cases = [(Q - 5, 0, fifth)]
        if k >= 4:
            sixth = r**3 * s**(Q - 15) + sum(
                s**(Q - 2**(j + 1) + 1) * sum(
                    coefficient**3 * x**(2**j - 5)
                    for coefficient, x in zip((a, b, c), complements)
                )
                for j in range(4, k)
            )
            cases.append((Q - 6, -1, sixth))
        for order, minimum, predicted in cases:
            leading = {}
            for i in range(order + 1):
                for j in range(order - i + 1):
                    remaining = order - i - j
                    valuation = digits[order] + digits[i] + digits[j] + digits[remaining]
                    valuation -= digits[i + j] + digits[j + remaining] + digits[remaining + i]
                    assert valuation >= minimum, (order, i, j, remaining)
                    if valuation == minimum:
                        leading[(i, j, remaining)] = 1
            assert ring(leading) == predicted, order


def check_residue_eliminations():
    parameters = PolynomialRing(GF(2), names=("h", "r"))
    h, r = parameters.gens()
    ring = PolynomialRing(parameters.fraction_field(), "A")
    A = ring.gen()
    f = A**3 + A**2 + h*A + r
    x = A**2 + A + h

    def trace(form):
        remainder = ring(form) % f
        return sum(remainder[j] for j in range(3))

    targets = (r**4 * (h**2 + r**2 + h + 1), r**5 * (r + 1))
    for exponent, target in zip((4, 5), targets):
        remainder = (A**exponent * (1 + A)) % f
        W, V, U = (remainder[j] for j in range(3))
        D = U + V + W
        product_difference = (
            U**2*x**2 + (U*V + V**2 + D**2)*x + U*V*r
            + U*W*(1 + A**2) + V*W*(1 + A) + W**2
        )
        required = r**2*(1 + h)**2 + r*trace(A**7) + r**2*trace(A**6) + r**4
        assert trace(A**8 * product_difference) + D**2*required == target, exponent


def check_frobenius_table():
    binary_polynomials = PolynomialRing(GF(2), "t")
    t = binary_polynomials.gen()
    parameters = binary_polynomials.fraction_field()
    ring = PolynomialRing(parameters, "A")
    A = ring.gen()
    h, r = t**2 + 1, t**2 + t + 1
    f = A**3 + A**2 + h*A + r
    image = (t*A**2 + h**2*A + h*r) / t**3
    columns = [((image * A**j) % f) for j in range(3)]
    multiplication = matrix(parameters, 3, 3, lambda i, j: columns[j][i])
    coefficients = multiplication.charpoly().list()
    assert coefficients == [(t**8 + t**4 + 1) / t**8, t**4 + 1, 1, 1]

    field = GF(8, "t", modulus=t**3 + t**2 + 1)
    t = field.gen()
    ring = PolynomialRing(field, "A")
    A = ring.gen()
    f = A**3 + A**2 + (t**2 + 1)*A + t**2 + t + 1
    expected = {
        1: A**2,
        4: t*A**2 + (t**2 + 1)*A + t**2 + t,
        7: (t + 1)*A**2 + (t**2 + 1)*A + t**2 + t + 1,
        9: A,
    }
    power = A
    for exponent in range(1, 10):
        power = power**2 % f
        if exponent in expected:
            assert power == expected[exponent], exponent

    power = A**2 + A + t**2 + 1
    running = ring.zero()
    for j in range(9):
        running += power
        if j == 6:
            assert (A**4 * running) % f == A**2 + (t**2 + 1)*A + t + 1
        power = power**2 % f
    assert running == 0


def check_infinite_subsequence():
    period = 16895838960
    primes_and_valuations = [(7, 1), (3, 2), (197, 1), (13, 1),
                             (3, 1), (17, 1), (11, 1), (3, 1)]
    for shift, (prime, valuation) in zip(range(9, 24, 2), primes_and_valuations):
        value = 2**25 - shift
        prime_power = prime**valuation
        modulus = prime_power * prime
        assert value > prime_power and value % prime_power == 0 and value % modulus != 0
        assert pow(2, period, modulus) == 1


if __name__ == "__main__":
    check_factorial_layers()
    check_residue_eliminations()
    check_frobenius_table()
    check_infinite_subsequence()
    print("PASS earlier dyadic layers, two eliminations, Frobenius table and infinite subsequence")
