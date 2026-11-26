"""Exact Legendre cross-check of the split Gauss–Manin jet thresholds."""

from fractions import Fraction


def multiply_polynomial(left, right):
    result = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] += a * b
    return result


def coefficient(series, index):
    if 0 <= index < len(series):
        return series[index]
    return 0


def legendre_jet(order):
    # The parameter is lambda = t - 1 and L = lambda(1-lambda).
    # In the frame (omega, nabla omega), the Picard–Fuchs equation
    # gives u = -w' + (L'/L)w and
    # 4 L^2 w'' - 4 L L' w' + (7 L + 4 (L')^2)w = 0.
    leading = [-2, 3, -1]
    derivative = [3, -2]
    second_order = [4 * x for x in multiply_polynomial(leading, leading)]
    first_order = [-4 * x for x in multiply_polynomial(leading, derivative)]
    derivative_square = multiply_polynomial(derivative, derivative)
    zeroth_order = [7 * x + 4 * y for x, y in zip(leading, derivative_square)]

    second = [Fraction(0), Fraction(-1)]
    for i in range(order):
        known = sum(
            second_order[k] * (i - k + 2) * (i - k + 1) * coefficient(second, i - k + 2)
            for k in range(1, len(second_order))
        )
        known += sum(
            a * (i - k + 1) * coefficient(second, i - k + 1)
            for k, a in enumerate(first_order)
        )
        known += sum(a * coefficient(second, i - k) for k, a in enumerate(zeroth_order))
        second.append(-known / (second_order[0] * (i + 2) * (i + 1)))

    first = []
    quotient = []
    for i in range(order + 1):
        numerator = -(i + 1) * leading[0] * second[i + 1]
        numerator -= sum(
            leading[k] * (i - k + 1) * coefficient(second, i - k + 1)
            for k in (1, 2)
        )
        numerator += sum(a * coefficient(second, i - k) for k, a in enumerate(derivative))
        numerator -= sum(leading[k] * coefficient(first, i - k) for k in (1, 2))
        first.append(numerator / leading[0])
        quotient.append(second[i] - sum(first[k] * quotient[i - k] for k in range(1, i + 1)))
    assert first[0] == 1 and second[0] == 0
    return first, second, quotient


def main():
    if not __debug__:
        raise SystemExit("The verifier requires active Python assertions.")
    primes = (3, 5, 7, 11, 13)
    order = max(2 * p**2 + p for p in primes if p % 4 == 3)
    first, second, quotient = legendre_jet(order)
    for p in primes:
        vector_order = next(
            i for i, (u, w) in enumerate(zip(first, second))
            if u.denominator % p == 0 or w.denominator % p == 0
        )
        projective_order = next(i for i, y in enumerate(quotient) if y.denominator % p == 0)
        expected = (p, p) if p % 4 == 1 else (p**2 + p, 2 * p**2 + p)
        assert (vector_order, projective_order) == expected
        assert first[vector_order].denominator % (p * p)
        assert second[vector_order].denominator % (p * p)
        assert quotient[projective_order].denominator % (p * p)
        print(f"p={p}: vector at {vector_order}, projective coordinate at {projective_order}")


if __name__ == "__main__":
    main()
