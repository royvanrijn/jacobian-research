"""Exact modular replay of the opposite-height specialization tests."""

from math import comb, factorial, prod

from sympy import Matrix, Poly, Symbol


Z = Symbol("Z")
Y = Symbol("Y")


def core_coefficients(rows, order):
    return [
        comb(order, k) * prod(factorial((order - k) * a + k * b) for a, b in rows)
        for k in range(order + 1)
    ]


def norm_polynomial(poly, value, prime):
    """Characteristic polynomial of multiplication by value modulo poly."""
    deg = poly.degree()
    columns = []
    coordinate = Poly(Z, Z, modulus=prime)
    current = value.rem(poly)
    for _ in range(deg):
        columns.append([int(current.nth(i)) % prime for i in range(deg)])
        current = (current * coordinate).rem(poly)
    norm = Matrix(columns).T.charpoly(Y).as_expr()
    return Poly(norm, Y, modulus=prime)


def check_specializations(rows, n, upper_order, r, s, prime):
    low = core_coefficients(rows, n)
    high = core_coefficients(rows, upper_order)
    a0 = -pow(upper_order, -1, prime)
    for h in range(1, upper_order + 1):
        b0 = h * a0
        least = int(h == upper_order)
        lower_values = []
        for k in range(least, n + 1):
            increase = prod((n * a0 + i) % prime for i in range(1, k + 1))
            decrease = prod((n * b0 + i) % prime for i in range(1, n - k + 1))
            lower_values.append(low[k] * pow(increase, r, prime) * pow(decrease, s, prime))
        p_h = Poly.from_list(lower_values[::-1], Z, modulus=prime)
        assert p_h.degree() == n - least and p_h.nth(0) != 0
        assert p_h.gcd(p_h.diff()).degree() == 0, h

        first_index = upper_order - h + 1
        high_values = [0] * first_index
        for k in range(first_index, upper_order + 1):
            increase = upper_order * factorial(k - 1)
            decrease = prod(i - h for i in range(1, upper_order - k + 1))
            high_values.append(high[k] * pow(increase, r, prime) * pow(decrease, s, prime))
        s_h = Poly.from_list(high_values[::-1], Z, modulus=prime)
        norm = norm_polynomial(p_h.monic(), s_h, prime)
        assert norm.degree() == p_h.degree()
        assert norm.gcd(norm.diff()).degree() == 0, h


if __name__ == "__main__":
    check_specializations(((2, 5), (7, 3)), 14, 15, 2, 2, 101)
    print("All 15 opposite-height specialization tests certified modulo 101.")
