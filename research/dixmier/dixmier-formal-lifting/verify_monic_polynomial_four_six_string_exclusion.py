#!/usr/bin/env python3
"""Replay the finite Ore coefficients for the last monic (4,6) weight."""

import sympy as sp

from verify_order_46_constant_skew import product, power


def main():
    x = sp.Symbol("x")
    U, V, Z, u0, v0, z0, a, b, c, e = sp.symbols("U V Z u0 v0 z0 a b c e")
    u, v, z = U * x + u0, V * x + v0, Z * x + z0
    p = {4: sp.Integer(1), 2: u, 1: U + v, 0: u**2 / 4 + V / 2 + z}
    root = {1: sp.Integer(1)}
    for i in range(1, 6):
        degree = 3 - i
        current = power(root, 4, degree)
        root[-i] = sp.expand((p.get(degree, 0) - current.get(degree, 0)) / 4)
    fourth_power = power(root, 4, -2)
    assert all(
        sp.expand(fourth_power.get(j, 0) - p.get(j, 0)) == 0
        for j in range(-2, 5)
    )

    q = {}
    for exponent, multiplier in ((6, 1), (5, a), (3, b), (2, c), (1, e)):
        for degree, coefficient in power(root, exponent).items():
            q[degree] = sp.expand(q.get(degree, 0) + multiplier * coefficient)
    pq, qp = product(p, q), product(q, p)
    degree_two = sp.expand(pq.get(2, 0) - qp.get(2, 0)).coeff(x, 2)
    degree_one = sp.expand(pq.get(1, 0) - qp.get(1, 0)).coeff(x, 2)
    assert degree_two == -sp.Rational(15, 32) * a * U**3
    assert degree_one.subs(a, 0) == sp.Rational(9, 4) * U * V**2
    print("Verified the exact Ore coefficients for weight (1, 1/2).")


if __name__ == "__main__":
    main()
