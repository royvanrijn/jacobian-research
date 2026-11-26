#!/usr/bin/env python3
"""Replay the formal Ore identities used in the polynomial (4,6) exclusion."""

from functools import cache

import sympy as sp

x = sp.Symbol("x")
u, v, z = (sp.Function(name)(x) for name in ("u", "v", "z"))
a, b, c, e, v0, z0 = sp.symbols("a b c e v0 z0")


@cache
def derivative(coefficient, order=1):
    return sp.diff(coefficient, x, order)


def product(left, right, cutoff=0):
    result = {}
    for i, left_coefficient in left.items():
        for j, right_coefficient in right.items():
            upper = i + j - cutoff
            if i >= 0:
                upper = min(i, upper)
            for h in range(upper + 1):
                degree = i + j - h
                contribution = (
                    sp.binomial(i, h)
                    * left_coefficient
                    * derivative(right_coefficient, h)
                )
                result[degree] = result.get(degree, 0) + contribution
    return {
        degree: expanded
        for degree, coefficient in result.items()
        if (expanded := sp.expand(coefficient)) != 0
    }


def power(operator, exponent, cutoff=0):
    result = {0: sp.Integer(1)}
    for factors_left in range(exponent - 1, -1, -1):
        result = product(result, operator, cutoff - factors_left)
    return result


def check_equal(actual, expected):
    assert sp.expand(actual - expected) == 0


def main():
    p = {
        4: sp.Integer(1),
        2: u,
        1: derivative(u) + v,
        0: u**2 / 4 + derivative(u, 2) / 2 + derivative(v) / 2 + z,
    }
    root = {1: sp.Integer(1)}
    for i in range(1, 6):
        degree = 3 - i
        current = power(root, 4, degree)
        root[-i] = sp.expand((p.get(degree, 0) - current.get(degree, 0)) / 4)
    assert power(root, 4, -2) == p

    q = {}
    for exponent, multiplier in ((6, 1), (5, a), (3, b), (2, c), (1, e)):
        for degree, coefficient in power(root, exponent).items():
            q[degree] = sp.expand(q.get(degree, 0) + multiplier * coefficient)
    pq, qp = product(p, q), product(q, p)
    r = [sp.expand(pq.get(t, 0) - qp.get(t, 0)) for t in range(11)]
    assert not any(r[3:])

    constant_v = {v: v0, **{derivative(v, j): 0 for j in range(1, 9)}}
    h = (
        (derivative(u, 2) / 2 - 3 * z) * v0
        - 2 * c * v0 - e * u
        + a * (
            -5 * u**3 / 32 - 5 * u * derivative(u, 2) / 16
            - 5 * derivative(u)**2 / 32 - 5 * u * z / 4
            - 5 * v0**2 / 8 - derivative(u, 4) / 16 - 5 * derivative(z, 2) / 4
        )
        + b * (-3 * u**2 / 8 - derivative(u, 2) / 4 - 3 * z)
    )
    check_equal(r[2].xreplace(constant_v), derivative(h))
    first_difference = (
        u * derivative(z, 3) + 3 * derivative(u) * derivative(z, 2) / 2
        + derivative(u, 2) * derivative(z) / 2 + derivative(z, 5) / 2
        - 3 * z * derivative(z) - (5 * a * v0 / 4 + 2 * c) * derivative(z)
        + 5 * a * v0 * (u * derivative(u) + derivative(u, 3)) / 16
        + 3 * v0 * (v0 + b) * derivative(u) / 4
    )
    check_equal((r[1] - derivative(r[2])).xreplace(constant_v), first_difference)

    final = {
        **constant_v,
        a: 0, b: 0, e: 0,
        z: derivative(u, 2) / 6 + z0,
        **{derivative(z, j): derivative(u, j + 2) / 6 for j in range(1, 9)},
    }
    check_equal(r[2].xreplace(final), 0)
    check_equal(
        r[1].xreplace(final),
        u * derivative(u, 5) / 6 + derivative(u) * derivative(u, 4) / 4
        + derivative(u, 7) / 12 - (c / 3 + z0 / 2) * derivative(u, 3)
        + 3 * v0**2 * derivative(u) / 4,
    )
    check_equal(
        (r[0] - derivative(r[1]) / 2).xreplace(final),
        v0 * derivative(u, 5) / 12 + v0 * (c + 3 * z0 / 2) * derivative(u) / 2,
    )
    print("Verified the (4,6) Ore coefficient identities exactly.")


if __name__ == "__main__":
    main()
