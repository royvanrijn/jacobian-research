#!/usr/bin/env python3
"""Exact algebra for the degree-nine conductor-fourteen cable exclusion."""

import sympy as sp


q = sp.Rational


def equal(left, right):
    assert sp.cancel(left - right) == 0, (left, right)


def cable_coefficients():
    z, r, s, T = sp.symbols("z r s T")
    dr, ds, dT = sp.symbols("dr ds dT")
    t = T - 2 * s**2 / r
    k = 2 * T / r - 3 * s**2 / r**2
    ell = -2 * s * T / r**2 + 2 * s**3 / r**3
    L = z**3 + k * z + ell
    M = r * z**2 + s * z + t
    alpha = r**3 / 8 - 3 * T**2 / 4
    A = r**2 * s
    X = r**2 * T - 3 * r * s**2
    U = L**2 + M
    V = L**3 + q(3, 2) * L * M + q(3, 8) * (r**2 * z + 2 * r * s)
    equal(M**2 - (r**2 * z + 2 * r * s) * L, T**2)

    reciprocal = sp.Symbol("v")
    tail = ((q(3, 8) * T**2 + alpha / 2) / L - M**3 / (16 * L**3))
    series = sp.series(tail.subs(z, 1 / reciprocal), reciprocal, 0, 6).removeO()
    equal(series, -3 * A * reciprocal**4 / 16 + X * reciprocal**5 / 16)

    def derivative(expression):
        return sum(
            expression.diff(variable) * rate
            for variable, rate in ((r, dr), (s, ds), (T, dT))
        )

    equal(
        U.diff(z) * derivative(V) - derivative(U) * V.diff(z),
        derivative(alpha) * (3 * z**2 + k) + q(9, 8) * derivative(A) * z - q(3, 8) * derivative(X),
    )


def pole_limit():
    w = sp.Symbol("w")
    L = w**3 + w
    M = q(3, 2) * (w**2 + q(1, 2))
    U = L**2 + M
    V = L**3 + q(3, 2) * L * M + q(27, 32) * w
    W = sp.expand(V**2 - U**3)
    assert sp.degree(V, w) == 9 and sp.degree(W, w) == 4
    equal(W.diff(w), -q(27, 512) * w * (16 * w**2 + 13))
    equal(V.diff(w).subs(w, 0), q(63, 32))
    equal(sp.rem(sp.cancel(U.diff(w) / (2 * w)), 16 * w**2 + 13, w), q(315, 256))


def zero_limit():
    w, mu, eps = sp.symbols("w mu eps")
    L = w**3 + (2 * mu - 3) * w + 2 * (1 - mu)
    N = w**2 + w + mu - 2
    U = L**2 + eps * N
    V = L**3 + q(3, 2) * eps * L * N + q(3, 8) * eps**2 * (w + 2)
    alpha = eps**3 / 8 - q(3, 4) * eps**2 * mu**2
    W = (N**3 - 3 * mu**2 * N - L**2) / 8 + eps * (q(9, 64) * (w + 2)**2 - N / 8)
    equal(N**2 - (w + 2) * L, mu**2)
    equal(V**2 - U**3 - alpha * U, eps**3 * W)

    base = {mu: 0, eps: 0}
    equal(W.subs(base), q(3, 8) * (w - 1)**3 * (w + 2)**2)
    assert sp.degree(V.subs(base), w) == 9
    common = sp.Poly(sp.gcd(V.subs(base).diff(w), W.subs(base).diff(w)), w).monic()
    equal(common.as_expr(), (w - 1)**2 * (w + 2))

    e, z, direction = sp.symbols("e z direction")
    zoom = {eps: e**3, w: 1 + e * z, mu: e * direction}
    scaled_u, scaled_w = (
        sp.expand(expression.subs(zoom, simultaneous=True)) for expression in (U, W)
    )
    for expression, order, limit in (
        (scaled_u, 4, 9 * z**4 + 3 * z),
        (scaled_w, 3, q(27, 8) * z**3 + q(81, 64)),
    ):
        for degree in range(order):
            equal(expression.coeff(e, degree), 0)
        equal(expression.coeff(e, order).subs(direction, 0), limit)

    equal(U.diff(w, 2).subs(base).subs(w, -2), 162)
    equal(V.diff(w, 2).subs(base).subs(w, -2), 0)
    equal(V.diff(w, 3).subs(base).subs(w, -2) / 6, 729)
    mid = -2 + 2 * mu / 3
    equal(sp.series(L.subs(w, mid), mu, 0, 2).removeO(), 0)
    equal(sp.series(N.subs(w, mid), mu, 0, 2).removeO(), -mu)
    equal(
        V.diff(w) - q(3, 2) * L * U.diff(w),
        q(3, 8) * eps * (4 * L.diff(w) * N + eps),
    )


def main():
    cable_coefficients()
    pole_limit()
    zero_limit()
    print("Verified the degree-nine cable identities and both finite-fibre limits.")


if __name__ == "__main__":
    main()
