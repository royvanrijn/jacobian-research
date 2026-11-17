#!/usr/bin/env python3
"""Verify the degree-five branch-scale constants in the H1/H2 obstruction."""

import sympy as sp


W, X, t = sp.symbols("W X t")


def scaled_limits(p: int, q: int):
    a = t**p
    b = 1 + t**q
    raw = W**2 * (W - 1) * (W - a) * (W - b)
    seed = sp.cancel(-raw / sp.diff(raw, W).subs(W, 1))

    zero_cluster = sp.cancel(seed.subs(W, t**p * X) / t ** (3 * p - q))
    one_cluster = sp.cancel(seed.subs(W, 1 + t**q * X) / t**q)
    return (
        sp.factor(sp.limit(zero_cluster, t, 0)),
        sp.factor(sp.limit(one_cluster, t, 0)),
    )


for valuations in ((1, 1), (2, 3)):
    zero_limit, one_limit = scaled_limits(*valuations)
    assert zero_limit == X**2 * (X - 1)
    assert one_limit == X * (X - 1)
    lambda_zero = zero_limit.subs(X, sp.Rational(2, 3))
    lambda_one = one_limit.subs(X, sp.Rational(1, 2))
    assert lambda_zero == -sp.Rational(4, 27)
    assert lambda_one == -sp.Rational(1, 4)
    assert sp.cancel(lambda_zero / lambda_one) == sp.Rational(16, 27)

print("H1_H2_COMPARISON_OBSTRUCTION_PASS")
