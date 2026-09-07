#!/usr/bin/env python3
"""Standard-library replay of the weighted-seed Gaussian bridge.

This checker deliberately does not import SymPy or project modules.  It
reconstructs the polynomial correction from sparse Fraction arithmetic,
checks its formal fixed branch and determinant, and evaluates bounded circular
Gaussian moments independently of the primary symbolic implementation.
"""

from __future__ import annotations

from fractions import Fraction as Q
from math import factorial


BPoly = dict[tuple[int, int], Q]
Series = dict[int, Q]


def clean(poly):
    return {exponent: coefficient for exponent, coefficient in poly.items() if coefficient}


def add(*polys):
    out = {}
    for poly in polys:
        for exponent, coefficient in poly.items():
            out[exponent] = out.get(exponent, Q(0)) + coefficient
    return clean(out)


def scale(poly, scalar):
    return clean({exponent: Q(scalar) * coefficient for exponent, coefficient in poly.items()})


def multiply(left, right, cutoff=None):
    out = {}
    for (az, ay), ac in left.items():
        for (bz, by), bc in right.items():
            exponent = az + bz, ay + by
            if cutoff is not None and (exponent[0] > cutoff or exponent[1] > cutoff):
                continue
            out[exponent] = out.get(exponent, Q(0)) + ac * bc
    return clean(out)


def power(poly, exponent, cutoff=None):
    out = {(0, 0): Q(1)}
    base = poly
    remaining = exponent
    while remaining:
        if remaining & 1:
            out = multiply(out, base, cutoff)
        base = multiply(base, base, cutoff)
        remaining //= 2
    return out


def derivative(poly, variable):
    out = {}
    for (z_degree, y_degree), coefficient in poly.items():
        degrees = [z_degree, y_degree]
        if degrees[variable] == 0:
            continue
        factor = degrees[variable]
        degrees[variable] -= 1
        out[tuple(degrees)] = coefficient * factor
    return clean(out)


def from_univariate(poly):
    return {(degree, 0): coefficient for degree, coefficient in poly.items() if coefficient}


def bridge(seed):
    h = add({(0, 0): Q(1)}, from_univariate(seed))
    hp = derivative(h, 0)
    zhp = multiply({(1, 0): Q(1)}, hp)
    D = add(h, scale(zhp, -1))
    D_minus_one = add(D, {(0, 0): Q(-1)})
    assert all(z_degree >= 1 for z_degree, _ in D_minus_one)
    R = {(z_degree - 1, y_degree): coefficient for (z_degree, y_degree), coefficient in D_minus_one.items()}
    one_plus_y = {(0, 0): Q(1), (0, 1): Q(1)}
    L = add(multiply(D, {(0, 1): Q(1)}), multiply({(1, 0): Q(1)}, R))
    B = add(
        multiply(multiply(h, R), one_plus_y),
        scale(multiply(multiply(h, hp), power(one_plus_y, 2)), -1),
    )
    phi2 = add(
        scale(multiply(multiply(h, R), one_plus_y), -1),
        multiply(L, B),
    )
    return h, D, R, L, B, h, phi2


def series_add(*series):
    out = {}
    for value in series:
        for degree, coefficient in value.items():
            out[degree] = out.get(degree, Q(0)) + coefficient
    return clean(out)


def series_multiply(left, right, cutoff):
    out = {}
    for a, ac in left.items():
        for b, bc in right.items():
            if a + b <= cutoff:
                out[a + b] = out.get(a + b, Q(0)) + ac * bc
    return clean(out)


def series_power(value, exponent, cutoff):
    out = {0: Q(1)}
    base = value
    remaining = exponent
    while remaining:
        if remaining & 1:
            out = series_multiply(out, base, cutoff)
        base = series_multiply(base, base, cutoff)
        remaining //= 2
    return out


def series_scale(value, scalar):
    return clean({degree: Q(scalar) * coefficient for degree, coefficient in value.items()})


def series_shift(value, amount, cutoff):
    return {degree + amount: coefficient for degree, coefficient in value.items() if degree + amount <= cutoff}


def evaluate(poly, z_value, y_value, cutoff):
    out = {}
    for (z_degree, y_degree), coefficient in poly.items():
        term = series_multiply(
            series_power(z_value, z_degree, cutoff),
            series_power(y_value, y_degree, cutoff),
            cutoff,
        )
        out = series_add(out, series_scale(term, coefficient))
    return out


def inverse_series(value, cutoff):
    constant = value.get(0, Q(0))
    assert constant
    out = {0: 1 / constant}
    for degree in range(1, cutoff + 1):
        convolution = sum(
            value.get(index, Q(0)) * out.get(degree - index, Q(0))
            for index in range(1, degree + 1)
        )
        out[degree] = -convolution / constant
    return clean(out)


def moment(phi1, phi2, observable, m):
    normalized = Q(0)
    for a in range(m + 1):
        b = m - a
        product = multiply(power(phi1, a, m), power(phi2, b, m), m)
        product = multiply(product, observable, m)
        normalized += product.get((a, b), Q(0))
    return factorial(m) * normalized


def main():
    seeds = {
        "canonical-cubic": {2: Q(1), 3: Q(-1)},
        "canonical-quartic": {3: Q(1), 4: Q(-1)},
        "split-quartic": {2: Q(3), 3: Q(-4), 4: Q(1)},
    }
    cutoff = 14
    for name, seed in seeds.items():
        h, D, R, L, B, phi1, phi2 = bridge(seed)

        # Solve g=u*h(g) by fixed iteration and set k=D(g)^(-1)-1.
        g = {}
        zero = {}
        for _ in range(cutoff + 1):
            g = series_shift(evaluate(h, g, zero, cutoff), 1, cutoff)
        D_of_g = evaluate(D, g, zero, cutoff)
        k = series_add(inverse_series(D_of_g, cutoff), {0: Q(-1)})

        fixed1 = series_add(series_shift(evaluate(phi1, g, k, cutoff), 1, cutoff), series_scale(g, -1))
        fixed2 = series_add(series_shift(evaluate(phi2, g, k, cutoff), 1, cutoff), series_scale(k, -1))
        assert not fixed1 and not fixed2

        # The triangular determinant is (1-u*h'(g))
        # (1-u*partial_y(phi2)(g,k)).
        first = series_add({0: Q(1)}, series_scale(series_shift(evaluate(derivative(phi1, 0), g, k, cutoff), 1, cutoff), -1))
        second = series_add({0: Q(1)}, series_scale(series_shift(evaluate(derivative(phi2, 1), g, k, cutoff), 1, cutoff), -1))
        determinant = series_multiply(first, second, cutoff)
        assert determinant == {0: Q(1)}

        for m in range(1, 11):
            assert moment(phi1, phi2, {(0, 0): Q(1)}, m) == 0
            mixed = moment(phi1, phi2, {(1, 0): Q(1)}, m)
            h_power = power(h, m, m)
            expected = factorial(m - 1) * h_power.get((m - 1, 0), Q(0))
            assert mixed == expected

        print(f"PASS (stdlib) weighted Gaussian bridge: {name} through order {cutoff}")

    print("PASS (stdlib) weighted Gaussian bridge: independent construction and Wick replay")


if __name__ == "__main__":
    main()
