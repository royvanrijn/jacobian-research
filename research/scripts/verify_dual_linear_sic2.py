#!/usr/bin/env python3
"""Exact regression for the dual-linear SIC(2) theorem.

The all-degree proof is in extended-geometry/DUAL_LINEAR_SIC2.md.  This
dependency-free checker verifies the second-contraction identity and
replays the normal form, pure moments, mixed cutoff, and Keller subcase.
"""

from __future__ import annotations

from math import prod


# Sparse polynomials in (w1,w2,x,y).
Exponent = tuple[int, int, int, int]
Polynomial = dict[Exponent, int]
ZERO: Exponent = (0, 0, 0, 0)


def clean(poly: Polynomial) -> Polynomial:
    return {exponent: coefficient for exponent, coefficient in poly.items() if coefficient}


def add(*polys: Polynomial) -> Polynomial:
    result: Polynomial = {}
    for poly in polys:
        for exponent, coefficient in poly.items():
            result[exponent] = result.get(exponent, 0) + coefficient
    return clean(result)


def scale(poly: Polynomial, scalar: int) -> Polynomial:
    return clean({exponent: scalar * coefficient for exponent, coefficient in poly.items()})


def multiply(left: Polynomial, right: Polynomial) -> Polynomial:
    result: Polynomial = {}
    for left_exponent, left_coefficient in left.items():
        for right_exponent, right_coefficient in right.items():
            exponent = tuple(a + b for a, b in zip(left_exponent, right_exponent))
            result[exponent] = result.get(exponent, 0) + left_coefficient * right_coefficient
    return clean(result)


def power(poly: Polynomial, exponent: int) -> Polynomial:
    result = {ZERO: 1}
    factor = poly
    while exponent:
        if exponent & 1:
            result = multiply(result, factor)
        factor = multiply(factor, factor)
        exponent //= 2
    return result


def derivative(poly: Polynomial, variable: int) -> Polynomial:
    result: Polynomial = {}
    for exponent, coefficient in poly.items():
        if exponent[variable] == 0:
            continue
        reduced = list(exponent)
        reduced[variable] -= 1
        reduced_tuple = tuple(reduced)
        result[reduced_tuple] = (
            result.get(reduced_tuple, 0) + coefficient * exponent[variable]
        )
    return clean(result)


def contraction(poly: Polynomial) -> Polynomial:
    """Apply E_2 and return a polynomial in x,y."""
    result: Polynomial = {}
    for (w1, w2, x, y), coefficient in poly.items():
        if x < w1 or y < w2:
            continue
        factor = prod(range(x - w1 + 1, x + 1)) * prod(range(y - w2 + 1, y + 1))
        exponent = (0, 0, x - w1, y - w2)
        result[exponent] = result.get(exponent, 0) + coefficient * factor
    return clean(result)


def univariate_at_linear(coefficients: dict[int, int], a: int, b: int) -> Polynomial:
    linear = clean({(0, 0, 1, 0): a, (0, 0, 0, 1): b})
    return add(
        *(
            scale(power(linear, degree), coefficient)
            for degree, coefficient in coefficients.items()
        )
    )


def coordinate_degree(poly: Polynomial) -> int:
    return max((x + y for (_, _, x, y) in poly), default=-1)


def dual_linear_polynomial(h1: Polynomial, h2: Polynomial) -> Polynomial:
    return add(
        multiply({(1, 0, 0, 0): 1}, h1),
        multiply({(0, 1, 0, 0): 1}, h2),
    )


def check_second_contraction_identity() -> None:
    h1 = {
        (0, 0, 2, 0): 1,
        (0, 0, 1, 1): -2,
        (0, 0, 0, 1): 3,
    }
    h2 = {
        (0, 0, 1, 0): -1,
        (0, 0, 0, 2): 2,
        (0, 0, 1, 1): 1,
    }
    p = dual_linear_polynomial(h1, h2)
    h1_x, h1_y = derivative(h1, 2), derivative(h1, 3)
    h2_x, h2_y = derivative(h2, 2), derivative(h2, 3)
    divergence = add(h1_x, h2_y)
    h_dot_gradient_divergence = add(
        multiply(h1, derivative(divergence, 2)),
        multiply(h2, derivative(divergence, 3)),
    )
    trace_square = add(
        power(h1_x, 2),
        scale(multiply(h1_y, h2_x), 2),
        power(h2_y, 2),
    )
    expected = add(
        scale(h_dot_gradient_divergence, 2),
        power(divergence, 2),
        trace_square,
    )
    assert contraction(power(p, 2)) == expected


def verify_normal_form(
    a: int,
    b: int,
    constant: tuple[int, int],
    coefficients: dict[int, int],
) -> None:
    f_of_l = univariate_at_linear(coefficients, a, b)
    constant_h1 = {(0, 0, 0, 0): constant[0]}
    constant_h2 = {(0, 0, 0, 0): constant[1]}
    h1 = add(constant_h1, scale(f_of_l, b))
    h2 = add(constant_h2, scale(f_of_l, -a))
    p = dual_linear_polynomial(h1, h2)

    assert add(derivative(h1, 2), derivative(h2, 3)) == {}
    jacobian_determinant = add(
        multiply(derivative(h1, 2), derivative(h2, 3)),
        scale(multiply(derivative(h1, 3), derivative(h2, 2)), -1),
    )
    assert jacobian_determinant == {}
    assert contraction(p) == {}
    assert contraction(power(p, 2)) == {}
    for exponent in range(1, 8):
        assert contraction(power(p, exponent)) == {}

    # Coordinate degree two; f has degree two, so the theorem gives m>8.
    g = {
        (0, 0, 0, 0): 3,
        (1, 0, 1, 1): -2,
        (0, 2, 0, 2): 1,
        (1, 1, 2, 0): 5,
    }
    assert coordinate_degree(g) == 2
    assert max(coefficients) == 2
    for exponent in (9, 10):
        assert contraction(multiply(g, power(p, exponent))) == {}


def verify_keller_subcase() -> None:
    a, b = 2, 3
    f_of_l = univariate_at_linear({2: 1, 3: -1}, a, b)
    h1 = scale(f_of_l, b)
    h2 = scale(f_of_l, -a)
    p = dual_linear_polynomial(h1, h2)
    g = {
        (0, 0, 2, 1): 1,
        (2, 0, 0, 2): -3,
    }
    assert coordinate_degree(g) == 3
    for exponent in (4, 5):
        assert contraction(multiply(g, power(p, exponent))) == {}


def main() -> None:
    check_second_contraction_identity()
    verify_normal_form(2, 3, (5, -4), {0: 1, 1: -2, 2: 1})
    verify_normal_form(1, 0, (2, 3), {0: -1, 1: 3, 2: 2})
    verify_normal_form(0, 1, (-3, 2), {0: 2, 1: 1, 2: -1})
    verify_keller_subcase()

    print("PASS dual-linear SIC(2): exact second-contraction identity")
    print("PASS dual-linear SIC(2): exact normal-form examples have zero first moments")
    print("PASS dual-linear SIC(2): all pure contractions vanish through order seven")
    print("PASS dual-linear SIC(2): mixed contractions obey the explicit cutoff")
    print("PASS dual-linear SIC(2): normalized Keller subcase has the sharper cutoff")


if __name__ == "__main__":
    main()
