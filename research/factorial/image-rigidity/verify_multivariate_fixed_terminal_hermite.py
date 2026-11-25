#!/usr/bin/env python3
"""Compute a rational radius for the linked fixed-terminal Hermite comparison."""

from fractions import Fraction as Q
from math import comb, factorial
from sys import argv

from sympy import Poly, Rational, symbols


X = symbols("x")


def hermite_coefficients(order, variance):
    return {
        order - 2*p: Q(factorial(order), factorial(p)*factorial(order - 2*p))
        * (variance/2)**p
        for p in range(order//2 + 1)
    }


def real_anchor(order, variance):
    return Poly(
        sum((-1)**((order - degree)//2) * coefficient * X**degree
            for degree, coefficient in hermite_coefficients(order, variance).items()),
        X,
        domain="QQ",
    )


def isolate_pair(orders):
    lower, upper = orders
    polynomials = (real_anchor(lower, Q(1)), real_anchor(upper, Q(lower, upper)))
    if polynomials[0].gcd(polynomials[1]).degree() != 0:
        raise ArithmeticError("the two Hermite anchors share a root")
    precision = 8
    while True:
        roots = []
        for order, poly in zip(orders, polynomials, strict=True):
            intervals = poly.intervals(eps=Rational(1, 2**precision))
            if len(intervals) != order or any(m != 1 for _, m in intervals):
                raise ArithmeticError(("the Hermite roots are not simple real roots", order))
            roots.extend((Q(lo), Q(hi), order) for (lo, hi), _ in intervals)
        roots.sort()
        gaps = [b[0] - a[1] for a, b in zip(roots, roots[1:])]
        if min(gaps) > 0:
            return roots, min(gaps)
        precision += 4


def coefficient_error(order, radius, lambda0, upper_order):
    bound = sum(
        2**(order - degree + 2) * comb(order, degree)
        * lambda0**(degree - order - 1) * radius**degree
        for degree in range(order + 1)
    )
    if order == upper_order:
        bound += sum(
            Q(p*factorial(order), 2**p * factorial(p) * factorial(order - 2*p) * order)
            * radius**(order - 2*p)
            for p in range(1, order//2 + 1)
        )
    return bound


def bound_at_circle(roots, radius):
    margins = []
    for i, (lo, hi, order) in enumerate(roots):
        margin = radius
        for k, (other_lo, other_hi, other_order) in enumerate(roots):
            if i != k and order == other_order:
                gap = max(other_lo - hi, lo - other_hi)
                margin *= gap - radius
        margins.append((order, margin))
    return margins


def certify(lower):
    if lower < 1:
        raise ValueError("the lower order must be positive")
    orders = lower, lower + 1
    upper = orders[1]
    roots, gap = isolate_pair(orders)
    circle = min(gap/8, Q(1, 8))
    radius = 2 + max(max(abs(a), abs(b)) for a, b, _ in roots)
    lambda0 = Q(1, 2*upper)
    errors = {j: coefficient_error(j, radius, lambda0, upper) for j in orders}
    margins = bound_at_circle(roots, circle)
    inverse_power = 1
    while True:
        lam = Q(1, 2**inverse_power)
        if (lam <= lambda0/2 and lam*radius <= 1 and lam**2 * radius**3/3 < circle
                and all(lam*errors[order] < margin for order, margin in margins)):
            break
        inverse_power += 1
    print(f"Orders {lower}, {upper} | {len(roots)} disjoint Hermite root intervals")
    print(f"Certified lambda = 2^-{inverse_power}, hence X_{lower} <= {lam**2}")
    print(f"Sufficient curvature bound: S <= {lower*lam**2 / (1 + lower*lam**2)}")


if __name__ == "__main__":
    certify(int(argv[1]) if len(argv) > 1 else 8)
