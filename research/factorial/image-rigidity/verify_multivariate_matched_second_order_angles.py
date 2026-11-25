#!/usr/bin/env python3
"""Check the second-order angular formula by exact finite polynomial expansion."""

from fractions import Fraction as Q
from math import comb, factorial

from sympy import Poly, QQ, Symbol


Z = Symbol("z")


def require(condition, context):
    if not condition:
        raise ArithmeticError(context)


def factor_logs(order, index, factors, max_power):
    logs = [Q(0)] * (max_power + 1)
    for reciprocal, orientation in factors:
        k = index if orientation == 1 else order - index
        for power in range(1, max_power + 1):
            chord = Q(sum(h**power for h in range(k))) - Q(
                k * sum(h**power for h in range(order)), order
            )
            logs[power] += Q((-1)**(power+1), power) * reciprocal**power * chord
    return logs


def exp_coefficients(logs):
    coeffs = [Q(1)]
    for degree in range(1, len(logs)):
        coeffs.append(sum(
            power * logs[power] * coeffs[degree-power]
            for power in range(1, degree+1)
        ) / degree)
    return coeffs


def direct_coefficients(order, factors):
    max_power = (order + 2) // 2
    terms = [dict() for _ in range(order + 3)]
    for k in range(order + 1):
        center = Q(2*k-order, 2)
        coeffs = exp_coefficients(factor_logs(order, k, factors, max_power))
        signed_binomial = (-1)**k * comb(order, k)
        for power, coeff in enumerate(coeffs):
            for degree in range(order + 3 - 2*power):
                total = terms[2*power + degree]
                contribution = signed_binomial * coeff * (-center)**degree / factorial(degree)
                total[(degree,)] = total.get((degree,), Q(0)) + contribution
    require(not any(any(row.values()) for row in terms[:order]), ("exponential-scale cancellation", order))
    return tuple(Poly.from_dict(row, Z, domain=QQ) for row in terms[order:])


def hermite(order, variance):
    return Poly.from_dict({
        (order - 2*power,): Q(factorial(order),
            2**power * factorial(power) * factorial(order - 2*power)) * variance**power
        for power in range(order // 2 + 1)
    }, Z, domain=QQ)


def verify_formula(order, factors):
    variance = sum(x for x, _ in factors)
    v2 = sum(x*x for x, _ in factors)
    v3 = sum(x**3 for x, _ in factors)
    signed2 = sum(sign*x*x for x, sign in factors)
    first, second, third = direct_coefficients(order, factors)
    anchor = hermite(order, variance)
    anchor_derivative = anchor.diff()
    d = signed2 / 6

    expected_second = d * anchor.diff((Z, 3))
    expected_third = (
        -variance * Q(order**2, 8) * anchor
        -Q(order - 1, 4) * v2 * anchor.diff((Z, 2))
        +v3 / 12 * anchor.diff((Z, 4))
        +d*d / 2 * anchor.diff((Z, 6))
        +Q(order, 24) * hermite(order + 2, variance)
    )
    require(first == anchor and second == expected_second and third == expected_third,
            ("formal coefficients", order, factors))

    # On an imaginary root, z**2 is minus the squared height.
    angular_model = Poly(
        Z * (
            -order * variance / 24 - Q(order-1, 4) * v2 / variance
            +v3 / (12 * variance**2) * (2*order-3 + Z**2 / variance)
            -d*d / variance**3 * (Q(7*order, 2)-5 + 2*Z**2 / variance)
        ), Z, domain=QQ
    )
    # Clear the denominators of the two root-displacement coefficients.
    # Division by the anchor tests the identity at all its roots simultaneously.
    numerator = (
        -2 * third * anchor_derivative**2
        +2 * second * second.diff() * anchor_derivative
        -second**2 * anchor_derivative.diff()
        -2 * angular_model * anchor_derivative**3
    )
    require(numerator.rem(anchor).is_zero, ("angular correction", order, factors))


def ordinary_factors(order, rows):
    return tuple(
        (Q(gap, order*lower + h), orientation)
        for lower, gap, orientation in rows
        for h in range(1, gap+1)
    )


def main():
    samples = (
        ((Q(1, 3), 1), (Q(1, 4), 1), (Q(1, 7), -1)),
        ((Q(1, 2), 1), (Q(1, 2), -1)),
        ((Q(1, 2), 1), (Q(1, 3), -1), (Q(1, 6), -1)),
        ((Q(4, 5), 1),),
    )
    ordinary = (
        ((7, 1, 1), (11, 3, -1)),
        ((12, 2, 1), (12, 2, -1)),
        ((19, 1, 1), (37, 1, -1), (63, 2, -1)),
    )
    tested = 0
    for order in range(2, 12):
        for factors in samples:
            verify_formula(order, factors)
            tested += 1
        for rows in ordinary:
            lower = ordinary_factors(order, rows)
            upper = ordinary_factors(order+1, rows)
            total = sum(x for x, _ in lower)
            for degree, factors in ((order, lower), (order+1, upper)):
                verify_formula(degree, tuple((x/total, sign) for x, sign in factors))
                tested += 1
    print(f"Verified {tested} exact second-order coefficient and root-reduction checks.")


if __name__ == "__main__":
    main()
