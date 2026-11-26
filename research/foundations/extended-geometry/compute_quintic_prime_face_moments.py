#!/usr/bin/env python3
"""Generate sparse exact moments for normalized stripped quintic charts."""

import argparse
from collections import defaultdict
from fractions import Fraction
from functools import cache
from math import comb, factorial


def chart(name):
    """Return (x exponent, y exponent, t exponent, coefficient variable)."""
    if name == "I":
        points = [(j, r) for j in range(6) for r in range(6 - j) if j + 2*r <= 7]
        positive = [(1, 0, 2, None)]
        fixed_negative = (3, 2)
        parameter = (1, 3)
    elif name == "II":
        points = [(j, r) for j in range(5) for r in range(5 - j) if 3*j + 2*r <= 9]
        positive = [(1, 0, 2, "q0"), (1, 0, 3, "q1"), (1, 0, 4, None)]
        fixed_negative = (3, 0)
        parameter = (1, 3)
    else:
        points = [(j, r) for j in range(6) for r in range(0, 6 - j, 2)]
        positive = [(1, 0, 2, "q0"), (1, 0, 4, None)]
        fixed_negative = (3, 2) if name == "even-D0" else (5, 0)
        parameter = (1, 4)
        if name == "even-D0":
            points.remove((5, 0))
    negative = []
    for j, r in points:
        variable = None if (j, r) == fixed_negative else f"p{j}{r}"
        if (j, r) == parameter:
            variable = "a"
        negative.append((0, j, r, variable))
    return positive + sorted(negative, key=lambda term: -term[1])


@cache
def monomial_moment(power, height):
    if height % 2:
        return Fraction(0)
    half_height = height // 2
    return Fraction(
        2**power * factorial(power) * factorial(height) * factorial(power + half_height),
        factorial(half_height) * factorial(2*power + height + 1),
    )


def moments(terms, order):
    """Calculate a moment with multinomial weights and angular pruning."""
    weights = [i - j for i, j, _, _ in terms]
    coeffs = defaultdict(Fraction)
    variables = []
    bounds = [(min(weights[k:]), max(weights[k:])) for k in range(len(terms))]

    def visit(k, remaining, x_power, y_power, height, multiplier):
        if not remaining:
            if x_power == y_power and height % 2 == 0:
                coeffs[tuple(variables)] += multiplier * monomial_moment(x_power, height)
            return
        if k == len(terms):
            return
        low, high = bounds[k]
        if not low * remaining <= y_power - x_power <= high * remaining:
            return
        i, j, r, variable = terms[k]
        for count in range(remaining + 1):
            if count and variable is not None:
                variables.append((variable, count))
            visit(k + 1, remaining - count, x_power + i*count, y_power + j*count,
                  height + r*count, multiplier * comb(remaining, count))
            if count and variable is not None:
                variables.pop()

    visit(0, order, 0, 0, 0, 1)
    return {power: value for power, value in coeffs.items() if value}


def literal_moment(terms, order):
    """Independent check by literal Laurent multiplication and integration."""
    phase = defaultdict(int)
    for i, j, r, variable in terms:
        for k in range(j + 1):
            phase[(i - j, r + 2*k, (variable,) if variable else ())] += Fraction(
                (-1)**k * comb(j, k), 2**j,
            )
    current = {(0, 0, ()): Fraction(1)}
    for _ in range(order):
        new = defaultdict(Fraction)
        for (i, r, names), a in current.items():
            for (j, s, others), b in phase.items():
                new[(i + j, r + s, tuple(sorted(names + others)))] += a*b
        current = {term: value for term, value in new.items() if value}
    result = defaultdict(Fraction)
    for (i, r, names), coefficient in current.items():
        if i == 0 and r % 2 == 0:
            result[tuple(sorted((name, names.count(name)) for name in set(names)))] += (
                coefficient / (r + 1)
            )
    return {power: value for power, value in result.items() if value}


def check():
    for name in ("I", "II", "even-D0", "even-BD"):
        terms = chart(name)
        for order in range(1, 5):
            actual = {tuple(sorted(power)): value for power, value in moments(terms, order).items()}
            assert actual == literal_moment(terms, order), (name, order)
    assert moments([(1, 0, 0, None), (0, 1, 0, None)], 2) == {(): Fraction(2, 3)}
    assert moments([(1, 0, 0, None)], 2) == {}
    # A triangular rational section proves that the finite-prefix elimination
    # ideal in Q[a] is zero. Restricting to these slices suffices for that claim.
    for name, pivots in (
        ("I", ["p00", "p10", "p20", "p30", "p40", "p50"]),
        ("II", ["p00", "p10", "p20", "p21"]),
    ):
        terms = []
        for i, j, r, variable in chart(name):
            if variable is None or variable == "a" or variable in pivots:
                terms.append((i, j, r, variable))
            elif name == "II" and variable == "p01":
                terms.append((i, j, r, None))
        for order, pivot in enumerate(pivots, 1):
            polynomial = moments(terms, order)
            pivot_terms = [power for power in polynomial if any(v == pivot for v, _ in power)]
            assert pivot_terms == [((pivot, 1),)], (name, order, pivot_terms)
            for power in polynomial:
                assert not {v for v, _ in power}.intersection(pivots[order:])


def format_polynomial(polynomial):
    return "+".join(
        "*".join([f"({value})"] + [f"{name}^{count}" for name, count in power])
        for power, value in polynomial.items()
    ) or "0"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("chart", choices=["I", "II", "even-D0", "even-BD"], nargs="?")
    parser.add_argument("order", type=int, nargs="?")
    args = parser.parse_args()
    if args.chart is None and args.order is None:
        if not __debug__:
            raise RuntimeError("run with assertions enabled")
        check()
        print("Independent rational moment checks passed for all four charts through order 4.")
        print("Triangular sections certify zero parameter elimination through I:6 and II:4.")
        return
    if args.chart is None or args.order is None or args.order < 1:
        parser.error("give a chart and a positive maximum order, or omit both to run the checks")
    terms = chart(args.chart)
    variables = sorted({variable for _, _, _, variable in terms if variable and variable != "a"})
    print(f"ring R=0,({','.join(variables + ['a'])}), (dp({len(variables)}),lp(1));")
    orders = []
    for order in range(1, args.order + 1):
        orders.append(f"m{order}")
        print(f"poly m{order}={format_polynomial(moments(terms, order))};", flush=True)
    print(f"ideal M={','.join(orders)};")


if __name__ == "__main__":
    main()
