#!/usr/bin/env python3
"""Exact arithmetic checks for the intermediate binary-quadratic sectors."""

from collections import Counter
from fractions import Fraction as Q
from math import factorial
import json


def digits(value):
    return value.bit_count()


def delta(residue, degree):
    return digits((residue-degree) % 8)-digits(residue)-digits((degree+7-residue)//8)


def pure(order, degree, scale):
    return degree*scale-digits(degree)+digits(order-degree)-digits(order)


def check_factorial_ratios(counts):
    for i in range(15):
        for ell in range(15):
            for s in range(15):
                k = i+ell+s
                denominator = factorial(k)*factorial(ell)*factorial(i)*factorial(s)
                numerator = factorial(2*i+ell)*factorial(2*s+ell)
                assert numerator % denominator == 0
                counts["integer_factorial_ratios"] += 1


def verify_dyadic_hierarchy(counts):
    for size in (4, 8, 16, 32, 64):
        for numerator in (5, 6, 7):
            r = Q(numerator, 4*size)
            for order in range(1, 3*size):
                expected = (order % size)*r-digits(order % size)
                values = [i*r+digits(order-i)-digits(order) for i in range(order+1)]
                assert values[order % size] == expected and values.count(min(values)) == 1
                assert min(values) == expected
                counts["axis_annulus_minima"] += 1
            for order in (2*size-1, 4*size-1, 6*size-1):
                selected = 2*size-1
                assert all(
                    pure(order, k, r) > pure(order, selected, r)
                    for k in range(order+1) if k != selected
                )
                counts["deeper_pure_minima"] += 1


def check_open_axis_sector(counts):
    for k in range(1, 55):
        assert all(
            (k-7)*r+6-2*digits(k) > 0
            for r in (Q(13, 50), Q(2, 5), Q(49, 100)) if k != 7
        )
        for ell in range(k+1):
            if (k, ell) != (1, 1):
                def margin(r):
                    return (k-4*ell+3)*r+2*ell-digits(ell)+delta(1, k)
                assert margin(Q(1, 4)) >= 0 and margin(Q(1, 2)) > 0
                counts["open_annulus_term_margins"] += 1


def check_high_axes(counts):
    half = Q(1, 2)
    endpoint = Q(2, 3)
    for k in range(1, 105):
        assert k*half-digits(k)+delta(0, k) > 0
        assert k*half+delta(4, k) >= 1
        if k != 3:
            assert all((k-3)*q+2-digits(k)+delta(4, k) > 0
                       for q in (Q(13, 50), Q(3, 10), half))
        counts["high_axis_order_bounds"] += 1
    for u, selected, value in ((1, 1, Q(-4, 3)), (3, 3, Q(-2)), (5, 1, Q(-4, 3))):
        for k in range(1, 65):
            for ell in range(k+1):
                if k == ell == selected:
                    continue
                assert k*half+(endpoint-half)*ell-digits(ell)+delta(u, k) > value
                counts["endpoint_neighbor_comparisons"] += 1

    for V in range(1, 34):
        j = 8*V+1
        for q in (Q(13, 50), Q(3, 10), Q(1, 3), Q(2, 5), Q(49, 100)):
            if V % 2 and q == Q(1, 3):
                order = j-1
                selected = 7
            else:
                order = j
                selected = 7 if V % 2 and q < Q(1, 3) else 1
            assert all(
                pure(order, k, q) > pure(order, selected, q)
                for k in range(order+1) if k != selected
            )
            counts["pure_cross_switches"] += 1


def run():
    counts = Counter()
    check_factorial_ratios(counts)
    verify_dyadic_hierarchy(counts)
    check_open_axis_sector(counts)
    check_high_axes(counts)
    print(json.dumps(counts, indent=2, sort_keys=True))


if __name__ == "__main__":
    run()
