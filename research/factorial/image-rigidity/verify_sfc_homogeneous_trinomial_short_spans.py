#!/usr/bin/env python3
"""Check a signed cubic remainder for homogeneous trinomial translates."""

import argparse
from math import comb, factorial, prod

from sympy import ZZ
from sympy.polys.rings import ring


def initial_remainder(middle, span):
    height_ring, left, right = ring("L H", ZZ)
    offsets = (0, middle, span)
    first = [
        prod(left + i for i in range(1, x + 1))
        * prod(right + i for i in range(1, span - x + 1))
        for x in offsets
    ]

    forms = []
    for order in (2, 3):
        coefficients = [height_ring.zero for _ in range(order + 1)]
        for a in range(order + 1):
            for b in range(order - a + 1):
                c = order - a - b
                x_offset = middle * b + span * c
                y_offset = order * span - x_offset
                multinomial = factorial(order) // (factorial(a) * factorial(b) * factorial(c))
                term = (
                    multinomial
                    * prod(order * left + i for i in range(1, x_offset + 1))
                    * prod(order * right + i for i in range(1, y_offset + 1))
                    * first[2] ** (a + b)
                    * (-1) ** c
                )
                # The first equation is eliminated by setting the column
                # coefficients to (first[2]*u, first[2]*v, -first[0]*u-first[1]*v).
                for k in range(c + 1):
                    coefficients[b + k] += term * comb(c, k) * first[0] ** (c - k) * first[1] ** k
        forms.append(coefficients)
    quadratic_left, quadratic_mixed, quadratic_right = forms[0]
    cubic_left, _, cubic_linear_left, cubic_right = forms[1]
    return (
        quadratic_left * quadratic_right * cubic_linear_left
        - quadratic_right**2 * cubic_left
        - quadratic_left * quadratic_mixed * cubic_right
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-span", type=int, default=12)
    args = parser.parse_args()
    if args.max_span < 2:
        parser.error("The maximum span must permit three distinct offsets.")

    checked = 0
    for span in range(2, args.max_span + 1):
        for middle in range(1, span // 2 + 1):
            remainder = initial_remainder(middle, span)
            if remainder.get((0, 0), 0) >= 0 or max(remainder.values()) >= 0:
                raise ArithmeticError(f"Sign test failed for offsets {(0, middle, span)}.")
            checked += 1
        print(f"Span {span}: all {span // 2} profiles passed.", flush=True)
    print(f"Verified {checked} profiles, each at all nonnegative integer heights.")


if __name__ == "__main__":
    main()
