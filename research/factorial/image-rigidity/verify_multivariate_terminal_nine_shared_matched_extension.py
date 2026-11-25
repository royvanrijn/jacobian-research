#!/usr/bin/env python3
"""Exact continuation of the matched light rays at terminal nine."""

from fractions import Fraction as Q
from math import comb

from verify_multivariate_terminal_eight_paired_bridge import radial
from verify_multivariate_terminal_nine_private_rows import outer_coefficient_box, outer_exponents
from verify_multivariate_terminal_nine_shared_matched_light import (
    ANNULUS, CAPS, ORDERS, RAYS as PRECEDING_RAYS, RAY_SCALE, coefficient_model,
)
from verify_multivariate_terminal_six_intermediate_cayley import interval
from verify_multivariate_terminal_six_private_angles import require


RAYS = {
    44: PRECEDING_RAYS[44],
    46: (238165, 104122, -302282, -409963, -720430, -803969, -953711),
    48: (196334, 81204, -329959, -433921, -732385, -812464, -955756),
}


def verify_band(order, left, right):
    outer = Q(left, 100), Q(right, 100)
    pending = [(outer, (Q(0), Q(1)), tuple(range(7)), 0, 0)]
    total = 0
    while pending:
        u, p, rays, mean_depth, skew_depth = pending.pop()
        model = coefficient_model(order, u, p)
        if model is None:
            continue
        unresolved = []
        for ray in rays:
            a, b = RAYS[left][ray], RAYS[right][ray]
            cosines = tuple(sorted(Q(a, RAY_SCALE) + Q(b-a, RAY_SCALE)*(100*v-left)/(right-left)
                                   for v in u))
            cells = radial(order, interval(*cosines), model) if model else None
            if cells is None:
                unresolved.append(ray)
            else:
                total += cells
        if not unresolved:
            continue
        require(max(mean_depth, skew_depth) < 16, ("unresolved ray", order, u, p, unresolved))
        if mean_depth < 2*skew_depth + 2:
            center = sum(u)/2
            pending.extend((((u[0], center), p, unresolved, mean_depth+1, skew_depth),
                            ((center, u[1]), p, unresolved, mean_depth+1, skew_depth)))
        else:
            center = sum(p)/2
            pending.extend(((u, (p[0], center), unresolved, mean_depth, skew_depth+1),
                            (u, (center, p[1]), unresolved, mean_depth, skew_depth+1)))
    return total


def main():
    small, large = ANNULUS
    for values in RAYS.values():
        boundaries = (RAY_SCALE, *values, -RAY_SCALE)
        require(len(values) == 7 and all(a > b for a, b in zip(boundaries, boundaries[1:])),
                "ordered rays")
    for order in ORDERS:
        raw = tuple((Q(comb(order, k)),)*2 for k in range(order+1))
        caps = outer_exponents(order, *CAPS)
        box = outer_coefficient_box(raw, (Q(0), Q(max(RAYS), 100)), caps)
        require(sum(cap*small**k for k, (_, cap) in enumerate(box) if k) < 1,
                ("inner radius", order))
        require(sum(cap*large**k for k, (_, cap) in enumerate(box[:-1])) < large**order,
                ("outer radius", order))
    knots = tuple(sorted(RAYS))
    for left, right in zip(knots, knots[1:]):
        print(left, right, tuple(verify_band(order, left, right) for order in ORDERS), flush=True)


if __name__ == "__main__":
    main()
