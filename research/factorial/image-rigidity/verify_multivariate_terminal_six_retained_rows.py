#!/usr/bin/env python3
"""Certify terminal-six angle separation by retaining the large rows."""

from fractions import Fraction as Q
from functools import lru_cache
from math import comb, factorial, prod

from flint import acb_poly, arb, ctx, fmpq

from verify_multivariate_terminal_six_joint_ray_bands import (
    BANDS, INNER, OUTER, coefficient_box, exponents, power_bounds, radial_replay,
)
from verify_multivariate_terminal_six_private_angles import multiply_interval, require


ORDERS = (5, 6)
MAX_CURVATURE = Q(4)
TAIL_CAP = Q(1, 8)
TAIL_EXPONENT_CAP = Q(6, 5)
EXCEPTIONAL_CORE = ((1, 1, -1),)
SECTOR_ORDERS = (6, 5, 6, 5, 6)
RAY_FRACTIONS = (Q(1, 2), Q(1, 3), Q(2, 3), Q(1, 5), Q(4, 5))


def ball(value):
    return arb(fmpq(value.numerator, value.denominator))


def rows_above(cap, budget):
    """Enumerate shared ordinary rows matching a strict Gauss cap and a weak budget."""
    for gap in range(1, 1 + int(budget * (1 + 1 / (5 * cap)))):
        for lower in range(1, 1 + int(Q(gap, 5) / cap - Q(1, 5))):
            curvature = Q(gap * gap, lower + gap)
            if curvature <= budget and Q(gap, 5 * lower + 1) > cap:
                yield lower, gap, curvature


def enumerate_cores():
    cores = {}
    for lower, gap, spent in rows_above(TAIL_CAP, MAX_CURVATURE):
        if spent >= MAX_CURVATURE:
            continue
        remainder = MAX_CURVATURE - spent
        candidates = tuple((a, d, sign, curvature)
                           for a, d, curvature in rows_above(TAIL_CAP, remainder)
                           for sign in (-1, 1))

        def collect(index, budget, core):
            if index == len(candidates):
                direct = tuple(sorted(core))
                reversed_columns = tuple(sorted((a, d, -sign) for a, d, sign in core))
                cores[min(direct, reversed_columns)] = budget
                return
            a, d, sign, curvature = candidates[index]
            for count in range(1 + int(budget / curvature)):
                collect(index + 1, budget - count * curvature,
                        core + ((a, d, sign),) * count)

        collect(0, remainder, ((lower, gap, 1),))
    return sorted(cores.items())


@lru_cache(maxsize=None)
def row_powers(order, a, d, sign):
    """The n-th powers of the endpoint-normalized multipliers of a single row."""
    first = factorial(order * a)
    last = factorial(order * (a + d))
    values = tuple(Q(factorial(order * a + d * k)**order, first**(order - k) * last**k)
                   for k in range(order + 1))
    return values if sign == 1 else values[::-1]


def core_powers(order, rows):
    return tuple(prod(column) for column in zip(*(row_powers(order, *row) for row in rows)))


def box_at(core_bounds, tail_bounds, tail):
    return {order: [multiply_interval(a, b) for a, b in zip(
                core_bounds[order], coefficient_box(order, tail, tail_bounds[order]), strict=True)]
            for order in ORDERS}


def nominal_cosines(core, value, tail_bounds):
    cosines = []
    for order in ORDERS:
        anchor = [arb(1)]
        for k, (bounds, power) in enumerate(zip(tail_bounds[order], core[order][1:-1],
                                                  strict=True), 1):
            exponent = Q(5 * k * (order - k), 4 * order)
            require(bounds[0] <= exponent <= bounds[1], ("nominal choice", order, k, bounds))
            anchor.append(comb(order, k) * ball(power)**ball(Q(1, order))
                          * ball(value)**ball(exponent))
        anchor.append(arb(1))
        try:
            roots = acb_poly(anchor).roots(maxprec=256)
        except ValueError:
            return None
        if len(roots) != order:
            return None
        upper = [root.real / abs(root) for root in roots if root.imag > 0]
        lower_count = sum(root.imag < 0 for root in roots)
        if len(upper) != order // 2 or lower_count != order // 2:
            return None
        cosines.extend((cosine, order) for cosine in upper)
    # Sorting midpoints proposes the order; interval comparisons certify it.
    cosines.sort(key=lambda record: float(record[0].mid()), reverse=True)
    if tuple(order for _, order in cosines) != SECTOR_ORDERS:
        return None
    for (left, _), (right, _) in zip(cosines, cosines[1:]):
        if not left > right:
            return None
    return [cosine for cosine, _ in cosines]


def ray_counts(angles, boxes):
    checked = 0
    for left, right in zip(angles, angles[1:]):
        for position in RAY_FRACTIONS:
            proposal = float(left.mid()) * float(position) + float(right.mid()) * float(1 - position)
            cosine = Q(round(100000 * proposal), 100000)
            if not left > ball(cosine) > right:
                continue
            current = 0
            for order in ORDERS:
                boxes_count = radial_replay(order, boxes[order], (cosine, cosine))
                if boxes_count is None:
                    break
                current += boxes_count
            else:
                checked += current
                break
        else:
            return None
    return checked


def verify_core(rows, remainder, tail_bounds):
    core = {order: core_powers(order, rows) for order in ORDERS}
    core_bounds = {order: [power_bounds(power, Q(1, order)) for power in core[order]]
                   for order in ORDERS}
    for order in ORDERS:
        maxima = [comb(order, k) * hi for k, (_, hi) in enumerate(core_bounds[order])]
        require(sum(x * INNER**k for k, x in enumerate(maxima) if k) < 1,
                ("inner radius", rows, order))
        require(sum(x * OUTER**k for k, x in enumerate(maxima[:-1])) < OUTER**order,
                ("outer radius", rows, order))
    bottom = power_bounds(Q(7, 11), TAIL_EXPONENT_CAP * remainder)[0]
    stack = [((bottom, Q(1)), 0)]
    accepted = checked = 0
    while stack:
        tail, depth = stack.pop()
        midpoint = sum(tail) / 2
        cosines = nominal_cosines(core, midpoint, tail_bounds)
        if cosines is not None:
            boxes = box_at(core_bounds, tail_bounds, tail)
            count = ray_counts(cosines, boxes)
            if count is not None:
                accepted += 1
                checked += count
                continue
        require(depth < 16, ("unresolved tail", rows, tail, cosines))
        stack.append(((tail[0], midpoint), depth + 1))
        stack.append(((midpoint, tail[1]), depth + 1))
    return accepted, checked


def verify():
    ctx.prec = 128
    constants = next(data for _, _, cap, data in BANDS if cap == TAIL_CAP)
    tail_bounds = exponents(TAIL_CAP, constants)
    cores = enumerate_cores()
    segments = radial = 0
    require(len(cores) == 823 and any(rows == EXCEPTIONAL_CORE for rows, _ in cores),
            "retained core enumeration")
    for rows, budget in cores:
        if rows == EXCEPTIONAL_CORE:
            continue
        intervals, cells = verify_core(rows, budget, tail_bounds)
        segments += intervals
        radial += cells
    print(f"Verified {len(cores) - 1} retained cores, {segments} tail intervals "
          f"and {radial} radial cells.")


if __name__ == "__main__":
    verify()
