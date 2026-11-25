#!/usr/bin/env python3
"""Exact Cayley projection bridge for terminal seven on 0.75 <= u <= 0.88."""

from fractions import Fraction as Q
from functools import lru_cache
from sys import argv

import verify_multivariate_terminal_seven_upper_bridge as bridge
from verify_multivariate_terminal_six_intermediate_cayley import (
    GRID, bisect, interval, projections, require,
)
from verify_multivariate_terminal_six_marked_row_tail import projected_interval


RAYS = {
    75: (125311, 99322, 70037, 51516, 28922),
    76: (124432, 98841, 69792, 51382, 28862),
    77: (123571, 98366, 69551, 51249, 28803),
    78: (122726, 97898, 69312, 51118, 28745),
    79: (121897, 97436, 69076, 50988, 28687),
    80: (121083, 96981, 68843, 50859, 28630),
    81: (120285, 96532, 68613, 50732, 28573),
    82: (119502, 96090, 68385, 50606, 28517),
    83: (118732, 95653, 68160, 50481, 28462),
    84: (117977, 95222, 67938, 50358, 28407),
    85: (117235, 94797, 67718, 50236, 28352),
    86: (116506, 94377, 67501, 50115, 28298),
    87: (115789, 93963, 67286, 49995, 28245),
    88: (115085, 93554, 67073, 49876, 28192),
}
RAY_SCALE = 100000
DISK = Q(19, 10)
ORIGINAL_CAP = bridge.reciprocal_cap


def reciprocal_cap(u):
    return min(ORIGINAL_CAP(u), max(Q(1, 13), (1-u)/(6*u-1)))


@lru_cache(maxsize=8192)
def disk_check(order, u):
    base, groups = bridge.cayley_model(order, u)
    total = []
    for degree, coefficient in base:
        for parameter, terms in groups:
            coefficient = coefficient + parameter*terms[degree][1]
        total.append(coefficient)
    return total[order].lo*DISK**order > sum(
        max(abs(c.lo), abs(c.hi))*DISK**degree for degree, c in enumerate(total[:-1])
    )


def ray_at(left, right, index, u):
    a, b = RAYS[left][index], RAYS[right][index]
    return Q(a, RAY_SCALE) + Q(b-a, RAY_SCALE)*(100*u-left)/(right-left)


def verify_ray(order, left, right, index):
    pending = [((Q(left, 100), Q(right, 100)), interval(Q(0), DISK), 0)]
    count = 0
    while pending:
        u, x, depth = pending.pop()
        if not disk_check(order, u):
            require(depth < 28, ("scaled disk", order, u))
            pending.extend((half, x, depth+1) for half in bridge.bisect_u(u))
            continue
        p = interval(1-u[1], 1-u[0])
        bounds = tuple(sorted(ray_at(left, right, index, value) for value in u))
        y = bridge.height(p, x, interval(*bounds))
        base, groups = bridge.cayley_model(order, u)
        for direction in projections(order, x, y):
            total = projected_interval(base, direction)
            for parameter, contribution in groups:
                total = total + parameter*projected_interval(contribution, direction)
            if total.lo > 0 or total.hi < 0:
                count += 1
                break
        else:
            require(depth < 34, ("unresolved Cayley ray", order, left, index, u, x.shown()))
            if 100*(u[1]-u[0]) > Q(x.hi-x.lo, GRID):
                pending.extend((half, x, depth+1) for half in bridge.bisect_u(u))
            else:
                pending.extend((u, half, depth+1) for half in bisect(x))
    return count


def main():
    require(Q(7, 12)**2 < Q(3, 4)**2, "unique marked row")
    require(Q(7, 12)*Q(8, 13) < Q(3, 4)**2, "marked residual cap")
    require(Q(12, 5)*Q(1, 3) < Q(17, 4), "marked curvature bound")
    bridge.verify_anchor()
    require(RAYS[88] == bridge.RAYS[88], "matching anchor rays")
    require(DISK**2/4 < 1, "left-half-plane disk")
    for row in RAYS.values():
        require(2*RAY_SCALE > row[0] > row[1] > row[2] > row[3] > row[4] > 0,
                "ray order")

    bridge.reciprocal_cap = reciprocal_cap
    bridge.cayley_model.cache_clear()
    disk_check.cache_clear()
    first, last = (75, 88) if len(argv) == 1 else tuple(map(int, argv[1:]))
    require(75 <= first < last <= 88, "replay band")
    total = 0
    for left in range(first, last):
        counts = [[verify_ray(j, left, left+1, i) for i in range(5)] for j in (6, 7)]
        total += sum(sum(row) for row in counts)
        print(left, left+1, counts, flush=True)
    print(f"Certified {total} Cayley projection cells.", flush=True)


if __name__ == "__main__":
    main()
