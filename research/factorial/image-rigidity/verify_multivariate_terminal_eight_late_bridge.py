#!/usr/bin/env python3
"""Exact shared terminal-eight projection bridge on .72 <= u <= .88."""

from fractions import Fraction as Q
from functools import lru_cache
from math import comb, factorial
from sys import argv

import verify_multivariate_terminal_eight_upper_bridge as bridge
from verify_multivariate_terminal_eight_private_linear import (
    ORDERS, coefficient_model, verify_skew_convexity,
)
from verify_multivariate_terminal_seven_private_rows import root_bounds
from verify_multivariate_terminal_six_intermediate_cayley import (
    GRID, I, ONE, ZERO, bisect, complex_mul, interval, modulus, projections, rational, require, scalar,
)
from verify_multivariate_terminal_six_joint_ray_bands import power_bounds
from verify_multivariate_terminal_six_marked_row_tail import projected_interval


ORIGINAL_CAP = bridge.reciprocal_cap
INNER = Q(1, 12)
XI_BOUND = (1-INNER**2)/(1+INNER**2)
RAY_SCALE = 100000
RAYS = {
    72: (131260, 106985, 79864, 62265, 41669, 27105),
    73: (130244, 106397, 79536, 62070, 41565, 27048),
    74: (129251, 105818, 79212, 61877, 41463, 26992),
    75: (128277, 105249, 78893, 61686, 41362, 26937),
    76: (127324, 104688, 78578, 61498, 41262, 26882),
    77: (126391, 104135, 78268, 61312, 41163, 26828),
    78: (125476, 103592, 77961, 61128, 41065, 26774),
    79: (124579, 103056, 77659, 60946, 40968, 26721),
    80: (123700, 102529, 77360, 60766, 40872, 26669),
    81: (122837, 102009, 77066, 60588, 40777, 26617),
    82: (121992, 101497, 76775, 60412, 40683, 26565),
    83: (121162, 100992, 76487, 60238, 40590, 26514),
    84: (120348, 100495, 76204, 60065, 40498, 26463),
    85: (119749, 100005, 75923, 59895, 40406, 26413),
    86: (119164, 99521, 75647, 59726, 40316, 26364),
    87: (118694, 99045, 75373, 59559, 40227, 26315),
    88: (118237, 98575, 75103, 59394, 40138, 26266),
}
MARKED_RAYS = {
    95: (-34250, -50995, -68791, -79732, -90283, -95776),
    96: (-35926, -52274, -69634, -80288, -90554, -95894),
    97: (-37600, -53551, -70477, -80843, -90824, -96013),
    98: (-39270, -54826, -71319, -81397, -91094, -96132),
    99: (-40935, -56098, -72160, -81951, -91365, -96250),
    100: (-42594, -57366, -73001, -82506, -91635, -96369),
}
MARKED_CENTERS = {
    95: (((-3302, 8493), (-7212, 6986), (-9809, 3899), (-10723, 0)), ((-2909, 8535), (-6514, 7464), (-9158, 5033), (-10565, 1772))),
    96: (((-3439, 8404), (-7291, 6907), (-9850, 3853), (-10750, 0)), ((-3052, 8448), (-6603, 7381), (-9209, 4975), (-10596, 1751))),
    97: (((-3575, 8311), (-7369, 6826), (-9891, 3807), (-10779, 0)), ((-3193, 8357), (-6692, 7297), (-9260, 4917), (-10628, 1731))),
    98: (((-3708, 8215), (-7448, 6743), (-9934, 3760), (-10809, 0)), ((-3331, 8262), (-6781, 7210), (-9313, 4858), (-10661, 1710))),
    99: (((-3839, 8114), (-7527, 6659), (-9978, 3713), (-10841, 0)), ((-3467, 8162), (-6869, 7122), (-9367, 4797), (-10697, 1688))),
    100: (((-3968, 8009), (-7606, 6572), (-10024, 3664), (-10876, 0)), ((-3601, 8059), (-6957, 7031), (-9422, 4736), (-10734, 1667))),
}

SECOND_RAYS = {
    94: (-51476, -64247, -77716, -85618, -93187, -97043),
    95: (-53421, -65699, -78639, -86218, -93474, -97169),
    96: (-55367, -67149, -79560, -86817, -93760, -97293),
    97: (-57311, -68598, -80479, -87414, -94045, -97418),
    98: (-59254, -70045, -81396, -88010, -94329, -97542),
    99: (-61192, -71489, -82313, -88605, -94613, -97666),
    100: (-63124, -72930, -83228, -89201, -94898, -97790),
}
SECOND_DISKS = (((-5992, 7136), (-8518, 5301), (-10033, 2810), (-10542, 0)), ((-5727, 7259), (-8100, 5747), (-9667, 3666), (-10457, 1258)))
SECOND_CENTERS = {
    94: (((-5036, 8065), (-8008, 6026), (-9787, 3201), (-10386, 0)), ((-4727, 8190), (-7514, 6521), (-9353, 4171), (-10280, 1433))),
    95: (((-5203, 7928), (-8093, 5913), (-9824, 3139), (-10406, 0)), ((-4901, 8054), (-7613, 6402), (-9402, 4092), (-10304, 1405))),
    96: (((-5366, 7785), (-8178, 5798), (-9862, 3076), (-10429, 0)), ((-5073, 7911), (-7711, 6279), (-9452, 4011), (-10330, 1377))),
    97: (((-5528, 7634), (-8263, 5680), (-9902, 3012), (-10453, 0)), ((-5242, 7760), (-7809, 6153), (-9503, 3928), (-10358, 1348))),
    98: (((-5686, 7476), (-8348, 5557), (-9943, 2946), (-10480, 0)), ((-5408, 7601), (-7906, 6022), (-9556, 3843), (-10388, 1319))),
    99: (((-5841, 7310), (-8433, 5431), (-9987, 2879), (-10510, 0)), ((-5570, 7434), (-8003, 5887), (-9610, 3756), (-10421, 1289))),
    100: (((-5992, 7136), (-8518, 5301), (-10033, 2810), (-10542, 0)), ((-5727, 7259), (-8100, 5747), (-9667, 3666), (-10457, 1258))),
}


MARKED_CASES = {
    1: (95, Q(3, 500), MARKED_RAYS, MARKED_CENTERS),
    2: (94, Q(1, 100), SECOND_RAYS, SECOND_CENTERS),
}


def reciprocal_cap(u, exclude_second):
    unit_cap = Q(1, 22) if exclude_second else Q(1, 15)
    return min(ORIGINAL_CAP(u), max(unit_cap, (1-u)/(7*u-1)))


@lru_cache(maxsize=8192)
def unmarked_model(order, u):
    base, groups = bridge.cayley_model(order, u)
    coarse = []
    for degree, coefficient in base:
        for parameter, terms in groups:
            coefficient = coefficient + parameter*terms[degree][1]
        coarse.append((degree, coefficient))
    sparse = tuple((parameter, tuple((d, c) for d, c in terms if c.lo or c.hi))
                   for parameter, terms in groups)
    p = interval(1-u[1], 1-u[0])
    lam = interval(power_bounds(1-u[1], Q(1, 2))[0], power_bounds(1-u[0], Q(1, 2))[1])
    return base, sparse, tuple(coarse), p, lam


def ray_at(rays, left, right, index, u):
    a, b = rays[left][index], rays[right][index]
    return Q(a, RAY_SCALE)+Q(b-a, RAY_SCALE)*(100*u-left)/(right-left)


def coordinates(p, lam, b, xi):
    y = I(0, b.hi)
    num = b*xi.square()
    for _ in range(3):
        z = b-num/(ONE+p*b*y)
        y = I(max(y.lo, z.lo), min(y.hi, z.hi))
    return xi/lam, y


def projection(order, x, y, base, groups, coarse):
    directions = tuple(projections(order, x, y))
    for direction in directions:
        bound = projected_interval(coarse, direction)
        if bound.lo > 0 or bound.hi < 0:
            return True
    for direction in directions:
        bound = projected_interval(base, direction)
        for parameter, group in groups:
            bound = bound+parameter*projected_interval(group, direction)
        if bound.lo > 0 or bound.hi < 0:
            return True
    return False


def verify_unmarked_ray(order, left, right, index):
    pending = [((Q(left, 100), Q(right, 100)), interval(Q(0), XI_BOUND), 0)]
    count = 0
    while pending:
        u, xi, depth = pending.pop()
        bounds = tuple(sorted(ray_at(RAYS, left, right, index, q) for q in u))
        base, groups, coarse, p, lam = unmarked_model(order, u)
        x, y = coordinates(p, lam, interval(*bounds), xi)
        if projection(order, x, y, base, groups, coarse):
            count += 1
            continue
        require(depth < 34, ("unresolved unmarked ray", order, left, index, u, xi.shown()))
        if 75*(u[1]-u[0]) > Q(xi.hi-xi.lo, GRID):
            pending.extend((half, xi, depth+1) for half in bridge.bisect_u(u))
        else:
            pending.extend((u, half, depth+1) for half in bisect(xi))
    return count


def marked_core(order, base):
    return tuple(root_bounds(Q(comb(order,k)**order*factorial(order*base+k)**order,
                               factorial(order*base)**(order-k)*factorial(order*(base+1))**k), order)
                 for k in range(order+1))


def marked_centre(table, order, left, right, index, v):
    first, last = (table[p][order-7][index] for p in (left, right))
    return tuple(Q(a, 10000)+Q(b-a, 10000)*(100*v-left)/(right-left)
                 for a, b in zip(first, last, strict=True))


def marked_disk_check(order, centre, model, radius):
    x, y = (rational(q) for q in centre)
    powers = [(ONE, ZERO)]
    for _ in range(order):
        powers.append(complex_mul(powers[-1], (x, y)))
    taylor = []
    for h in range(order+1):
        real = imag = ZERO
        for parameter, vector in model:
            term_real = term_imag = ZERO
            for k in range(h, order+1):
                factor = scalar(comb(k, h))*vector[k]
                a, b = powers[k-h]
                term_real = term_real+factor*a
                term_imag = term_imag+factor*b
            real = real+parameter*term_real
            imag = imag+parameter*term_imag
        taylor.append(modulus((real, imag)))
    error = taylor[0]
    for k in range(2, order+1):
        error = error+taylor[k]*rational(radius**k)
    return (rational(radius)*taylor[1]).lo > error.hi


def marked_sector(centre, rays, sector, radius):
    x, y = map(rational, centre)
    for i, cosine in enumerate(rays):
        sine = (ONE-cosine.square()).sqrt()
        side = cosine*y-sine*x
        if i < sector and side.lo <= rational(radius).hi:
            return False
        if i >= sector and side.hi >= rational(-radius).lo:
            return False
    return True


def verify_marked_band(order, left, right, core, rays, table, radius):
    pending = [((Q(left, 100), Q(right, 100)), 0)]
    sectors = (1, 3, 5) if order == 7 else (0, 2, 4, 6)
    count = 0
    while pending:
        v, depth = pending.pop()
        c = (v[0]**-2-1)/6
        model = coefficient_model(order, core, v, (c, 7*c/(8-c)))
        midpoint = sum(v)/2
        centres = [marked_centre(table, order, left, right, k, midpoint)
                   for k in range(len(table[left][order-7]))]
        ray_box = [interval(*sorted(ray_at(rays, left, right, i, end) for end in v))
                for i in range(6)]
        full_centres = centres+[(x, -y) for x, y in centres if y]
        require(len(full_centres) == order, "marked disk count")
        disjoint = all((a-x)**2+(b-y)**2 > 4*radius**2
                       for i, (x, y) in enumerate(full_centres) for a, b in full_centres[:i])
        real_disk = order == 8 or (centres[-1][1] == 0 and centres[-1][0] < -radius)
        sectors_certified = all(y > radius and marked_sector((x, y), ray_box, sector, radius)
                                for (x, y), sector in zip(centres[:len(sectors)], sectors, strict=True))
        if disjoint and real_disk and sectors_certified:
            if all(marked_disk_check(order, centre, model, radius) for centre in centres):
                count += 1
                continue
        require(depth < 14, ("unresolved marked disk", order, left, v))
        pending.extend((half, depth+1) for half in bridge.bisect_u(v))
    return count


def verify_marked(base):
    require(Q(4, 7) < Q(18, 25), "unique (1,2) row")
    require(Q(18, 25)**2 > Q(4, 7)*Q(19, 20)**2, "(1,2) residual cutoff")
    require(Q(5, 7) < Q(4, 5), "unique (2,3) row")
    require(Q(4, 5)**2 > Q(5, 7)*Q(47, 50)**2, "(2,3) residual cutoff")
    lower, radius, rays, table = MARKED_CASES[base]
    for row in rays.values():
        values = (RAY_SCALE, *row, -RAY_SCALE)
        require(all(a > b for a, b in zip(values, values[1:])), "marked ray order")
    total = 0
    for order in ORDERS:
        core = marked_core(order, base)
        total += sum(verify_marked_band(order, left, left+1, core, rays, table, radius)
                     for left in range(lower, 100))
    print(f"Marked ({base},{base+1}): {total} root-disk intervals.", flush=True)


def verify_unmarked(first, last):
    bridge.verify_anchor()
    require(RAYS[88] == bridge.RAYS[88], "matching reference rays")
    require((1+INNER)**8 < 2, "physical annulus")
    for left, right in zip(sorted(RAYS), sorted(RAYS)[1:]):
        require((1-Q(left, 100))*Q(max(RAYS[left][0],RAYS[right][0]),RAY_SCALE)**2 < 1,
                ("negative ray cosines", left, right))
    for row in RAYS.values():
        values = (*row, 0)
        require(all(a > b for a, b in zip(values, values[1:])), "unmarked ray order")
    require(72 <= first < last <= 88, "unmarked replay band")
    total = 0
    for left in range(first, last):
        exclude_second = 80 <= left < 85
        bridge.reciprocal_cap = lambda u: reciprocal_cap(u, exclude_second)
        bridge.cayley_model.cache_clear()
        unmarked_model.cache_clear()
        counts = [[verify_unmarked_ray(j, left, left+1, index) for index in range(6)]
                  for j in ORDERS]
        total += sum(sum(row) for row in counts)
        print(left, left+1, counts, flush=True)
    print(f"Certified {total} unmarked Cayley cells.", flush=True)


def main():
    verify_skew_convexity()
    if len(argv) == 1 or argv[1] == "marked":
        verify_marked(1)
        verify_marked(2)
    if len(argv) == 1:
        verify_unmarked(72, 88)
    elif argv[1] != "marked":
        verify_unmarked(*map(int, argv[1:]))


if __name__ == "__main__":
    main()
