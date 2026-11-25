#!/usr/bin/env python3
"""Exact paired ray certificates for terminal-eight marked singleton rows."""

from fractions import Fraction as Q
from math import comb, factorial
from sys import argv

from verify_multivariate_terminal_eight_private_linear import (
    ORDERS, coefficient_model, possible_rows, reciprocal_caps, verify_skew_convexity,
)
from verify_multivariate_terminal_eight_private_rows import verify_low_anchor
from verify_multivariate_terminal_seven_paired_light import radial
from verify_multivariate_terminal_seven_private_rows import root_bounds
from verify_multivariate_terminal_six_intermediate_cayley import interval, require


RECIPROCAL_CAP = Q(1, 8)
INNER = Q(1, 10)
MARKED_GAPS = {1: 2, 2: 2, 3: 3}
CASES = {
    1: ({
        9: (89035, 59693, 23500, -14003, -54490, -79472),
        10: (88701, 59138, 22755, -14624, -54895, -79669),
        15: (86824, 56271, 19083, -17619, -56804, -80588),
        20: (84652, 53289, 15510, -20458, -58561, -81424),
        25: (82233, 50222, 12031, -23165, -60198, -82196),
        30: (79601, 47090, 8642, -25760, -61739, -82918),
        35: (76784, 43905, 5335, -28257, -63198, -83598),
        40: (73804, 40678, 2103, -30670, -64588, -84242),
        45: (70679, 37414, -1060, -33007, -65918, -84856),
        50: (67425, 34120, -4159, -35276, -67194, -85443),
        55: (64054, 30800, -7199, -37485, -68424, -86006),
        60: (60579, 27457, -10185, -39638, -69612, -86549),
        65: (57010, 24096, -13122, -41742, -70763, -87074),
        70: (53357, 20718, -16013, -43801, -71881, -87582),
        75: (49629, 17327, -18862, -45820, -72969, -88075),
        80: (45835, 13925, -21673, -47801, -74031, -88556),
        85: (41985, 10515, -24450, -49750, -75069, -89025),
        90: (38087, 7099, -27195, -51670, -76087, -89484),
        95: (34152, 3680, -29913, -53565, -77087, -89934),
        100: (30188, 262, -32607, -55438, -78073, -90378),
    }, (
        ((2369, 8416), (-3888, 9280), (-9072, 5880), (-11057, 0)),
        ((2869, 8127), (-2622, 9515), (-7694, 7367), (-10688, 2753)),
    )),
    2: ({
        9: (87936, 57897, 21115, -15982, -55781, -80095),
        10: (87425, 57118, 20119, -16793, -56298, -80343),
        15: (84550, 53098, 15245, -20681, -58719, -81496),
        20: (81229, 48925, 10542, -24332, -60925, -82535),
        25: (77543, 44643, 5993, -27790, -62965, -83487),
        30: (73550, 40278, 1584, -31089, -64873, -84371),
        35: (69294, 35846, -2701, -34250, -66670, -85199),
        40: (64807, 31360, -6875, -37292, -68375, -85981),
        45: (60118, 26829, -10947, -40230, -69998, -86722),
        50: (55249, 22259, -14927, -43074, -71552, -87428),
        55: (50218, 17656, -18824, -45835, -73043, -88104),
        60: (45041, 13025, -22643, -48521, -74478, -88752),
        65: (39733, 8369, -26392, -51138, -75863, -89376),
        70: (34305, 3691, -30075, -53693, -77204, -89979),
        75: (28771, -1004, -33698, -56191, -78504, -90561),
        80: (23141, -5714, -37265, -58638, -79768, -91126),
        85: (17428, -10436, -40782, -61038, -80999, -91675),
        90: (11644, -15166, -44252, -63397, -82202, -92210),
        95: (5803, -19900, -47682, -65719, -83380, -92734),
        100: (-77, -24632, -51076, -68011, -84538, -93248),
    }, (
        ((-249, 9035), (-5591, 8349), (-9405, 4897), (-10784, 0)),
        ((233, 8938), (-4606, 8767), (-8429, 6239), (-10532, 2248)),
    )),
    3: ({
        9: (89012, 59636, 23415, -14079, -54547, -79498),
        10: (88674, 59074, 22660, -14708, -54958, -79697),
        15: (86779, 56173, 18941, -17743, -56895, -80629),
        20: (84589, 53156, 15320, -20620, -58679, -81477),
        25: (82150, 50055, 11794, -23365, -60343, -82262),
        30: (79498, 46888, 8356, -25998, -61910, -82995),
        35: (76660, 43668, 5001, -28533, -63394, -83686),
        40: (73658, 40404, 1721, -30983, -64808, -84341),
        45: (70509, 37104, -1490, -33356, -66161, -84965),
        50: (67230, 33773, -4637, -35662, -67461, -85562),
        55: (63831, 30414, -7726, -37905, -68712, -86135),
        60: (60325, 27032, -10760, -40093, -69921, -86687),
        65: (56722, 23629, -13744, -42231, -71092, -87220),
        70: (53028, 20208, -16681, -44321, -72228, -87735),
        75: (49254, 16771, -19575, -46369, -73332, -88235),
        80: (45405, 13320, -22429, -48378, -74407, -88721),
        85: (41489, 9857, -25245, -50351, -75456, -89194),
        90: (37513, 6384, -28027, -52291, -76482, -89656),
        95: (33483, 2903, -30777, -54201, -77486, -90107),
        100: (29406, -584, -33498, -56084, -78471, -90549),
    }, (
        ((2416, 8824), (-3957, 9227), (-8866, 5682), (-10689, 0)),
        ((2937, 8576), (-2734, 9528), (-7585, 7153), (-10344, 2634)),
    )),
}


def marked_core(order, base):
    gap = MARKED_GAPS[base]
    first = factorial(order*base)
    last = factorial(order*(base+gap))
    return tuple(root_bounds(Q(comb(order,k)**order*factorial(order*base+gap*k)**order,
                               first**(order-k)*last**k), order)
                 for k in range(order+1))


def tail_caps(rows, lower):
    c7, c8 = reciprocal_caps(rows, lower)
    c7 = min(c7, RECIPROCAL_CAP)
    return c7, min(c8, 7*c7/(8-c7))


def ray_at(rays, left, right, index, v):
    first, last = Q(rays[left][index], 100000), Q(rays[right][index], 100000)
    return first+(last-first)*(100*v-left)/(right-left)


def verify_ray(order, rays, left, right, index, rows, core):
    pending = [((Q(left, 100), Q(right, 100)), 0)]
    count = 0
    while pending:
        v, depth = pending.pop()
        caps = tail_caps(rows, v[0])
        model = coefficient_model(order, core, v, caps)
        bounds = tuple(sorted(ray_at(rays, left, right, index, value) for value in v))
        cells = radial(order, interval(*bounds), model)
        if cells is not None:
            count += cells
            continue
        require(depth < 15, ("unresolved marked ray", order, index, v, caps))
        midpoint = sum(v)/2
        pending.extend((((v[0], midpoint), depth+1), ((midpoint, v[1]), depth+1)))
    return count


def verify_case(base, rays, disks, rows):
    for ray in rays.values():
        cosines = (100000, *ray, -100000)
        require(all(a > b for a, b in zip(cosines, cosines[1:])), "ordered marked rays")
    require(min(rays) == 9 and max(rays) == 100, "residual outer-mean coverage")
    cores = {}
    for order, centers in zip(ORDERS, disks, strict=True):
        core = marked_core(order, base)
        require(sum(c[1]*INNER**k for k, c in enumerate(core) if k) < 1,
                ("inner radius", base, order))
        require(sum(c[1]*INNER**(order-k) for k, c in enumerate(core[:-1])) < 1,
                ("outer radius", base, order))
        verify_low_anchor(core, centers, tuple(Q(x, 100000) for x in rays[100]), order)
        cores[order] = core
    knots = sorted(rays)
    total = 0
    for left, right in zip(knots, knots[1:]):
        count = sum(verify_ray(j, rays, left, right, i, rows, cores[j])
                    for j in ORDERS for i in range(6))
        total += count
        print(base, left, right, count, flush=True)
    return total


def main():
    verify_skew_convexity()
    require(Q(9, 100)**4 < Q(7, 13)**15, "tail outer-mean cutoff")
    rows = tuple(possible_rows())
    for base, (rays, disks) in CASES.items():
        if len(argv) > 1 and base not in map(int, argv[1:]):
            continue
        count = verify_case(base, rays, disks, rows)
        print(f"Marked ({base},{base+MARKED_GAPS[base]}): {count} radial cells.", flush=True)


if __name__ == "__main__":
    main()
