#!/usr/bin/env python3
"""Exact terminal-eight angle separation for total private degree at least two."""

from fractions import Fraction as Q
from math import factorial

from verify_multivariate_terminal_eight_high_curvature import COSINES
from verify_multivariate_terminal_seven_private_rows import (
    DISK_RADIUS,
    coefficient_box,
    private_core,
    root_bounds,
    verify_annulus,
    verify_disk,
)
from verify_multivariate_terminal_six_joint_ray_bands import radial_replay
from verify_multivariate_terminal_six_private_angles import (
    coefficient_bounds, exponential_lower, multiply_interval, require,
)


ORDERS = (7, 8)
LOSS = ((Q(1), Q(3, 2)), (Q(9, 4), Q(9, 5)),
        (Q(3), Q(201, 100)), (Q(13, 4), Q(221, 100)))
HIGH_CASES = (
    ((1, 1, 1, 1), (), Q(0)),
    ((3,), (), Q(0)),
    ((1, 1, 1), (1,), Q(0)),
    ((1, 1), (1, 1), Q(0)),
    ((1, 1, 1), (), Q(5, 4)),
    ((1, 2), (), Q(3, 4)),
    ((1, 1), (), Q(13, 4)),
    ((2,), (), Q(11, 4)),
    ((1,), (1,), Q(3)),
    ((1, 1), (1,), Q(1)),
    ((2,), (1,), Q(1, 2)),
)
CASES = (
    ((1, 1, 1), (), (
        (8678, 7130, 1838, -1538, -5698, -8436),
        (8751, 7230, 1990, -1409, -5618, -8404),
        (8812, 7315, 2123, -1185, -5590, -8426),
        (8864, 7388, 2240, -1084, -5527, -8452),
        (8925, 7476, 2387, -1066, -5446, -8421),
        (8961, 7530, 2481, -984, -5393, -8454),
    ), (
        ((7734, 4863), (929, 9938), (-7217, 7718), (-10771, 0)),
        ((8072, 4224), (2635, 9476), (-4842, 9217), (-10033, 3773)),
    )),
    ((1, 2), (), (
        (8791, 7279, 2045, -1253, -5595, -8443),
        (8845, 7356, 2172, -1143, -5569, -8416),
        (8890, 7422, 2283, -1045, -5507, -8444),
        (8928, 7478, 2381, -960, -5453, -8475),
    ), (
        ((7915, 4693), (1149, 9904), (-7097, 7773), (-10722, 0)),
        ((8223, 4092), (2824, 9421), (-4704, 9249), (-9969, 3801)),
    )),
    ((1, 1), (), (
        (7570, 5388, 430, -2619, -6412, -8535),
        (7833, 5563, 717, -2264, -6304, -8471),
        (8037, 5586, 1000, -2033, -6176, -8414),
        (8218, 6032, 1251, -1827, -6021, -8408),
        (8408, 6297, 1569, -1773, -5830, -8391),
        (8539, 6572, 1774, -1605, -5729, -8397),
        (8646, 6709, 1952, -1565, -5641, -8411),
        (8732, 6932, 2074, -1437, -5604, -8379),
        (8802, 7134, 2209, -1324, -5535, -8403),
        (8880, 7151, 2363, -1277, -5435, -8416),
        (8917, 7319, 2462, -1193, -5426, -8396),
        (8952, 7369, 2546, -1121, -5379, -8378),
        (8994, 7439, 2635, -1129, -5309, -8405),
        (9014, 7573, 2693, -1079, -5321, -8392),
    ), (
        ((6152, 5839), (-534, 10033), (-8106, 7460), (-11346, 0)),
        ((6576, 5168), (1168, 9731), (-5900, 9039), (-10701, 3639)),
    )),
    ((2,), (), (
        (8012, 5749, 872, -2138, -6210, -8444),
        (8192, 5740, 1141, -1916, -6082, -8388),
        (8345, 6278, 1378, -1720, -5968, -8384),
        (8466, 6554, 1587, -1546, -5867, -8388),
        (8613, 6642, 1851, -1539, -5736, -8382),
        (8707, 6876, 1988, -1398, -5650, -8347),
        (8783, 7088, 2135, -1275, -5575, -8366),
        (8839, 7173, 2262, -1167, -5508, -8391),
        (8889, 7244, 2372, -1072, -5449, -8368),
        (8943, 7338, 2492, -1053, -5365, -8389),
        (8974, 7385, 2572, -984, -5364, -8371),
        (9000, 7528, 2640, -924, -5325, -8411),
    ), (
        ((6590, 5396), (-45, 10031), (-7909, 7641), (-11329, 0)),
        ((6940, 4759), (1622, 9652), (-5629, 9179), (-10641, 3730)),
    )),
    ((1,), (1,), (
        (7388, 5840, 935, -2692, -6036, -8592),
        (7703, 6156, 1155, -2343, -5955, -8470),
        (7966, 6421, 1361, -2194, -5877, -8439),
        (8169, 6641, 1550, -1953, -5803, -8409),
        (8397, 6797, 1800, -1764, -5743, -8368),
        (8539, 6955, 1966, -1743, -5672, -8389),
        (8659, 7085, 2113, -1629, -5607, -8363),
        (8742, 7192, 2209, -1636, -5548, -8391),
        (8814, 7386, 2323, -1547, -5494, -8370),
        (8885, 7395, 2483, -1530, -5416, -8392),
        (8928, 7451, 2534, -1462, -5372, -8375),
        (8962, 7601, 2607, -1515, -5334, -8360),
        (8995, 7638, 2669, -1464, -5300, -8402),
    ), (
        ((7071, 7071), (0, 10000), (-7071, 7071), (-10000, 0)),
        ((7660, 6428), (1736, 9848), (-5000, 8660), (-9397, 3420)),
    )),
    ((1, 1), (1,), (
        (8693, 7422, 2175, -1575, -5478, -8552),
        (8771, 7516, 2284, -1491, -5429, -8483),
        (8842, 7594, 2380, -1524, -5384, -8519),
        (8893, 7660, 2466, -1456, -5388, -8504),
        (8936, 7714, 2542, -1394, -5350, -8491),
    ), (
        ((8232, 5147), (1235, 9943), (-6774, 7587), (-10216, 0)),
        ((8612, 4469), (2977, 9524), (-4440, 9105), (-9502, 3701)),
    )),
    ((2,), (1,), (
        (8818, 7579, 2371, -1306, -5382, -8517),
        (8874, 7647, 2457, -1236, -5342, -8503),
        (8919, 7704, 2534, -1285, -5305, -8490),
    ), (
        ((8327, 4876), (1446, 9909), (-6717, 7691), (-10273, 0)),
        ((8669, 4246), (3153, 9450), (-4341, 9178), (-9528, 3756)),
    )),
)


def verify_row_loss():
    count = 0
    for limit, factor in LOSS:
        require(7*(factor-1)**2/factor > limit, ("base cutoff", limit))
        require(Q(7**2, 6+7) > limit, ("gap cutoff", limit))
        for base in range(1, 7):
            for gap in range(1, 7):
                strength = Q(gap**2, base+gap)
                if strength > limit or Q(gap, base) <= factor-1:
                    continue
                if limit == Q(9, 4) and strength == limit:
                    continue
                for order in ORDERS:
                    first = factorial(order*base)
                    last = factorial(order*(base+gap))
                    for k in range(1, order):
                        power = Q(first**(order-k)*last**k, factorial(order*base+k*gap)**order)
                        exponent = factor*strength*k*(order-k)/2
                        require(power <= exponential_lower(exponent),
                                ("row loss", limit, base, gap, order, k))
                        count += 1
    return count


def verify_binomial_anchor(cosines):
    knots = (Q(1), *cosines, Q(-1))
    require(all(a > b for a, b in zip(knots, knots[1:])), "binomial ray order")
    for values, signs in (
        ((8*x**3-4*x**2-4*x+1 for x in knots), (1, 1, -1, -1, 1, 1, -1, -1)),
        ((8*x**4-8*x**2+1 for x in knots), (1, -1, -1, 1, 1, -1, -1, 1)),
    ):
        require(all(value*sign > 0 for value, sign in zip(values, signs, strict=True)),
                "binomial anchor sector")


def verify_high():
    count = 0
    for plus, minus, strength in HIGH_CASES:
        cosines = COSINES
        if plus == (3,):
            cosines = (*COSINES[:2], Q(7, 25), *COSINES[3:])
        verify_binomial_anchor(cosines)
        for order in ORDERS:
            caps = coefficient_bounds(order, plus, strength, minus)
            box = [(Q(1), Q(1))] + [(Q(0), c) for c in caps[1:-1]] + [(Q(1), Q(1))]
            verify_annulus(box, order)
            for cosine in cosines:
                cells = radial_replay(order, box, (cosine, cosine))
                require(cells is not None, ("unresolved high ray", plus, minus, strength, order, cosine))
                count += cells
    return count


def verify_low_anchor(box, centers, cosines, order):
    upper = tuple((Q(x, 10000), Q(y, 10000)) for x, y in centers if y > 0)
    require(len(upper) == order//2, ("upper anchor disks", order))
    disks = upper + tuple((x, -y) for x, y in upper)
    for _, y in upper:
        require(y > DISK_RADIUS, ("upper disk", order))
    if order == 7:
        x, y = centers[-1]
        require(y == 0 and Q(x, 10000) < -DISK_RADIUS, "negative real disk")
        disks += ((Q(x, 10000), Q(0)),)
    require(len(disks) == order == len(centers)+len(upper), "anchor disk count")
    for i, (x, y) in enumerate(disks):
        verify_disk(box, (x, y))
        for a, b in disks[:i]:
            require((x-a)**2+(y-b)**2 > (2*DISK_RADIUS)**2, "disjoint anchor disks")

    sectors = range(1, 6, 2) if order == 7 else range(0, 7, 2)
    for (x, y), sector in zip(upper, sectors, strict=True):
        for j, c in enumerate(cosines):
            sine = root_bounds(1-c*c, 2)
            lo, hi = multiply_interval(sine, (-x, -x))
            if j < sector:
                require(c*y+lo > DISK_RADIUS, ("anchor ray side", sector, j, order))
            else:
                require(c*y+hi < -DISK_RADIUS, ("anchor ray side", sector, j, order))


def ray_at(knots, index, ray, value):
    a, b = knots[index][ray], knots[index+1][ray]
    return Q(a, 10000) + Q(b-a, 10000)*(4*value-index)


def verify_low_band(order, core, knots, index, ray):
    first, last = Q(index, 4), Q(index+1, 4)
    factor = next(f for limit, f in LOSS if last <= limit)
    pending = [((first, last), 0)]
    count = 0
    while pending:
        strength, depth = pending.pop()
        cosines = tuple(sorted(ray_at(knots, index, ray, s) for s in strength))
        cells = radial_replay(order, coefficient_box(order, core, factor, strength), cosines)
        if cells is not None:
            count += cells
            continue
        require(depth < 12, ("unresolved low ray", order, index, ray, strength))
        middle = sum(strength)/2
        pending.extend((((strength[0], middle), depth+1), ((middle, strength[1]), depth+1)))
    return count


def verify_low_case(plus, minus, knots, centers):
    require((plus, minus, Q(len(knots)-1, 4)) in HIGH_CASES,
            ("upper shared-curvature cutoff", plus, minus))
    for values in knots:
        cosines = (10000, *values, -10000)
        require(all(a > b for a, b in zip(cosines, cosines[1:])), ("ray order", plus, minus))
    total = 0
    for order, disks in zip(ORDERS, centers, strict=True):
        core = private_core(order, plus, minus)
        verify_annulus(core, order)
        verify_low_anchor(core, disks, tuple(Q(c, 10000) for c in knots[0]), order)
        total += sum(verify_low_band(order, core, knots, index, ray)
                     for index in range(len(knots)-1) for ray in range(6))
    return total


def main():
    print(f"Checked {verify_row_loss()} retained row-index inequalities.", flush=True)
    print(f"Checked {verify_high()} high-band radial cells.", flush=True)
    total = 0
    for plus, minus, knots, centers in CASES:
        cells = verify_low_case(plus, minus, knots, centers)
        total += cells
        print(f"{plus!r}, {minus!r}, through shared curvature {Q(len(knots)-1, 4)}: "
              f"{cells} radial cells.", flush=True)
    print(f"Verified {total} low-band radial cells.", flush=True)


if __name__ == "__main__":
    main()
