#!/usr/bin/env python3
"""Exact rational certification of terminal seven with private degree at least two."""

from fractions import Fraction as Q
from functools import lru_cache
from math import comb, factorial

from verify_multivariate_terminal_seven_high_curvature import COSINES, verify_anchor
from verify_multivariate_terminal_six_joint_ray_bands import INNER, OUTER, radial_replay
from verify_multivariate_terminal_six_private_angles import (
    coefficient_bounds, exponential_lower, multiply_interval, require,
)
from verify_multivariate_terminal_six_two_private_linear import exp_negative_bounds


ORDERS = (6, 7)
LOSS = (Q(3, 2), Q(7, 4), Q(2))
SCALE = 10**8
DISK_RADIUS = Q(1, 1000)
HIGH_CASES = (
    ((1, 1, 1, 1), (), Q(0)),
    ((3,), (), Q(0)),
    ((1, 1, 1), (), Q(1)),
    ((1, 2), (), Q(1, 2)),
    ((1, 1), (), Q(3)),
    ((1,), (1,), Q(3)),
    ((1, 1), (1,), Q(2)),
    ((2,), (1,), Q(2)),
)
CASES = (
    ((1, 1, 1), (), (
        (8194, 5931, -221, -4101, -7986),
        (8294, 6060, -66, -4146, -7943),
        (8378, 5980, 136, -4051, -7905),
        (8450, 6267, 257, -4124, -7870),
        (8511, 6350, 365, -4050, -7839),
    ), (
        ((7214, 5685), (-1396, 10068), (-9545, 4861)),
        ((7734, 4863), (929, 9938), (-7217, 7718), (-10771, 0)),
    )),
    ((1, 2), (), (
        (8360, 6139, 72, -3939, -7926),
        (8433, 6048, 202, -4004, -7889),
        (8496, 6325, 317, -3923, -7855),
    ), (
        ((7445, 5473), (-1152, 10079), (-9466, 4920)),
        ((7915, 4693), (1149, 9904), (-7097, 7773), (-10722, 0)),
    )),
    ((1, 1), (), (
        (6843, 3311, -1552, -4637, -8330),
        (7145, 3443, -1310, -4422, -8256),
        (7404, 3534, -1042, -4369, -8189),
        (7626, 3589, -865, -4343, -8129),
        (7864, 5421, -573, -4288, -8110),
        (8028, 5611, -373, -4152, -8059),
        (8165, 5774, -198, -4184, -8014),
        (8278, 5912, -42, -4234, -7973),
        (8402, 5914, 184, -4082, -7913),
        (8477, 6011, 303, -4163, -7942),
        (8537, 6092, 406, -4093, -7914),
        (8587, 5974, 496, -4033, -7825),
        (8628, 6032, 573, -3980, -7803),
    ), (
        ((5517, 6649), (-2783, 9950), (-10203, 4624)),
        ((6152, 5839), (-534, 10033), (-8106, 7460), (-11346, 0)),
    )),
    ((2,), (), (
        (7387, 3671, -1073, -4293, -8219),
        (7605, 3729, -878, -4109, -8153),
        (7793, 3757, -652, -4091, -8095),
        (7953, 5691, -452, -4097, -8098),
        (8139, 5719, -270, -4078, -8035),
        (8257, 5865, -105, -4119, -7991),
        (8356, 5988, 103, -4020, -7952),
        (8424, 5904, 231, -4091, -7918),
        (8522, 6052, 353, -3966, -7866),
        (8576, 6126, 450, -3899, -7839),
        (8620, 6003, 534, -3841, -7815),
        (8656, 6425, 607, -3790, -7793),
        (8686, 6468, 669, -4247, -7775),
    ), (
        ((6057, 6184), (-2292, 10060), (-10122, 4785)),
        ((6590, 5396), (-45, 10031), (-7909, 7641), (-11329, 0)),
    )),
    ((1,), (1,), (
        (6699, 4157, -1113, -5052, -8094),
        (7041, 5021, -905, -4793, -8050),
        (7342, 5314, -772, -4675, -8064),
        (7591, 5381, -590, -4564, -8024),
        (7853, 5484, -358, -4420, -7971),
        (8028, 5670, -195, -4318, -7933),
        (8171, 5827, -49, -4224, -7897),
        (8288, 5959, 83, -4139, -7865),
        (8413, 5956, 277, -4011, -7877),
        (8488, 6232, 381, -3942, -7787),
        (8548, 6308, 472, -4207, -7762),
        (8597, 6371, 552, -4156, -7740),
        (8637, 6423, 684, -4110, -7720),
    ), (
        ((6235, 7818), (-2225, 9749), (-9010, 4339)),
        ((7071, 7071), (0, 10000), (-7071, 7071), (-10000, 0)),
    )),
    ((1, 1), (1,), (
        (8207, 6392, 196, -4498, -7816),
        (8312, 6329, 245, -4429, -7788),
        (8400, 6432, 347, -4365, -7762),
        (8473, 6521, 440, -4308, -7738),
        (8556, 6623, 553, -4236, -7708),
        (8607, 6689, 631, -4186, -7686),
        (8649, 6745, 700, -4306, -7667),
        (8684, 6792, 760, -4267, -7650),
        (8713, 6649, 813, -4233, -7634),
    ), (
        ((7649, 6004), (-1107, 10010), (-9037, 4753)),
        ((8232, 5147), (1235, 9943), (-6774, 7587), (-10216, 0)),
    )),
    ((2,), (1,), (
        (8387, 6425, 352, -4358, -7758),
        (8463, 6514, 443, -4302, -7735),
        (8526, 6591, 525, -4250, -7713),
        (8578, 6656, 599, -4367, -7694),
        (8637, 6548, 689, -4310, -7669),
        (8674, 6598, 751, -4271, -7651),
        (8705, 6639, 805, -4237, -7636),
        (8730, 6675, 852, -4206, -7622),
        (8750, 6704, 894, -4179, -7609),
    ), (
        ((7807, 5683), (-877, 10044), (-9054, 4852)),
        ((8327, 4876), (1446, 9909), (-6717, 7691), (-10273, 0)),
    )),
)


def verify_row_loss():
    count = 0
    for limit, factor in enumerate(LOSS, 1):
        require(7 * (factor - 1)**2 / factor > limit, ("base cutoff", limit))
        require(Q(7**2, 6 + 7) > limit, ("gap cutoff", limit))
        for base in range(1, 7):
            for gap in range(1, 7):
                strength = Q(gap**2, base + gap)
                if strength > limit or Q(gap, base) <= factor - 1:
                    continue
                for n in ORDERS:
                    const, lead = factorial(n * base), factorial(n * (base + gap))
                    for k in range(1, n):
                        reciprocal = Q(const**(n-k) * lead**k, factorial(n*base+k*gap)**n)
                        require(reciprocal <= exponential_lower(factor*strength*k*(n-k)/2),
                                ("row loss", base, gap, n, k, limit))
                        count += 1
    return count


def rounded(interval):
    lo, hi = (endpoint * SCALE for endpoint in interval)
    return Q(lo.numerator // lo.denominator, SCALE), Q(-((-hi.numerator) // hi.denominator), SCALE)


@lru_cache(maxsize=None)
def root_bounds(value, order):
    require(value >= 0, ("root bound", value))
    if value in (0, 1):
        return value, value
    target = value.numerator * SCALE**order
    denominator = value.denominator
    lo, hi = 0, SCALE * (1 + comb(order, order // 2))
    require(hi**order * denominator >= target, ("root bracket", value, order))
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if mid**order * denominator < target:
            lo = mid
        else:
            hi = mid
    return Q(lo, SCALE), Q(hi, SCALE)


def private_core(order, plus, minus):
    values = []
    for k in range(order + 1):
        power = Q(comb(order, k)**order)
        for index, gaps in ((k, plus), (order-k, minus)):
            for d in gaps:
                power *= Q(factorial(d*index)**order, factorial(d*order)**index)
        values.append(root_bounds(power, order))
    return values


def verify_annulus(box, order):
    require(sum(hi*INNER**k for k, (_, hi) in enumerate(box) if k) < 1,
            ("inner radius", order))
    require(sum(hi*OUTER**k for k, (_, hi) in enumerate(box[:-1])) < OUTER**order,
            ("outer radius", order))


def verify_high():
    verify_anchor()
    count = 0
    for plus, minus, strength in HIGH_CASES:
        for order in ORDERS:
            maxima = coefficient_bounds(order, plus, strength, minus)
            box = [(Q(1), Q(1))] + [(Q(0), x) for x in maxima[1:-1]] + [(Q(1), Q(1))]
            verify_annulus(box, order)
            for cosine in COSINES:
                result = radial_replay(order, box, (cosine, cosine))
                require(result is not None, ("high ray", plus, minus, strength, cosine, order))
                count += result
    return count


def complex_powers(x, y, order):
    powers = [(Q(1), Q(0))]
    for _ in range(order):
        u, v = powers[-1]
        powers.append((u*x-v*y, u*y+v*x))
    return powers


def interval_sum(terms):
    return tuple(sum(side) for side in zip(*terms, strict=True))


def verify_disk(box, center):
    order = len(box)-1
    powers = complex_powers(*center, order)
    local = []
    for m in range(order+1):
        coordinates = [
            interval_sum(multiply_interval(box[k], (comb(k, m)*power[h],)*2)
                         for k, power in enumerate(powers[:order+1-m], m))
            for h in (0, 1)
        ]
        local.append(coordinates)

    def max_abs(coordinates):
        return sum(max(abs(a), abs(b)) for a, b in coordinates)

    derivative = max(max(0, a, -b) for a, b in local[1])
    error = max_abs(local[0]) + sum(max_abs(v)*DISK_RADIUS**m for m, v in enumerate(local[2:], 2))
    require(derivative * DISK_RADIUS > error, ("Rouche disk", order, center))


def verify_low_anchor(box, centers, cosines, order):
    upper = tuple((Q(x, 10000), Q(y, 10000)) for x, y in centers[:3])
    disks = upper + tuple((x, -y) for x, y in upper)
    for _, y in upper:
        require(y > DISK_RADIUS, ("upper disk", order))
    if order == 7:
        x, y = centers[3]
        require(y == 0 and Q(x, 10000) < -DISK_RADIUS, "negative real disk")
        disks += ((Q(x, 10000), Q(0)),)
    require(len(disks) == order, "anchor disk count")
    for i, (x, y) in enumerate(disks):
        verify_disk(box, (x, y))
        for a, b in disks[:i]:
            require((x-a)**2+(y-b)**2 > (2*DISK_RADIUS)**2, "disjoint anchor disks")
    sectors = (1, 3, 5) if order == 6 else (0, 2, 4)
    for (x, y), sector in zip(upper, sectors, strict=True):
        for j, c in enumerate(cosines):
            sine = root_bounds(1-c*c, 2)
            lo, hi = multiply_interval(sine, (-x, -x))
            if j < sector:
                require(c*y + lo > DISK_RADIUS, ("anchor ray side", sector, j, order))
            else:
                require(c*y + hi < -DISK_RADIUS, ("anchor ray side", sector, j, order))


def ray_at(knots, index, ray, value):
    left, right = knots[index][ray], knots[index+1][ray]
    return Q(left, 10000) + Q(right-left, 10000) * (4*value-index)


def coefficient_box(order, core, factor, strength):
    coefficients = []
    for k, interval in enumerate(core):
        weight = Q(k*(order-k), 2*order)
        if weight == 0:
            coefficients.append(interval)
            continue
        lower = exp_negative_bounds(factor*strength[1]*weight)[0]
        upper = exp_negative_bounds(strength[0]*weight)[1]
        coefficients.append(rounded(multiply_interval(interval, (lower, upper))))
    return coefficients


def verify_low_band(order, core, knots, index, ray):
    factor = LOSS[index//4]
    pending = [((Q(index, 4), Q(index+1, 4)), 0)]
    count = 0
    while pending:
        strength, depth = pending.pop()
        cosines = tuple(sorted(ray_at(knots, index, ray, s) for s in strength))
        cell_count = radial_replay(order, coefficient_box(order, core, factor, strength), cosines)
        if cell_count is not None:
            count += cell_count
            continue
        require(depth < 12, ("unresolved low band", order, index, ray, strength))
        middle = sum(strength)/2
        pending.extend((((strength[0], middle), depth+1), ((middle, strength[1]), depth+1)))
    return count


def verify_low_case(plus, minus, knots, centers):
    for rays in knots:
        require(10000 > rays[0] > rays[1] > rays[2] > rays[3] > rays[4] > -10000,
                ("ray order", plus, minus))
    total = 0
    for order, disks in zip(ORDERS, centers, strict=True):
        core = private_core(order, plus, minus)
        verify_annulus(core, order)
        verify_low_anchor(core, disks, tuple(Q(c, 10000) for c in knots[0]), order)
        total += sum(verify_low_band(order, core, knots, index, ray)
                     for index in range(len(knots)-1) for ray in range(5))
    return total


def main():
    print(f"Checked {verify_row_loss()} retained row-index inequalities.", flush=True)
    print(f"Checked {verify_high()} high-band radial cells.", flush=True)
    for plus, minus, knots, centers in CASES:
        cells = verify_low_case(plus, minus, knots, centers)
        print(f"{plus!r}, {minus!r}, through shared curvature {Q(len(knots)-1, 4)}: "
              f"{cells} radial cells.", flush=True)


if __name__ == "__main__":
    main()
