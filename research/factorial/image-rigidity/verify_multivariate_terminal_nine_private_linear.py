#!/usr/bin/env python3
"""Exact paired outer-mean replay for one private linear row at terminal nine."""

from fractions import Fraction as Q
from functools import lru_cache

from verify_multivariate_terminal_nine_high_curvature import anchor_value
from verify_multivariate_terminal_nine_private_rows import (
    ORDERS, RAY_SCALE, REFERENCE_RATIO, outer_coefficient_box, outer_exponents,
    possible_rows, reciprocal_caps, verify_skew_convexity,
)
from verify_multivariate_terminal_seven_paired_light import radial
from verify_multivariate_terminal_seven_private_linear import skew_weights
from verify_multivariate_terminal_seven_private_rows import private_core, verify_annulus
from verify_multivariate_terminal_six_intermediate_cayley import (
    GRID, I, ONE, ZERO, interval, log, rational, scalar,
)
from verify_multivariate_terminal_six_joint_ray_bands import log_bounds, power_bounds, radial_replay
from verify_multivariate_terminal_six_private_angles import coefficient_bounds, multiply_interval, require


LIMIT = Q(21, 4)
HIGH_RAYS = (92414, 74840, 41696, 9543, -30137, -63050, -87083)
PAIRED = (24, 88)

RAYS = {
    0: (93248, 74751, 44698, 11450, -27144, -60668, -85173),
    2: (92934, 74988, 42934, 10706, -29149, -62175, -86264),
    4: (92618, 75120, 41551, 9941, -30259, -62863, -86835),
    6: (92210, 74746, 40140, 9749, -31511, -63136, -87368),
    8: (91813, 74330, 38955, 8994, -32336, -63404, -87599),
    10: (91309, 73861, 37649, 8820, -33506, -63635, -88088),
    12: (90751, 73354, 36216, 8108, -34435, -63848, -88555),
    14: (90092, 72815, 34913, 7410, -35334, -64046, -88870),
    16: (89405, 71788, 33486, 7264, -36403, -63891, -89173),
    18: (88643, 71194, 32195, 6587, -37247, -64072, -89591),
    20: (87995, 69225, 31148, 5346, -37953, -64298, -89517),
    22: (87216, 68594, 29872, 4691, -38954, -64135, -89795),
    24: (86336, 66983, 28999, 2472, -39932, -66585, -88819),
    26: (85516, 65341, 27741, 1326, -40700, -66724, -88978),
    28: (84733, 64168, 26492, -7033, -41638, -67174, -89132),
    30: (83855, 62981, 25251, -3483, -42372, -67614, -89405),
    32: (82951, 61780, 24019, -4055, -43091, -68043, -89551),
    34: (82021, 61069, 22796, -5127, -43798, -68161, -89692),
    36: (81067, 59847, 21581, -6185, -44491, -68577, -89948),
    38: (79970, 59119, 20243, -7726, -45174, -68985, -90083),
    39: (79411, 58500, 19642, -8243, -45333, -69186, -90149),
    40: (78846, 58383, 19044, -8756, -45668, -69092, -90328),
    41: (78507, 57357, 18397, -9351, -45938, -69372, -90199),
    42: (77935, 56736, 17801, -9861, -46273, -69282, -90265),
    44: (76853, 55992, 16617, -10387, -46930, -69678, -90506),
    46: (75751, 55240, 15441, -11390, -47577, -70068, -90631),
    48: (74630, 54482, 14271, -11908, -48216, -70453, -90862),
    50: (73491, 54219, 13109, -12896, -48847, -70554, -91089),
    52: (72336, 53945, 11824, -13408, -49471, -70932, -91206),
    54: (71214, 53166, 10675, -14384, -50087, -71304, -91424),
    56: (69992, 51891, 9659, -14897, -50704, -71409, -91539),
    58: (68795, 48098, 8522, -15406, -51470, -71776, -91650),
    60: (67536, 47314, 7518, -16367, -52227, -71878, -91761),
    62: (66317, 47026, 6391, -17322, -52817, -71981, -91870),
    64: (65087, 45232, 5397, -17829, -53559, -72345, -91978),
    66: (63850, 44437, 4278, -18777, -54293, -72451, -92086),
    68: (62606, 43141, 3290, -19284, -54867, -72812, -92192),
    70: (61356, 41845, 2179, -20227, -55590, -72922, -92390),
    72: (60101, 40550, 1197, -20738, -56156, -73281, -92494),
    73: (59417, 39903, 643, -21208, -56438, -73216, -92545),
    74: (58788, 39751, 91, -21677, -56720, -73396, -92597),
    75: (58158, 39103, -461, -22145, -57001, -73576, -92649),
    76: (57639, 38497, -1144, -22635, -57467, -73545, -92712),
    78: (56374, 37696, -2367, -23150, -58023, -73667, -92814),
    80: (55109, 36403, -3586, -24082, -58577, -74028, -92915),
    82: (53845, 35114, -4803, -25012, -59130, -74156, -93101),
    84: (52583, 33829, -6017, -25940, -59682, -74519, -93201),
    86: (51393, 32555, -7229, -26870, -60380, -74660, -93303),
    88: (50199, 30304, -8438, -27795, -60929, -74800, -93403),
    90: (48890, 29037, -10008, -27937, -62290, -74047, -93344),
    92: (47646, 28259, -11091, -28868, -62696, -74432, -93446),
    94: (46347, 27001, -12294, -29799, -63238, -75036, -93626),
    96: (45117, 26233, -13375, -31109, -63651, -75418, -93727),
    98: (43897, 24989, -14455, -32036, -64068, -76014, -93904),
    100: (42687, 24234, -15535, -32962, -64490, -76811, -94151),
}


def verify_binomial_anchor(values):
    sequence = (RAY_SCALE, *values, -RAY_SCALE)
    require(len(values) == 7 and all(a > b for a, b in zip(sequence, sequence[1:])),
            "ordered binomial rays")
    for order in ORDERS:
        for k, value in enumerate(sequence):
            roots = (k + (order == 9)) // 2
            require((-1)**roots * anchor_value(order, Q(value, RAY_SCALE)) > 0,
                    ("binomial anchor", order, k))


def verify_high():
    verify_binomial_anchor(HIGH_RAYS)
    count = 0
    for order in ORDERS:
        caps = coefficient_bounds(order, (1,), LIMIT)
        box = [(Q(1), Q(1))] + [(Q(0), cap) for cap in caps[1:-1]] + [(Q(1), Q(1))]
        verify_annulus(box, order)
        for value in HIGH_RAYS:
            cosine = Q(value, RAY_SCALE)
            cells = radial_replay(order, box, (cosine, cosine))
            require(cells is not None, ("high ray", order, value))
            count += cells
    return count


@lru_cache(maxsize=None)
def shape(order):
    m = order - 1
    distances = tuple(min(k, order-k) for k in range(order+1))
    gaussian = tuple(rational(Q(k*(order-k), m)) for k in range(order+1))
    skew_rows = tuple(skew_weights(order, k) for k in distances)
    points = tuple(Q(1)-Q(2*h,m) for h in range(order))
    skew = []
    for j in (1, 2, 3):
        moment = tuple(sum(w * p**(2*j) for w, p in zip(row, points, strict=True))
                       for row in skew_rows)
        skew.append(tuple(rational(-a if 2*k <= order else a) for k, a in enumerate(moment)))
    defect = []
    for j in (1, 2):
        defect.append(tuple(rational(sum(p-p**(2*j+1) for p in points[:r]))
                            for r in distances))
    skew_sizes = tuple(rational(sum(abs(w) for w in row)) for row in skew_rows)
    return gaussian, tuple(skew), tuple(defect), skew_sizes, distances


@lru_cache(maxsize=None)
def ratio_bounds(cap):
    lo, hi = log_bounds(1+7*cap)
    lower = log_bounds(1+Q(64, 9)*cap)[0]/hi
    upper = 1+log_bounds(1/(1-cap/9))[1]/lo
    require(lower <= REFERENCE_RATIO <= upper, ("ratio", cap))
    return interval(lower, upper)


def paired_box(order, core, means, caps):
    midpoint = sum(means)/2
    ref_ratio = Q(1) if order == 8 else REFERENCE_RATIO
    reference = tuple(interval(*multiply_interval(c, power_bounds(midpoint, ref_ratio*Q(k*(order-k), order-1))))
                      if 0 < k < order else ONE for k, c in enumerate(core))
    outer_log = -log(interval(*means))
    outer_mean = rational(ref_ratio)*log(interval(means[0]/midpoint, means[1]/midpoint))
    if order == 9:
        ratio = ratio_bounds(caps[0])
        outer_mean += (rational(ref_ratio)-ratio)*outer_log
        outer_log = ratio*outer_log
    m = order - 1
    cap = caps[order-8]
    T = rational(m*cap/(2+m*cap))
    gaussian, skew, defect, skew_sizes, distances = shape(order)
    pieces = [(outer_mean, gaussian)]
    for j, vector in enumerate(skew, 1):
        parameter = T**(2*j-1)*outer_log*scalar(1, 2*j)
        pieces.append((parameter.symmetric_error(), vector))
    for j, vector in enumerate(defect, 1):
        parameter = T**(2*j)*outer_log*scalar(1, 2*j+1)
        pieces.append((parameter.positive_error(), vector))
    skew_tail = T**7*outer_log/(scalar(8)*(ONE-T.square()))
    defect_tail = T**6*outer_log/(scalar(7)*(ONE-T.square()))
    errors = []
    for k, c in enumerate(reference):
        tail = skew_sizes[k]*skew_tail+scalar(distances[k])*defect_tail
        size = tail
        for parameter, vector in pieces:
            term = parameter*vector[k]
            size += I(0, max(-term.lo, term.hi))
        if size.hi >= GRID:
            return None
        error = scalar(1, 6)*size**3/(ONE-size)+tail*(ONE+size)
        errors.append(c*(ONE+error.symmetric_error()))

    products = []
    for i, (first, a) in enumerate(pieces):
        for j, (second, b) in enumerate(pieces[:i+1]):
            parameter = scalar(1, 2)*first.square() if i == j else first*second
            products.append((parameter, tuple(v*w for v, w in zip(a, b, strict=True))))
    return ((ONE, tuple(errors)),) + tuple(
        (parameter, tuple(c*v for c, v in zip(reference, vector, strict=True)))
        for parameter, vector in (*pieces, *products)
    )


def verify_knot(order, left, right, rows, core):
    pending = [((Q(left, 100), Q(right, 100)), tuple(range(7)), 0)]
    total = 0
    may_pair = PAIRED[0] <= left < right <= PAIRED[1]
    while pending:
        means, rays, depth = pending.pop()
        caps = reciprocal_caps(rows, means[0])
        independent = outer_coefficient_box(core, means, outer_exponents(order, *caps))
        paired = None
        failed = []
        for ray in rays:
            a, b = RAYS[left][ray], RAYS[right][ray]
            cosines = tuple(sorted(Q(a, RAY_SCALE)+Q(b-a, RAY_SCALE)*(100*u-left)/(right-left)
                                   for u in means))
            cells = radial_replay(order, independent, cosines)
            if cells is None and may_pair:
                if paired is None:
                    paired = paired_box(order, core, means, caps)
                if paired is not None:
                    cells = radial(order, interval(*cosines), paired)
            if cells is None:
                failed.append(ray)
            else:
                total += cells
        if failed:
            require(depth < 13, ("low ray", order, means, failed))
            middle = sum(means)/2
            pending.extend((((means[0], middle), failed, depth+1),
                            ((middle, means[1]), failed, depth+1)))
    return total


def main():
    print(f"Verified {verify_skew_convexity()} shared skew inequalities and "
          f"{verify_high()} high-strength radial cells.", flush=True)
    verify_binomial_anchor(RAYS[0])
    for values in RAYS.values():
        sequence = (RAY_SCALE, *values, -RAY_SCALE)
        require(len(values) == 7 and all(a > b for a, b in zip(sequence, sequence[1:])),
                "ordered continuation rays")
    rows = tuple(possible_rows(LIMIT))
    require(all(0 < x9 <= x8 <= 1 for x8, x9, _ in rows), "reciprocal cap range")
    cores = {order: private_core(order, (1,), ()) for order in ORDERS}
    for order, core in cores.items():
        verify_annulus(core, order)
        require(sum(hi*Q(1, 10)**k for k, (_, hi) in enumerate(core) if k) < 1,
                ("paired inner radius", order))
        require(sum(hi*Q(10)**k for k, (_, hi) in enumerate(core[:-1])) < 10**order,
                ("paired outer radius", order))
    total = 0
    knots = tuple(sorted(RAYS))
    for left, right in zip(knots, knots[1:]):
        total += sum(verify_knot(order, left, right, rows, cores[order]) for order in ORDERS)
    print(f"Verified {total} low-strength radial cells with {len(rows)} retained rows.", flush=True)


if __name__ == "__main__":
    main()
