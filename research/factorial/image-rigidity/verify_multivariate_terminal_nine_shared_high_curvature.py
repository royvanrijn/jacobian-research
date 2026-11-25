#!/usr/bin/env python3
"""Certify fully shared terminal-nine separation by exact moving rays."""

from fractions import Fraction as Q
from functools import lru_cache
from math import comb, prod

from verify_multivariate_terminal_nine_high_curvature import (
    CURVATURE as EXISTING_CURVATURE,
    anchor_value,
)
from verify_multivariate_terminal_six_joint_factor_links import check_nonnegative_on_interval, x
from verify_multivariate_terminal_eight_shared_high_curvature import (
    skew_weights,
    weighted_log_bounds,
)
from verify_multivariate_terminal_six_joint_ray_bands import (
    INNER,
    OUTER,
    coefficient_box,
    log_bounds,
    radial_replay,
    rounded_interval,
)
from verify_multivariate_terminal_six_private_angles import (
    coefficient_bounds,
    exponential_lower,
    require,
)


CURVATURE = Q(53, 10)
FLOOR_CAP = Q(1, 8)
GAUSSIAN_RATIO = Q(64, 63)
ATTENUATION = Q(1, 153)

# Integer numerators of ray cosines, by outer-mean knot in thousandths.
RAYS = {
    0: (93243, 74896, 44698, 11263, -27144, -60668, -85173),
    1: (93175, 74816, 44152, 11082, -27875, -61289, -85539),
    2: (93125, 75023, 43836, 10899, -28185, -61652, -85796),
    3: (93074, 74937, 43515, 10715, -28493, -62011, -85945),
    4: (93032, 75140, 43350, 10531, -28661, -62118, -86093),
    5: (92989, 75051, 43104, 10722, -28966, -62223, -86239),
    6: (92933, 74959, 42857, 10537, -29269, -62574, -86383),
    7: (92887, 75157, 42607, 10352, -29432, -62676, -86525),
    8: (92828, 75063, 42356, 10167, -29733, -62777, -86667),
    9: (92813, 74968, 42264, 10358, -29755, -62877, -86706),
    10: (92763, 75163, 42011, 10173, -30053, -62975, -86845),
    11: (92711, 75065, 41757, 9989, -30211, -63073, -86983),
    12: (92658, 74965, 41582, 9804, -30368, -63169, -87119),
    13: (92604, 75157, 41326, 9995, -30524, -63265, -87156),
    14: (92549, 75056, 41069, 9811, -30817, -63360, -87290),
    15: (92492, 74953, 40892, 9627, -30971, -63454, -87423),
    16: (92434, 75142, 40634, 9443, -31124, -63547, -87555),
    17: (92423, 75038, 40701, 9634, -31138, -63398, -87397),
    18: (92374, 74932, 40523, 9451, -31289, -63490, -87528),
    19: (92313, 74825, 40263, 9642, -31440, -63581, -87562),
    20: (92251, 75012, 40085, 9459, -31589, -63672, -87691),
    21: (92199, 74904, 39824, 9277, -31738, -63762, -87820),
    22: (92134, 74794, 39645, 9095, -32022, -63851, -87853),
    23: (92081, 74979, 39465, 9286, -32169, -63940, -87979),
    24: (92014, 74868, 39203, 9104, -32315, -64028, -88011),
    25: (91945, 74756, 39023, 8923, -32460, -64116, -88136),
    26: (91889, 74940, 38843, 9114, -32604, -64202, -88260),
    27: (91818, 74826, 38662, 8934, -32748, -64289, -88291),
    28: (91760, 74712, 38399, 8754, -32891, -64138, -88414),
    29: (91700, 74597, 38218, 8945, -33033, -64223, -88444),
    30: (91719, 74481, 38370, 8765, -33039, -64072, -88290),
    31: (91658, 74662, 38189, 8586, -33180, -64157, -88320),
    32: (91583, 74545, 38008, 8777, -33321, -64242, -88442),
    33: (91520, 74427, 37744, 8598, -33460, -64326, -88471),
    34: (91456, 74607, 37563, 8420, -33599, -64410, -88591),
    35: (91378, 74488, 37382, 8611, -33738, -64258, -88620),
    36: (91312, 74368, 37201, 8433, -33875, -64341, -88738),
    37: (91245, 74247, 37020, 8255, -34012, -64424, -88767),
    38: (91178, 74426, 36839, 8447, -34149, -64506, -88795),
    39: (91110, 74305, 36574, 8269, -34285, -64588, -88911),
    40: (91026, 74182, 36393, 8092, -34420, -64436, -88939),
    41: (90956, 74059, 36212, 8284, -34555, -64518, -89054),
    42: (90885, 74237, 36030, 8107, -34689, -64598, -89081),
    43: (90814, 74113, 35849, 7931, -34822, -64679, -89195),
    44: (90741, 73988, 35668, 8123, -34955, -64528, -89222),
    45: (90668, 74165, 35487, 7947, -35088, -64608, -89335),
    46: (90593, 74039, 35306, 7772, -35219, -64687, -89361),
    47: (90518, 73913, 35040, 7963, -35351, -64766, -89386),
    48: (90427, 73786, 34859, 7788, -35481, -64845, -89498),
    49: (90472, 73658, 35100, 7614, -35480, -64694, -89179),
    50: (90395, 73530, 34919, 7439, -35610, -64543, -89205),
    51: (90318, 73706, 34739, 7631, -35739, -64622, -89231),
    52: (90255, 73577, 34558, 7457, -35868, -64700, -89256),
    53: (90176, 73448, 34378, 7284, -35997, -64778, -89368),
    54: (90096, 73317, 34197, 7475, -36125, -64627, -89393),
    55: (90000, 73493, 34017, 7302, -36252, -64705, -89418),
    56: (89903, 73362, 33836, 7129, -36379, -64782, -89527),
    57: (89820, 73230, 33656, 7320, -36506, -64859, -89552),
    58: (89721, 73098, 33476, 7148, -36632, -64708, -89576),
    59: (89638, 73273, 33295, 6976, -36757, -64785, -89685),
    60: (89537, 73140, 33115, 6804, -36883, -64861, -89709),
    61: (89435, 73007, 32935, 6995, -37007, -64937, -89733),
    62: (89349, 72874, 32756, 6824, -37131, -64787, -89839),
    63: (89245, 73048, 32576, 6653, -37255, -64863, -89863),
    64: (89158, 72914, 32396, 6844, -37379, -64939, -89886),
    65: (89053, 72779, 32216, 6673, -37502, -65014, -89992),
    66: (88963, 72643, 31952, 6502, -37624, -64864, -90015),
    67: (88856, 72817, 31772, 6693, -37746, -64939, -90038),
    68: (88749, 72681, 31593, 6523, -37868, -65014, -90142),
    69: (88657, 72545, 31414, 6353, -37989, -64864, -90164),
    70: (88548, 72408, 31320, 6544, -37982, -64939, -90106),
    71: (88455, 72581, 31141, 6374, -37975, -64789, -90128),
    72: (88536, 72133, 31302, 5845, -38351, -64864, -89743),
    73: (88426, 71995, 31124, 6036, -38471, -64938, -89766),
    74: (88332, 71856, 30945, 5867, -38590, -65012, -89789),
    75: (88237, 71717, 30767, 5699, -38709, -65085, -89812),
    76: (88124, 71891, 30588, 5889, -38827, -65158, -89916),
    77: (88028, 71751, 30410, 5721, -38946, -65010, -89939),
    78: (87931, 68404, 30232, 5553, -39064, -65083, -89962),
    79: (87816, 68586, 30054, 5386, -39055, -64934, -89984),
    80: (87718, 68440, 29961, 5576, -39046, -65007, -89926),
    81: (87619, 68294, 29783, 5409, -39037, -64859, -89867),
    82: (87520, 68147, 29606, 5598, -39028, -64711, -89809),
    83: (87401, 68000, 29513, 5431, -39019, -64563, -89831),
    84: (87319, 67853, 29421, 5265, -39010, -64636, -89773),
    85: (87218, 67705, 29414, 5454, -39126, -64488, -89795),
    86: (87117, 67887, 29322, 5288, -39243, -64561, -89818),
    87: (87033, 67739, 29315, 5477, -39233, -64634, -89840),
    88: (86930, 67591, 29309, 5665, -39349, -64707, -89862),
    89: (86846, 67442, 29217, 5499, -39464, -64779, -89804),
    90: (86761, 67293, 29210, 5333, -39579, -64851, -89826),
    91: (86695, 67474, 29204, 5521, -39694, -64704, -89768),
    92: (86610, 67325, 29197, 5356, -39808, -64776, -89710),
    93: (86543, 67175, 29190, 5190, -39922, -64848, -89652),
    94: (86476, 67025, 29099, 5378, -40036, -64920, -89674),
    95: (86409, 66875, 29093, 5213, -40149, -64991, -89616),
    96: (86362, 67056, 29086, 5048, -40262, -65062, -89558),
    97: (86314, 66906, 28995, 4883, -40375, -65133, -89580),
    98: (86266, 66755, 28904, 5070, -40364, -65204, -89603),
    99: (86217, 66604, 28898, 4906, -40476, -65274, -89545),
}


def possible_rows(floor):
    """Return each row above the given reciprocal floor in the curvature band."""
    gap_bound = int(EXISTING_CURVATURE * (1 + 1 / (8 * floor)))
    for gap in range(1, gap_bound + 1):
        lower_bound = int(gap / (8 * floor) - Q(1, 8))
        for lower in range(1, lower_bound + 1):
            curvature = Q(gap**2, lower + gap)
            if Q(gap, 8 * lower + 1) <= floor or curvature >= EXISTING_CURVATURE:
                continue
            outer_square = prod(Q(8 * lower + h, 8 * lower + 7 * gap + h)
                                for h in range(1, gap + 1))
            yield lower, gap, curvature, outer_square


def profile_outer_square(core_square, core_curvature):
    # The individual loss is 7/16 for the log mean, or 7/8 for its square.
    other_loss = exponential_lower(Q(7, 8) * max(CURVATURE - core_curvature, 0))
    return core_square / other_loss


def reciprocal_caps(rows, outer_lower):
    lower = FLOOR_CAP
    upper = 8 * FLOOR_CAP / (9 - FLOOR_CAP)
    for x8, x9, u_squared_cap in rows:
        if outer_lower**2 <= u_squared_cap:
            lower = max(lower, x8)
            upper = max(upper, x9)
    return lower, upper


def verify_skew_convexity():
    for order in (8, 9):
        for k in range(1, (order + 1) // 2):
            first = sum(w * m / (1 + m*x) for m, w in enumerate(skew_weights(order, k)))
            derivative_scale = 1 + (order - 1)*x
            second = derivative_scale * (derivative_scale * first).diff(x)
            check_nonnegative_on_interval(second.subs(x, x/2))


@lru_cache(maxsize=None)
def exponent_bounds(cap8, cap9):
    denominator_lower, denominator_upper = log_bounds(1 + 7*cap8)
    rho_lower = log_bounds(1 + Q(64, 9)*cap8)[0] / denominator_upper
    rho_upper = 1 + log_bounds(1 / (1-cap8/9))[1] / denominator_lower
    require(rho_lower <= GAUSSIAN_RATIO <= rho_upper, "Gaussian ratio")
    result = {}
    for order, cap, ratio in ((8, cap8, (1, 1)), (9, cap9, (rho_lower, rho_upper))):
        local_outer = log_bounds(1 + (order-1)*cap)[1] / 2
        skew_outer = log_bounds(1 + (order-1)*cap)[0] / 2
        bounds = []
        for k in range(1, order):
            distance = min(k, order-k)
            mean_lower = sum(log_bounds((1+(order-1-m)*cap)/(1+m*cap))[0]
                             for m in range(distance)) / (2 * local_outer)
            skew_upper = weighted_log_bounds(skew_weights(order, distance), cap)[1] / skew_outer
            mean_upper = Q(k * (order-k), order-1)
            bounds.append(rounded_interval(ratio[0] * (mean_lower-skew_upper),
                                           ratio[1] * (mean_upper+skew_upper)))
            require(0 < bounds[-1][0] <= mean_upper * (1 if order == 8 else GAUSSIAN_RATIO)
                    <= bounds[-1][1], "Gaussian anchor envelope")
        result[order] = tuple(bounds)
    return result


def boxes(u, bounds, total_caps):
    put = {
        order: [(lo, min(hi, cap)) for (lo, hi), cap in
                zip(coefficient_box(order, u, bounds[order]), total_caps[order], strict=True)]
        for order in (8, 9)
    }
    lower = [hi / comb(8, k) for k, (_, hi) in enumerate(put[8])]
    for k in range(1, 9):
        lo, hi = put[9][k]
        left = lower[k] / exponential_lower(CURVATURE * k*k * ATTENUATION)
        right = lower[k-1] / exponential_lower(CURVATURE * (9-k)**2 * ATTENUATION)
        put[9][k] = lo, min(hi, comb(9, k) * min(left, right))
    return put


def ray_at(index, left, right, u):
    first, last = Q(RAYS[left][index], 100000), Q(RAYS[right][index], 100000)
    return first + (last-first) * (1000*u-left) / (right-left)


def verify_knot(order, index, left, right, rows, total_caps):
    stack = [((Q(left, 1000), Q(right, 1000)), 0)]
    count = 0
    while stack:
        u, depth = stack.pop()
        caps = reciprocal_caps(rows, u[0])
        coefficients = boxes(u, exponent_bounds(*caps), total_caps)
        if any(lo > hi for array in coefficients.values() for lo, hi in array):
            continue
        cosines = tuple(sorted(ray_at(index, left, right, end) for end in u))
        cells = radial_replay(order, coefficients[order], cosines)
        if cells is not None:
            count += cells
            continue
        require(depth < 16, ("unresolved ray", order, index, u, caps))
        midpoint = sum(u) / 2
        stack.extend((((u[0], midpoint), depth+1), ((midpoint, u[1]), depth+1)))
    return count

def verify_anchor():
    for values in RAYS.values():
        require(len(values) == 7, "number of rays")
        cosines = (100000, *values, -100000)
        require(all(left > right for left, right in zip(cosines, cosines[1:])),
                "ray order")
    knots = (Q(1), *(Q(value, 100000) for value in RAYS[0]), Q(-1))
    for order in (8, 9):
        for index, cosine in enumerate(knots):
            preceding = (index + (order == 9)) // 2
            require((-1) ** preceding * anchor_value(order, cosine) > 0,
                    ("anchor sign", order, index))


def verify():
    verify_skew_convexity()
    verify_anchor()
    knots = sorted(RAYS)
    require(Q(1, 144) > ATTENUATION, "Gaussian endpoint attenuation")
    require(exponential_lower(Q(7, 16) * CURVATURE) * Q(knots[-1], 1000) > 1,
            "outer mean cutoff")
    rows = tuple((Q(d, 8*a+1), Q(d, 9*a+1), profile_outer_square(v2, s))
                 for a, d, s, v2 in possible_rows(FLOOR_CAP))
    require(all(x8 <= 1 and x9 <= 1 for x8, x9, _ in rows), "convexity range")
    total_caps = {order: coefficient_bounds(order, (), CURVATURE)
                  for order in (8, 9)}
    for order, caps in total_caps.items():
        require(sum(cap * INNER**k for k, cap in enumerate(caps) if k) < 1,
                ("small circle", order))
        require(sum(cap * OUTER**k for k, cap in enumerate(caps[:-1])) < OUTER**order,
                ("large circle", order))
    total = 0
    for left, right in zip(knots, knots[1:]):
        cells = sum(verify_knot(order, ray, left, right, rows, total_caps)
                    for order in (8, 9) for ray in range(7))
        total += cells
        print(f"[{left}/1000,{right}/1000]: {cells} radial cells.", flush=True)
    print(f"Verified {total} radial cells.", flush=True)


if __name__ == "__main__":
    verify()
