#!/usr/bin/env python3
"""Exact terminal-seven angular certificate for one private linear row."""

from fractions import Fraction as Q
from functools import lru_cache
from math import prod

from verify_multivariate_terminal_seven_high_curvature import COSINES
from verify_multivariate_terminal_seven_private_rows import (
    private_core,
    verify_annulus,
)
from verify_multivariate_terminal_six_joint_ray_bands import (
    log_bounds,
    power_bounds,
    radial_replay,
    rounded_interval,
)
from verify_multivariate_terminal_six_private_angles import (
    coefficient_bounds,
    multiply_interval,
    require,
)


ORDERS = (6, 7)
STRENGTH_LIMIT = Q(11, 2)
RECIPROCAL_FLOOR = Q(1, 40)
REFERENCE_RATIO = Q(36, 35)
RAYS = {
    0: (88559, 72445, 12016, -41150, -74961),
    5: (87040, 67796, 7022, -42443, -77627),
    10: (85216, 65323, 3548, -38855, -79081),
    15: (83046, 61261, -202, -39964, -80818),
    20: (80492, 58400, -3792, -40992, -81584),
    25: (77823, 55424, -6797, -41966, -82292),
    30: (74698, 52356, -9686, -41878, -83339),
    35: (71670, 45078, -12905, -41831, -83948),
    40: (68356, 35053, -15596, -42793, -84526),
    45: (65387, 33249, -17373, -43746, -84728),
    50: (62015, 31394, -20338, -45600, -85267),
    55: (58557, 25515, -23227, -46527, -86113),
    60: (54852, 22307, -25652, -48310, -86605),
    65: (51262, 20373, -28031, -50052, -87084),
    70: (47629, 14614, -30750, -51762, -87551),
    75: (43968, 12653, -33052, -53445, -87725),
    80: (40487, 10649, -34241, -56580, -87913),
    85: (37017, 9789, -36175, -58887, -88381),
    90: (33573, -9725, -38117, -61134, -88593),
    95: (29967, -11498, -40075, -62672, -89071),
    100: (26412, -14416, -42384, -64845, -89552),
}


def possible_rows():
    largest_gap = int(STRENGTH_LIMIT * (1 + 1 / (6 * RECIPROCAL_FLOOR)))
    for gap in range(1, largest_gap + 1):
        largest_base = int((gap / RECIPROCAL_FLOOR - 1) / 6)
        for base in range(1, largest_base + 1):
            x6 = Q(gap, 6*base+1)
            if x6 <= RECIPROCAL_FLOOR or Q(gap**2, base+gap) > STRENGTH_LIMIT:
                continue
            mean_squared = prod(Q(6*base+h, 6*base+5*gap+h) for h in range(1, gap+1))
            yield x6, Q(gap, 7*base+1), mean_squared


def reciprocal_caps(rows, outer_lower):
    c6 = c7 = RECIPROCAL_FLOOR
    bound = outer_lower**2
    for x6, x7, mean_squared in rows:
        if bound <= mean_squared:
            c6, c7 = max(c6, x6), max(c7, x7)
    return c6, c7


def skew_weights(order, distance):
    return tuple((Q(distance <= h < order-distance) - Q(order-2*distance, order))/2
                 for h in range(order))


def multiply_polynomials(left, right):
    result = [Q(0)] * (len(left)+len(right)-1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i+j] += a*b
    return result


def verify_skew_convexity():
    for order in ORDERS:
        m = order-1
        for distance in range(1, (order+1)//2):
            weights = skew_weights(order, distance)
            require(sum(h*w for h, w in enumerate(weights)) == 0, ("skew derivative", order, distance))
            numerator = [Q(0)] * (2*m-3)
            for h in range(1, m):
                term = [weights[h]*h*(m-h)]
                for l in range(1, m):
                    if h != l:
                        term = multiply_polynomials(term, (1, 2*l, l*l))
                numerator = [a+b for a, b in zip(numerator, term, strict=True)]

            # For 0 <= x <= 1 one may replace a negative c*x**i
            # by its lower bound c*x**(i-2).
            for i in range(len(numerator)-1, 1, -1):
                if numerator[i] < 0:
                    numerator[i-2] += numerator[i]
                    numerator[i] = 0
            require(min(numerator) >= 0, ("skew convexity", order, distance))


def weighted_log_upper(weights, cap):
    upper = Q(0)
    for h, weight in enumerate(weights):
        if h and weight:
            lo, hi = log_bounds(1+h*cap)
            upper += max(weight*lo, weight*hi)
    return upper


@lru_cache(maxsize=None)
def exponent_bounds(cap6, cap7):
    denominator_lo, denominator_hi = log_bounds(1+5*cap6)
    ratio_lo = log_bounds(1+Q(36, 7)*cap6)[0] / denominator_hi
    ratio_hi = 1 + log_bounds(1/(1-cap6/7))[1] / denominator_lo
    require(ratio_lo <= REFERENCE_RATIO <= ratio_hi, "outer ratio enclosure")
    result = {}
    for order, cap, ratio in ((6, cap6, (1, 1)), (7, cap7, (ratio_lo, ratio_hi))):
        m = order-1
        log_lo, log_hi = log_bounds(1+m*cap)
        bounds = []
        reference_ratio = 1 if order == 6 else REFERENCE_RATIO
        for k in range(1, order):
            distance = min(k, order-k)
            mean_lo = sum(log_bounds((1+(m-h)*cap)/(1+h*cap))[0]
                          for h in range(distance)) / log_hi
            mean_hi = Q(k*(order-k), m)
            skew = 2 * weighted_log_upper(skew_weights(order, distance), cap) / log_lo
            lower, upper = rounded_interval(ratio[0]*(mean_lo-skew), ratio[1]*(mean_hi+skew))
            require(0 < lower <= mean_hi*reference_ratio <= upper,
                    ("canonical exponent", order, k, cap, lower, upper))
            bounds.append((lower, upper))
        result[order] = tuple(bounds)
    return result


def box_at(order, core, u, exponents):
    result = [core[0]]
    for coefficient, (a, b) in zip(core[1:-1], exponents, strict=True):
        multiplier = (power_bounds(u[0], b)[0], power_bounds(u[1], a)[1])
        result.append(multiply_interval(coefficient, multiplier))
    return result + [core[-1]]


def ray_at(index, left, right, u):
    first, last = Q(RAYS[left][index], 100000), Q(RAYS[right][index], 100000)
    return first + (last-first)*(100*u-left)/(right-left)


def verify_anchor(cosines):
    knots = (Q(1), *cosines, Q(-1))
    for values, signs in (
        ((8*c**3-6*c for c in knots), (1, 1, -1, -1, 1, 1, -1)),
        ((8*c**3-4*c**2-4*c+1 for c in knots), (1, -1, -1, 1, 1, -1, -1)),
    ):
        require(all(value*sign > 0 for value, sign in zip(values, signs, strict=True)),
                "anchor sector signs")


def verify_high():
    verify_anchor(COSINES)
    count = 0
    for order in ORDERS:
        caps = coefficient_bounds(order, (1,), STRENGTH_LIMIT)
        box = [(Q(1), Q(1))] + [(Q(0), x) for x in caps[1:-1]] + [(Q(1), Q(1))]
        verify_annulus(box, order)
        for cosine in COSINES:
            result = radial_replay(order, box, (cosine, cosine))
            require(result is not None, ("high-strength ray", order, cosine))
            count += result
    return count


def verify_knot(order, index, left, right, rows, core):
    pending = [((Q(left, 100), Q(right, 100)), 0)]
    count = 0
    while pending:
        u, depth = pending.pop()
        cap = reciprocal_caps(rows, u[0])
        box = box_at(order, core, u, exponent_bounds(*cap)[order])
        cosines = tuple(sorted(ray_at(index, left, right, value) for value in u))
        cells = radial_replay(order, box, cosines)
        if cells is not None:
            count += cells
            continue
        require(depth < 15, ("unresolved ray", order, index, u, cap))
        midpoint = sum(u)/2
        pending.extend((((u[0], midpoint), depth+1), ((midpoint, u[1]), depth+1)))
    return count


def main():
    verify_skew_convexity()
    print(f"Certified {verify_high()} radial cells with shared curvature at least 11/2.", flush=True)
    for values in RAYS.values():
        cosines = (100000, *values, -100000)
        require(all(a > b for a, b in zip(cosines, cosines[1:])), "ray order")
    verify_anchor(tuple(Q(x, 100000) for x in RAYS[0]))

    rows = tuple(possible_rows())
    require(all(0 < x7 <= x6 <= 1 for x6, x7, _ in rows), "skew convexity range")
    cores = {j: private_core(j, (1,), ()) for j in ORDERS}
    for order, core in cores.items():
        verify_annulus(core, order)
    knots = sorted(RAYS)
    total = 0
    for left, right in zip(knots, knots[1:]):
        cells = sum(verify_knot(j, i, left, right, rows, cores[j])
                    for j in ORDERS for i in range(5))
        total += cells
        print(f"[{left}/100,{right}/100]: {cells} radial cells.", flush=True)
    print(f"Verified {total} low-strength radial cells using {len(rows)} candidate rows.",
          flush=True)


if __name__ == "__main__":
    main()
