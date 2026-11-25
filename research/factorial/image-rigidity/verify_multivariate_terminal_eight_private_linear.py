#!/usr/bin/env python3
"""Exact terminal-eight ray exclusions for one private linear row."""

from fractions import Fraction as Q
from functools import lru_cache
from math import prod
from sys import argv

from verify_multivariate_terminal_eight_high_curvature import COSINES
from verify_multivariate_terminal_eight_private_rows import verify_binomial_anchor
from verify_multivariate_terminal_seven_paired_light import radial
from verify_multivariate_terminal_seven_private_linear import (
    multiply_polynomials, skew_weights, weighted_log_upper,
)
from verify_multivariate_terminal_seven_private_rows import private_core, verify_annulus
from verify_multivariate_terminal_six_intermediate_cayley import (
    GRID, I, ONE, ZERO, interval, log, rational, scalar,
)
from verify_multivariate_terminal_six_joint_ray_bands import (
    log_bounds, power_bounds, radial_replay, rounded_interval,
)
from verify_multivariate_terminal_six_private_angles import (
    coefficient_bounds, multiply_interval, require,
)


ORDERS = (7, 8)
STRENGTH_LIMIT = Q(21, 4)
RECIPROCAL_FLOOR = Q(1, 40)
REFERENCE_RATIO = Q(49, 48)
PAIRED_BANDS = ((30, 85),)
RAYS = {
    0: (91375, 80130, 30876, -21974, -49846, -84301),
    5: (90289, 77566, 26807, -10027, -53469, -84567),
    10: (88817, 75576, 23347, -10758, -55376, -83686),
    15: (87050, 72220, 19646, -12501, -56661, -84359),
    20: (85001, 69732, 16694, -15211, -58240, -84486),
    25: (82547, 65958, 13510, -16759, -59709, -84600),
    30: (80043, 62034, 10414, -18252, -61087, -85157),
    34: (77815, 56303, 7996, -20207, -62135, -85143),
    35: (77265, 55714, 7399, -19703, -62391, -85248),
    38: (75580, 52788, 5306, -21139, -63142, -85134),
    40: (74428, 51585, 4139, -22082, -63632, -85342),
    42: (73469, 52654, 3303, -23961, -64113, -85545),
    43: (72878, 52043, 2410, -24417, -64350, -85237),
    45: (71680, 50815, 1270, -25321, -64820, -85440),
    50: (68494, 44309, -1851, -27539, -65963, -85545),
    55: (65321, 40056, -5220, -29702, -67405, -85659),
    57: (64028, 38799, -6303, -30553, -67833, -85855),
    60: (62194, 33566, -8219, -31817, -68795, -85783),
    63: (60213, 32794, -9813, -32223, -69414, -85720),
    65: (58881, 32636, -11169, -33058, -69822, -85918),
    67: (57542, 31369, -12216, -33888, -70542, -85768),
    70: (55522, 27289, -13778, -35125, -71137, -85728),
    72: (54304, 22784, -14813, -36744, -71836, -85930),
    75: (52270, 21995, -16359, -37951, -72417, -85904),
    76: (51729, 20302, -15996, -39134, -72312, -86654),
    78: (50370, 19061, -17027, -39924, -72700, -86848),
    80: (49011, 16763, -18055, -39949, -73086, -87041),
    85: (45619, 11599, -20334, -42673, -73766, -87521),
    90: (42099, 10621, -22615, -45345, -74733, -87999),
    95: (38914, 8591, -24903, -47970, -75703, -88758),
    100: (35634, 7561, -27471, -49882, -76425, -89228),
}


def possible_rows():
    largest_gap = int(STRENGTH_LIMIT * (1+1/(7*RECIPROCAL_FLOOR)))
    for gap in range(1, largest_gap+1):
        largest_base = int((gap/RECIPROCAL_FLOOR-1)/7)
        for base in range(1, largest_base+1):
            x7 = Q(gap, 7*base+1)
            if x7 <= RECIPROCAL_FLOOR or Q(gap**2, base+gap) > STRENGTH_LIMIT:
                continue
            mean_square = prod(Q(7*base+h, 7*base+6*gap+h) for h in range(1, gap+1))
            yield x7, Q(gap, 8*base+1), mean_square


def reciprocal_caps(rows, lower_mean):
    c7 = c8 = RECIPROCAL_FLOOR
    lower_square = lower_mean**2
    for x7, x8, mean_square in rows:
        if lower_square <= mean_square:
            c7, c8 = max(c7, x7), max(c8, x8)
    return c7, c8


def link_weights(order):
    link_coefficients = ((1, -1, 0), (1, -3, 5)) if order == 7 else ((8, -7, 0), (3, -7, 7))
    base = tuple(skew_weights(order, k) for k in (1, 2, 3))
    return tuple(tuple(sum(a*b for a, b in zip(coefficients, column, strict=True))
                       for column in zip(*base, strict=True))
                 for coefficients in link_coefficients)


def verify_skew_convexity():
    for order in ORDERS:
        m = order-1
        weights = tuple(skew_weights(order, k) for k in (1, 2, 3)) + link_weights(order)
        for row in weights:
            require(sum(h*w for h, w in enumerate(row)) == 0, "skew linear term")
            numerator = [Q(0)] * (2*m-3)
            for h in range(1, m):
                term = [row[h]*h*(m-h)]
                for a in range(1, m):
                    if a != h:
                        term = multiply_polynomials(term, (1, 2*a, a*a))
                numerator = [a+b for a, b in zip(numerator, term, strict=True)]
            for i in range(len(numerator)-1, 1, -1):
                if numerator[i] < 0:
                    numerator[i-2] += numerator[i]
                    numerator[i] = 0
            require(min(numerator) >= 0, ("skew convexity", order, row))


def relative_outer_ratio(cap7):
    d_lo, d_hi = log_bounds(1+6*cap7)
    rho_lo = log_bounds(1+Q(49, 8)*cap7)[0]/d_hi
    rho_hi = 1 + log_bounds(1/(1-cap7/8))[1]/d_lo
    require(rho_lo <= REFERENCE_RATIO <= rho_hi, "canonical outer ratio")
    return rho_lo, rho_hi


@lru_cache(maxsize=None)
def gaussian_exponents(order):
    return tuple(Q(k*(order-k), order-1) for k in range(order+1))


@lru_cache(maxsize=None)
def exponent_bounds(cap7, cap8):
    rho = relative_outer_ratio(cap7)
    result = {}
    for order, cap, ratio in ((7, cap7, (1, 1)), (8, cap8, rho)):
        m = order-1
        d_lo, d_hi = log_bounds(1+m*cap)
        bounds = []
        reference = 1 if order == 7 else REFERENCE_RATIO
        for k, gaussian in enumerate(gaussian_exponents(order)[1:-1], 1):
            r = min(k, order-k)
            mean_lo = sum(log_bounds((1+(m-h)*cap)/(1+h*cap))[0] for h in range(r))/d_hi
            skew = 2*weighted_log_upper(skew_weights(order, r), cap)/d_lo
            lo, hi = rounded_interval(ratio[0]*(mean_lo-skew), ratio[1]*(gaussian+skew))
            require(0 < lo <= reference*gaussian <= hi, ("canonical exponent", order, k, cap))
            bounds.append((lo, hi))
        result[order] = tuple(bounds)
    return result


def coefficient_box(core, u, exponents):
    result = [ONE]
    for c, (lo, hi) in zip(core[1:-1], exponents, strict=True):
        bounds = (power_bounds(u[0], hi)[0], power_bounds(u[1], lo)[1])
        result.append(interval(*multiply_interval(c, bounds)))
    return ((ONE, tuple(result + [ONE])),)


def mean_defect(order, a, cap):
    m = order-1
    lo, hi = log_bounds(1+m*cap)
    numerator = (1+Q(m, 2)*a*cap + m*cap/2)/(1-Q(m, 2)*a*cap + m*cap/2)
    n_lo, n_hi = log_bounds(numerator)
    return a-n_hi/lo, a-n_lo/hi


def reflected_vector(order, values):
    middle = (0,) if order % 2 == 0 else ()
    return (0, *values, *middle, *(-x for x in reversed(values)), 0)


@lru_cache(maxsize=None)
def grouped_envelope(order, cap):
    denom = log_bounds(1+(order-1)*cap)[0]
    skew = [2*weighted_log_upper(weights, cap)/denom
            for weights in (skew_weights(order, 3), *link_weights(order))]
    if order == 7:
        directions = ((Q(5, 2), Q(5, 2), 1), (Q(3, 2), Q(1, 2), 0),
                      (Q(-1, 2), Q(-1, 2), 0))
        h = mean_defect(7, Q(2, 3), cap)
        j = Q(4, 5)*h[1] - mean_defect(7, Q(1, 3), cap)[0]
        defects = (((0, 0, 1, Q(9, 5), Q(9, 5), 1, 0, 0), h[1]),
                   ((0, 0, 0, -1, -1, 0, 0, 0), j))
    else:
        directions = ((Q(7, 5), Q(8, 5), 1), (Q(1, 5), Q(3, 35), 0),
                      (Q(-1, 5), Q(-8, 35), 0))
        h = mean_defect(8, Q(5, 7), cap)
        j = h[1] - mean_defect(8, Q(3, 7), cap)[0]
        k = Q(2, 5)*h[1] - mean_defect(8, Q(1, 7), cap)[0]
        defects = (((0, 0, 1, 2, Q(12, 5), 2, 1, 0, 0), h[1]),
                   ((0, 0, 0, -1, -1, -1, 0, 0, 0), j),
                   ((0, 0, 0, 0, -1, 0, 0, 0, 0), k))
    groups = tuple((interval(-bound, bound), reflected_vector(order, values))
                   for bound, values in zip(skew, directions, strict=True))
    groups += tuple((interval(Q(0), bound), vector) for vector, bound in defects)
    return tuple((p, tuple(rational(v) for v in vector)) for p, vector in groups)


def coefficient_model(order, core, u, caps):
    ratio = interval(*relative_outer_ratio(caps[0]))
    reference_ratio = Q(1) if order == 7 else REFERENCE_RATIO
    midpoint = sum(u)/2
    gaussian = gaussian_exponents(order)
    reference = tuple(interval(*multiply_interval(c, power_bounds(midpoint, e*reference_ratio)))
                      if e else ONE for c, e in zip(core, gaussian, strict=True))
    mean = log(interval(u[0]/midpoint, u[1]/midpoint))
    outer_log = -log(interval(*u))
    if order == 8:
        mean = rational(REFERENCE_RATIO)*mean+(rational(REFERENCE_RATIO)-ratio)*outer_log
        outer_log = ratio*outer_log
    pieces = tuple((p*outer_log, vector)
                   for p, vector in grouped_envelope(order, caps[order-7]))
    pieces += ((mean, tuple(rational(e) for e in gaussian)),)
    errors = []
    for k, c in enumerate(reference):
        radius = ZERO
        for parameter, vector in pieces:
            term = parameter*vector[k]
            radius = radius+I(0, max(-term.lo, term.hi))
        require(radius.hi < GRID, ("exponential radius", order, u))
        error = scalar(1, 6)*radius**3/(ONE-radius)
        errors.append(c*(ONE+error.symmetric_error()))
    products = []
    for i, (first, a) in enumerate(pieces):
        for j, (second, b) in enumerate(pieces[:i+1]):
            p = scalar(1, 2)*first.square() if i == j else first*second
            products.append((p, tuple(v*w for v, w in zip(a, b, strict=True))))
    return ((ONE, tuple(errors)),) + tuple(
        (p, tuple(c*v for c, v in zip(reference, vector, strict=True)))
        for p, vector in (*pieces, *products)
    )


def verify_high():
    verify_binomial_anchor(COSINES)
    count = 0
    for order in ORDERS:
        caps = coefficient_bounds(order, (1,), STRENGTH_LIMIT)
        box = [(Q(1), Q(1))] + [(Q(0), c) for c in caps[1:-1]] + [(Q(1), Q(1))]
        verify_annulus(box, order)
        for cosine in COSINES:
            cells = radial_replay(order, box, (cosine, cosine))
            require(cells is not None, ("high-strength ray", order, cosine))
            count += cells
    return count


def ray_at(left, right, index, u):
    first, last = Q(RAYS[left][index], 100000), Q(RAYS[right][index], 100000)
    return first+(last-first)*(100*u-left)/(right-left)


def verify_knot(order, index, left, right, rows, core):
    paired = any(a <= left < right <= b for a, b in PAIRED_BANDS)
    pending = [((Q(left, 100), Q(right, 100)), 0)]
    count = 0
    while pending:
        u, depth = pending.pop()
        caps = reciprocal_caps(rows, u[0])
        if paired:
            model = coefficient_model(order, core, u, caps)
        else:
            model = coefficient_box(core, u, exponent_bounds(*caps)[order])
        bounds = tuple(sorted(ray_at(left, right, index, value) for value in u))
        cells = radial(order, interval(*bounds), model)
        if cells is not None:
            count += cells
            continue
        require(depth < 16, ("unresolved ray", order, index, u, caps, paired))
        midpoint = sum(u)/2
        pending.extend((((u[0], midpoint), depth+1), ((midpoint, u[1]), depth+1)))
    return count


def main():
    verify_skew_convexity()
    print(f"Certified {verify_high()} high-strength radial cells.", flush=True)
    verify_binomial_anchor(tuple(Q(v, 100000) for v in RAYS[0]))
    for rays in RAYS.values():
        bounds = (100000, *rays, -100000)
        require(all(a > b for a, b in zip(bounds, bounds[1:])), "ordered rays")
    rows = tuple(possible_rows())
    require(all(0 < x8 <= x7 <= 1 for x7, x8, _ in rows), "skew convexity range")
    cores = {j: private_core(j, (1,), ()) for j in ORDERS}
    for order, core in cores.items():
        r = Q(1, 10)
        require(sum(c[1]*r**k for k, c in enumerate(core) if k) < 1, ("inner radius", order))
        require(sum(c[1]*r**(order-k) for k, c in enumerate(core[:-1])) < 1,
                ("outer radius", order))
    lo, hi = (0, 100) if len(argv) == 1 else tuple(map(int, argv[1:]))
    knots = tuple(sorted(RAYS))
    require(lo in RAYS and hi in RAYS and lo < hi, "replay band")
    total = 0
    for left, right in zip(knots, knots[1:]):
        if not lo <= left < right <= hi:
            continue
        count = sum(verify_knot(j, i, left, right, rows, cores[j])
                    for j in ORDERS for i in range(6))
        print(left, right, count, flush=True)
        total += count
    print(f"Certified {total} low-strength cells using {len(rows)} candidate rows.", flush=True)


if __name__ == "__main__":
    main()
