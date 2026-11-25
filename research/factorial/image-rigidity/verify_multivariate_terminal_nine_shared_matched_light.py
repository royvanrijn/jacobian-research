#!/usr/bin/env python3
"""Exact matched-skew ray proof for fully shared terminal nine."""

from fractions import Fraction as Q
from functools import lru_cache
from math import comb, factorial

from verify_multivariate_terminal_eight_paired_bridge import radial
from verify_multivariate_terminal_nine_private_linear import ratio_bounds, shape
from verify_multivariate_terminal_nine_private_rows import (
    ORDERS, REFERENCE_RATIO, outer_coefficient_box, outer_exponents,
    verify_private_anchor, verify_skew_convexity,
)
from verify_multivariate_terminal_six_intermediate_cayley import (
    GRID, I, ONE, ZERO, interval, log, rational, scalar,
)
from verify_multivariate_terminal_six_joint_ray_bands import power_bounds, log_bounds
from verify_multivariate_terminal_six_private_angles import require


CAPS = (Q(1, 8), Q(8, 71))
RAY_SCALE = 10**6
ETA = Q(1, 206)
SKEW_TRANSFER_CENTER = REFERENCE_RATIO**2 * (1 + ETA**2)
SKEW_TRANSFER_WIDTH = 2 * REFERENCE_RATIO**2 * ETA
ANNULUS = Q(1, 12), Q(12)

RAYS = {
    9: (868038, 653209, 283037, 30643, -407679, -655760, -888401),
    10: (856139, 637745, 265507, -60638, -418686, -662616, -890750),
    12: (833940, 605986, 226298, -90725, -439830, -675703, -895211),
    14: (808687, 573332, 192214, -101819, -459958, -677362, -903434),
    16: (780347, 559070, 158669, -112775, -479207, -689488, -907238),
    18: (751689, 486507, 121202, -123686, -503484, -690959, -910851),
    20: (720155, 471647, 88719, -151560, -521104, -702405, -917817),
    22: (688728, 456322, 52388, -162243, -538115, -713410, -920983),
    24: (656122, 459564, 20885, -189262, -554577, -724018, -924022),
    26: (622429, 424116, -14365, -199943, -570541, -725381, -926947),
    28: (587732, 407151, -40881, -226300, -590966, -735619, -932724),
    30: (552101, 370909, -71190, -252183, -605890, -745549, -935338),
    32: (515599, 352865, -101144, -277622, -620430, -747120, -937868),
    34: (478283, 278884, -130762, -288433, -634612, -756802, -940322),
    36: (440203, 260203, -160061, -313477, -652724, -766234, -942705),
    38: (399743, 222732, -189057, -338154, -666111, -775433, -945022),
    40: (360271, 202776, -217764, -349453, -679207, -777388, -947277),
    42: (320163, 164541, -246195, -373873, -695825, -786441, -949475),
    44: (279451, 143173, -274364, -397983, -708248, -795298, -951618),
}

ANCHORS = {8: ((8545, 5195), (2182, 9759), (-4930, 8700), (-9397, 3419)), 9: ((8831, 4693), (3460, 9382), (-3069, 9518), (-8120, 5837), (-10000, 0))}


def exp_small(value):
    size = I(0, max(abs(value.lo), abs(value.hi)))
    require(size.hi < GRID, ("small exponential", value.shown()))
    term = result = ONE
    for j in range(1, 13):
        term = term * value * scalar(1, j)
        result += term
    tail = size**13 * scalar(1, factorial(13)) / (ONE-size)
    return result + tail.symmetric_error()


@lru_cache(maxsize=None)
def constants(order):
    cap = CAPS[order-8]
    m = order - 1
    T = rational(m*cap/(2+m*cap))
    ls = log_bounds(1+m*cap)
    s = interval(ls[0]/2, ls[1]/2)
    require(T.square().hi < GRID//2, "moment quotient monotonicity range")
    skew_bound = T.square() / (scalar(2)*s)
    return T, s, skew_bound


@lru_cache(maxsize=None)
def leading_skew_direction(order):
    _, skew, _, _, _ = shape(order)
    T, _, _ = constants(order)
    weights = tuple(T**(2*l-2)*scalar(1, 2*l) for l in (2, 3))
    second, third = weights
    direction = tuple(a + second*b + third*c for a, b, c in zip(*skew, strict=True))
    return direction, weights


@lru_cache(maxsize=8192)
def coefficient_model(order, means, signed):
    midpoint = sum(means)/2
    skew_midpoint = sum(signed)/2
    ref_ratio = Q(1) if order == 8 else REFERENCE_RATIO
    T, s, skew_bound = constants(order)
    Pmax = constants(8)[2]
    U8 = -log(interval(*means))
    leading = interval(*signed)*Pmax*U8
    nominal = rational(skew_midpoint)*Pmax*(-log(rational(midpoint)))
    outer_mean = rational(ref_ratio)*log(interval(means[0]/midpoint, means[1]/midpoint))
    strength = U8
    if order == 9:
        rho = ratio_bounds(CAPS[0])
        outer_mean += (rational(ref_ratio)-rho)*U8
        strength = rho*U8
        transfer = rational(SKEW_TRANSFER_CENTER)
        leading = transfer*leading + (rational(SKEW_TRANSFER_WIDTH)*Pmax*U8).symmetric_error()
        nominal = transfer*nominal
    gaussian, skew, mean, skew_sizes, distances = shape(order)
    direction, weights = leading_skew_direction(order)
    reference = []
    for k, v in enumerate(direction):
        if k == 0 or k == order:
            reference.append(ONE)
            continue
        exponent = ref_ratio*Q(k*(order-k), order-1)
        value = scalar(comb(order,k))*interval(*power_bounds(midpoint, exponent))
        reference.append(value*exp_small(nominal*v))
    pieces = [(outer_mean, gaussian), (leading-nominal, direction)]
    for vector, weight in zip(skew[1:], weights, strict=True):
        bound = weight*skew_bound*strength
        pieces.append((bound.symmetric_error(), vector))
    lower = I(leading.lo, leading.lo) if leading.lo > 0 else ZERO
    upper_strength = I(strength.hi, strength.hi)
    for j, vector in enumerate(mean,1):
        lo = (scalar(2)*lower/upper_strength)**(2*j)*upper_strength*scalar(1, 2*j+1)
        hi = (T**(2*j+1))*strength/(scalar(2*j+1)*s)
        if lo.lo > hi.hi:
            # The necessary lower and upper bounds are disjoint.
            return None
        pieces.append((I(lo.lo, hi.hi), vector))
    skew_tail = T**7*strength/(scalar(8)*(ONE-T.square()))
    defect_tail = T**6*strength/(scalar(7)*(ONE-T.square()))
    errors = []
    for k, c in enumerate(reference):
        tail = skew_sizes[k]*skew_tail+scalar(distances[k])*defect_tail
        size = tail
        for parameter, vector in pieces:
            term = parameter*vector[k]
            size += I(0, max(-term.lo, term.hi))
        if size.hi >= GRID:
            # This rectangle needs subdivision before using the exponential bound.
            return ()
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


def verify_band(order, left, right):
    means = Q(left, 100), Q(right, 100)
    stack = [(means, (Q(0), Q(1)), tuple(range(7)), 0, 0)]
    total = 0
    while stack:
        u, p, rays, du, dp = stack.pop()
        model = coefficient_model(order, u, p)
        if model is None:
            continue
        bad = []
        for ray in rays:
            a, b = RAYS[left][ray], RAYS[right][ray]
            cosines = tuple(sorted(Q(a, RAY_SCALE)+Q(b-a,RAY_SCALE)*(100*v-left)/(right-left)
                                   for v in u))
            cells = radial(order, interval(*cosines), model) if model else None
            if cells is None:
                bad.append(ray)
            else:
                total += cells
        if not bad:
            continue
        require(max(du, dp) < 16, ("unresolved shared ray", order, left, right, u, p, bad))
        if du < 2*dp+2:
            center = sum(u)/2
            stack.extend((((u[0], center), p, bad, du+1, dp),
                          ((center, u[1]), p, bad, du+1, dp)))
        else:
            center = sum(p)/2
            stack.extend(((u, (p[0], center), bad, du, dp+1),
                          (u, (center, p[1]), bad, du, dp+1)))
    return total


def verify_setup():
    verify_skew_convexity()
    left, right = Q(min(RAYS),100), Q(max(RAYS),100)
    for values in RAYS.values():
        rays = (RAY_SCALE, *values, -RAY_SCALE)
        require(len(values) == 7 and all(a > b for a,b in zip(rays,rays[1:])), "ordered rays")
    small, large = ANNULUS
    for order in ORDERS:
        raw = tuple((Q(comb(order,k)),)*2 for k in range(order+1))
        exponents = outer_exponents(order, *CAPS)
        box = outer_coefficient_box(raw, (Q(0), right), exponents)
        require(sum(cap*small**k for k,(_,cap) in enumerate(box) if k) < 1, ("inner radius",order))
        require(sum(cap*large**k for k,(_,cap) in enumerate(box[:-1])) < large**order,
                ("outer radius",order))
        ratio = Q(1) if order == 8 else REFERENCE_RATIO
        anchor = [(Q(1),Q(1))]
        for k in range(1,order):
            lo,hi = power_bounds(left, ratio*Q(k*(order-k),order-1))
            anchor.append((comb(order,k)*lo, comb(order,k)*hi))
        anchor.append((Q(1),Q(1)))
        verify_private_anchor(order, anchor, ANCHORS[order],
                              tuple(Q(x,RAY_SCALE) for x in RAYS[min(RAYS)]))


def main():
    verify_setup()
    total = 0
    knots = tuple(sorted(RAYS))
    for left, right in zip(knots, knots[1:]):
        counts = tuple(verify_band(order, left, right) for order in ORDERS)
        total += sum(counts)
        print(left, right, counts, flush=True)
    print(f"Verified {total} radial cells.")


if __name__=="__main__":
    main()
