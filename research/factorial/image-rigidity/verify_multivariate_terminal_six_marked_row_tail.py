#!/usr/bin/env python3
"""Exact marked-row continuation for a residual order-five cap of 1/8."""

from fractions import Fraction as Q
from functools import lru_cache
from sys import argv

from verify_multivariate_terminal_six_intermediate_cayley import (
    GRID, I, ONE, bisect, exception_contributions, interval, linear_form,
    projections, rational, require, scalar,
)
from verify_multivariate_terminal_six_marked_row_extension import (
    RAYS as EARLIER_RAYS, ray_domain, verify_extension_constants,
)


# Residual deficit in thousandths followed by four physical ray cosines.
RAYS = {
    161: EARLIER_RAYS[161],
    180: (-233, -488, -756, -891),
    200: (-202, -467, -746, -887),
    220: (-171, -445, -735, -883),
    240: (-140, -424, -724, -879),
    250: (-124, -413, -719, -876),
    265: (-101, -397, -711, -872),
    275: (-86, -386, -705, -870),
    285: (-70, -375, -700, -867),
    295: (-55, -364, -694, -864),
    310: (-32, -348, -686, -860),
    325: (-9, -331, -677, -856),
    340: (13, -315, -669, -852),
    360: (44, -293, -657, -847),
    380: (74, -271, -645, -842),
    400: (103, -250, -633, -836),
    420: (133, -228, -621, -830),
    440: (163, -206, -609, -824),
    460: (192, -184, -596, -817),
    480: (221, -162, -584, -808),
    500: (250, -140, -571, -802),
    520: (278, -118, -558, -796),
    540: (307, -97, -545, -789),
    560: (335, -75, -532, -783),
    580: (363, -53, -518, -776),
    600: (390, -31, -504, -770),
    620: (418, -9, -490, -763),
    640: (445, 13, -476, -756),
    660: (471, 35, -462, -749),
    680: (498, 57, -447, -742),
    700: (524, 79, -432, -734),
    720: (549, 101, -417, -727),
    740: (574, 123, -402, -719),
    760: (599, 145, -386, -711),
    780: (623, 167, -370, -703),
    800: (647, 189, -353, -695),
    820: (670, 211, -336, -686),
    840: (693, 233, -319, -677),
    860: (715, 255, -301, -668),
}

CAPS = (
    (161, 250, Q(1, 11)),
    (250, 275, Q(11, 100)),
    (275, 295, Q(117, 1000)),
    (295, 460, Q(1, 8)),
    (460, 520, Q(1, 8)),
    (520, 580, Q(1, 8)),
    (580, 640, Q(1, 8)),
    (640, 700, Q(1, 8)),
    (700, 760, Q(1, 8)),
    (760, 820, Q(1, 8)),
    (820, 860, Q(1, 8)),
)


def angular_imaginary_coordinate(b, x):
    square = x.square()
    lower = max(scalar(1, 15).lo, (b*(ONE-square)).lo)
    y = I(lower,b.hi)
    for _ in range(3):
        update=b-b*square/(ONE+b*y)
        y=I(max(y.lo,update.lo),min(y.hi,update.hi))
    return y


def small_ray_coordinates(left, right, index, p, x):
    first, last = RAYS[left][index], RAYS[right][index]
    cosine = scalar(first, 1000) + scalar(last-first, right-left)*(p-scalar(left, 1000))
    b = ((ONE + cosine) / (ONE - cosine)).sqrt()
    return (angular_imaginary_coordinate(b, x),)


def ray_coordinates(left, right, index, p, x):
    first, last = RAYS[left][index], RAYS[right][index]
    cosine = scalar(first, 1000) + scalar(last-first, right-left)*(p-scalar(left, 1000))
    sine2 = ONE-cosine.square()
    discriminant = ONE/sine2-x.square()
    if discriminant.hi < 0:
        return ()
    height = I(max(0,discriminant.lo),discriminant.hi).sqrt()
    center = cosine/sine2.sqrt()
    return tuple(I(max(value.lo,scalar(1,15).lo),value.hi)
                 for value in (center+height,center-height)
                 if value.hi >= scalar(1,15).lo)


def projected_interval(terms, directions):
    lower = upper = 0
    for degree, coefficient in terms:
        value = directions[degree]
        a, b, c, d = coefficient.lo, coefficient.hi, value.lo, value.hi
        products = a*c, a*d, b*c, b*d
        lower += min(products)
        upper += max(products)
    return I(lower//GRID,-(-upper//GRID))


def verify_ray(order, index, left, right, cap):
    rho, cap, s_cap = ray_domain(order, cap)
    if left < 460:
        coordinates = small_ray_coordinates
        dot = linear_form
    else:
        coordinates = ray_coordinates
        dot = projected_interval
    p = interval(Q(left, 1000), Q(right, 1000))
    cos_upper = Q(max(RAYS[left][index], RAYS[right][index]), 1000)
    if cos_upper <= Q(1,5):
        limit = 63/(65-16*cos_upper)
        radial = interval(-limit, limit)
    else:
        limit = (ONE/(ONE-rational(cos_upper).square()).sqrt()).hi
        radial = I(-limit, limit)
    initial = p, rho, radial
    widths = [max(value.hi-value.lo, 1) for value in initial]

    @lru_cache(maxsize=8192)
    def coefficients(plo, phi, rlo, rhi):
        return exception_contributions(order, I(plo, phi), I(rlo, rhi), cap, s_cap)

    stack = [(initial, 0)]
    count = 0
    while stack:
        box, depth = stack.pop()
        p, r, x = box
        count += 1
        nominal, families = coefficients(p.lo, p.hi, r.lo, r.hi)
        for y in coordinates(left, right, index, p, x):
            for directions in projections(order, x, y):
                total = dot(nominal, directions)
                for parameter, contribution in families:
                    total = total + parameter*dot(contribution, directions)
                if total.lo > 0 or total.hi < 0:
                    break
            else:
                break
        else:
            continue
        require(depth < 39, (
            "unresolved extension ray", order, index, left, right,
            [value.shown() for value in box],
        ))
        weights = (120, 1 if depth >= 16 else 0, 32)
        k = max(range(3),
                key=lambda i: (box[i].hi-box[i].lo)*weights[i]*GRID//widths[i])
        for half in bisect(box[k]):
            stack.append((box[:k]+(half,)+box[k+1:], depth+1))
    return count


def verify_tail_constants():
    verify_extension_constants()
    require(all(all(800 >= first > second >= -900
                    for first, second in zip(rays, rays[1:]))
                for rays in RAYS.values()), "ordered rays")
    require(Q(16)**2*(1-Q(9, 10)**2) > (Q(1, 15)*(65+16*Q(9, 10)))**2,
            "uniform positive Cayley imaginary coordinate")
    require(Q(1, 5)*(1+Q(1, 8)**2) < 2*Q(1, 8), "radial monotonicity")
    require(Q(149, 1000)**5 < Q(7, 11)**21, "residual deficit when curvature is below four")
    for _, right, cap in CAPS:
        if cap < Q(1, 8):
            u = 1-Q(right, 1000)
            require(Q(1, 11) <= cap and (1-u)/(5*u-1) <= cap, "integer residual cap")
        rho, _, _ = ray_domain(6, cap)
        require(rho.lo <= scalar(25, 24).lo <= scalar(25, 24).hi <= rho.hi,
                "common Gaussian anchor")


def main():
    verify_tail_constants()
    for low, high, cap in CAPS:
        if len(argv) > 1 and low not in map(int, argv[1:]):
            continue
        knots = [knot for knot in RAYS if low <= knot <= high]
        for left, right in zip(knots, knots[1:]):
            sizes = [verify_ray(order, ray, left, right, cap)
                     for order in (5, 6) for ray in range(4)]
            print(f"Marked residual deficit [{left}, {right}]/1000: {sizes}", flush=True)


if __name__ == "__main__":
    main()
