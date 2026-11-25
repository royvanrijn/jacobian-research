#!/usr/bin/env python3
"""Exact ray exclusions for a shared (1,2) row and a residual deficit up to .161."""

from fractions import Fraction as Q
from functools import lru_cache

from verify_multivariate_terminal_six_intermediate_cayley import (
    EXCEPTION_RAYS, GRID, I, ONE, bisect, exception_contributions,
    imaginary_coordinate, interval, linear_form, log, projections, rational,
    require, scalar, verify_scalar_constants,
)


# Residual deficit in thousandths, followed by the four physical ray cosines.
RAYS = {
    71: EXCEPTION_RAYS[71],
    83: (-386, -579, -832, -891),
    94: (-369, -579, -814, -900),
    105: (-351, -568, -797, -901),
    115: (-335, -558, -791, -901),
    125: (-320, -547, -785, -900),
    135: (-304, -536, -780, -899),
    143: (-291, -527, -775, -898),
    152: (-277, -516, -771, -897),
    161: (-263, -507, -766, -895),
}


# Each residual cap covers all actual profiles in its deficit band.
CAPS = (
    (71, 83, Q(1, 26)),
    (83, 105, Q(1, 21)),
    (105, 143, Q(1, 16)),
    (143, 161, Q(1, 11)),
)


def ray_parameter(left, right, index, p):
    first, last = RAYS[left][index], RAYS[right][index]
    cosine = scalar(first, 1000) + scalar(last-first, right-left)*(p-scalar(left, 1000))
    return ((ONE + cosine) / (ONE - cosine)).sqrt()


def ray_domain(order, cap):
    t = rational(cap)
    log_base = log(ONE + scalar(4)*t)
    if order == 5:
        return ONE, t, log_base*scalar(1, 2)
    lower = log(ONE + scalar(25, 6)*t) / log_base
    upper = ONE - log(ONE - t*scalar(1, 6)) / log_base
    upper_cap = scalar(5)*t/(scalar(6)-t)
    return I(lower.lo, upper.hi), upper_cap, scalar(1, 2)*log(ONE+scalar(5)*upper_cap)


def verify_ray(order, index, left, right, cap):
    rho, cap, s_cap = ray_domain(order, cap)
    p = interval(Q(left, 1000), Q(right, 1000))
    initial = p, rho, interval(Q(-21, 23), Q(21, 23))
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
        y = imaginary_coordinate(ONE, ray_parameter(left, right, index, p), x)
        for directions in projections(order, x, y):
            total = linear_form(nominal, directions)
            for parameter, contribution in families:
                total = total + parameter*linear_form(contribution, directions)
            if total.lo > 0 or total.hi < 0:
                break
        else:
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


def verify_extension_constants():
    verify_scalar_constants()
    require(all(all(-250 >= first > second > -1000
                    for first, second in zip(rays, rays[1:]))
                for rays in RAYS.values()), "ordered left-half-plane rays")
    require(Q(13, 20)**2 / Q(3, 5) > Q(839, 1000)**2, "target deficit")
    for _, right, cap in CAPS:
        u = 1-Q(right, 1000)
        minimum = (1/cap-1)/5
        require(minimum.denominator == 1 and minimum >= 2, "integral minimum row")
        require(u**2 > (5*(minimum-1)+1)/(5*(minimum-1)+5),
                "exclude smaller unit-gap rows")
        require((1-u)/(5*u-1) < cap, "larger-gap row cap")
        rho, _, _ = ray_domain(6, cap)
        require(rho.lo <= scalar(25, 24).lo <= scalar(25, 24).hi <= rho.hi,
                "common Gaussian anchor")


def main():
    verify_extension_constants()
    for low, high, cap in CAPS:
        knots = [knot for knot in RAYS if low <= knot <= high]
        for left, right in zip(knots, knots[1:]):
            sizes = [verify_ray(order, ray, left, right, cap)
                     for order in (5, 6) for ray in range(4)]
            print(f"Marked residual deficit [{left}, {right}]/1000: {sizes}", flush=True)


if __name__ == "__main__":
    main()
