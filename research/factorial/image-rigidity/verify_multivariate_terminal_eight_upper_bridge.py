#!/usr/bin/env python3
"""Exact Cayley projections for fully shared terminal eight, .88 <= u <= .95."""

from fractions import Fraction as Q
from functools import lru_cache
from math import comb
from sys import argv

from verify_multivariate_terminal_eight_private_linear import (
    ORDERS, REFERENCE_RATIO, gaussian_exponents, grouped_envelope,
    relative_outer_ratio, verify_skew_convexity,
)
from verify_multivariate_terminal_seven_upper_bridge import (
    bisect_u, cayley_matrix, height,
)
from verify_multivariate_terminal_six_intermediate_cayley import (
    GRID, I, ONE, ZERO, bisect, interval, log, projections, rational, require, scalar,
)
from verify_multivariate_terminal_six_joint_ray_bands import power_bounds
from verify_multivariate_terminal_six_marked_row_tail import projected_interval


RAY_SCALE = 100000
DISK = Q(2)
RAYS = {
    88: (118237, 98575, 75103, 59394, 40138, 26266),
    89: (117494, 98111, 74836, 59231, 40050, 26218),
    90: (116764, 97654, 74572, 59069, 39963, 26170),
    91: (116047, 97204, 74311, 58909, 39877, 26122),
    92: (115342, 96759, 74053, 58750, 39792, 26076),
    93: (114649, 96320, 73799, 58593, 39300, 26029),
    94: (113967, 95887, 73547, 58438, 38500, 25983),
    95: (113297, 95459, 73298, 58284, 38000, 25937),
}


def reciprocal_cap(u):
    lower = 1
    while u*u > Q(7*lower+1, 7*lower+7):
        lower += 1
    return min((u**-2-1)/6, max(Q(1, 7*lower+1), (1-u)/(7*u-1)))


def nominal(order, u):
    ratio = 1 if order == 7 else REFERENCE_RATIO
    return tuple(scalar(comb(order,k))*interval(*power_bounds(u, ratio*exponent))
                 if exponent else ONE
                 for k, exponent in enumerate(gaussian_exponents(order)))


def log_pieces(order, u):
    c7 = reciprocal_cap(u[0])
    cap = c7 if order == 7 else 7*c7/(8-c7)
    rho = interval(*relative_outer_ratio(c7))
    midpoint = sum(u)/2
    outer_log = -log(interval(*u))
    mean = log(interval(u[0]/midpoint, u[1]/midpoint))
    if order == 8:
        ratio = rational(REFERENCE_RATIO)
        mean = ratio*mean+(ratio-rho)*outer_log
        outer_log = rho*outer_log
    pieces = tuple((parameter*outer_log, vector)
                   for parameter, vector in grouped_envelope(order, cap))
    pieces += ((mean, tuple(rational(e) for e in gaussian_exponents(order))),)
    return nominal(order, midpoint), pieces


@lru_cache(maxsize=8192)
def cayley_model(order, u):
    p = (1-u[1], 1-u[0])
    lam = interval(power_bounds(p[0], Q(1, 2))[0], power_bounds(p[1], Q(1, 2))[1])
    reference, pieces = log_pieces(order, u)
    matrix = cayley_matrix(order)
    scales = tuple(lam**degree for degree in range(order+1))

    def transform(vector):
        result = []
        for degree, weights in enumerate(matrix):
            value = ZERO
            for k in range((order+1)//2):
                opposite = order-k
                pair = scalar(weights[k])*vector[k]+scalar(weights[opposite])*vector[opposite]
                value = value+reference[k]*pair
            if order%2 == 0:
                k = order//2
                value = value+reference[k]*scalar(weights[k])*vector[k]
            result.append((degree, value*scales[degree]))
        return tuple(result)

    base = list(transform((ONE,)*(order+1)))
    errors = []
    for k, coefficient in enumerate(reference):
        radius = ZERO
        for parameter, vector in pieces:
            term = parameter*vector[k]
            radius = radius+I(0, max(-term.lo, term.hi))
        require(radius.hi < GRID, ("exponential radius", order, u))
        error = scalar(1, 6)*radius**3/(ONE-radius)
        errors.append(coefficient*error.symmetric_error())
    for degree, weights in enumerate(matrix):
        error = ZERO
        for a, b in zip(errors, weights, strict=True):
            error = error+a*scalar(b)
        base[degree] = degree, base[degree][1]+error*scales[degree]

    products = []
    for i, (first, v) in enumerate(pieces):
        for j, (second, w) in enumerate(pieces[:i+1]):
            parameter = first*second if i != j else scalar(1, 2)*first.square()
            products.append((parameter, tuple(a*b for a, b in zip(v, w, strict=True))))
    groups = tuple((parameter, transform(vector))
                   for parameter, vector in (*pieces, *products))
    return tuple(base), groups


@lru_cache(maxsize=8192)
def disk_check(order, u):
    base, groups = cayley_model(order, u)
    total = []
    for degree, coefficient in base:
        for parameter, terms in groups:
            coefficient = coefficient+parameter*terms[degree][1]
        total.append(coefficient)
    return total[order].lo*DISK**order > sum(
        max(abs(c.lo), abs(c.hi))*DISK**degree for degree, c in enumerate(total[:-1])
    )


def ray_at(left, right, index, u):
    a, b = RAYS[left][index], RAYS[right][index]
    return Q(a, RAY_SCALE)+Q(b-a, RAY_SCALE)*(100*u-left)/(right-left)


def verify_ray(order, left, right, index):
    pending = [((Q(left, 100), Q(right, 100)), interval(Q(0), DISK), 0)]
    count = 0
    while pending:
        u, x, depth = pending.pop()
        if not disk_check(order, u):
            require(depth < 28, ("scaled root disk", order, u))
            pending.extend((half, x, depth+1) for half in bisect_u(u))
            continue
        p = interval(1-u[1], 1-u[0])
        values = tuple(sorted(ray_at(left, right, index, t) for t in u))
        y = height(p, x, interval(*values))
        base, groups = cayley_model(order, u)
        for direction in projections(order, x, y):
            total = projected_interval(base, direction)
            for parameter, group in groups:
                total = total+parameter*projected_interval(group, direction)
            if total.lo > 0 or total.hi < 0:
                count += 1
                break
        else:
            require(depth < 34, ("unresolved Cayley ray", order, left, index, u, x.shown()))
            if 120*(u[1]-u[0]) > Q(x.hi-x.lo, GRID):
                pending.extend((half, x, depth+1) for half in bisect_u(u))
            else:
                pending.extend((u, half, depth+1) for half in bisect(x))
    return count


def verify_anchor():
    for rays in RAYS.values():
        endpoints = (2*RAY_SCALE, *rays, 0)
        require(all(a > b for a, b in zip(endpoints, endpoints[1:])), "ray order")
    u = Q(22, 25)
    require(DISK**2*(1-u) < 1, "Cayley disk in the left half-plane")
    a,b,c = nominal(7,u)[1:4]
    d,e,f,g = nominal(8,u)[1:5]
    abscissae = (scalar(2),) + tuple(
        rational(2*((1-u)*v*v-1)/((1-u)*v*v+1))
        for v in (Q(raw, RAY_SCALE) for raw in RAYS[88])
    ) + (scalar(-2),)

    def septic(x):
        return x**3-x.square()-scalar(2)*x+ONE + a*(x.square()-x-ONE) + b*(x-ONE)+c

    def octic(x):
        return (x**4-scalar(4)*x.square()+scalar(2)+d*(x**3-scalar(3)*x)
                + e*(x.square()-scalar(2)) + f*x+g)

    for polynomial, signs in (
        (septic, (1, 1, -1, -1, 1, 1, -1, -1)),
        (octic, (1, -1, -1, 1, 1, -1, -1, 1)),
    ):
        for x, sign in zip(abscissae, signs, strict=True):
            require((scalar(sign)*polynomial(x)).lo > 0,
                    ("Gaussian anchor sector", polynomial.__name__, x.shown()))


def main():
    verify_skew_convexity()
    verify_anchor()
    first, last = (88, 95) if len(argv) == 1 else tuple(map(int, argv[1:]))
    require(88 <= first < last <= 95, "replay band")
    total = 0
    for left in range(first, last):
        counts = [[verify_ray(j, left, left+1, i) for i in range(6)] for j in ORDERS]
        total += sum(sum(row) for row in counts)
        print(left, left+1, counts, flush=True)
    print(f"Certified {total} Cayley projection cells.", flush=True)


if __name__ == "__main__":
    main()
