#!/usr/bin/env python3
"""Exact Cayley projection bridge for terminal seven on 0.88 <= u <= 0.95."""

from fractions import Fraction as Q
from functools import lru_cache
from math import comb
from sys import argv

from verify_multivariate_terminal_seven_paired_light import (
    ANCHOR_RATIO, envelope, exponents, nominal,
)
from verify_multivariate_terminal_six_intermediate_cayley import (
    GRID, I, ONE, ZERO, bisect, interval, log, projections, rational, require, scalar,
)
from verify_multivariate_terminal_six_joint_ray_bands import power_bounds
from verify_multivariate_terminal_six_marked_row_tail import projected_interval


RAYS = {
    88: (115085, 93554, 67073, 49876, 28192),
    89: (114393, 93151, 66863, 49759, 28140),
    90: (113713, 92752, 66655, 49642, 28087),
    91: (113044, 92359, 66449, 49527, 28036),
    92: (112386, 91970, 66246, 49412, 27985),
    93: (111739, 91587, 66045, 49299, 27934),
    94: (111102, 91208, 65846, 49187, 27884),
    95: (110476, 90834, 65648, 49076, 27834),
}
RAY_SCALE = 100000


def reciprocal_cap(u):
    lower = 1
    while u*u > Q(6*lower+1, 6*lower+6):
        lower += 1
    row_cap = max(Q(1, 6*lower+1), (1-u)/(6*u-1))
    return min((u**-2-1)/5, row_cap)


@lru_cache(maxsize=None)
def log_envelope(cap, order):
    return envelope(cap, order)


def log_pieces(order, u):
    rho, variables = log_envelope(reciprocal_cap(u[0]), order)
    midpoint = sum(u) / 2
    upper_log = -log(interval(*u))
    mean = log(interval(u[0]/midpoint, u[1]/midpoint))
    if order == 7:
        mean = ANCHOR_RATIO*mean + (ANCHOR_RATIO-rho)*upper_log
        upper_log = rho*upper_log
    pieces = tuple((parameter*upper_log, vector) for vector, parameter in variables)
    pieces += ((mean, tuple(rational(e) for e in exponents(order))),)
    return nominal(order, midpoint), pieces


@lru_cache(maxsize=None)
def cayley_matrix(order):
    return tuple(tuple(
        sum((-1)**(k-i)*comb(k, i)*comb(order-k, degree-i)
            for i in range(max(0, degree-(order-k)), min(k, degree)+1))
        for k in range(order+1)
    ) for degree in range(order+1))


@lru_cache(maxsize=8192)
def cayley_model(order, u):
    p = (1-u[1], 1-u[0])
    lam = interval(power_bounds(p[0], Q(1, 2))[0], power_bounds(p[1], Q(1, 2))[1])
    reference, pieces = log_pieces(order, u)
    matrix = cayley_matrix(order)
    scales = tuple(lam**degree for degree in range(order+1))

    def transform(vector):
        result = []
        for degree, row in enumerate(matrix):
            value = ZERO
            for k in range((order+1)//2):
                reflected = order-k
                pair = scalar(row[k])*vector[k] + scalar(row[reflected])*vector[reflected]
                value = value + reference[k]*pair
            if order % 2 == 0:
                k = order//2
                value = value + reference[k]*scalar(row[k])*vector[k]
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
    for degree, row in enumerate(matrix):
        error = ZERO
        for a, b in zip(errors, row, strict=True):
            error = error + a*scalar(b)
        base[degree] = degree, base[degree][1] + error*scales[degree]

    squares = []
    for i, (first, v) in enumerate(pieces):
        for j, (second, w) in enumerate(pieces[:i+1]):
            parameter = first*second if i != j else scalar(1, 2)*first.square()
            squares.append((parameter, tuple(a*b for a, b in zip(v, w, strict=True))))
    model = tuple((parameter, transform(vector)) for parameter, vector in (*pieces, *squares))
    return tuple(base), model


@lru_cache(maxsize=8192)
def disk_check(order, u):
    base, model = cayley_model(order, u)
    total = []
    for degree, nominal_bound in base:
        coefficient = nominal_bound
        for parameter, group in model:
            coefficient = coefficient + parameter*group[degree][1]
        total.append(coefficient)
    return total[order].lo*2**order > sum(
        max(abs(coefficient.lo), abs(coefficient.hi))*2**degree
        for degree, coefficient in enumerate(total[:-1])
    )


def height(p, x, b):
    num = p*b*x.square()
    y = I(0, b.hi)
    for _ in range(3):
        update = b - num/(ONE+p*b*y)
        y = I(max(y.lo, update.lo), min(y.hi, update.hi))
    return y


def ray_at(left, right, index, u):
    a, b = RAYS[left][index], RAYS[right][index]
    return Q(a, RAY_SCALE) + Q(b-a, RAY_SCALE) * (100*u-left)/(right-left)


def bisect_u(u):
    middle = sum(u)/2
    return (u[0], middle), (middle, u[1])


def verify_ray(order, left, right, index):
    pending = [((Q(left, 100), Q(right, 100)), interval(Q(0), Q(2)), 0)]
    count = 0
    while pending:
        u, x, depth = pending.pop()
        if not disk_check(order, u):
            require(depth < 28, ("scaled root disk", order, u))
            pending.extend((half, x, depth+1) for half in bisect_u(u))
            continue
        p = interval(1-u[1], 1-u[0])
        bounds = tuple(sorted(ray_at(left, right, index, value) for value in u))
        y = height(p, x, interval(*bounds))
        base, groups = cayley_model(order, u)
        for direction in projections(order, x, y):
            total = projected_interval(base, direction)
            for parameter, contribution in groups:
                total = total + parameter*projected_interval(contribution, direction)
            if total.lo > 0 or total.hi < 0:
                count += 1
                break
        else:
            require(depth < 32, ("unresolved Cayley ray", order, left, index, u, x.shown()))
            if 120*(u[1]-u[0]) > Q(x.hi-x.lo, GRID):
                pending.extend((half, x, depth+1) for half in bisect_u(u))
            else:
                pending.extend((u, half, depth+1) for half in bisect(x))
    return count


def verify_anchor():
    for rays in RAYS.values():
        require(2*RAY_SCALE > rays[0] > rays[1] > rays[2] > rays[3] > rays[4] > 0,
                "ray order")
    u = Q(22, 25)
    require(4*(1-u) < 1, "Cayley disk remains in the left half-plane")
    a, b, c = nominal(6, u)[1:4]
    d, e, f = nominal(7, u)[1:4]
    abscissae = (scalar(2),) + tuple(
        rational(2*((1-u)*v*v-1)/((1-u)*v*v+1))
        for v in (Q(raw, RAY_SCALE) for raw in RAYS[88])
    ) + (scalar(-2),)

    def sextic(x):
        return x**3-scalar(3)*x + a*(x.square()-scalar(2)) + b*x + c

    def septic(x):
        return x**3-x.square()-scalar(2)*x+ONE + d*(x.square()-x-ONE) + e*(x-ONE)+f

    for polynomial, signs in (
        (sextic, (1, 1, -1, -1, 1, 1, -1)),
        (septic, (1, -1, -1, 1, 1, -1, -1)),
    ):
        for x, sign in zip(abscissae, signs, strict=True):
            require((scalar(sign)*polynomial(x)).lo > 0,
                    ("Gaussian anchor sector", polynomial.__name__, x.shown()))


def main():
    verify_anchor()
    leftmost, rightmost = (88, 95) if len(argv) == 1 else tuple(map(int, argv[1:]))
    require(88 <= leftmost < rightmost <= 95, "replay band")
    total = 0
    for left in range(leftmost, rightmost):
        counts = [[verify_ray(order, left, left+1, i) for i in range(5)]
                  for order in (6, 7)]
        total += sum(sum(row) for row in counts)
        print(left, left+1, counts, flush=True)
    print(f"Certified {total} Cayley projection cells.", flush=True)


if __name__ == "__main__":
    main()
