#!/usr/bin/env python3
"""Replay terminal-seven shared profiles above the middle-band reciprocal caps."""

from fractions import Fraction as Q
from functools import lru_cache
from math import comb, prod

from flint import ctx

from verify_multivariate_terminal_seven_paired_light import (
    coefficient_model,
    envelope,
    radial,
)
from verify_multivariate_terminal_seven_retained_rows import (
    MAX_CURVATURE,
    ORDERS,
    RAY_FRACTIONS,
    ball,
    boxes,
    core_powers,
    nominal_cosines,
    power_bounds,
    rays,
    rows_above,
)
from verify_multivariate_terminal_seven_shared_high_curvature import exponent_bounds
from verify_multivariate_terminal_six_intermediate_cayley import interval, rational
from verify_multivariate_terminal_six_joint_ray_bands import INNER, OUTER, log_bounds
from verify_multivariate_terminal_six_private_angles import require


BANDS = (
    (Q(1, 16), Q(62, 100), Q(72, 100), 9),
    (Q(7, 200), Q(72, 100), Q(75, 100), 49),
)
PAIRED_RAY_FRACTIONS = tuple(dict.fromkeys((*RAY_FRACTIONS, *(Q(k, 16) for k in range(1, 16)))))


def row_outer_square(a, d):
    return prod(Q(6*a+h, 6*a+5*d+h) for h in range(1, d+1))


def enumerate_cores(cap, outer_lower):
    budget = min(MAX_CURVATURE, Q(12, 5)*log_bounds(1/outer_lower)[1])
    minimum = outer_lower**2
    rows = sorted(
        ((a, d, orientation), curvature, row_outer_square(a, d))
        for a, d, curvature in rows_above(cap, budget)
        if Q(d, 6*a+1) <= Q(1, 8) and row_outer_square(a, d) >= minimum
        for orientation in (-1, 1)
    )
    cores = []

    def collect(index, left, square, selected):
        if square < minimum:
            return
        if index == len(rows):
            reflected = tuple(sorted((a, d, -sign) for a, d, sign in selected))
            if selected and selected <= reflected:
                cores.append((selected, left))
            return

        row, curvature, contribution = rows[index]
        for multiplicity in range(1+int(left/curvature)):
            collect(index+1, left-multiplicity*curvature,
                    square*contribution**multiplicity, selected+(row,)*multiplicity)

    collect(0, budget, Q(1), ())
    return budget, cores


def annular_dominance(j, core_bounds):
    maxima = [comb(j, k)*hi for k, (_, hi) in enumerate(core_bounds)]
    return (
        sum(a*INNER**k for k, a in enumerate(maxima) if k) < 1
        and sum(a*OUTER**k for k, a in enumerate(maxima[:-1])) < OUTER**j
    )


@lru_cache(maxsize=None)
def paired_envelope(cap, order):
    return envelope(cap, order)


def paired_ray_replay(angles, core_arrays, tail, cap):
    models = {}
    for order in ORDERS:
        pieces = coefficient_model(order, tail, paired_envelope(cap, order))
        retained = core_arrays[order]
        models[order] = tuple(
            (parameter, tuple(a*b for a, b in zip(vector, retained, strict=True)))
            for parameter, vector in pieces
        )

    count = 0
    for left, right in zip(angles, angles[1:]):
        for position in PAIRED_RAY_FRACTIONS:
            proposal = float(left.mid())*float(position)+float(right.mid())*float(1-position)
            cosine = Q(round(100000*proposal), 100000)
            if not left > ball(cosine) > right:
                continue
            dyadic = rational(cosine)
            cells6 = radial(6, dyadic, models[6])
            if cells6 is None:
                continue
            cells7 = radial(7, dyadic, models[7])
            if cells7 is not None:
                count += cells6+cells7
                break
        else:
            return None
    return count


def verify_core(rows, retained_curvature, budget, cap, outer_lower, outer_upper):
    powers = {j: core_powers(j, rows) for j in ORDERS}
    coefficients = {
        j: [power_bounds(value, Q(1, j)) for value in powers[j]] for j in ORDERS
    }
    dyadic_cores = {j: tuple(interval(lo, hi) for lo, hi in coefficients[j]) for j in ORDERS}
    dominance = {j: annular_dominance(j, coefficients[j]) for j in ORDERS}

    core_lo, core_hi = power_bounds(powers[6][1]*powers[6][-2], Q(1, 12))
    lower = max(outer_lower/core_hi, power_bounds(Q(7, 12), Q(13, 12)*budget)[0])
    upper = min(Q(1), outer_upper/core_lo)
    if lower > upper:
        return 0, 0

    attenuation = {
        k: power_bounds(Q(90, 91), retained_curvature*k*k)[1] for k in range(1, 7)
    }
    pending = [((lower, upper), 0)]
    leaves = radial_cells = 0
    while pending:
        tail, depth = pending.pop()
        midpoint = sum(tail)/2
        proposed = 2000*(tail[0]**-2-1)/5
        cap6 = min(cap, Q(max(1, -(-proposed.numerator//proposed.denominator)), 2000))
        shape = exponent_bounds(cap6, 6*cap6/(7-cap6))
        cosines = nominal_cosines(powers, midpoint, shape)
        if cosines is not None:
            enclosure = boxes(coefficients, shape, tail, attenuation)
            annulus = all(
                dominance[j] or min(enclosure[j][1][0], enclosure[j][-2][0]) > j*Q(3, 8)
                for j in ORDERS
            )
            if annulus:
                count = rays(cosines, enclosure)
                if count is None and depth >= 4:
                    try:
                        count = paired_ray_replay(cosines, dyadic_cores, tail, cap6)
                    except ArithmeticError:
                        # The current tail interval can be too wide for the exponential remainder.
                        count = None
                if count is not None:
                    leaves += 1
                    radial_cells += count
                    continue

        require(depth < 15, ("unresolved retained array", rows, tail))
        pending.append(((tail[0], midpoint), depth+1))
        pending.append(((midpoint, tail[1]), depth+1))

    return leaves, radial_cells


def main():
    ctx.prec = 128
    require(Q(91, 90)**720 > Q(12, 7)**13, "tail attenuation from outer loss")
    for cap, lower, upper, expected_count in BANDS:
        budget, cores = enumerate_cores(cap, lower)
        require(len(cores) == expected_count, ("retained lists", lower, upper))
        leaves = radial_cells = 0
        for rows, left in cores:
            tail_count, radial_count = verify_core(rows, budget-left, left, cap, lower, upper)
            leaves += tail_count
            radial_cells += radial_count
        print(f"Verified [{lower},{upper}]: {len(cores)} retained lists, "
              f"{leaves} tail intervals and {radial_cells} radial cells.", flush=True)


if __name__ == "__main__":
    main()
