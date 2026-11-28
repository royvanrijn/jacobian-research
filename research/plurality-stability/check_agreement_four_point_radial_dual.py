#!/usr/bin/env python3
"""Exact rational interval check for the variable-slope radial dual."""

from fractions import Fraction as Q
from functools import lru_cache


ALPHA = Q(29, 200)
BETA = Q(3, 40)
TARGET = Q(997, 1000)
HALF = Q(1, 2)
RADIUS_MIN = Q(1, 10)
RADIUS_MAX = Q(33, 10)
SCALE = 10**25


def rounded_box(lo, hi):
    """Round rational endpoints outward onto a fixed rational grid."""
    floor = lo.numerator * SCALE // lo.denominator
    ceil = -((-hi.numerator * SCALE) // hi.denominator)
    return Q(floor, SCALE), Q(ceil, SCALE)


def atan_box(x):
    term = x
    total = Q(0)
    for index in range(16):
        total += term / (2 * index + 1)
        term *= -(x * x)
    return total, total + term / 33


def pi_box():
    small_lo, small_hi = atan_box(Q(1, 239))
    big_lo, big_hi = atan_box(Q(1, 5))
    return rounded_box(16 * big_lo - 4 * small_hi, 16 * big_hi - 4 * small_lo)


def gaussian_box(r):
    """Bounds for exp(-r*r/2) and its integral from zero to r."""
    assert 0 <= r <= RADIUS_MAX
    x = r * r / 2
    term = Q(1)
    exp_sum = Q(0)
    integral_sum = Q(0)
    for index in range(40):
        exp_sum += term
        integral_sum += term / (2 * index + 1)
        term *= -x / (index + 1)
    # Term 39 is negative, and every term in the remaining tail decreases.
    return (
        rounded_box(exp_sum, exp_sum + term),
        rounded_box(r * integral_sum, r * (integral_sum + term / 81)),
    )


PI_LO, PI_HI = pi_box()
INVA_LO, INVA_HI = 3 / (2 * PI_HI), 3 / (2 * PI_LO)
CUTOFF_LO, CUTOFF_HI = INVA_LO / ALPHA, INVA_HI / ALPHA
assert CUTOFF_HI < RADIUS_MAX
assert RADIUS_MIN < CUTOFF_LO


def primitive_box(r_lo, r_hi):
    """Enclose the untruncated angular primitive at r in [r_lo, r_hi]."""
    (exp_lo, _), (_, integral_hi) = gaussian_box(r_hi)
    (_, exp_hi), (integral_lo, _) = gaussian_box(r_lo)
    return rounded_box(
        INVA_LO * integral_lo - ALPHA * (1 - exp_lo),
        INVA_HI * integral_hi - ALPHA * (1 - exp_hi),
    )


Z_LO, Z_HI = primitive_box(CUTOFF_LO, CUTOFF_HI)
assert 0 < Z_LO <= Z_HI < 1


@lru_cache(None)
def endpoint(r):
    (s_lo, s_hi), _ = gaussian_box(r)
    if r >= CUTOFF_HI:
        return s_lo, s_hi, Z_HI, Q(0)
    k_lo, k_hi = primitive_box(min(r, CUTOFF_LO), min(r, CUTOFF_HI))
    j_hi = max(Q(0), Z_HI - k_lo) / s_lo
    return s_lo, s_hi, k_hi, j_hi


def interval_upper(r_lo, r_hi):
    _, s_hi, _, j_hi = endpoint(r_lo)
    s_lo, _, k_hi, _ = endpoint(r_hi)
    c_lo, c_hi = 1 - s_hi, 1 - s_lo
    d_hi = max(abs(c_lo - HALF), abs(c_hi - HALF))
    radial_hi = min(c_hi, s_hi, HALF) * (
        2 * ALPHA + BETA * (-Q(2, 5) + Q(32, 5) * d_hi)
    )
    if c_hi <= HALF:
        # S Z + (2c-1) J = [c*c Z + (1-2c) K]/S, with K = integral_0^c f.
        angular_hi = 2 * (c_hi * c_hi * Z_HI + (1 - 2 * c_lo) * k_hi) / s_lo
    else:
        angular_hi = 2 * (s_hi * Z_HI + (2 * c_hi - 1) * j_hi)
    return (angular_hi + radial_hi) / (s_lo * r_lo)


def main():
    # Below r=1/10 use c <= r*r/2, S >= 1-r*r/2, K <= r/a,
    # and (1-2c)/S**2 <= 1. Above r=33/10 the angular tail is zero.
    s_min = 1 - RADIUS_MIN**2 / 2
    small_upper = (
        2 * INVA_HI
        + Z_HI * RADIUS_MIN**3 / (2 * s_min**2)
        + RADIUS_MIN / (4 * s_min)
    )
    large_upper = (2 * Z_HI + HALF) / RADIUS_MAX
    assert small_upper < TARGET and large_upper < TARGET

    pending = [(RADIUS_MIN, RADIUS_MAX, 0)]
    maximum = max(small_upper, large_upper)
    accepted = 0
    max_depth = 0
    while pending:
        lo, hi, depth = pending.pop()
        bound = interval_upper(lo, hi)
        if bound < TARGET:
            maximum = max(maximum, bound)
            accepted += 1
            max_depth = max(max_depth, depth)
            continue
        assert depth < 20, (lo, hi, bound)
        mid = (lo + hi) / 2
        pending.extend(((lo, mid, depth + 1), (mid, hi, depth + 1)))
    printed_upper = 1 + (maximum.numerator * 10**8) // maximum.denominator
    print(f"Certified {accepted} intervals, maximum subdivision depth {max_depth}.")
    integer, decimals = divmod(printed_upper, 10**8)
    print(f"Global ratio upper bound < {integer}.{decimals:08d}.")


if __name__ == "__main__":
    main()
