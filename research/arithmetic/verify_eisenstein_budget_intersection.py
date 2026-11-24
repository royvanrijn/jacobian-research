#!/usr/bin/env python3
"""Exact exponent and shape enclosures for the Eisenstein sieve budgets."""

from bisect import bisect_right
from functools import lru_cache
from fractions import Fraction
from math import comb, factorial

from verify_eisenstein_factor_moments import (
    factor_constant as earlier_factor_constant,
    linear_interval,
    show_interval,
)


SCALE = 10**15
GRID = 192
HALF = Fraction(1, 2)
BETA = Fraction(4267, 1000)


def logarithm(x):
    u = (x - 1) / (x + 1)
    terms = 64
    lower = 2 * sum(u**(2*k + 1) / (2*k + 1) for k in range(terms))
    tail = 2*u**(2*terms + 1) / ((2*terms + 1)*(1-u*u))
    return lower, lower + tail


def cube_simplex_fraction(dimension, cutoff):
    """Volume fraction of the unit cube below sum of coordinates <= cutoff."""
    return Fraction(
        sum((-1)**i * comb(dimension, i) * (cutoff - i)**dimension
            for i in range(cutoff+1)),
        factorial(dimension),
    )


def truncated_power(weights, degree):
    """Coefficients of the indicated power through degree len(weights)-1."""
    result = weights
    for _ in range(degree-1):
        result = [
            sum(weights[i] * result[k-i] for i in range(k+1))
            for k in range(len(weights))
        ]
    return result


def factor_constant(eta, rank):
    if rank == 1:
        return Fraction(1), Fraction(1)
    if rank * eta >= 1:
        return Fraction(0), Fraction(0)
    if rank == 2:
        return logarithm((1-eta)/eta)

    dimension = rank - 1
    step = (1-rank*eta)/GRID
    midpoints = []
    for index in range(GRID):
        mass = step / (eta + (index + HALF)*step)
        midpoints.append((SCALE*mass.numerator)//mass.denominator)
    assert min(midpoints) > 0
    coefficients = truncated_power(midpoints, dimension)

    full = sum(coefficients[:GRID-dimension+1])
    partial = sum(
        coefficients[GRID-slack]*cube_simplex_fraction(dimension, slack)
        for slack in range(1, dimension)
    )
    # Convexity bounds a full cell by midpoint and trapezoid.
    # A boundary cell has known volume and endpoint bounds instead.
    relative_half_step = step/(2*eta)
    rounding_upper = (1+Fraction(1,min(midpoints)))**dimension
    norm = SCALE**dimension * factorial(dimension)
    lower = (full + (1-relative_half_step)**dimension * partial)/norm
    upper = rounding_upper * (
        (1+relative_half_step**2)**dimension * full
        + (1+relative_half_step)**dimension * partial
    ) / norm
    return lower, upper


@lru_cache(maxsize=None)
def rank_constants(eta):
    return tuple(
        factor_constant(eta, rank)
        for rank in range(1, (eta.denominator-1)//eta.numerator + 1)
    )


@lru_cache(maxsize=None)
def total_mass(eta):
    constants = rank_constants(eta)
    return linear_interval(constants, [HALF**i for i in range(len(constants))])


def capped_mass(eta, cap, step):
    """Lower/upper enclosure for the finite union of low-product bin pairs."""
    nodes = [eta]
    for i in range(1, int(HALF/step)+1):
        if eta < i*step <= cap-eta:
            nodes.append(i*step)
    totals = [total_mass(node) for node in nodes]
    lo = [Fraction(0)]
    for _, hi in totals[1:]:
        lo.append(max(lo[-1], totals[0][0]-hi))
    hi = [totals[0][1]-lower for lower, _ in totals]
    for j in range(len(hi)-2, -1, -1):
        hi[j] = min(hi[j], hi[j+1])
    hi[0] = Fraction(0)

    partner = [
        max(0, bisect_right(nodes, cap-node)-1)
        for node in nodes[1:]
    ] + [0]

    def abel(cdf):
        return sum(cdf[i]*(cdf[partner[i-1]]-cdf[partner[i]])
                   for i in range(1, len(nodes)))
    return abel(lo), abel(hi)


def missed_mass_lower(eta, cap, step):
    """Charge only rectangles whose every composite-minimum pair is missed."""
    nodes = [eta] + [i*step for i in range(1, int(HALF/step)+1) if i*step>eta]
    totals = [total_mass(node) for node in nodes]
    differences = [
        max(Fraction(0), totals[i][0]-totals[i+1][1])
        for i in range(len(nodes)-1)
    ]
    return sum(
        first*second
        for i, first in enumerate(differences)
        for j, second in enumerate(differences)
        if nodes[i]+nodes[j] >= cap
    )


def common_ceiling_margin(level, step):
    eta = level / BETA
    lo, hi = total_mass(eta)
    inner = capped_mass(eta, level, step)[0]
    outer = capped_mass(eta, level + 2*step, step)[1]
    return 2/hi - 1 + inner/hi**2, 2/lo - 1 + outer/lo**2


def main():
    for rank in (2, 3, 4):
        previous = earlier_factor_constant(Fraction(1, 5), rank)
        current = factor_constant(Fraction(1, 5), rank)
        assert max(previous[0], current[0]) <= min(previous[1], current[1])

    for eta, sign in ((Fraction(164,1000), -1), (Fraction(165,1000), 1)):
        lower, upper = total_mass(eta)
        signed = 2/upper-1, 2/lower-1
        assert signed[1] < 0 if sign < 0 else signed[0] > 0
        print("single mark", eta, show_interval(signed))

    eta = Fraction(29,250)
    cap = Fraction(99,200)
    floor, ceiling = total_mass(eta)
    marked = capped_mass(eta, cap, Fraction(1,200))
    normalized = marked[0]/ceiling**2, marked[1]/floor**2
    diff = 2/ceiling-1+normalized[0], 2/floor-1+normalized[1]
    assert diff[0] > 0
    composite = 1-1/floor, 1-1/ceiling
    safe_common_error = diff[0]/(1+2*composite[1]+normalized[1])
    assert safe_common_error > Fraction(19,10000)
    print("cross mark", eta, cap, "R-normalized composite", show_interval(composite),
          "part recovered", show_interval(normalized), "margin", show_interval(diff))
    print("cross safe common fractional error", show_interval((safe_common_error, safe_common_error)))
    band_low, band_high = total_mass(cap)
    band = (band_low-1)/ceiling, (band_high-1)/floor
    assert band[0] > diff[1]/2
    print("semiprime fraction with both factors above a^0.495", show_interval(band))
    print("separate exponents", BETA*eta, cap, "product exponent", BETA*eta+cap)

    step = Fraction(1, 2000)
    for level, sign in ((Fraction(487, 1000), -1), (Fraction(489, 1000), 1)):
        margin = common_ceiling_margin(level, step)
        assert margin[1] < 0 if sign < 0 else margin[0] > 0
        print("best full-cap margin at common ceiling", level, show_interval(margin))

    missed = missed_mass_lower(Fraction(2,25), HALF, Fraction(1,40))
    assert missed > 1
    assert HALF/(BETA+2) < Fraction(2,25)
    print("missed unnormalized mass lower bound", show_interval((missed, missed)))
    missed_linear = missed_mass_lower(Fraction(1,8), Fraction(17,50), Fraction(1,20))
    assert missed_linear > 1
    print("missed at eta=1/8 and tau=0.34", show_interval((missed_linear, missed_linear)))

    direct_eta = 1/(1+BETA)
    constants = rank_constants(direct_eta)
    proper = linear_interval(constants, [1, *[-i*HALF**i for i in range(1,len(constants))]])
    low, high = total_mass(direct_eta)
    assert proper[0] > 0
    print("separate proper-divisor intersection", direct_eta, "level", 1-direct_eta,
          "margin", show_interval((proper[0]/high, proper[1]/low)))


if __name__ == "__main__":
    main()
