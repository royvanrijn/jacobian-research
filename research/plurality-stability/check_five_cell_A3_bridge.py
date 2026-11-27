#!/usr/bin/env python3
"""Certify the spherical chamber map and the radial five-cell comparison."""

from functools import cache
import json
import platform

import flint
from flint import arb, ctx, fmpq, fmpq_mat

from check_five_cell_weighted_gap import covariance_key, orthant_probability, raw_moment
from check_tetrahedral_standard_gap import require

ctx.prec = 192
HALF_PI = arb.pi() / 2
BINS = 256
DYADIC_BITS = 42
SLABS = (
    "41/50", "33/40", "83/100", "167/200", "21/25",
    "17/20", "43/50", "87/100", "22/25", "89/100",
    "9/10", "23/25", "47/50", "19/20",
)

SCHUR_CAP = fmpq(990793, 10**6)
GAP_MARGIN = fmpq(1, 200)


def interval(lo, hi):
    return arb(lo).union(arb(hi))

def angular(a, b, t):
    sine, cosine = t.sin_cos()
    d = a * cosine + b * sine
    dprime = -a * sine + b * cosine
    r = 1 + d*d
    return 1 - d/r.sqrt(), -dprime/(r*r.sqrt())

def primitive(a, b, t):
    sine, cosine = t.sin_cos()
    denom = (1 + a*a + b*b).sqrt()
    return t - ((a*sine - b*cosine)/denom).asin() - (b/denom).asin()

def verify_map():
    a, b = (arb(3)/5).sqrt(), 2/arb(5).sqrt()
    target = 1/arb(2).sqrt()
    area = primitive(a,b,HALF_PI)
    area_target = primitive(target,target,HALF_PI)
    area_ratio = area_target/area
    require(area > 0 and area_ratio > 1, "Invalid area comparison")
    D = 2**DYADIC_BITS

    def bracket_inverse(i):
        if i in (0,BINS):
            s = fmpq(i,BINS)
            return s,s
        phi = HALF_PI*arb(fmpq(i,BINS))
        aim = primitive(a,b,phi)*area_ratio
        lo,hi = 0,D
        while hi-lo > 1:
            mid = (lo+hi)//2
            psi = HALF_PI*arb(fmpq(mid,D))
            value = primitive(target,target,psi)
            if value < aim:
                lo = mid
            elif value > aim:
                hi = mid
            else:
                raise ArithmeticError(f"Unable to bracket inverse at endpoint {i}")
        return fmpq(lo,D),fmpq(hi,D)

    inverse_intervals = [bracket_inverse(i) for i in range(BINS+1)]
    parameters = []
    c_global = arb(1)
    for i in range(BINS):
        phi = HALF_PI*interval(fmpq(i,BINS),fmpq(i+1,BINS))
        psi = HALF_PI*interval(inverse_intervals[i][0],inverse_intervals[i+1][1])
        W, Wprime = angular(a,b,phi)
        Ws, Wsprime = angular(target,target,psi)
        require(0 < W < 1 and 0 < Ws < 1, f"Invalid angular weight {i}")
        c = Ws/W
        cprime = (Wsprime*area_ratio/c - c*Wprime)/W
        require(c > 1, f"The radial angle does not expand in box {i}")
        x_lo = c.lower()
        x_hi = (c*(2-W)/(2-Ws)).upper()
        B2 = abs(cprime).upper()**2/(c*(2-W)*(2-Ws))
        def metric_cap(x):
            return (x-1)*(area_ratio*area_ratio/x-1)
        require(B2 < metric_cap(x_lo) and B2 < metric_cap(x_hi),
                f"Metric expansion failed in box {i}")
        parameters.append((c,W,Ws,phi,psi))
        c_global = max(c_global,c.upper())

    derivative = (area_ratio/c_global).lower()
    require(derivative*2/arb(3).sqrt() > 1, "Smallest gap ratio failed")
    suffix_max = arb(1)
    last_gap = arb(2)
    for i,(c,W,Ws,phi,psi) in reversed(list(enumerate(parameters))):
        suffix_max = max(suffix_max,c.upper())
        # The suffix bound is valid including the endpoint where both cosines vanish.
        cosine_lower = min(arb(1),(area_ratio/suffix_max).lower())
        if i != BINS-1 and phi.cos() > 0:
            cosine_lower = max(cosine_lower,(psi.cos()/phi.cos()).lower())
        radial_lower = (c*(2-Ws)/(2-W)).sqrt()
        bound = (cosine_lower*radial_lower).lower()
        require(bound > 1, f"Largest gap ratio failed in box {i}")
        last_gap = min(last_gap,bound)
    return {"angular_bins":BINS, "inverse_dyadic_bits":DYADIC_BITS,
            "area_ratio":area_ratio.str(14), "azimuth_derivative_lower":derivative.str(14),
            "last_gap_ratio_lower":last_gap.str(14)}

def radial_polynomial(inverse):
    # R = |x|² = wᵀ C_w⁻¹ w. The output is (R−3)²/6.
    poly = {(0,0,0): fmpq(-3)}
    for i in range(3):
        for j in range(3):
            power = tuple((k==i)+(k==j) for k in range(3))
            poly[power] = poly.get(power,fmpq(0)) + inverse[i,j]
    square = {}
    for i,a in poly.items():
        for j,b in poly.items():
            power = tuple(x+y for x,y in zip(i,j))
            square[power] = square.get(power,fmpq(0)) + a*b/6
    return square

def verify_radial():
    base = fmpq_mat([[2,-1,0],[-1,2,-1],[0,-1,fmpq(3,2)]])
    inverse = base.inv()
    polynomial = radial_polynomial(inverse)
    p_full = arb(1)/8+3*(arb(1)/3).asin()/(4*arb.pi())
    prob = orthant_probability(covariance_key(base))
    require((6*prob).overlaps(p_full),"Ordered cone does not have the expected Gaussian mass")

    @cache
    def tilted_radial(rate):
        precision = fmpq_mat([[inverse[i,j]+(2*rate if i==2 and j==2 else 0)
                              for j in range(3)] for i in range(3)])
        covariance = precision.inv()
        log_free = arb(covariance.det()/base.det()).sqrt()/prob
        key = covariance_key(covariance)
        return log_free*sum((arb(c)*raw_moment(key,exponent) for exponent,c in polynomial.items()),arb(0))

    require(tilted_radial(fmpq(0)).contains(1),"Untilted radial second moment must be one")
    L = 2/(5*p_full)
    require(L > fmpq(3,2), "Row cap does not give monotone prefactor")
    q = (-15/(16*arb.pi())).exp()/4
    side = L*q
    def rate(rho):
        return fmpq(17,15)*rho*rho/(1-rho*rho)
    def offset(rho):
        return 1+arb(rho)+arb(rho*rho)-L

    rows = []
    for first,last in zip(SLABS, SLABS[1:]):
        a,c = fmpq(first),fmpq(last)
        require(a < c and offset(a) > 0, "Invalid radial slab")
        X = tilted_radial(rate(c))
        loss = X-tilted_radial(rate(a)+rate(c))
        require(loss > 0, "The positive quadratic secant correction was not certified")
        factor = side/(offset(c)+side)
        power = arb(a*a)/offset(a)*(1-factor*X-factor*factor*loss)
        require(0 < power < arb(SCHUR_CAP), f"Radial Schur bound failed on {a} to {c}: {power}")
        rows.append({"left":str(a), "right":str(c), "schur_upper":power.str(14)})
    require(offset(fmpq(SLABS[0]))*(1-arb(SCHUR_CAP)) > arb(GAP_MARGIN),
            "The bridge does not give the claimed gap margin")
    require(1+arb(fmpq(SLABS[-1]))-L > arb(GAP_MARGIN),
            "The bare row-cap interval does not give the claimed gap margin")
    return {"row_cap":L.str(14), "schur_cap":str(SCHUR_CAP),
            "scaled_gap_margin":str(GAP_MARGIN), "radial_slabs":rows}

if __name__ == "__main__":
    print(json.dumps({"python": platform.python_version(), "python_flint": flint.__version__,
                      "precision_bits":ctx.prec, "map":verify_map(),
                      "radial":verify_radial()}, indent=2))

