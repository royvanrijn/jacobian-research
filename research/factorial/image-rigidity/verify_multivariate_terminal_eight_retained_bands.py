#!/usr/bin/env python3
"""Exact paired replay above the light-row caps at terminal eight."""

from fractions import Fraction as Q
from functools import lru_cache
from math import prod

from flint import ctx

from verify_multivariate_terminal_eight_paired_light import coefficient_model, envelope, radial
from verify_multivariate_terminal_eight_retained_rows import (
    MAX_CURVATURE, ORDERS, RAY_FRACTIONS, ball, boxes, core_powers,
    nominal_cosines, power_bounds, rays, rows_above,
)
from verify_multivariate_terminal_eight_shared_high_curvature import exponent_bounds
from verify_multivariate_terminal_six_intermediate_cayley import interval, rational
from verify_multivariate_terminal_six_joint_ray_bands import log_bounds
from verify_multivariate_terminal_six_private_angles import require

PAIRED_RAY_FRACTIONS = tuple(dict.fromkeys((*RAY_FRACTIONS, *(Q(k,16) for k in range(1,16)))))

def row_outer_square(a,d):
    return prod(Q(7*a+h,7*a+6*d+h) for h in range(1,d+1))

def enumerate_cores(cap,outer_lower):
    budget = min(MAX_CURVATURE,Q(7,3)*log_bounds(1/outer_lower)[1])
    minimum = outer_lower**2
    rows = sorted(
        ((a,d,orientation),curvature,row_outer_square(a,d))
        for a,d,curvature in rows_above(cap,budget)
        if Q(d,7*a+1)<=Q(1,8) and row_outer_square(a,d)>=minimum
        for orientation in (-1,1)
    )
    cores = []
    def collect(index,left,square,selected):
        if square < minimum:
            return
        if index == len(rows):
            reflected = tuple(sorted((a,d,-sign) for a,d,sign in selected))
            if selected and selected <= reflected:
                cores.append((selected,left))
            return
        row,curvature,contribution = rows[index]
        for n in range(1+int(left/curvature)):
            collect(index+1,left-n*curvature,square*contribution**n,selected+(row,)*n)
    collect(0,budget,Q(1),())
    return budget,cores

@lru_cache(maxsize=None)
def paired_envelope(cap,order):
    return envelope(cap,order)

def paired_ray_replay(angles,core_arrays,tail,cap):
    models={}
    lower_bound=paired_envelope(cap,7)[1][0][1]
    for j in ORDERS:
        pieces=coefficient_model(j,tail,(Q(-1),Q(1)),paired_envelope(cap,j),lower_bound)
        retained=core_arrays[j]
        models[j]=tuple((parameter,tuple(a*b for a,b in zip(vector,retained,strict=True)))
                        for parameter,vector in pieces)
    count=0
    for left,right in zip(angles,angles[1:]):
        for position in PAIRED_RAY_FRACTIONS:
            proposal = float(left.mid())*float(position)+float(right.mid())*float(1-position)
            cosine=Q(round(100000*proposal),100000)
            if not left>ball(cosine)>right:
                continue
            j7=radial(7,rational(cosine),models[7])
            if j7 is None:
                continue
            j8=radial(8,rational(cosine),models[8])
            if j8 is not None:
                count+=j7+j8
                break
        else:
            return None
    return count

def verify_core(rows,retained_curvature,budget,cap,outer_lower,outer_upper):
    powers={j:core_powers(j,rows) for j in ORDERS}
    coefficients={j:[power_bounds(value,Q(1,j)) for value in powers[j]] for j in ORDERS}
    dyadic_cores={j:tuple(interval(lo,hi) for lo,hi in coefficients[j]) for j in ORDERS}
    core_lo,core_hi=power_bounds(powers[7][1]*powers[7][-2],Q(1,14))
    lower=max(outer_lower/core_hi,power_bounds(Q(7,13),budget)[0])
    upper=min(Q(1),outer_upper/core_lo)
    if lower>upper:
        return 0,0
    attenuation={k:power_bounds(Q(119,120),retained_curvature*k*k)[1] for k in range(1,8)}
    pending=[((lower,upper),0)]
    leaves=cells=0
    while pending:
        tail,depth=pending.pop()
        midpoint=sum(tail)/2
        proposed=2000*(tail[0]**-2-1)/6
        cap7=min(cap,Q(max(1,-(-proposed.numerator//proposed.denominator)),2000))
        shape=exponent_bounds(cap7,7*cap7/(8-cap7))
        cosines=nominal_cosines(powers,midpoint,shape)
        if cosines is not None:
            enclosure=boxes(coefficients,shape,tail,attenuation)
            count = None
            if all(min(enclosure[j][1][0],enclosure[j][-2][0]) > j*Q(7,16) for j in ORDERS):
                count=rays(cosines,enclosure)
            if count is None and depth>=3:
                try:
                    count=paired_ray_replay(cosines,dyadic_cores,tail,cap7)
                except ArithmeticError:
                    count=None
            if count is not None:
                leaves+=1
                cells+=count
                continue
        require(depth<14,("unresolved retained array",rows,tuple(map(float,tail))))
        pending.extend((((tail[0],midpoint),depth+1),((midpoint,tail[1]),depth+1)))
    return leaves,cells

BANDS = (
    (Q(2, 25), Q(13, 25), Q(3, 5), (5, 55, 4031)),
    (Q(13, 200), Q(3, 5), Q(33, 50), (12, 94, 7323)),
    (Q(9, 200), Q(33, 50), Q(18, 25), (27, 220, 16711)),
)


def main():
    ctx.prec = 128
    require(Q(120, 119)**75 > Q(13, 7), "tail attenuation from outer loss")
    require(Q(13, 12)**8 < 2, "paired annulus")
    for cap, low, high, expected in BANDS:
        budget, cores = enumerate_cores(cap, low)
        leaves = cells = 0
        for rows, left in cores:
            row_leaves, row_cells = verify_core(rows, budget-left, left, cap, low, high)
            leaves += row_leaves
            cells += row_cells
        totals = len(cores), leaves, cells
        require(totals == expected, ("retained replay counts", low, high, totals))
        print(f"{low} <= u <= {high}: {len(cores)} lists, {leaves} tail intervals, "
              f"and {cells} radial cells.", flush=True)


if __name__ == "__main__":
    main()
