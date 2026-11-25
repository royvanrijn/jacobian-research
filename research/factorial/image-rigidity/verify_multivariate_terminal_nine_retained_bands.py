#!/usr/bin/env python3
"""Exact terminal-nine pairs with retained rows in three reduced cap bands."""

from fractions import Fraction as Q
from math import comb, prod
from sys import argv

from flint import ctx
from verify_multivariate_terminal_nine_retained_rows import (
    ORDERS,CURVATURE_CAP,rows_above,core_powers,power_bounds,ball,POSITIONS,
    exponents,boxes,nominal_cosines,rays,require,paired_box,paired_radial,rational,
)
from verify_multivariate_terminal_six_joint_ray_bands import log_bounds

def outer_square(a,d):
    return prod(Q(8*a+h,8*a+7*d+h) for h in range(1,d+1))

def enumerate_band(cap,lower):
    bud=min(CURVATURE_CAP,Q(16,7)*log_bounds(1/lower)[1])
    minimum=lower**2
    rows=sorted(((a,d,o),s,outer_square(a,d))
                for a,d,s in rows_above(cap,bud)
                if Q(d,8*a+1)<=Q(1,8) and outer_square(a,d)>=minimum
                for o in (-1,1))
    result=[]
    def collect(index,left,square,selected):
        if square<minimum:
            return
        if index==len(rows):
            reversed_=tuple(sorted((a,d,-o) for a,d,o in selected))
            if selected and selected<=reversed_:
                result.append((selected,left))
            return
        row,s,v2=rows[index]
        for m in range(1+int(left/s)):
            collect(index+1,left-m*s,
                    square*v2**m,selected+(row,)*m)
    collect(0,bud,Q(1),())
    return bud,result

def paired_rays(angles,models):
    if any(v is None for v in models.values()):
        return None
    count=0
    for left,right in zip(angles,angles[1:]):
        for p in POSITIONS[:12]:
            guess=float(left.mid())*float(p)+float(right.mid())*float(1-p)
            c=Q(round(100000*guess),100000)
            if not left>ball(c)>right:
                continue
            low=paired_radial(8,rational(c),models[8])
            if low is None:
                continue
            high=paired_radial(9,rational(c),models[9])
            if high is not None:
                count+=low+high
                break
        else:
            return None
    return count

def verify_band_core(rows,budget,curvature,cap,outer_lower,outer_upper,refine):
    powers={j:core_powers(j,rows) for j in ORDERS}
    core_bounds={j:[power_bounds(value,Q(1,j)) for value in powers[j]] for j in ORDERS}
    full_cores={j:tuple((comb(j,k)*lo,comb(j,k)*hi)
                       for k,(lo,hi) in enumerate(core_bounds[j])) for j in ORDERS}
    core_lo,core_hi=power_bounds(powers[8][1]*powers[8][-2],Q(1,16))
    lo=max(outer_lower/core_hi,power_bounds(Q(8,15),budget)[0])
    hi=min(Q(1),outer_upper/core_lo)
    if lo>hi:
        return 0,0
    attenuation={k:power_bounds(Q(152,153),curvature*k*k)[1] for k in range(1,9)}
    pending=[((lo,hi),0)]
    leaves=cells=0
    while pending:
        tail,depth=pending.pop()
        midpoint=sum(tail)/2
        effective=((tail[0]**-2-1)/7)*2000
        t8=min(cap,Q(max(1,-(-effective.numerator//effective.denominator)),2000))
        if refine and Q(1,8)<tail[0]<1:
            v=tail[0]
            unit_gap=max(Q(1),(8*v*v-1)/(8*(1-v*v)),(1/t8-1)/8)
            a0=-(-unit_gap.numerator//unit_gap.denominator)
            t8=min(t8,max(Q(1,8*a0+1),(1-v)/(8*v-1)))
        shape=exponents(t8)
        angles=nominal_cosines(powers,midpoint,shape)
        if angles is not None:
            values=boxes(core_bounds,shape,tail,attenuation)
            count=None
            if all(min(values[j][1][0],values[j][-2][0])>j*Q(j-1,16)
                   and max(values[j][1][1],values[j][-2][1])<8
                   for j in ORDERS):
                count=rays(angles,values)
            if count is None and depth>=2 and all(
                    min(values[j][1][0],values[j][-2][0])>j*Q(j-1,20)
                    for j in ORDERS):
                models={j:paired_box(j,full_cores[j],tail,(t8,8*t8/(9-t8))) for j in ORDERS}
                count=paired_rays(angles,models)
            if count is not None:
                leaves+=1
                cells+=count
                continue
        require(depth<11,("unresolved tail",rows,tuple(map(float,tail)),angles))
        pending.extend((((tail[0],midpoint), depth+1),((midpoint,tail[1]),depth+1)))
    return leaves,cells


BANDS = (
    (Q(3,40), Q(31,50), Q(7,10), True),
    (Q(11,200), Q(7,10), Q(19,25), True),
    (Q(7,200), Q(3,4), Q(17,20), False),
)


def main():
    ctx.prec = 128
    bands = BANDS if len(argv) == 1 else (BANDS[int(argv[1])],)
    for cap, lower, upper, refine in bands:
        budget, lists = enumerate_band(cap, lower)
        leaves = cells = 0
        for rows, remainder in lists:
            a, b = verify_band_core(rows, remainder, budget-remainder, cap, lower, upper, refine)
            leaves += a
            cells += b
        print(f"{lower} <= u <= {upper}: {len(lists)} lists, "
              f"{leaves} outer cells and {cells} radial cells.", flush=True)


if __name__ == "__main__":
    main()
