#!/usr/bin/env python3
"""Exact linked ray separation for small terminal-nine outer means."""

from fractions import Fraction as Q
from functools import lru_cache
from math import comb

from flint import acb_poly, arb, ctx
import sympy as sym

from verify_multivariate_terminal_six_joint_factor_links import check_nonnegative_on_interval, x
from verify_multivariate_terminal_six_joint_ray_bands import (
    INNER, OUTER, coefficient_box, log_bounds, radial_replay, rounded_interval, power_bounds,
)
from verify_multivariate_terminal_six_private_angles import exponential_lower, multiply_interval, require
from verify_multivariate_terminal_six_retained_rows import RAY_FRACTIONS, ball

LOWER, UPPER = 8, 9
ORDERS = LOWER, UPPER
RATIO = Q(64,63)
CAP = Q(1,8)
TOP = Q(19,200)
RAY_POSITIONS = tuple(dict.fromkeys((*RAY_FRACTIONS, *(Q(k,32) for k in range(1,32)))))


def skew_weights(order,k):
    return tuple((Q(k<=m<order-k)-Q(order-2*k,order))/2 for m in range(order))


def verify_convexity():
    for j in ORDERS:
        for k in range(1,(j+1)//2):
            first=sum(w*m/(1+m*x) for m,w in enumerate(skew_weights(j,k)))
            r=1+(j-1)*x
            check_nonnegative_on_interval((r*sym.diff(r*first,x)).subs(x,x/2))


def exponent(j, k):
    return Q(k*(j-k),j-1)*(RATIO if j==UPPER else 1)


@lru_cache(maxsize=None)
def exponents(cap):
    denom_lo,denom_hi=log_bounds(1+7*cap)
    vals={
        LOWER:(cap,Q(1),Q(1)),
        UPPER:(8*cap/(9-cap),log_bounds(1+Q(64,9)*cap)[0]/denom_hi,
              1+log_bounds(1/(1-cap/9))[1]/denom_lo)
    }
    boxes={}
    for j,(c,ratio_lo,ratio_hi) in vals.items():
        scale_lo,scale_hi=log_bounds(1+(j-1)*c)
        factors=[]
        for k in range(1,j):
            d=min(k,j-k)
            lo=sum(log_bounds((1+(j-1-m)*c)/(1+m*c))[0] for m in range(d))/scale_hi
            skew=2*sum(multiply_interval((w,w),log_bounds(1+m*c))[1]
                       for m,w in enumerate(skew_weights(j,d)) if m and w)/scale_lo
            hi=Q(k*(j-k),j-1)
            enclosure=rounded_interval(ratio_lo*(lo-skew),ratio_hi*(hi+skew))
            require(0<enclosure[0]<=exponent(j,k)<=enclosure[1],("anchor power",j,k,enclosure))
            factors.append(enclosure)
        boxes[j]=tuple(factors)
    return boxes


def anchor_cosines(u):
    roots=[]
    for j in ORDERS:
        coeff=[arb(1)]+[comb(j,k)*ball(u)**ball(exponent(j,k)) for k in range(1,j)]+[arb(1)]
        roots_j=acb_poly(coeff).roots(maxprec=256)
        upper=[z.real/abs(z) for z in roots_j if z.imag>0]
        require(len(upper)==j//2 and sum(z.imag<0 for z in roots_j)==j//2,("anchor root signs",j,u))
        roots.extend((c,j) for c in upper)
    roots.sort(key=lambda p:float(p[0].mid()),reverse=True)
    require(tuple(j for _,j in roots)==(9,8)*4,("anchor sector pattern",u))
    require(all(a>b for (a,_),(b,_) in zip(roots,roots[1:])),"certified angular order")
    return tuple(c for c,_ in roots)


def ray_choices(left,right,box):
    for p in RAY_POSITIONS:
        proposal=Q(round(10**5*(float(left.mid())*float(p)+float(right.mid())*float(1-p))),10**5)
        if not left>ball(proposal)>right:
            continue
        low=radial_replay(8,box[8],(proposal,proposal))
        if low is None:
            continue
        high=radial_replay(9,box[9],(proposal,proposal))
        if high is not None:
            return low+high
    return None


def joint_boxes(lower,upper,shape):
    box={j:coefficient_box(j,(lower,upper),shape[j]) for j in ORDERS}
    for j in ORDERS:
        for k in range(1,j):
            lo,hi=box[j][k]
            total=comb(j,k)*power_bounds(upper,Q(25*k*(j-k),32*j))[1]
            box[j][k]=lo,min(hi,total)
    for k in range(1,9):
        lo,hi=box[9][k]
        left=box[8][k][1]/comb(8,k)*power_bounds(upper,Q(k*k,96))[1]
        right=box[8][k-1][1]/comb(8,k-1)*power_bounds(upper,Q((9-k)**2,96))[1]
        box[9][k]=lo,min(hi,comb(9,k)*left,comb(9,k)*right)
    return box

def verify(left,right,shape):
    pending=[(left,right,0)]
    count=radial_total=0
    while pending:
        lower,upper,depth=pending.pop()
        mid=(lower+upper)/2
        box=joint_boxes(lower,upper,shape)
        angles=anchor_cosines(mid)
        cells=0
        for a,b in zip(angles,angles[1:]):
            used=ray_choices(a,b,box)
            if used is None:
                require(depth<12,("unresolved ray", lower,upper,float(a.mid()),float(b.mid())))
                pending.extend(((lower,mid,depth+1),(mid,upper,depth+1)))
                break
            cells+=used
        else:
            count+=1
            radial_total+=cells
    print(left, right, count, radial_total, flush=True)
    return count, radial_total


if __name__=="__main__":
    ctx.prec=128
    verify_convexity()
    require(Q(153,152)**96>Q(15,8),"tail attenuation")
    require(exponential_lower(Q(16,25))>Q(15,8),"total curvature")
    shape=exponents(CAP)
    for j in ORDERS:
        full=coefficient_box(j,(Q(0),TOP),shape[j])
        require(sum(hi*INNER**k for k,(_,hi) in enumerate(full) if k)<1,"inner radius")
        require(sum(hi*OUTER**k for k,(_,hi) in enumerate(full[:-1]))<OUTER**j,"outer radius")
    steps=int(200*TOP)
    require(Q(steps,200)==TOP,"endpoint partition")
    u_cells = radial_cells = 0
    for k in range(steps):
        cells, radials = verify(Q(k,200),Q(k+1,200),shape)
        u_cells += cells
        radial_cells += radials
    print(f"Verified {u_cells} outer cells and {radial_cells} radial cells.")
