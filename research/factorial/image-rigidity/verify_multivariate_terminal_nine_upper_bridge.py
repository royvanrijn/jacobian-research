#!/usr/bin/env python3
"""Exact terminal-nine Cayley projections on the shared upper moderate band."""

from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations_with_replacement
from math import comb
from sys import argv

from verify_multivariate_terminal_nine_private_linear import ratio_bounds, shape
from verify_multivariate_terminal_nine_upper_marked_rows import precise_power
from verify_multivariate_terminal_seven_upper_bridge import bisect_u, cayley_matrix, height
from verify_multivariate_terminal_six_intermediate_cayley import (
    GRID, I, ONE, ZERO, bisect, interval, log, projections, rational, require, scalar,
)
from verify_multivariate_terminal_six_joint_ray_bands import power_bounds
from verify_multivariate_terminal_six_marked_row_tail import projected_interval


ORDERS = (8, 9)
REFERENCE_RATIO = Q(64, 63)
SCALE = 100000
DISK = Q(2)
FINE_PRECISION_START = Q(47, 50)
REPLAY_BANDS = ((880, 895), (900, 905), (910, 915))
RAYS = {
    880: (120667, 106504, 81134, 71533, 48649, 39002, 18433),
    885: (120273, 106199, 80975, 71407, 48588, 38958, 18416),
    890: (119883, 105896, 80817, 71281, 48527, 38914, 18399),
    895: (119341, 105595, 80660, 71156, 48467, 38870, 18382),
    900: (118959, 105297, 80503, 71032, 48406, 38826, 18366),
    905: (118580, 105002, 80348, 70909, 48346, 38782, 18349),
    910: (118205, 104709, 80194, 70786, 48287, 37910, 18333),
    915: (117683, 104418, 80040, 70664, 48227, 37868, 18316),
    920: (117315, 104130, 79888, 70543, 48168, 37826, 18300),
    925: (117099, 103844, 79736, 70422, 48110, 38610, 18284),
    930: (116736, 103561, 79586, 70302, 48051, 38568, 18267),
    935: (116377, 103279, 79436, 70183, 47993, 38525, 18251),
    940: (116165, 103000, 79288, 70064, 47935, 38483, 18235),
    945: (115811, 104158, 79140, 69947, 47877, 38441, 18219),
    950: (115460, 103876, 78993, 69829, 47820, 38400, 18203),
}


def reciprocal_cap(u):
    lower = 1
    while u*u > Q(8*lower+1, 8*lower+8):
        lower += 1
    return min(Q(1,41),(u**-2-1)/7, max(Q(1, 8*lower+1), (1-u)/(8*u-1)))


def nominal(order, u):
    ratio = 1 if order == 8 else REFERENCE_RATIO
    def power(exponent):
        if u >= FINE_PRECISION_START:
            return precise_power(u, exponent)
        return interval(*power_bounds(u, exponent))
    return tuple(scalar(comb(order,k))*power(ratio*Q(k*(order-k),order-1))
                 if 0 < k < order else ONE for k in range(order+1))


@lru_cache(maxsize=64)
def vectors(order):
    gauss, skew, defect, skew_sizes, distances = shape(order)
    pure = (gauss, *skew, *defect)
    products = tuple(tuple(a*b for a,b in zip(pure[i],pure[j],strict=True))
                     for i,j in combinations_with_replacement(range(len(pure)),2))
    matrix = cayley_matrix(order)
    half = order//2
    def pair(vector):
        ans = []
        for weights in matrix:
            ans.append(tuple(scalar(weights[k])*vector[k]+scalar(weights[-1-k])*vector[-1-k]
                             for k in range((order+1)//2))
                       + ((scalar(weights[half])*vector[half],) if order%2 == 0 else ()))
        return tuple(ans)
    base = pair((ONE,)*(order+1))
    return pure, skew_sizes, distances, base, tuple(map(pair,(*pure,*products)))


@lru_cache(maxsize=8192)
def model(order, u):
    c8 = reciprocal_cap(u[0])
    cap = c8 if order == 8 else 8*c8/(9-c8)
    midpoint = sum(u)/2
    p = interval(1-u[1],1-u[0])
    lam = p.sqrt()
    scale = tuple(lam**degree for degree in range(order+1))
    reference = nominal(order, midpoint)
    outer_log = -log(interval(*u))
    mean = log(interval(u[0]/midpoint,u[1]/midpoint))
    if order == 9:
        ratio = ratio_bounds(c8)
        flat = rational(REFERENCE_RATIO)
        mean = flat*mean + (flat-ratio)*outer_log
        outer_log = ratio*outer_log
    m = order-1
    T = rational(m*cap/(2+m*cap))
    parameters = (mean,
                  *((T**(2*k-1)*outer_log*scalar(1,2*k)).symmetric_error() for k in range(1,4)),
                  *((T**(2*k)*outer_log*scalar(1,2*k+1)).positive_error() for k in range(1,3)))
    pure, sizes, distances, base_vectors, grouped = vectors(order)
    skew_tail = T**7*outer_log/(scalar(8)*(ONE-T.square()))
    defect_tail = T**6*outer_log/(scalar(7)*(ONE-T.square()))

    errors=[]
    for k, coefficient in enumerate(reference):
        tail = sizes[k]*skew_tail+scalar(distances[k])*defect_tail
        radius = tail
        for parameter, vector in zip(parameters,pure,strict=True):
            value = parameter*vector[k]
            radius += I(0,max(-value.lo,value.hi))
        require(radius.hi < GRID, ("exponential radius",order,u))
        error = scalar(1,6)*radius**3/(ONE-radius)+tail*(ONE+radius)
        errors.append(coefficient*error.symmetric_error())

    def convert(group, degree):
        return projected_interval(tuple(enumerate(group)),reference)*scale[degree]

    base = tuple((d,convert(group,d)+projected_interval(tuple(enumerate(errors)),
                 tuple(map(scalar,cayley_matrix(order)[d])))*scale[d])
                 for d,group in enumerate(base_vectors))
    products = tuple((scalar(1,2)*parameters[i].square() if i==j else parameters[i]*parameters[j])
                     for i,j in combinations_with_replacement(range(len(parameters)),2))
    dense = tuple(tuple(convert(v,d) for d,v in enumerate(terms)) for terms in grouped)
    groups = tuple((s,tuple((d,c) for d,c in enumerate(terms) if c.lo or c.hi))
                  for s,terms in zip((*parameters,*products),dense,strict=True))
    coarse = tuple((d,c+sum_intervals(s*v[d] for s,v in zip((*parameters,*products),dense,strict=True)))
                   for d,c in base)
    return base, groups, coarse, p


def sum_intervals(terms):
    total=ZERO
    for term in terms:
        total += term
    return total


@lru_cache(maxsize=8192)
def disk_check(order,u):
    _,_,coarse,_=model(order,u)
    return coarse[-1][1].lo*DISK**order > sum(
        max(abs(c.lo),abs(c.hi))*DISK**d for d,c in coarse[:-1])


def ray_at(left,right,index,u):
    a,b=RAYS[left][index],RAYS[right][index]
    return Q(a,SCALE)+Q(b-a,SCALE)*(1000*u-left)/(right-left)


def projected_pair(terms,re,im):
    return projected_interval(terms,re), projected_interval(terms,im)


def axis_encloses(rectangle, a, b):
    real, imag = rectangle
    lo = a*(real.lo if a>=0 else real.hi) + b*(imag.lo if b>=0 else imag.hi)
    hi = a*(real.hi if a>=0 else real.lo) + b*(imag.hi if b>=0 else imag.lo)
    return lo,hi


def zonotope_sign(order,x,y,base,groups,coarse,fine):
    re=[ONE]; im=[ZERO]
    for _ in range(order):
        re.append(re[-1]*x-im[-1]*y)
        im.append(re[-2]*y+im[-1]*x)
    rect=projected_pair(coarse,re,im)
    if any(a.lo>0 or a.hi<0 for a in rect):
        return True
    rect=projected_pair(base,re,im)
    vectors=tuple((a,projected_pair(v,re,im)) for a,v in groups)
    candidates=[(1,0),(0,1)]
    if fine:
        candidates.extend((v.lo+v.hi,-w.lo-w.hi) for w,v in zip(re,im,strict=True))
    candidates.extend((pair[1].lo+pair[1].hi,-pair[0].lo-pair[0].hi)
                      for parameter,pair in vectors
                      if parameter.lo!=parameter.hi)
    for a,b in candidates:
        if (a,b)==(0,0):
            continue
        lo,hi=axis_encloses(rect,a,b)
        lo*=GRID; hi*=GRID
        for parameter,pair in vectors:
            v,w=axis_encloses(pair,a,b)
            s,t=parameter.lo,parameter.hi
            terms=v*s,v*t,w*s,w*t
            lo+=min(terms); hi+=max(terms)
        if lo>0 or hi<0:
            return True
    return False


def projected_sign(order,x,y,base,groups,coarse,fine):
    if zonotope_sign(order,x,y,base,groups,coarse,fine):
        return True
    directions=tuple(projections(order,x,y))
    for direction in directions[1:]:
        total=projected_interval(coarse,direction)
        if total.lo>0 or total.hi<0:
            return True
    return False


def verify_ray(order,left,right,index):
    pending=[((Q(left,1000),Q(right,1000)),interval(Q(0),DISK),0)]
    accepted=0
    while pending:
        u,x,depth=pending.pop()
        if not disk_check(order,u):
            require(depth<24,("disk failure",order,u))
            pending.extend((half,x,depth+1) for half in bisect_u(u))
            continue
        base,groups,coarse,p=model(order,u)
        levels=tuple(sorted(ray_at(left,right,index,t) for t in u))
        y=height(p,x,interval(*levels))
        if projected_sign(order,x,y,base,groups,coarse,Q(left,1000)>=FINE_PRECISION_START):
            accepted+=1
            continue
        require(depth<36,("ray failure",order,left,right,index,u,x.shown()))
        if 120*(u[1]-u[0])>Q(x.hi-x.lo,GRID):
            pending.extend((half,x,depth+1) for half in bisect_u(u))
        else:
            pending.extend((u,half,depth+1) for half in bisect(x))
    return accepted


def reciprocal_quotient(order,u,z):
    coeffs=nominal(order,u)
    if order%2:
        coeffs=tuple(sum_intervals(scalar((-1)**(k-l))*coeffs[l] for l in range(k+1))
                     for k in range(order))
    m=len(coeffs)//2
    s,t=scalar(2),z
    value=coeffs[m]
    for k in range(1,m+1):
        value+=coeffs[m-k]*t
        s,t=t,z*t-s
    return value


def verify_anchor(anchor):
    for row in RAYS.values():
        vertices=(2*SCALE,*row,0)
        require(all(a>b for a,b in zip(vertices,vertices[1:])), "ray order")
    u=Q(anchor,1000)
    require(DISK**2*(1-u)<1,"left half-plane")
    positions=(scalar(2),)+tuple(rational(2*((1-u)*b*b-1)/((1-u)*b*b+1))
                               for b in (Q(k,SCALE) for k in RAYS[anchor]))+(scalar(-2),)
    for order in ORDERS:
        last=switches=0
        for i,z in enumerate(positions):
            val=reciprocal_quotient(order,u,z)
            new=1 if val.lo>0 else -1 if val.hi<0 else 0
            require(new != 0,("anchor sign",order,i,val.shown()))
            if i:
                require((new!=last)==(order==9 if i%2 else order==8),("anchor separation",order,i,val.shown()))
                switches+=new!=last
            last=new
        require(switches==4,"anchor root count")


def main():
    keys=tuple(RAYS)
    bands = REPLAY_BANDS if len(argv)==1 else (tuple(map(int,argv[1:3])),)
    total = 0
    for first,last in bands:
        require(first in keys and last in keys and first<last,"replay band")
        verify_anchor(first)
        for left,right in zip(keys,keys[1:]):
            if not first<=left<right<=last:
                continue
            for order in ORDERS:
                for i in range(7):
                    if len(argv)>3 and (order,i)!=tuple(map(int,argv[3:])):
                        continue
                    n=verify_ray(order,left,right,i)
                    total+=n
                    print(left,right,order,i,n,flush=True)
    print("Certified",total,"Cayley cells",flush=True)


if __name__=="__main__":
    main()
