#!/usr/bin/env python3
"""Exact terminal-nine Cayley exclusions between outer means .76 and .88."""

from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations_with_replacement
from math import comb
from sys import argv

from verify_multivariate_terminal_nine_private_linear import ratio_bounds, shape
from verify_multivariate_terminal_nine_upper_bridge import projected_sign
from verify_multivariate_terminal_nine_upper_marked_rows import precise_power
from verify_multivariate_terminal_seven_upper_bridge import bisect_u, cayley_matrix, height
from verify_multivariate_terminal_six_intermediate_cayley import (
    GRID, I, ONE, ZERO, bisect, interval, log, rational, require, scalar,
)
from verify_multivariate_terminal_six_joint_ray_bands import power_bounds
from verify_multivariate_terminal_six_marked_row_tail import projected_interval


ORDERS = (8, 9)
REFERENCE_RATIO = Q(64, 63)
SCALE = 100000
DISK = Q(2)
FAST_BAND_ORDERS = {(830, 9), (870, 9)}
RAYS = {
    760: (131154, 112993, 85173, 74548, 50198, 39624, 18816),
    770: (130165, 112259, 85019, 74256, 50060, 39525, 18778),
    780: (129196, 111540, 84650, 73968, 49923, 39428, 18740),
    790: (128121, 110833, 84287, 73683, 49787, 39331, 18703),
    800: (127194, 110140, 83928, 74185, 49653, 39236, 18667),
    810: (126285, 109460, 83574, 73903, 49521, 39141, 18630),
    820: (125274, 109930, 83226, 73624, 49389, 39048, 18595),
    830: (124403, 109259, 82881, 73349, 49260, 38955, 18559),
    840: (123432, 108601, 82542, 73078, 49131, 38864, 18524),
    850: (123049, 107954, 82206, 72053, 49004, 38773, 18489),
    860: (122109, 107318, 81876, 71793, 48879, 38683, 18454),
    870: (121188, 106694, 81549, 72284, 48754, 38594, 18420),
    880: (120394, 106080, 81227, 72026, 48631, 38506, 18386),
}
EXPECTED_CELLS = {
    760: ((52150, 860, 14685, 2799, 4100, 881, 135),
          (82344, 4403, 14206, 11055, 6725, 1157, 240)),
    770: ((55955, 910, 15258, 2888, 4701, 921, 140),
          (88609, 4501, 19788, 11461, 7537, 1157, 228)),
    780: ((71486, 941, 15392, 3006, 4862, 968, 140),
          (88172, 4565, 21930, 11830, 7700, 1156, 228)),
    790: ((99790, 978, 15501, 3240, 5350, 974, 140),
          (92537, 4658, 24826, 12484, 8322, 1156, 239)),
    800: ((116476, 1053, 16797, 3495, 5836, 974, 167),
          (108958, 4898, 27226, 13176, 8736, 1261, 276)),
    810: ((155513, 1028, 17339, 3561, 6110, 974, 188),
          (172154, 4906, 29723, 13656, 9525, 1302, 304)),
    820: ((253317, 1071, 18106, 3638, 6699, 1012, 188),
          (224959, 5848, 32082, 14255, 10091, 1785, 329)),
    830: ((481968, 1077, 19385, 3692, 6808, 1479, 238),
          (48931, 2514, 7284, 2124, 1736, 387, 100)),
    840: ((77477, 915, 16014, 3298, 4894, 970, 140),
          (71806, 4091, 24607, 12458, 7836, 1136, 226)),
    850: ((62883, 946, 17137, 3248, 5740, 970, 141),
          (96137, 4322, 27204, 13035, 8442, 1160, 226)),
    860: ((102189, 965, 17456, 3497, 5907, 970, 174),
          (120454, 4475, 29655, 13587, 9080, 1250, 258)),
    870: ((132956, 1056, 19032, 3635, 6374, 987, 188),
          (41739, 2627, 6431, 2109, 1375, 328, 68)),
}


def reciprocal_cap(u):
    return Q(2,49) if u<Q(39,50) else Q(1,25) if u<Q(21,25) else Q(1,33)


def nominal(order, u, precise=False):
    ratio = 1 if order == 8 else REFERENCE_RATIO

    def power(exponent):
        if precise:
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
def model(order, u, fast):
    c8 = reciprocal_cap(u[0])
    cap = c8 if order == 8 else 8*c8/(9-c8)
    midpoint = sum(u)/2
    p = interval(1-u[1],1-u[0])
    lam = p.sqrt()
    scale = tuple(lam**degree for degree in range(order+1))
    reference = nominal(order, midpoint, fast)
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
def disk_check(order,u,fast):
    _,_,coarse,_=model(order,u,fast)
    return coarse[-1][1].lo*DISK**order > sum(
        max(abs(c.lo),abs(c.hi))*DISK**d for d,c in coarse[:-1])


def ray_at(left,right,index,u):
    a,b=RAYS[left][index],RAYS[right][index]
    return Q(a,SCALE)+Q(b-a,SCALE)*(1000*u-left)/(right-left)


def verify_ray(order,left,right,index):
    pending=[((Q(left,1000),Q(right,1000)),interval(Q(0),DISK),0)]
    accepted=0
    fast = (left, order) in FAST_BAND_ORDERS
    split_scale = 24 if fast else 120
    while pending:
        u,x,depth=pending.pop()
        if not disk_check(order,u,fast):
            require(depth<24,("disk failure",order,u))
            pending.extend((half,x,depth+1) for half in bisect_u(u))
            continue
        base,groups,coarse,p=model(order,u,fast)
        levels=tuple(sorted(ray_at(left,right,index,t) for t in u))
        y=height(p,x,interval(*levels))
        if projected_sign(order,x,y,base,groups,coarse,fast):
            accepted+=1
            continue
        require(depth<36,("ray failure",order,left,right,index,u,x.shown()))
        if split_scale*(u[1]-u[0])>Q(x.hi-x.lo,GRID):
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


def verify_cap_thresholds():
    low, mid, high = Q(19,25), Q(39,50), Q(21,25)
    require(Q(9,16) < low*low, "first unit row excluded")
    require(Q(51,65)**3 < low*low and Q(41*42,55*56) < low*low,
            "nonunit cap below .78")
    require(Q(26,33)**3 < mid*mid and Q(49*50,63*64) < mid*mid,
            "nonunit cap below .84")
    require((1-high)/(8*high-1) < Q(1,33), "nonunit cap from .84")


def main():
    keys = tuple(RAYS)
    first, last = (keys[0], keys[-1]) if len(argv) == 1 else tuple(map(int,argv[1:3]))
    require(first in keys and last in keys and first < last, "replay band")
    orders = ORDERS if len(argv) < 4 else (int(argv[3]),)
    require(all(order in ORDERS for order in orders), "replay order")
    verify_cap_thresholds()
    verify_anchor(first)
    total = 0
    for left,right in zip(keys,keys[1:]):
        if first <= left < right <= last:
            for order in orders:
                cells = tuple(verify_ray(order,left,right,i) for i in range(7))
                require(cells == EXPECTED_CELLS[left][order-8],
                        ("replay cell counts", order, left, right, cells))
                total += sum(cells)
                print(order,left,right,cells,flush=True)
    print(f"Certified {total} Cayley projection cells.",flush=True)


if __name__=="__main__":
    main()
