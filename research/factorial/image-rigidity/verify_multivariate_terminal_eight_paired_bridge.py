#!/usr/bin/env python3
"""Exact terminal-eight paired continuation beyond u=.48."""
from fractions import Fraction as Q
from functools import lru_cache
from math import comb
from sys import argv

from verify_multivariate_terminal_eight_paired_light import (
    RAYS as PRECEDING_RAYS, bernstein_transform, check_anchor,
    coefficient_model, envelope, halves, signed, validate_skew,
)
from verify_multivariate_terminal_six_diffuse_bands import bernstein_coefficients
from verify_multivariate_terminal_six_intermediate_cayley import GRID,I,interval,require


RAYS={
    480: PRECEDING_RAYS[480],
    500: (129, -85, -449, -589, -842, -902),
    520: (88, -120, -473, -598, -853, -903),
    540: (48, -139, -490, -626, -854, -919),
    560: (5, -175, -517, -644, -861, -920),
    580: (-37, -212, -541, -662, -872, -921),
    600: (-80, -233, -567, -671, -881, -922),
    620: (-124, -270, -590, -689, -885, -929),
    630: (-146, -289, -602, -698, -890, -930),
    640: (-169, -308, -613, -707, -895, -931),
    660: (-213, -345, -636, -724, -901, -935),
    680: (-256, -370, -659, -741, -905, -942),
    690: (-279, -389, -670, -750, -909, -943),
    700: (-302, -409, -681, -759, -914, -944),
}

CAPS=((480,520,Q(1,8)),(520,560,Q(2,25)),(560,600,Q(2,25)),(600,630,Q(13,200)),
      (630,660,Q(13,200)),(660,690,Q(9,200)),(690,700,Q(9,200)))

@lru_cache(maxsize=None)
def sine_polynomials(order):
    result=[(0,),(1,)]
    for _ in range(1,order):
        a=(0,)+tuple(2*c for c in result[-1])
        b=result[-2]+(0,)*(len(a)-len(result[-2]))
        result.append(tuple(v-w for v,w in zip(a,b,strict=True)))
    return tuple(result)


@lru_cache(maxsize=16384)
def sine_intervals(order,cos_lo,cos_hi):
    lo,hi=Q(cos_lo,GRID),Q(cos_hi,GRID)
    out=[]
    for p in sine_polynomials(order):
        values=bernstein_coefficients(p,lo,hi)
        out.append(interval(min(values),max(values)))
    return tuple(out)


def projected(order,cosines,model):
    sines=sine_intervals(order,cosines.lo,cosines.hi)
    transform=bernstein_transform(order)
    result=[]
    for h in range(order+1):
        directions=tuple(sines[k-h] if k>=h else -sines[h-k] for k in range(order+1))
        groups=[]
        for _,coefficients in model:
            product=tuple(a*b for a,b in zip(coefficients,directions,strict=True))
            bounds=[]
            for row in transform:
                lo=hi=0
                for a,b in zip(row,product,strict=True):
                    vals=(a.lo*b.lo,a.lo*b.hi,a.hi*b.lo,a.hi*b.hi)
                    lo+=min(vals)
                    hi+=max(vals)
                bounds.append(I(lo//GRID,-(-hi//GRID)))
            groups.append(tuple(bounds))
        result.append(tuple(groups))
    return tuple(result)


@lru_cache(maxsize=None)
def binomial_weights(order):
    return tuple(comb(order,k) for k in range(order+1))


def midpoint_signed(projection,model,order):
    weights=binomial_weights(order)
    lo=hi=0
    for controls,(parameter,_) in zip(projection,model,strict=True):
        lower=sum(a.lo*b for a,b in zip(controls,weights,strict=True))
        upper=sum(a.hi*b for a,b in zip(controls,weights,strict=True))
        products=(parameter.lo*lower,parameter.lo*upper,parameter.hi*lower,parameter.hi*upper)
        lo+=min(products)
        hi+=max(products)
    return hi<0 or lo>0


def radial(order,cosines,model):
    stack=[(projected(order,cosines,model),0)]
    used=0
    while stack:
        polys,depth=stack.pop()
        if any(signed(candidate,model) for candidate in polys):
            used+=1
            continue
        if depth>=14 or not any(midpoint_signed(candidate,model,order) for candidate in polys):
            return None
        first,second=[],[]
        for polys_at_rotation in polys:
            children=tuple(halves(group) for group in polys_at_rotation)
            fl,fr=zip(*children,strict=True)
            first.append(fl)
            second.append(fr)
        stack.extend(((tuple(first),depth+1),(tuple(second),depth+1)))
    return used


def ray_cosines(left,right,index,u):
    a,b=RAYS[left][index],RAYS[right][index]
    values=tuple(Q(a,1000)+Q(b-a,right-left)*(v-Q(left,1000)) for v in u)
    return interval(min(values),max(values))


def verify(order,index,left,right,env,lower_bound):
    u0,u1=Q(left,1000),Q(right,1000)
    stack=[((u0,u1),(Q(0),Q(1)),0,0)]
    count=0
    while stack:
        u,normalized,du,ds=stack.pop()
        coefficients=coefficient_model(order,u,normalized,env,lower_bound)
        if coefficients is None:
            continue
        used=radial(order,ray_cosines(left,right,index,u),coefficients)
        if used is not None:
            count+=used
            continue
        require(max(du,ds)<15,("unresolved ray",order,index,left,right,u,normalized))
        if du<2*ds+3:
            mid=sum(u)/2
            stack.extend((((u[0],mid),normalized,du+1,ds),((mid,u[1]),normalized,du+1,ds)))
        else:
            mid=sum(normalized)/2
            stack.extend(((u,(normalized[0],mid),du,ds+1),(u,(mid,normalized[1]),du,ds+1)))
    return count


def main():
    validate_skew()
    check_anchor()
    for rays in RAYS.values():
        require(all(1000>a>b>-1000 for a,b in zip(rays,rays[1:])),"strict ray order")
    selection=tuple(map(int,argv[1:])) if len(argv)>1 else ()
    for a,b,cap in CAPS:
        if selection and a not in selection:
            continue
        env={j:envelope(cap,j) for j in (7,8)}
        limit=env[7][1][0][1].hi
        lower_bound=I(limit,limit)
        knots=[k for k in RAYS if a<=k<=b]
        for left,right in zip(knots,knots[1:]):
            counts=[verify(j,i,left,right,env[j],lower_bound) for j in (7,8) for i in range(6)]
            print(left,right,cap,counts,flush=True)


if __name__=="__main__":
    main()
