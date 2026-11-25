#!/usr/bin/env python3
"""Exact terminal-eight separation for light profiles with small outer mean."""
from fractions import Fraction as Q
from verify_multivariate_terminal_eight_paired_light import RAYS as PAIRED_RAYS, envelope, exponents, validate_skew
from verify_multivariate_terminal_six_intermediate_cayley import GRID,ONE,ZERO,rational,require
from verify_multivariate_terminal_six_joint_ray_bands import (
    INNER, OUTER, coefficient_box, rounded_interval, radial_replay,
)


RAYS={
    0:(913,794,309,-207,-498,-843),
    20:(904,789,270,-113,-532,-845),
    40:(890,769,234,-132,-552,-847),
    60:(872,735,195,-149,-575,-843),
    80:(851,698,160,-166,-592,-845),
    100:(828,684,132,-207,-607,-856),
}

def bounds(order,cap):
    rho,pieces=envelope(cap,order)
    if order==7:
        rho=ONE
    answer=[]
    for k,e in enumerate(exponents(order)[1:-1],1):
        delta=ZERO
        for vector,value in pieces:
            delta=delta+vector[k]*value
        exponent=(rational(e)-delta)*rho
        lo,hi=rounded_interval(Q(exponent.lo,GRID),Q(exponent.hi,GRID))
        require(0<lo and lo<=e*(Q(1) if order==7 else Q(49,48))<=hi,"exponent anchor")
        answer.append((lo,hi))
    return tuple(answer)


def verify_anchor():
    require(RAYS[100]==PAIRED_RAYS[100],"matching paired endpoint rays")
    for row in RAYS.values():
        require(all(1000>a>b>-1000 for a,b in zip(row,row[1:])),"strict ray order")
    points=(Q(1),)+tuple(Q(x,1000) for x in RAYS[0])+(Q(-1),)
    for form,signs in (
        (lambda x:8*x**3-4*x*x-4*x+1,(1,1,-1,-1,1,1,-1,-1)),
        (lambda x:8*x**4-8*x*x+1,(1,-1,-1,1,1,-1,-1,1)),
    ):
        require(all(form(x)*sign>0 for x,sign in zip(points,signs,strict=True)),
                "binomial sector signs")


def cosines(left,right,index,u):
    a,b=RAYS[left][index],RAYS[right][index]
    values=tuple(Q(a,1000)+Q(b-a,right-left)*(v-Q(left,1000)) for v in u)
    return min(values),max(values)


def verify_ray(left,right,index,order,exp):
    stack=[((Q(left,1000),Q(right,1000)),0)]
    count=0
    while stack:
        u,depth=stack.pop()
        cells=radial_replay(order,coefficient_box(order,u,exp),cosines(left,right,index,u))
        if cells is not None:
            count+=cells
            continue
        midpoint=sum(u)/2
        require(depth<15,("unresolved small-outer ray",left,right,index,order,u))
        stack.extend((((u[0],midpoint),depth+1),((midpoint,u[1]),depth+1)))
    return count


def main():
    validate_skew()
    verify_anchor()
    exponents={j:bounds(j,Q(1,8)) for j in (7,8)}
    for j in (7,8):
        full=coefficient_box(j,(Q(0),Q(1,10)),exponents[j])
        require(sum(b*INNER**k for k,(_,b) in enumerate(full) if k)<1,"inner radius")
        require(sum(b*OUTER**k for k,(_,b) in enumerate(full[:-1]))<OUTER**j,"outer radius")
    knots=tuple(RAYS)
    for left,right in zip(knots,knots[1:]):
        sizes=[verify_ray(left,right,i,j,exponents[j]) for j in (7,8) for i in range(6)]
        print(left,right,sizes,flush=True)


if __name__=="__main__":
    main()
