#!/usr/bin/env python3
"""Exact terminal-seven separation for light shared rows and small outer mean."""
from fractions import Fraction as Q
from verify_multivariate_terminal_seven_private_linear import exponent_bounds, verify_skew_convexity
from verify_multivariate_terminal_six_joint_ray_bands import (
    INNER, OUTER, coefficient_box, radial_replay,
)
from verify_multivariate_terminal_six_private_angles import require


RAYS = {
    0: (886, 572, 117, -371, -752),
    20: (869, 574, 72, -398, -771),
    40: (850, 582, 36, -407, -785),
    60: (833, 580, 5, -412, -799),
    80: (811, 583, -26, -416, -806),
    100: (789, 583, -59, -419, -813),
}


def verify_anchor():
    for row in RAYS.values():
        require(all(1000>a>b>-1000 for a,b in zip(row,row[1:])),"ordered rays")
    points=(Q(1),)+tuple(Q(x,1000) for x in RAYS[0])+(Q(-1),)
    for form,signs in (
        (lambda x: 8*x**3-6*x, (1,1,-1,-1,1,1,-1)),
        (lambda x: 8*x**3-4*x*x-4*x+1,(1,-1,-1,1,1,-1,-1)),
    ):
        require(all(form(x)*sign>0 for x,sign in zip(points,signs,strict=True)),
                "pure binomial sector signs")


def ray_cosines(left,right,index,u):
    a,b=RAYS[left][index],RAYS[right][index]
    cos=tuple(Q(a,1000)+Q(b-a,right-left)*(v-Q(left,1000)) for v in u)
    return min(cos),max(cos)


def verify_ray(left,right,index,order,bounds):
    stack=[((Q(left,1000),Q(right,1000)),0)]
    count=0
    while stack:
        u,depth=stack.pop()
        cells=radial_replay(order,coefficient_box(order,u,bounds),ray_cosines(left,right,index,u))
        if cells is not None:
            count+=cells
            continue
        midpoint=sum(u)/2
        require(depth<15,("unresolved small-outer ray",left,right,index,order,u))
        stack.extend((((u[0],midpoint),depth+1),((midpoint,u[1]),depth+1)))
    return count


def main():
    verify_skew_convexity()
    verify_anchor()
    cap=Q(1,8)
    bounds=exponent_bounds(cap,6*cap/(7-cap))
    for order in (6,7):
        full=coefficient_box(order,(Q(0),Q(1,10)),bounds[order])
        require(sum(b*INNER**k for k,(_,b) in enumerate(full) if k)<1,"inner radius")
        require(sum(b*OUTER**k for k,(_,b) in enumerate(full[:-1]))<OUTER**order,"outer radius")
    knots=tuple(RAYS)
    for left,right in zip(knots,knots[1:]):
        counts=[verify_ray(left,right,i,j,bounds[j]) for j in (6,7) for i in range(5)]
        print(left,right,counts,flush=True)


if __name__=="__main__":
    main()
