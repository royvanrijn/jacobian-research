#!/usr/bin/env python3
"""Exact light shared terminal-eight continuation on .70 <= u <= .72."""
from fractions import Fraction as Q

from verify_multivariate_terminal_eight_paired_bridge import RAYS, radial
from verify_multivariate_terminal_eight_paired_light import (
    check_anchor, coefficient_model, envelope, validate_skew,
)
from verify_multivariate_terminal_six_intermediate_cayley import I, interval, require


LEFT, RIGHT = Q(7,10), Q(18,25)
START=RAYS[700]
FINISH=(-347,-447,-703,-770,-920,-946)
CAP=Q(9,200)


def ray_cosines(index,u):
    a,b=START[index],FINISH[index]
    v=tuple(Q(a,1000)+Q(b-a,20)*(endpoint-LEFT) for endpoint in u)
    return interval(min(v),max(v))


def verify(order,index,env,lower_bound):
    stack=[((LEFT,RIGHT),(Q(0),Q(1)),0,0)]
    count=0
    while stack:
        u,normalized,du,ds=stack.pop()
        coefficients=coefficient_model(order,u,normalized,env,lower_bound)
        if coefficients is None:
            continue
        used=radial(order,ray_cosines(index,u),coefficients)
        if used is not None:
            count+=used
            continue
        require(max(du,ds)<15,("unresolved ray",order,index,u,normalized))
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
    for row in (START,FINISH):
        require(all(1000>a>b>-1000 for a,b in zip(row,row[1:])),"strict ray order")
    env={j:envelope(CAP,j) for j in (7,8)}
    upper=env[7][1][0][1].hi
    lower_bound=I(upper,upper)
    for j in (7,8):
        for i in range(6):
            print(j,i,verify(j,i,env[j],lower_bound),flush=True)


if __name__=="__main__":
    main()
