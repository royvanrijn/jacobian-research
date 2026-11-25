#!/usr/bin/env python3
"""Exact scaled Cayley replay at the fully shared terminal-eight flat endpoint."""
from fractions import Fraction as Q
from functools import lru_cache
from math import comb
from sys import argv

from verify_multivariate_terminal_six_intermediate_cayley import (
    GRID, I, ONE, ZERO, bisect, interval, projections, rational, require, scalar,
)
from verify_multivariate_terminal_six_marked_row_tail import projected_interval
from verify_multivariate_terminal_seven_near_flat import (
    coefficients as septic_coefficients, height, uniform, positive,
    validate_constants as validate_septic,
)

MAX_LAM=Q(9,40)
R0,R1=Q(1019,1000),Q(128,125)
RAYS=((Q(1101,1000),Q(13,20)),(Q(921,1000),Q(2,5)),
      (Q(719,1000),Q(1,4)),(Q(14,25),Q(1,5)),
      (Q(389,1000),Q(7,100)),(Q(6,25),Q(1,20)))
EXPONENTS=(Q(0),Q(1),Q(12,7),Q(15,7),Q(16,7))


def cayley(k,i):
    return Q(comb(8,i),256)*sum(
        (-1)**(i-l)*comb(i,l)*comb(8-i,k-l)
        for l in range(max(0,k-8+i),min(k,i)+1)
    )


@lru_cache(maxsize=None)
def parity_weights(k):
    if k%2==0:
        return tuple((2 if i<4 else 1)*cayley(k,i) for i in range(5))
    return tuple(2*cayley(k,i) if 1<=i<=3 else Q(0) for i in range(5))


@lru_cache(maxsize=None)
def powers_series(a):
    coeff=[Q(1)]
    for k in range(1,10):
        coeff.append(coeff[-1]*(k-1-a)/k)
    return tuple(coeff)


@lru_cache(maxsize=None)
def series(weights):
    terms=[powers_series(e) for e in EXPONENTS]
    coefficients=[sum(w*a[k] for w,a in zip(weights,terms,strict=True)) for k in range(9)]
    first=next((k for k,c in enumerate(coefficients) if c),None)
    if first is None:
        require(not any(weights),"identically zero Cayley family")
        return None
    tail=sum(abs(w*a[9]) for w,a in zip(weights,terms,strict=True))
    return first,tuple(rational(c) for c in coefficients[first:]),rational(tail)


def scaled_family(q,lam,r,q_power,vector,even):
    terms=[]
    for k in range(0 if even else 1,9,2):
        weights=tuple(a*b for a,b in zip(parity_weights(k),vector,strict=True))
        expansion=series(weights)
        if expansion is None:
            continue
        offset,coefficients,tail=expansion
        total=coefficients[-1]
        for coeff in reversed(coefficients[:-1]):
            total=coeff+q*total
        remainder=(tail*q**(9-offset)/(ONE-q)).symmetric_error()
        degree=2*(q_power+offset)-(8-k)
        require(degree>=0,("scaled Cayley divisibility",k,q_power,vector))
        terms.append((k,(total+remainder)*r**(q_power+offset)*lam**degree))
    return tuple(terms)


def coefficients(order,lam,r):
    if order==7:
        return septic_coefficients(7,lam,ONE)
    q=r*lam.square()
    nominal=dict(scaled_family(q,lam,r,0,(1,1,1,1,1),True))
    for k in range(9):
        even=k%2==0
        cap=(Q(0),Q(28,1000),Q(27,1000),Q(42,1000),Q(58,1000)) if even else (
            Q(0),Q(1,1000),Q(57,1000),Q(70,1000),Q(0)
        )
        weight=sum(abs(a)*b for a,b in zip(parity_weights(k),cap,strict=True))
        raw_power=6 if even else 5
        error=uniform(weight)*r**raw_power*lam**(2*raw_power-(8-k))
        nominal[k]=nominal.get(k,ZERO)+error
    families=(
        (positive(Q(14,100)),scaled_family(q,lam,r,3,(0,0,1,2,Q(12,5)),True)),
        (interval(Q(-31,1000),Q(0)),scaled_family(q,lam,r,5,(0,0,0,1,0),True)),
        (interval(Q(-5,100),Q(0)),scaled_family(q,lam,r,5,(0,0,0,0,1),True)),
        (positive(Q(244,1000)**2),scaled_family(q,lam,r,4,(0,Q(49,50),Q(32,25),Q(1,2),0),True)),
        (uniform(Q(244,1000)),scaled_family(q,lam,r,2,(0,Q(7,5),Q(8,5),1,0),False)),
        (uniform(Q(399,1000)),scaled_family(q,lam,r,4,(0,Q(1,5),Q(3,35),0,0),False)),
        (uniform(Q(56,1000)),scaled_family(q,lam,r,6,(0,0,Q(-1,7),0,0),False)),
    )
    return tuple(nominal.items()),families


def verify_ray(order,index):
    initial=(interval(Q(0),MAX_LAM),interval(R0,R1) if order==8 else ONE,interval(Q(0),Q(2)))
    widths=tuple(max(1,v.hi-v.lo) for v in initial)

    @lru_cache(maxsize=16000)
    def shapes(alo,ahi,rlo,rhi):
        return coefficients(order,I(alo,ahi),I(rlo,rhi))

    base,families=shapes(initial[0].lo,initial[0].hi,initial[1].lo,initial[1].hi)
    enclosure=dict(base)
    for uncertainty,vector in families:
        for k,a in vector:
            enclosure[k]=enclosure[k]+uncertainty*a
    lead=enclosure[order].lo*2**order
    others=sum(max(abs(a.lo),abs(a.hi))*2**k for k,a in enclosure.items() if k!=order)
    require(lead>others,("scaled root disk",order))
    stack=[(initial,0)]
    count=0
    while stack:
        (lam,r,x),depth=stack.pop()
        count+=1
        base,families=shapes(lam.lo,lam.hi,r.lo,r.hi)
        y=height(lam,x,RAYS[index])
        for direction in projections(order,x,y):
            total=projected_interval(base,direction)
            for parameter,contribution in families:
                total=total+parameter*projected_interval(contribution,direction)
            if total.lo>0 or total.hi<0:
                break
        else:
            require(depth<38,("unresolved flat ray",order,index,[v.shown() for v in (lam,r,x)]))
            box=(lam,r,x)
            weights=(100,1 if depth>=14 else 0,50)
            k=max(range(3),key=lambda j: weights[j]*(box[j].hi-box[j].lo)*GRID//widths[j])
            for half in bisect(box[k]):
                stack.append((box[:k]+(half,)+box[k+1:],depth+1))
    return count


def anchor():
    def seven(v):
        return ((v-Q(7,4))*v+Q(35,48))*v-Q(35,576)
    def eight(v):
        return (((v-2)*v+Q(15,14))*v-Q(15,98))*v+Q(15,5488)
    vals7=((111,112),(466,467),(1172,1173))
    vals8=((20,21),(191,192),(560,561),(1226,1227))
    intervals=[]
    for polynomial,values in ((seven,vals7),(eight,vals8)):
        enclosures=[]
        for a,b in values:
            lo,hi=Q(a,1000),Q(b,1000)
            require(polynomial(lo)*polynomial(hi)<0,"anchor Hermite bracket")
            enclosures.append((lo,hi) if polynomial is seven else (R0*lo,R1*hi))
        intervals.append(enclosures)
    a,b=intervals
    sequence=(b[3],a[2],b[2],a[1],b[1],a[0],b[0])
    for (high,low),(ray,_) in zip(zip(sequence,sequence[1:]),RAYS,strict=True):
        require(high[0]>ray**2>low[1],"anchor alternating sectors")



def validate_constants():
    validate_septic()
    q=Q(11,200)
    u=1-MAX_LAM**2
    require(R1*MAX_LAM**2<q and 2*MAX_LAM<1,"flat disk transform")
    require((u**-2-1)/6<Q(37,2000) and Q(49,6)/(8-Q(37,2000))<R1,
            "upper deficit ratio")
    require((u*(1+u)/96)**2>(R0-1)**2*(1+(1-u*u)/48)**3,
            "lower deficit ratio")
    x=Q(312,1000)
    g,L,K=Q(244,1000),Q(399,1000),Q(56,1000)
    require((2-q)/(7*(1-q)**2)<x and Q(5,2)*x*x<g and 42*x**4<L and 60*x**6<K,
            "octic reciprocal and linked skew bounds")
    s=q/(1-q)
    h=Q(40,343)/(1-Q(25,49)*s*s)/(1-q)**3
    J3=Q(120,343)*Q(16,49)/(5*(1-Q(25,49)*s*s)*(1-Q(9,49)*s*s))/(1-q)**5
    J4=J3+Q(48,343)*Q(24,49)/(5*(1-Q(25,49)*s*s)*(1-Q(1,49)*s*s))/(1-q)**5
    require(h<Q(14,100) and J3<Q(31,1000) and J4<Q(5,100),"octic mean defect bounds")
    gamma_bases=(Q(7,5)*g,Q(8,5)*g,g)
    skew_caps=(Q(343,1000),Q(391,1000),g)
    shift_caps=(L/5,Q(3,35)*L+K*q*q/7,Q(0))
    shape_caps=(Q(0),Q(14,100),Q(28,100),Q(336,1000))
    even_caps=(Q(28,1000),Q(27,1000),Q(42,1000),Q(58,1000))
    odd_caps=(Q(1,1000),Q(57,1000),Q(70,1000))
    z=Q(1,100)
    require(max(skew_caps)*q*q<z and
            Q(1,2)/(1-z*z/12)<Q(501,1000) and
            Q(1,24)/(1-z*z/30)<Q(1,20) and Q(1,6)/(1-z*z/20)<Q(1,5),
            "hyperbolic Taylor bounds")
    for k,mean in enumerate(shape_caps):
        shape=mean/(1-mean*q**3)
        exp_remainder=mean*shape/2
        if k==3:
            require(exp_remainder<even_caps[k],"central exponential error")
            continue
        A=skew_caps[k]
        require(gamma_bases[k]+shift_caps[k]*q*q<=A,"individual skew magnitude")
        square_error=shift_caps[k]*(A+gamma_bases[k])/2
        cosh_error=A**4*q*q/20+shape*Q(501,1000)*A*A*q
        require(exp_remainder+square_error+cosh_error<even_caps[k],"even error cap")
        require(shape*A+A**3*q/5<odd_caps[k],"odd error cap")
    for p in (Q(0),MAX_LAM**2):
        b=tuple(a+d*p for a,d in RAYS)
        require(all(a>c>0 for a,c in zip(b,b[1:])),"scaled ray order")


def main():
    validate_constants()
    anchor()
    for order in (7,8):
        if len(argv)>1 and order not in map(int,argv[1:]):
            continue
        for index in range(6):
            print(order,index,verify_ray(order,index),flush=True)


if __name__=="__main__":
    main()
