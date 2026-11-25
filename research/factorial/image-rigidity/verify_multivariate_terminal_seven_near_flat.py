#!/usr/bin/env python3
"""Exact joint Cayley ray exclusions at the terminal-seven flat endpoint."""
from fractions import Fraction as Q
from functools import lru_cache
from sys import argv
from verify_multivariate_terminal_six_intermediate_cayley import (
    GRID, I, ONE, ZERO, bisect, interval, projections, rational, require, scalar,
)
from verify_multivariate_terminal_six_marked_row_tail import projected_interval
from verify_multivariate_terminal_six_near_flat_extension import (
    check_constants as validate_sextic,
    coefficients as six_coeffs,
)

MAX_LAM=Q(9,40)
R0,R1=Q(513,500),Q(129,125)
RAYS=((Q(1073,1000),Q(3,4)),(Q(7,8),Q(1,3)),(Q(129,200),Q(1,5)),
      (Q(47,100),Q(1,10)),(Q(267,1000),Q(1,20)))


def positive(q):
    return interval(Q(0),q)


def uniform(q):
    return interval(-q,q)


def coefficients(order,lam,r):
    p=lam.square()
    if order==6:
        exact=six_coeffs(6,(Q(lam.lo,GRID),Q(lam.hi,GRID)),(Q(1),Q(1)))
        return tuple(enumerate(interval(a,b) for a,b in exact)),()
    q=r*p
    v=scalar(5,81)+positive(Q(22,1000))*q
    D=positive(Q(3,20))
    Egeo=interval(Q(-1,10),Q(1,10000))
    E=interval(Q(-1,10),Q(91,10000))
    A=positive(Q(56,1000))
    B=positive(Q(57,1000))
    q2,q3,q4=q.square(),q**3,q**4
    r2,r3,r4=r.square(),r**3,r**4
    k7=ONE-scalar(7,4)*q+scalar(35,48)*q2+(
        (scalar(42)*v+scalar(168)*D)*q3+(scalar(14)*A+scalar(42)*B+scalar(70)*E)*q4
    )*scalar(1,128)
    k5=r*(scalar(7,4)-scalar(35,24)*q+(
        (scalar(42)*v-scalar(336)*D)*q2+(scalar(126)*A+scalar(42)*B-scalar(210)*E)*q3
    )*scalar(1,128))
    k3=r2*(scalar(35,48)+(
        (-scalar(210)*v+scalar(168)*D)*q+(-scalar(70)*A-scalar(210)*B+scalar(210)*E)*q2
    )*scalar(1,128))
    k1=r3*(scalar(126)*v+(-scalar(70)*A+scalar(126)*B-scalar(70)*E)*q)*scalar(1,128)
    # Nominal paired means, with all cancellations applied before division by q.
    w=ONE-q
    y=ONE-scalar(5,3)*q+scalar(5,9)*q2+(v+D)*q3
    z=ONE-scalar(2)*q+q2+scalar(9,5)*D*q3+Egeo*q4
    wz=scalar(1)-q-scalar(9,5)*D*q2-Egeo*q3
    wy=scalar(2,3)-scalar(5,9)*q-(v+D)*q2
    difference=scalar(1,3)+(-scalar(3)*v+scalar(3,5)*D)*q+scalar(2)*Egeo*q2
    h6=r2*lam**3*(scalar(175)*w+scalar(315)*y+scalar(70)*z)*scalar(1,128)
    h4=r2*lam*(scalar(175)*w-scalar(525)*y-scalar(210)*z)*scalar(1,128)
    h2=-r3*lam*scalar(105,128)*(wy+scalar(2)*wz)
    h0=-r4*lam*scalar(35,128)*difference
    central=((1,k1),(3,k3),(5,k5),(7,k7),
             (0,uniform(Q(4,10000))*r4*lam),(2,uniform(Q(45,10000))*r3*lam),
             (4,uniform(Q(2,100000))*r2*lam),(6,uniform(Q(3,10000))*r2*lam**3))
    return central,((uniform(Q(134,1000)),((0,h0),(2,h2),(4,h4),(6,h6))),)


def height(lam,x,ray):
    p=lam.square()
    b=rational(ray[0])+rational(ray[1])*p
    num=p*b*x.square()
    y=I(0,b.hi)
    for _ in range(3):
        update=b-num/(ONE+p*b*y)
        y=I(max(y.lo,update.lo),min(y.hi,update.hi))
    return y


def verify_ray(order,index):
    initial=(interval(Q(0),MAX_LAM),interval(R0,R1) if order==7 else ONE,interval(Q(0),Q(2)))
    widths=tuple(max(1,v.hi-v.lo) for v in initial)

    @lru_cache(maxsize=16000)
    def shapes(alo,ahi,rlo,rhi):
        return coefficients(order,I(alo,ahi),I(rlo,rhi))

    nominal,family=shapes(initial[0].lo,initial[0].hi,initial[1].lo,initial[1].hi)
    magnitude=0
    upper=dict(nominal)
    for k in range(order):
        a=upper[k]
        coeff=max(abs(a.lo),abs(a.hi))
        for parm,values in family:
            b=(dict(values).get(k,ZERO))*parm
            coeff+=max(abs(b.lo),abs(b.hi))
        magnitude+=coeff*2**k
    require(upper[order].lo*2**order>magnitude,("scaled root disk",order))
    stack=[(initial,0)]
    count=0
    while stack:
        (lam,r,x),depth=stack.pop()
        count+=1
        nominal,family=shapes(lam.lo,lam.hi,r.lo,r.hi)
        y=height(lam,x,RAYS[index])
        for direction in projections(order,x,y):
            total=projected_interval(nominal,direction)
            for parameter,contribution in family:
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
    def six(v):
        return ((v-Q(3,2))*v+Q(9,20))*v-Q(3,200)
    def seven(v):
        return ((v-Q(7,4))*v+Q(35,48))*v-Q(35,576)
    vals6=((38,39),(356,357),(1105,1106))
    vals7=((111,112),(466,467),(1172,1173))
    intervals=[]
    for polynomial,values in ((six,vals6),(seven,vals7)):
        enclosures=[]
        for a,b in values:
            lo,hi=Q(a,1000),Q(b,1000)
            require(polynomial(lo)*polynomial(hi)<0,"anchor Hermite bracket")
            enclosures.append((lo,hi) if polynomial is six else (R0*lo,R1*hi))
        intervals.append(enclosures)
    a,b=intervals
    sequence=(b[2],a[2],b[1],a[1],b[0],a[0])
    for (high,low),(ray,_) in zip(zip(sequence,sequence[1:]),RAYS,strict=True):
        require(high[0]>ray**2>low[1],"anchor alternating sectors")



def validate_constants():
    validate_sextic()
    q = Q(11,200)
    u = 1 - MAX_LAM**2
    require(R1 * MAX_LAM**2 < q, "septic deficit cap")
    require((u**-2-1)/5 < Q(11,500)
            and Q(36,5)/(7-Q(11,500)) < R1, "upper deficit ratio")
    require((u*(1+u)/70)**2 > (R0-1)**2*(1+(1-u*u)/35)**3,
            "lower deficit ratio")
    x = Q(365,1000)
    A,B,C = Q(334,1000),Q(335,1000),Q(134,1000)
    L,K = Q(54,1000),Q(21,1000)
    require((2-q)/(6*(1-q)**2)<x, "septic reciprocal cap")
    require(Q(5,2)*x*x<A and x*x<C and 3*x**4<L and Q(60,7)*x**6<K
            and A+L*q*q<B, "septic log-skew bounds")
    require(Q(501,1000)*A*A<Q(56,1000) and Q(501,1000)*B*B<Q(57,1000)
            and Q(501,1000)*C*C<Q(9,1000), "even hyperbolic corrections")
    b = Q(10,81)/(1-Q(4,9)*(q/(1-q))**2)/(1-q)**3
    D = Q(3,20)
    require(b/(1-b*q**3)<D, "septic paired mean defect")
    require(Q(5,243)/(1-q)<Q(22,1000), "pure fifth-third power tail")
    require(Q(1,3)**3 < Q(7,20)**3*(1-q)**2 and
            Q(18,25)**3<Q(4,5)**3*(1-q)**4, "central Taylor derivatives")
    require(Q(9,5)*Q(7,20)*D+q/(25*(1-q)**6)<Q(1,10)
            and Q(4,5)*D*D*q*q<Q(1,10000), "linked central mean")
    require(max(A,B,C)*q*q < Q(1,100), "small hyperbolic arguments")
    for weights,left,right,divisor,cap in (
        ((70,126,70),168*q*q,98*q**4,q**4,Q(3,10000)),
        ((70,210,210),70*q**3,105*q**4,q**4,Q(2,100000)),
        ((126,42,210),189*q,63*q**3,q**3,Q(45,10000)),
        ((14,42,70),14*q,14*q*q,q*q,Q(4,10000)),
    ):
        cubic = sum(w*a**3 for w,a in zip(weights,(A,B,C),strict=True))/5
        require((left*L+right*K+divisor*cubic)/128<cap, "septic odd-part remainder")
    for p in (Q(0),MAX_LAM**2):
        b = tuple(start+slope*p for start,slope in RAYS)
        require(all(x>y>0 for x,y in zip(b,b[1:])), "strict ray order")


def main():
    validate_constants()
    anchor()
    for order in (6,7):
        if len(argv)>1 and order not in map(int,argv[1:]):
            continue
        for index in range(5):
            print(order,index,verify_ray(order,index),flush=True)


if __name__=="__main__":
    main()
