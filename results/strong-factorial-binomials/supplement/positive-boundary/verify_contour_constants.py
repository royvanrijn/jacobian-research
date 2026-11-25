#!/usr/bin/env python3
"""Rational checks for the gap-six contours and finite-domain cutoffs.

The contour checks retain the inequalities from the gap-six source.
Only standard-library arithmetic is needed to replay this proof route.
"""
from fractions import Fraction as F
from math import factorial

if not __debug__:
    raise SystemExit("Run without -O so that all proof assertions execute.")


def exp_bounds(x):
    x = F(x)
    assert abs(x) < 122
    # For x>=0, bound the tail after the degree-120 Taylor sum by a geometric series.
    if x<0:
        lo,hi=exp_bounds(-x)
        return 1/hi,1/lo
    lo=sum(x**j/factorial(j) for j in range(121))
    hi=lo+(x**121/factorial(121))/(1-x/122)
    return lo,hi
def check_contour_constants():
    jmin,jmax=F(34,5),F(36,5)
    kmin,kmax=F(99,100),F(101,100)
    dj,dk=F(1,5),F(1,100)
    rj,rk=F(3,4),F(1,5)
    cj,ck=F(7,200),F(1,4)
    bj,bk=F(1,30),F(5,6)
    aj,ak=F(3,25),F(9,10)
    amaxj,amaxk=F(1,5),F(11,10)
    xmaxj,xminj=F(9),F(28,5)
    assert dj*jmin**6/(6+dj)>3180
    assert kmax**6*(6+dk)/1000<dk
    jvar=dj*(6*jmax-7)/(7*jmin**2)
    kvar=7*dk/kmin**2
    assert F(1,7)-jvar>aj and F(1,7)+jvar<amaxj
    assert 1-kvar>ak and 1+kvar<amaxk
    def log_tail(r):
        return r*r/(2*(1-r))
    assert dj/6+7*log_tail(dj/7)+log_tail(dj/6)<F(1,25)
    assert exp_bounds(F(681,350))[1]<7
    assert log_tail(dk)+ (kmax**6/1000)/(1-kmax**6/1000)<F(1,500)
    def j_env(x):
        return x*(3200+x**6)*exp_bounds(-x)[1]
    assert max(j_env(F(3)),j_env(xminj))<705
    assert 705*(1+dj*dj/18)**7<exp_bounds(F(263,40))[0]
    assert exp_bounds(F(329,50))[0]>720
    # Taylor derivative bounds for log terms at the two saddles.
    rho=3200/(jmin-rj)**6
    y=rho/(1-rho)
    assert (14+336*y+756*y*y+432*y**3)/(jmin-rj)**3<F(1,5)
    rho=(kmax+rk)**6/1000
    y=rho/(1-rho)
    assert (2+120*y+540*y*y+432*y**3)/(kmin-rk)**3<5
    assert aj/2-bj*rj>=cj and ak/2-bk*rk>=ck
    def curvature_log(terms,rho,r):
        a,b=terms; t=rho/(1-rho)
        return (1-r*r)/(1+r*r)**2-a*t-b*t*t
    assert (7*curvature_log((0,0),F(0),dj/xminj)
            - (42*3200/(xminj**6-3200)
               + 36*(3200/(xminj**6-3200))**2))>F(3,2)
    assert curvature_log((30,36),F(151,100)**6/1000,2*dk)>F(1,2)
    assert F(7,9)+F(6,9)*3200/(9**6-3200)-1<-F(1,5)
    # Envelope comparisons for the short contour and the outer/inner domains.
    def k_env(x):
        return x*(1+x**6/1000)*exp_bounds(-x)[1]
    assert max(k_env(F(3,2)),k_env(F(3)))*(1+dk*dk/F(9,2))**7<F(17,50)
    assert (F(1,2)+dk)*(1+(F(1,2)+dk)**6/1000)*exp_bounds(-F(1,2))[1]<F(17,50)
    assert F(17,50)<exp_bounds(-F(501,500)-F(1,225))[0]
    assert F(9,25)<exp_bounds(-F(501,500))[0]
    assert 1100*exp_bounds(F(-1))[1]+3**7*exp_bounds(F(-3))[1]<520
    assert 3*exp_bounds(F(-3))[1]+F(4**7,3080)*exp_bounds(F(-4))[1]<F(8,25)
    assert 4*exp_bounds(F(-4))[1]+F(7**7,3080)*exp_bounds(F(-7))[1]<F(8,25)
    def local_bound(s,r,b,c):
        return b/((s-r)*c*c)+1/(s*(s-r)*c)+1/(s*c*r)
    assert local_bound(jmin,rj,bj,cj)<11
    assert local_bound(kmin,rk,bk,ck)<43
    assert 3*factorial(4)*200**4<10000**3
    assert 5*factorial(4)*225**4<10000**3
    assert 6/amaxj/jmax**2>F(3,4)**2
    assert 6/amaxk/kmax**2>4
    assert F(13)/(F(3,4))<25 and F(44)/2<25
    assert F(1,3102)-1/(jmin**6-3100)>F(1,3300)
    assert kmax**6<2
    for s,v,a in [(jmin,6+dj,aj),(kmin,dk,ak)]:
        vprime=v*v/(6*s**7*a)
        assert vprime*(1/s+(7/s**2+14/s**3)/(2*a))<F(1,50)
    assert F(1,50)+F(1,300)<F(1,20)
    assert F(10000,3300)>F(1,10)

def log_bounds(x):
    """Enclose log(x) after reducing to [1, 2]."""
    def reduced(z):
        u = (z - 1) / (z + 1)
        power = u
        lo = F(0)
        for j in range(24):
            lo += 2 * power / (2 * j + 1)
            power *= u * u
        return lo, lo + 2 * power / (49 * (1 - u * u))

    assert x >= 1
    x = F(x)
    scale = 0
    while x >= 2:
        x /= 2
        scale += 1
    lo, hi = reduced(x)
    two_lo, two_hi = reduced(F(2))
    return lo + scale * two_lo, hi + scale * two_hi


def check_cutoffs():
    def log_upper(x):
        return log_bounds(x)[1]

    assert 2 * log_upper(199) + 3 < 14
    assert 2 * log_upper(9999) + 2 < 21
    assert F(75, 4) - 2 * log_upper(F(512 * 25, 3)) - 2 > F(3, 100)
    assert F(400, 21) - 2 * log_upper(5000) - 2 > F(1, 100)
    assert F(39**2, 40) - 2 * log_upper(128 * 40) - 2 > 18

    def count(bounds, ratio_filter):
        return sum(
            max(0, upper - 25)
            for a, d, upper in bounds
            if not ratio_filter or d * d < 21 * (a + d)
        )

    families = [
        [(1, 6, 10000)],
        [(1, d, 5000) for d in range(10, 20)],
        [
            (a, d, max(4000, 128 * (a + d)))
            for d in range(6, 25)
            for a in range(2, d // 3 + 1)
        ],
    ]
    assert sum(count(family, False) for family in families) == 338071
    assert [count(family, True) for family in families] == [9975, 49750, 266421]


if __name__ == "__main__":
    check_contour_constants()
    check_cutoffs()
    print("Gap-six contour constants and all positive-boundary cutoffs passed.")
