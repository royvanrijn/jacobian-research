#!/usr/bin/env python3
"""Exact-rational endpoint certificate for (T42)–(T43) in the ordinary binomial proof."""
from fractions import Fraction as Q

class Interval:
    """Exact rational interval."""

    def __init__(self, lo, hi=None):
        if isinstance(lo, Interval):
            self.lo, self.hi = lo.lo, lo.hi
        else:
            self.lo = Q(lo)
            self.hi = self.lo if hi is None else Q(hi)
    def __add__(self, other):
        other = Interval(other)
        return Interval(self.lo + other.lo, self.hi + other.hi)
    __radd__ = __add__
    def __neg__(self):
        return Interval(-self.hi, -self.lo)
    def __sub__(self, other):
        return self + -Interval(other)
    def __rsub__(self, other):
        return Interval(other) + -self
    def __mul__(self, other):
        other = Interval(other)
        products = [a*b for a in [self.lo, self.hi]
                          for b in [other.lo, other.hi]]
        return Interval(min(products), max(products))
    __rmul__ = __mul__
    def __truediv__(self, other):
        other = Interval(other)
        assert other.lo * other.hi > 0
        return self * Interval(1 / other.hi, 1 / other.lo)
    def __rtruediv__(self, other):
        return Interval(other) / self

def log_base(r):
    assert 1 <= r <= 2
    z = (r - 1) / (r + 1)
    lower = 2 * sum(z ** (2*j+1) / (2*j+1) for j in range(14))
    return Interval(lower, lower + 2*z**29/(29*(1-z*z)))

LOG2 = log_base(Q(2))
LOG4 = 2*LOG2
PI = Interval("3.14159", "3.14160")

def log_rational(r):
    r = Q(r)
    assert r > 0
    exponent = 0
    while r < 1:
        r *= 2
        exponent -= 1
    while r > 2:
        r /= 2
        exponent += 1
    return log_base(r) + exponent*LOG2

def log(r):
    r = Interval(r)
    return Interval(log_rational(r.lo).lo, log_rational(r.hi).hi)

def root(r, n):
    r = Q(r)
    assert r > 0
    scale = 10**13
    target = (r.numerator*scale**n) // r.denominator
    lo, hi = 0, 1
    while hi**n <= target:
        hi *= 2
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if mid**n <= target:
            lo = mid
        else:
            hi = mid
    return Interval(Q(lo, scale), Q(lo+1, scale))

ROWS = [
    (3150, 3295, 282, 140, Q("2.504"), Q("2.885"), Q(".023")),
    (3296, 3467, 288, 143, Q("2.506"), Q("2.886"), Q(".146")),
    (3468, 4150, 292, 145, Q("2.514"), Q("2.889"), Q(".31")),
    (4151, 5000, 320, 159, Q("2.522"), Q("2.892"), Q(".44")),
    (5001, 6000, 354, 176, Q("2.531"), Q("2.896"), Q(".86")),
    (6001, 10000, 456, 227, Q("2.58"), Q("3.01"), Q(".59")),
    (10001, 19999, 468, 233, Q("2.58"), Q("2.92"), Q("2.37")),
]
SMALL = [(2, Q(1,2), Q(".2")), (3, Q(2,3), Q(".435")),
         (9, Q(".75"), Q(".68")), (19, Q(".9"), Q(".925"))]

def bounds(N, B, D, vp, digit, gain=None):
    U, w = log(N), log(B)
    v = 2*U+w
    z = LOG4+v/3
    T = U+w-log(D)
    P = 4*root(N*N*B, 3)
    xi = T+Q(".337")
    R = (z-Q(".49")+(P/N)*(digit+2*(v+Q(".001"))/64)+Q(".71")*v/64
         +Q("2.004")*D*xi/P+2*vp*D*xi/N)
    if gain is not None:
        f, cost = gain
        L = U+f*T-cost
    else:
        A = Q(D, N)/B
        L = (U-Q(".001")+(1-Q(1,D))*(T-Q(".6932"))
             +(2*log(2*D/PI)-Q("1.01"))/D-Q("1.159")*A
             -(Q(".641")*(Q("1.7032")+Q(".128")*D*A)+Q(".883"))*root(A,2))
    return L-R

def main():
    assert 354*(Q(6000*(5*6000-11),6)+2) < 10726904659
    assert Q(5*176*5000**2,6)+353 < 3842610773
    assert Q(5*159*4150**2,6)+319 < 2300942549
    assert Q(5*145*3467**2,6)+291 < 1453168141
    assert Q(5*143*3295**2,6)+287 < 1294268491
    assert (4*root(3000**2,3)).lo > 830
    assert 1-Q(1,3000)-Q(2*467,3000**2) > Q(".9995")
    assert Q(1)/(Q(".9995")-Q(1,830)) < Q("1.002")
    assert Q("1.002")*(Q(2,3)+Q(1,24)) < Q(".71")
    assert LOG2.hi < Q(".6932")
    assert Q(1)/Q(".97") < Q("1.031")
    assert (log(Q(2)/Q(".97"))).hi < Q(".7237")
    assert (2*Q("1.006")/PI).hi < Q(".641")
    assert (2*LOG2/PI*2).hi < Q(".883")
    assert (Q(".883")*root(Q(".03"),2)).hi < Q(".153")
    assert (Q(".641")*root(Q(".03"),2)).hi < Q(".112")
    assert (PI.lo/18-(PI.hi/18)**3/6)**2 > Q(".03")
    loss2 = (1+Q(".01")+Q(".00384")*2-log(Q(".97"))-log(Q("1.94")))/2+Q(".001")
    loss3 = (Q(2,3)*Q(".7237")-(2*log(6/PI)-Q("1.02152"))/3
             +2*LOG2*Q(".03")*Q(8,9)+Q(".001"))
    loss49 = Q(".75")*Q(".7237")+2*LOG2*Q(".03")*Q(80,27)+Q(".001")
    assert (2*log(8/PI)).lo > Q("1.045")
    loss1019 = (Q(".7237")+Q(".1083")+Q(".112")*(Q("1.083")+Q(".6932"))
                +Q(".153")+Q(".001")-2*log(38/PI)/19)
    assert loss2.hi < Q(".2") and loss3.hi < Q(".435")
    assert loss49.hi < Q(".68") and loss1019.hi < Q(".925")
    for N, Nmax, G, D, vp, digit, margin in ROWS:
        assert D == G//2-1
        assert N >= 1001
        assert (4*root(N*N,3)).lo > max(Q(5*G,2),830)
        U = log(N)
        xi = U-log(D)+Q(".337")
        assert xi.lo > 3
        # Convexity uniformly in w and N. This also handles each smaller D.
        curv = 16*root(N,3)*(digit+2*(2*U+Q(".001")+6)/64) - 2*Q("1.002")*D*(6-xi)
        assert curv.lo > 0
        # The digit ratio is increasing in z. Check its maximizing corner.
        zmax = LOG4 + (2*log(Nmax)+log(Q(5*D,6)))/3
        k = 3*LOG4-Q(".001")
        assert k.lo > 3*Q("1.2762")
        amax = 3-k/zmax
        assert amax.hi < vp
        assert (amax*(1+Q("1.2762")/zmax)).hi < digit
        endpoints = [Q(1), Q(5*D,6)]
        vals = [bounds(N,B,D,vp,digit) for B in endpoints]
        assert min(x.lo for x in vals) > margin
        print("large d", N, [float(x.lo) for x in vals])
        for small_D, f, cost in SMALL:
            vals = [bounds(N,B,small_D,vp,digit,(f,cost))
                    for B in [Q(1),Q(5*small_D,6)]]
            assert min(x.lo for x in vals) > Q("3.8")
            print("small d", N, small_D, [float(x.lo) for x in vals])

if __name__ == "__main__":
    main()

