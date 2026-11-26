# A one-sided phase bound for private exponent totals one and two

Let the two private exponent totals of an ordinary two-monomial
profile be one and two, up to exchange of the columns. Let
\[
                   S_{\rm sh}=\tau+\sum_{\rm shared}
                              \frac{(b_i-a_i)^2}{\max(a_i,b_i)}
\]
include any common normalized Gaussian strength \(\tau\ge0\).

**Theorem.** Every adjacent pair with terminal \(N\ge14\) has
disjoint root arguments under either of these shared reserve bounds.

- \(S_{\rm sh}\ge19/50\) if the exponent total two comes from one
  private row of exponent two.
- \(S_{\rm sh}\ge79/100\) if it comes from two private unit rows.

The [complete terminal-thirteen theorem](MULTIVARIATE_TERMINAL_THIRTEEN_ALL_PROFILES.md)
gives the all-degree ordinary result. The matched curvature statistic
is \(\sum 2u/(1+u)\) over the shared slots, where
\(u=(n/T'-(n-1)/T)^{-1}\). The bounds allow arbitrary orientations
on the shared rows.

The [existing shared reserve proof](MULTIVARIATE_PRIVATE_ONE_TWO_SHARED_ANGLE_RESERVE.md)
gives the conclusion for \(S_{\rm sh}\ge83/80\), so assume
\(S_{\rm sh}\le83/80\). This proof uses the exact shape of the
right half to give a one-sided sign throughout its root annulus.

## Monotone weighted coefficients in the right half

Orient the private side of total two toward index \(k\). After
reserving \(k!\) and \((m-k)!\), the right endpoint contributions
from the remaining private factor are
\[
               f_{m,m-k}^{(a)}=\frac{(a)_{m-k}}{(a)_m^{1-k/m}},
                          \qquad a=\tfrac12\ \text{or}\ 1.
\]
We work at \(m\ge13\) and radii \(R_m=e^{7/(10m)}\).
For \(U_m(z)=1+\sum_{k\ge1}u_k z^k\) the right half, including
half weight on an even middle term, we show
\[
                    k u_k R_m^k\ge(k+1)u_{k+1}R_m^{k+1}.       \tag{1}
\]
This also holds at all smaller radii.

For any shared row of gap \(d\ge1\) and smaller endpoint \(\ell\),
\(d/\ell\le1\) in the present range. If \(d>\ell\ge1\), then
\(d\ge2\) and its curvature is at least \(d^2/(2d-1)\ge4/3\).
Let \(S_i=d^2/(\ell+d)\) and \(\varphi(x)\) be its
endpoint-normalized log-factorial contribution in right-half index
\(x\in[0,m]\). The trigamma integral gives
\[
 \varphi''(x)\ge\frac{d^2}{m(\ell+d)+1},\qquad
             |\varphi'''(x)|\le\frac{d^3}{m^2\ell^2}.
\]
These follow from lower and upper integral comparisons for the
positive trigamma and second polygamma series respectively.
Since \(\varphi'(m/2)\) is the difference of the unnormalized
midpoint derivative and its mean, Taylor expansion and averaging
about the midpoint give
\[
                  |\varphi'(m/2)|\le\frac{d^3}{24\ell^2}
                                                  \le S_i/12 .
\]
Therefore \(\varphi'(x)\le S_i/12\) throughout the right half.
In addition
\(\varphi'(2)\le S_i/12-(m/2-2)d^2/[m(\ell+d)+1]<0\),
since \(m\ge13\) and \(\ell+d\ge2\). The Gaussian has a
nonpositive derivative on the same half and can be omitted
in these upper bounds. Shared factors therefore multiply the
ratio in (1) by at most one when \(k=1\), and by at most
\(e^{83/960}\) everywhere else.

The private factor in the successive weighted ratio, before
the possible even middle halving, is
\[
             \frac{k+1}{k}\frac{R_m(a)_m^{1/m}}{m+a-k-1}.
\]
For \(k+1\le m/2\), this increases with \(a\), since
\[
    \frac1m\sum_{h=0}^{m-1}\frac1{a+h}
      \ge \frac1{a+(m-1)/2}\ge\frac1{m+a-k-1}.
\]
Use \(a=1\) and put \(G_m=R_m(m!)^{1/m}/m\), which decreases
by the fixed-index descending cap comparison in the preceding
shared reserve proof. At \(k=1\) the ratio is at most
\(26G_{13}/12<1\). For \(m=13,14\) and \(k\ge2\), the
factor \((1+1/k)m/(m-k)\), after halving at an even middle
term, is at most \(39/20\). For \(m\ge15\) and \(k\ge2\) it is
at most \(450/221\). In fact its log is convex on
\(2\le k\le m/2-1\), so the maximum is at an endpoint.
The two endpoint values are \(3m/[2(m-2)]\) and
\(2m^2/(m^2-4)\), both decreasing in \(m\).
The three scalar bounds used here, including
\((39/20)G_{13}e^{83/960}<1\) and
\((450/221)G_{15}e^{83/960}<1\), are checked below.
This proves (1).

## The right phase sign and the actual root annulus

For \(0<\theta<\pi\) the classical Fejér sine and Young cosine
inequalities give
\(\sum_{k=1}^j\sin(k\theta)/k>0\) and
\(\sum_{k=1}^j\cos(k\theta)/k>-1\), respectively.
The weak endpoint versions follow by continuity.
[Brown and Yin, equations (1.1) and (1.2)](https://acta.hu/download.phtml?id=2604)

Apply Abel summation with the nonincreasing weights
\(k u_k r^k\), for \(0\le r\le R_m\), and with a final
weight zero. The resulting inequalities are
\[
       \Im U_m(r e^{i\theta})\ge0,\qquad
               \Re U_m(r e^{i\theta})\ge1-u_1 r.                \tag{2}
\]
In particular the right half is zero-free throughout the disk
as soon as \(u_1 R_m<1\). This real part bound is sharper than
the bound obtained by summing all the coefficient magnitudes.

Use the caps \(B_{H,k}\) from the
[shared reserve proof](MULTIVARIATE_PRIVATE_ONE_TWO_SHARED_ANGLE_RESERVE.md#the-exact-private-residual-and-its-two-envelopes)
with \(s=19/50\) if \(a=1/2\), and \(s=79/100\) if \(a=1\).
Write \(A_H=\sum B_{H,k}\).
The exact scalar bounds below give \(A_L<1/2\),
\(B_{U,1}<1/2\), and
\[
     \frac{1+A_L}{1-B_{U,1}}<e^{7/10},\qquad
                         \frac{1+A_U}{1-A_L}<e^{7/10}.
\]
On the outer circle use the first inequality and (2) to dominate
the left half by \(z^m U_m(1/z)\). For the reversal use the
second inequality and the elementary lower bound on the left
half. Rouché puts all roots in \(R_m^{-1}<|z|<R_m\).

## Comparing the full phases

For \(a=1/2\) use upper slopes \((.262,.671)\) for
\(H=L,U\), with lower slope \(.054\) for \(L\).
For \(a=1\) use upper slopes \((.312,.613)\),
with lower slope \(.066\) for \(L\).
Write the upper slopes as \(\alpha_H\) and the left lower
slope as \(\beta\). The signed weighted scalar bounds are
\[
   \sum_k(k-\alpha_H)B_{H,k}<\alpha_H,\qquad
 B_{L,1}+D(\beta)\sum_{k\ge2}(k+\beta)B_{L,k}<1,\qquad
 D(\beta)=\frac{2}{9\beta(1-(22\beta/7)^2/6)} .
\]
They give \(-\beta\theta<\arg L_m<\alpha_L\theta\)
on the ray at all radii in the disk, by the signed sine
bounds of the shared reserve proof. The right half has
\(0\le\arg U_m<\alpha_U\theta\) by (2) and the upper signed
sine bound. For the latter the sinc comparison remains
valid for \(\alpha_U\le7/10\): the seventh-degree sine Taylor
lower bound divided by its argument is greater than
\(1/\pi\) on \([0,7\pi/10]\), as checked at \(11/5\);
that polynomial decreases over this interval. The rotated
right argument lies in \((-7\pi/10,\pi/2)\) by (2), and
the rotated left arguments also lie in \((-\pi,\pi)\).

Thus the full correction to \(m\theta\) has interval
width below \((\alpha_L+\alpha_U+\beta)\theta<\theta\)
at the root radius. Subtracting the consecutive phase
equations excludes all shared arguments in \(0<\theta\le\pi\).
Positivity and conjugation give the rest. This proves the
claimed improvement for \(S_{\rm sh}\le83/80\); the existing
shared reserve theorem covers the larger values.

## Exact scalar replay

~~~python
from fractions import Fraction as F
from math import factorial, prod

def rising(a, k):
    return prod(a+j for j in range(k))

def lo_exp(x):
    return sum(x**j/factorial(j) for j in range(21))

def hi_exp(x):
    return lo_exp(x)+x**21/factorial(21)/(1-x/22)

def f_upper(a, m, j):
    den = 10**8
    target = rising(a,j)**m * den**m
    weight = rising(a,m)**j
    lo,hi = -1,den
    while hi-lo > 1:
        mid = (hi+lo)//2
        if mid**m*weight >= target:
            hi = mid
        else:
            lo = mid
    assert hi**m*weight >= target
    return F(hi,den)

def cap(a,s,m,k,j):
    loss = s*k*(1-F(k,m))/2-F(7*k,10*m)
    assert loss > 0
    return f_upper(a,m,j)/lo_exp(loss)

def check(left,right,angles,beta,t0=F(0),t1=F(0)):
    A = [sum(head)+t0 for head in (left,right)]
    B = [sum(k*v for k,v in enumerate(head,1))+t1 for head in (left,right)]
    assert A[0] < F(1,2) and right[0] < F(1,2)
    assert (1+A[0])/(1-right[0]) < lo_exp(F(7,10))
    assert (1+A[1])/(1-A[0]) < lo_exp(F(7,10))
    for a,mass,slope in zip(angles,A,B):
        assert slope-a*mass < a
    D = F(2,9)/(beta*(1-(F(22,7)*beta)**2/6))
    assert left[0]+D*(B[0]-left[0]+beta*(A[0]-left[0])) < 1

g13 = f_upper(F(1),13,12)*hi_exp(F(7,130))
g15 = f_upper(F(1),15,14)*hi_exp(F(7,150))
assert F(26,12)*g13 < 1
assert F(39,20)*g13*hi_exp(F(83,960)) < 1
assert F(450,221)*g15*hi_exp(F(83,960)) < 1
assert max(F(3,2)*F(13,11),F(6,5)*F(13,8),
           F(3,2)*F(14,12),F(6,5)*F(14,9),F(7,12)*F(14,8)) == F(39,20)
assert F(3,2)*F(15,13) < F(450,221)
x = F(11,5)
assert 1-x*x/6+x**4/120-x**6/5040 > F(8,25)
assert F(25,8)/3-1-F(8,243) > 0
for a,s,angles,beta in [
    (F(1,2),F(19,50),(F(262,1000),F(671,1000)),F(54,1000)),
    (F(1),F(79,100),(F(312,1000),F(613,1000)),F(66,1000)),
]:
    assert sum(angles)+beta < 1
    assert F(7,13)/lo_exp(s/4) < F(1,2)
    assert cap(a,s,20,10,10) < F(1,2)**10
    for m in range(13,20):
        heads = [
            [cap(a,s,m,k,m-k if side else k)/(2 if 2*k == m else 1)
             for k in range(1,m//2+1)] for side in (0,1)
        ]
        check(*heads,angles,beta)
    q = F(1,2)
    t0 = q**10/(1-q)
    t1 = t0*(10+q/(1-q))
    heads = [
        [cap(a,s,20,k,20-k if side else k) for k in range(1,10)]
        for side in (0,1)
    ]
    check(*heads,angles,beta,t0,t1)
print("All one-sided private totals one/two phase checks passed.")
~~~
