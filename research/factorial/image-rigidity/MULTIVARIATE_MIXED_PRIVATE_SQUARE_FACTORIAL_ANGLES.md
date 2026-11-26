# All-order angles with private exponent totals at least one and three

An ordinary two-monomial profile has **disjoint consecutive root
arguments in every degree** if its two private exponent totals are
at least one and three, respectively. Shared rows may have either
orientation and no size restriction applies to the remaining rows.
A common Gaussian coefficient strength \(\tau\ge0\) is allowed.
The proof uses only a one-sided bound on each normalized array.

In the [matched tanh grid](MULTIVARIATE_MATCHED_TANH_INTEGER_MULTIPLICITY.md#recovering-the-gauss-grid)
the private slots are exactly those with
\((n/T'-(n-1)/T)^{-1}=1\). Each private exponent contributes
that many slots in its orientation. The criterion can therefore
be tested directly from the exact matched rows.

## Square factorial bounds for both halves

Orient the side with private degree at least three toward index \(k\).
In degree \(m\) put \(f_{m,k}=k!/(m!)^{k/m}\). After normalizing the
endpoint coefficients to one by positive scalings, the coefficients
satisfy
\[
           0<c_{m,k}\le f_{m,k}^2 .
\]
Indeed a private row of exponent \(d\) contributes at most \(f_{m,k}^d\)
in the chosen orientation. The consecutive ratios of \((dk)!/(k!)^d\)
are \(\prod_{h=1}^d(d-(d-h)/(k+1))\) and hence nondecreasing.
Subtracting the endpoint chord of the log of this sequence proves
the bound. The other orientation similarly contributes at most
\(f_{m,m-k}^d\). One factor of each orientation cancels the binomial
coefficient. Every remaining normalized row factor is at most one
by log convexity, as is the Gaussian factor. This gives the bound.

Take \(R_m=e^{7/(10m)}\). Split
\(P_m(z)=L_m(z)+z^m U_m(1/z)\), giving each half weight \(1/2\)
at the even middle coefficient. The coefficient caps at radius
\(R_m\) for the two halves are
\[
             b_{L,k}=R_m^k f_{m,k}^2,\qquad
             b_{U,k}=R_m^k f_{m,m-k}^2
             \quad (1\le k\le \lfloor m/2\rfloor),
\]
each halved at the even middle index. The second cap dominates
the first. In fact for \(\phi(x)=\log\Gamma(1+x)-x\log(m!)/m\),
the function \(\phi(m-x)-\phi(x)\) is concave on \([0,m/2]\)
because the trigamma decreases. Its endpoint values vanish.

For fixed \(k\) the unhalved caps decrease with \(m\ge k\).
For \(f_{m,k}\) use the increasing geometric mean \((m!)^{1/m}\).
The logarithmic derivative of \(f_{m,m-k}\) at real \(m\ge k\) is
\[
 \frac{k}{m^2}\int_0^m s\psi'(1+s)\,ds
                  -\sum_{h=0}^{k-1}\frac1{m-h}<0,
\]
using \(\psi'(1+s)<1/s\) for \(s>0\) by the decreasing sum
and integral comparison. The factor \(R_m^k\) also decreases.

For \(m\ge14\) use the unhalved bounds at \(m=14\) for the first
three coefficients in either half. For \(k\ge4\), the uniform
cap in either half is \(3^{-k}\). To prove this first take
\(m=2h\), \(h\ge7\). The logarithm of the unhalved right cap is
convex in \(k\) and vanishes at zero, so its ratio to \(k\)
is bounded by its value at the middle divided by \(h\).
There the unhalved cap is \(e^{7/20}/\binom{2h}{h}<3^{-h}\),
as holds at \(h=7\) and the consecutive central-binomial ratio
is at least \(3\) thereafter. An odd degree \(2h+1\) uses the
unhalved right caps from \(2h\) by fixed-index monotonicity.

## Signed sectors at the actual radii

Use \((a_L,a_U)=(71/200,43/100)\) and
\((\beta_L,\beta_U)=(2/25,3/25)\). For each half \(H=L,U\)
the caps just given satisfy
\[
 \sum_k b_{U,k}<\tfrac12,\qquad
 \sum_k(k-a_H)b_{H,k}<a_H,\qquad
 b_{H,1}+\frac{10}{43\beta_H}\sum_{k\ge2}(k+\beta_H)b_{H,k}<1 .
                                                               \tag{1}
\]
At \(m=1\) the halves are one and these sums vanish.
The exact arithmetic below verifies (1) using scalar caps for
\(2\le m\le13\) and the uniform bound above for larger \(m\).

On radius \(R_m\), compare \(L_m(z)\) with \(U_m(1/\bar z)\) on
the same ray. Their coefficient difference in absolute value is
bounded by \(\sum_k b_{U,k}\), because each of the two terms at
index \(k\) lies in \([0,b_{U,k}]\) after applying its radial
factor. Also \(|U_m(1/z)|\ge1-\sum_k b_{U,k}\). Thus the ratio
of the two half moduli is less than \(2<R_m^m\). The halves
are nonzero throughout \(|z|\le R_m\) by the same cap.
Rouché in both orientations confines all roots to
\(R_m^{-1}<|z|<R_m\).

For \(H=L,U\), \(0\le t\le R_m\) and \(0<\theta\le\pi\), the
last two bounds in (1) imply
\[
                   -\beta_H\theta<\arg H(t e^{i\theta})<a_H\theta.
                                                               \tag{2}
\]
The real part of \(H\) is positive. For the upper bound rotate
by \(-a_H\theta\) and use
\(\lvert\sin((k-a_H)\theta)\rvert
 \le [(k-a_H)/a_H]\sin(a_H\theta)\).
This follows by decreasing sinc on \([0,\pi]\) and
\(\lvert\sin x\rvert/x\le1/\pi\) past \(\pi\).
The sum in the second inequality of (1) is exactly the
upper bound on the nonconstant part after division by
\(\sin(a_H\theta)\).

For the lower bound rotate by \(\beta_H\theta\).
The negative part of \(\sin((1+\beta_H)\theta)\) is at
most \(\sin(\beta_H\theta)\) by \(\theta\le\pi\).
For any \(k\ge2\), the ratio of the negative part of
\(\sin((k+\beta_H)\theta)\) to \(\sin(\beta_H\theta)\)
is at most \(10(k+\beta_H)/(43\beta_H)\).
Indeed \(-\sin x\le (2/9)x\) for \(x\ge\pi\).
For \(x\ge9/2\) use \(-\sin x\le1\). For
\(\pi\le x\le9/2\), write \(y=3\pi/2-x\) and
\(-\sin x=\cos y\le1-3y^2/8\) by Taylor's inequality
and \(y^2\le\pi^2/4<3\). Completing the square gives
\[
      \frac{2x}9-1+\frac38 y^2
          =\frac38(y-8/27)^2+\frac{\pi}3-1-\frac8{243}>0.
\]
Also \(\sin(\beta_H\theta)\ge(43/45)\beta_H\theta\)
by \(\beta_H\le3/25\) and the cubic sine bound.
These bounds prove the claimed ratios and the lower
inequality of (2). Both rotated arguments lie in
\((-\pi,\pi)\), so the imaginary signs determine the branches.

At a root on a ray of angle \(0<\theta\le\pi\), the
phase correction to \(m\theta\) is the sum of the two
half arguments evaluated at the root radius and its
reciprocal on that ray. By (2) all these corrections
lie in \((-\theta/5,157\theta/200)\), an interval of
length \(197\theta/200<\theta\). Subtracting the two
consecutive phase equations would leave a value
strictly between zero and \(2\theta\) equal to a
multiple of \(2\pi\). Positive coefficients exclude
angle zero, and reflection covers the lower half-plane.

## Exact scalar replay

~~~python
from fractions import Fraction as F
from math import comb, factorial

def exp_upper(x):
    return sum(x**j/factorial(j) for j in range(13)) + x**13/factorial(13)/(1-x/14)

def root_cap(m, j):
    den = 10**6
    target = factorial(j)**(2*m) * den**m
    weight = factorial(m)**(2*j)
    lo, hi = -1, den
    while hi-lo > 1:
        mid = (lo+hi)//2
        if mid**m * weight >= target:
            hi = mid
        else:
            lo = mid
    assert hi**m * weight >= target
    return F(hi, den)

def cap(m, k, j):
    return root_cap(m, j) * exp_upper(F(7*k, 10*m))

def check(mass, head1, slope, higher, a, beta):
    assert mass < F(1, 2)
    assert slope-a*mass < a
    assert head1+F(10, 43)/beta*(higher+beta*(mass-head1)) < 1

A = (F(71, 200), F(43, 100))
B = (F(2, 25), F(3, 25))
assert sum(A)+sum(B) == F(197, 200)
assert 1+F(7, 10)+F(7, 10)**2/2+F(7, 10)**3/6 > 2
assert F(25, 8)/3-1-F(8, 243) > 0
assert 1-(F(22, 7)*max(B))**2/6 > F(43, 45)

for m in range(2, 14):
    for side, (a, beta) in enumerate(zip(A, B)):
        head = [
            cap(m, k, m-k if side else k)/(2 if 2*k == m else 1)
            for k in range(1, m//2+1)
        ]
        check(sum(head), head[0], sum(k*c for k,c in enumerate(head, 1)),
              sum(k*c for k,c in enumerate(head, 1) if k >= 2), a, beta)

assert exp_upper(F(7, 20))*3**7 < comb(14, 7)
q = F(1, 3)
tail = q**4/(1-q)
tail_slope = tail*(4+q/(1-q))
head = [cap(14, k, 14-k) for k in range(1, 4)]
for a, beta in zip(A, B):
    check(sum(head)+tail, head[0],
          sum(k*c for k,c in enumerate(head, 1))+tail_slope,
          sum(k*c for k,c in enumerate(head, 1) if k >= 2)+tail_slope, a, beta)
print("All asymmetric square factorial half-phase scalar bounds passed.")
~~~
