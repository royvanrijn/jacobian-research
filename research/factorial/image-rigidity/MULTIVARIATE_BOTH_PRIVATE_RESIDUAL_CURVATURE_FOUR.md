# Both-private factorial interlacing at residual curvature four

Consider ordinary two-monomial supports with a private variable
in each monomial. The remaining exponent rows may have either
orientation. Define
\[
                    S_*=\sum_{\max(a_i,b_i)>0}
                         \frac{(b_i-a_i)^2}{\max(a_i,b_i)}-2 .
\]

**Theorem.** If \(S_*\ge4\), every consecutive pair of factorial
moment polynomials has strictly interlacing root arguments,
listed from the positive real axis. In particular the ordinary
profile satisfies SFC(2) in every degree over every
characteristic-zero field. This is total curvature at least six.
When the only private factors are one unit on each side, it is
shared curvature at least four.

The argument strengthens the
[residual-curvature-five theorem](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_CURVATURE_INTERLACING.md).
Keeping a few terms of the exact geometric resolvent makes
its first-order contribution cancel in the phase comparison.
The remaining terms are controlled by a narrower annular bound
and geometric tails.

## Endpoint normalization and the annulus

As in the
[both-private angle-box theorem](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_ANGLE_BOXES.md),
reserve a factorial on each side. For a private exponent \(d\)
the positive residual factor comes from
\((dj)!/j!=d^{dj}\prod_{h=1}^{d-1}(h/d)_j\).
The [residual curvature bound](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_CURVATURE_INTERLACING.md#a-common-radial-bound)
then gives the following endpoint-normalized polynomial and
coefficient estimate:
\[
                  P_m(z)=\sum_{k=0}^m c_{m,k}z^k,\qquad
                  c_{m,0}=c_{m,m}=1,\qquad
             0<c_{m,k}\le e^{-2k(m-k)/m}.                      \tag{1}
\]
Its coefficients are the geometric moments of a positive law
under the left endpoint normalization and likewise under
reversal.

For \(m\ge4\) we show
\[
                         e^{-2/(3m)}<|z|<e^{2/(3m)}
                                                               \tag{2}
\]
at every root. On the outer circle pair the interior indices
\(k,m-k\). Their combined radius multiplier is at most
\(1+e^{2/3}<3\), and the even middle term takes at most half
that multiplier. At \(m=4,5,6\), the corresponding sums of
the paired envelopes, with the indicated half weight, are
at most
\[
                 e^{-3/2}+\tfrac12 e^{-2},\qquad
                 e^{-8/5}+e^{-12/5},\qquad
                 e^{-5/3}+e^{-8/3}+\tfrac12 e^{-3},
\]
respectively. For \(m\ge7\), bound the first three coefficients
using \(m\ge7\), and the rest using \(k\le m/2\). The sum is
at most
\[
             e^{-12/7}+e^{-20/7}+e^{-24/7}
                                      +\frac{e^{-4}}{1-e^{-1}}.
\]
All four bounds are below \(3/10\). Thus the interior terms
have modulus below \(9/10<e^{2/3}-1\).
Rouché with \(1+z^m\), and then the same argument after
reversal, proves (2).

## The geometric expansion at either endpoint

At an upper-half-plane root \(z=r e^{i\theta}\) use the
[two endpoint resolvents](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_ANGLE_BOXES.md#the-geometric-resolvent-without-reflection-symmetry).
For either one let \(s\) be the endpoint ratio including
the appropriate radius \(r\) or \(1/r\), and set
\[
        c=\cos\theta,\quad q=\sin\theta,\quad
        D_\theta(s)=1-2cs+s^2,\qquad
        A=\mathbb E[1/D_\theta(s)],\quad
        B=\mathbb E[s/D_\theta(s)] .
\]
By (1) and (2), the endpoint moments \(u_j=\mathbb E s^j\)
satisfy
\[
                 u_j\le a_m(j):=\exp\left(-2j+
                                         \frac{2j^2+2j/3}{m}\right).
                                                               \tag{3}
\]
Here one endpoint uses \(r^j c_{m,j}\), and the other
\(r^{-j}c_{m,m-j}\).
We need the range
\[
                         \frac{\pi}{m-1}\le\theta\le\frac\pi2 .
                                                               \tag{4}
\]
The inequality \(1/D_\theta(s)\ge2-D_\theta(s)\) gives
\(A\ge1+2cu_1-u_2\). For use at small \(m\), the two
[square bounds on the endpoint](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_CURVATURE_INTERLACING.md#controlling-both-endpoint-phases)
give \(B\le u_1+u_2/[2(1-c)]\). Thus
\[
 (1+2c)B-A+1\le
                  u_1+\left(1+\frac{1+2c}{2(1-c)}\right)u_2 .
                                                               \tag{5}
\]

For any integer \(3\le L\le m\), expand
\((1-se^{i\theta})^{-1}\) through exponent \(L-1\) with its
exact geometric remainder. Taking the imaginary part and
using \(|1-se^{i\theta}|\ge q\) gives
\[
         B\le u_1+2c u_2+
                   \sum_{j=3}^{L-1} j u_j+\frac{u_L}{q^2}.
\]
Here \(\sin(j\theta)/q\le j\) for the remaining finite
terms. Therefore
\[
 (1+2c)B-A+1
 \le u_1+\bigl(1+2c(1+2c)\bigr)u_2
                   +\sum_{j=3}^{L-1}3j u_j+
                                         \frac{1+2c}{q^2}u_L .
                                                               \tag{6}
\]
For the longest ranges use the elementary
[resolvent inequality](MULTIVARIATE_PRIVATE_UNIT_MIRRORED_RESOLVENT_INTERLACING.md#bounding-the-resolvent-with-endpoint-moments)
\((1+2c)/q^2\le3/\theta^2\le(m-1)^2/3\).
It holds on (4).

## Scalar phase bounds

For \(4\le m\le9\), substitute (3) and the rational
cosine caps into (5) or (6), taking \(L=3\) in the
last two rows:
\[
\begin{array}{c|c|c|c}
m & \text{upper bound on }c & \text{inequality}
          & \text{upper bound on }(1+2c)B-A+1\\ \hline
4 &1/2&(5)&21/25\\
5 &71/100&(5)&17/20\\
6 &81/100&(5)&9/10\\
7 &13/15&(5)&49/50\\
8 &901/1000&(6)&199/200\\
9 &37/40&(6)&24/25
\end{array}
\]
The trigonometric bounds follow from the elementary
values of \(\cos(\pi/3),\cos(\pi/4),\cos(\pi/5),
\cos(\pi/6),\cos(\pi/8)\) and, for \(m=8\),
\(1-x^2/2+x^4/24<901/1000\) at \(x=333/742<\pi/7\).
All exponential comparisons, including the annulus
bounds, follow from the positive Taylor lower sum
for \(e^x\) through degree twelve.

For \(10\le m\le11\), take \(L=4\) in (6) and use
the bound \((m-1)^2/3\) on the tail factor. The first
three terms decrease with \(m\).
The tail bound \((m-1)^2 a_m(4)/3\) also decreases,
since \(6m^2<104(m-1)\) in this range.
At ten the sum
\[
                       a_{10}(1)+7a_{10}(2)+9a_{10}(3)
                                                 +27a_{10}(4)
\]
is below \(24/25\).

For \(m\ge12\), let \(L=\lfloor m/2\rfloor\) in (6).
The first five bounds
\(u_1,7u_2,9u_3,12u_4,15u_5\) are respectively below
\(17/100,7/25,3/25,3/40,3/50\), using (3) at twelve.
For \(6\le j\le L\), (3) also gives \(u_j\le e^{1/3-j}\).
The resulting sum of terms with \(6\le j<L\) is bounded
by \(3e^{1/3}\sum_{j\ge6}j e^{-j}<11/100\).
The last term is bounded by
\((4/3)L^2 e^{1/3-L}\le48 e^{-17/3}<17/100\).
This sequence decreases for \(L\ge6\).
The seven rational bounds sum to \(197/200<1\).
Thus \((1+2c)B<A\) throughout (4), for both endpoints.

## Consecutive separation

The inequality \((1+2c)B<A\) gives
\(\arg(1-(B/A)e^{-i\theta})<\theta/2\), by comparing
its tangent to \(\tan(\theta/2)=q/(1+c)\).
At the upper-order root with label \(k\), the exact
phase equation is
\[
                   m\theta=(2k-1)\pi
                              +\phi_{\rm left}+\phi_{\rm right},
\]
with the two endpoint phases just bounded. If that root
were in (4) and had angle at least the lower left
endpoint \((2k-1)\pi/(m-1)\), the strict bound on the
phase sum would give the opposite inequality.
If its angle exceeds \(\pi/2\), its own angle box
already places it below the lower left endpoint.
The other side of the interlacing follows from the
angle boxes directly. The pairs with upper order two
and three also follow from their boxes and the odd
negative real root. This proves strict consecutive
interlacing in every degree. The ordinary rational
Bézout identity transfers coprimality to every
characteristic-zero field.

## Reproducing the scalar bounds

The following rational comparisons certify the
scalar tables and tails. The Machin formula, which
follows by the tangent addition formula, gives
\(\pi>333/106\) from a four-term alternating lower
bound on \(\arctan(1/5)\) and the upper bound
\(\arctan(1/239)<1/239\).

~~~python
from fractions import Fraction as F
from math import factorial

def lo(x):
    return sum(F(x)**j/factorial(j) for j in range(13))

def hi(x):
    x = F(x)
    return lo(x) + x**13/factorial(13)/(1-x/14)

def neg(x):
    return 1/lo(x)

def a(m, j):
    return neg(2*j - F(2*j*j + F(2*j, 3), m))

assert hi(F(2,3)) < 2 and lo(F(2,3)) > F(19,10)
annuli = [
    neg(F(3,2)) + neg(2)/2,
    neg(F(8,5)) + neg(F(12,5)),
    neg(F(5,3)) + neg(F(8,3)) + neg(3)/2,
    neg(F(12,7)) + neg(F(20,7)) + neg(F(24,7)) + neg(4)/(1-neg(1))
]
assert max(annuli) < F(3,10)
x = F(1,5)
assert 4*(4*(x-x**3/3+x**5/5-x**7/7)-F(1,239)) > F(333,106)
x = F(333,742)
assert 1-x*x/2+x**4/24 < F(901,1000)
assert F(71,100)**2 > F(1,2)
assert (4*F(81,100)-1)**2 > 5
assert F(13,15)**2 > F(3,4)
assert F(17,12)**2 > 2 and F(37,40)**2 > (2+F(17,12))/4

small = [
    a(m,1)+(1+(1+2*c)/(2*(1-c)))*a(m,2)
    for m,c in [(4,F(1,2)),(5,F(71,100)),(6,F(81,100)),(7,F(13,15))]
]
small += [
    a(m,1)+(1+2*c*(1+2*c))*a(m,2)+(1+2*c)/(1-c*c)*a(m,3)
    for m,c in [(8,F(901,1000)),(9,F(37,40))]
]
caps = [F(21,25), F(17,20), F(9,10), F(49,50), F(199,200), F(24,25)]
assert all(x < y for x,y in zip(small,caps))
assert a(10,1)+7*a(10,2)+9*a(10,3)+27*a(10,4) < F(24,25)

first = [a(12,1), 7*a(12,2), 9*a(12,3), 12*a(12,4), 15*a(12,5)]
caps = [F(17,100), F(7,25), F(3,25), F(3,40), F(3,50)]
assert all(x < y for x,y in zip(first,caps))
q = neg(1)
assert 3*hi(F(1,3))*q**6*(6-5*q)/(1-q)**2 < F(11,100)
assert 48*neg(F(17,3)) < F(17,100)
~~~
