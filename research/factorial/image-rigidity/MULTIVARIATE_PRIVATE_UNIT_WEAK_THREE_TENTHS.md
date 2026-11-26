# All degree physical coprimality for private unit curvature at most three tenths

Take an ordinary two-monomial profile with exactly one private
unit row in each orientation and arbitrary shared rows. For the
shared rows write \(d_i=b_i-a_i\), \(\ell_i=\min(a_i,b_i)\),
and set
\[
           V=\sum_i d_i^2/\ell_i,\qquad \Delta=\sum_i d_i .
\]

**Theorem.** If \(V\le3/10\), the ordinary consecutive physical
factorial moments are coprime in every degree, over every
characteristic-zero field. The curvature bound permits arbitrary
shared orientations and signed total degree difference.
No eventual cutoff or fixed-terminal verifier is used.

This follows from the exact adjacent positive-moment link and
a comparison of normalized root radii. Weak curvature bounds
each row ratio \(|d_i|/\ell_i\) by \(V\). This improves the
exponential moment estimate for the linked logarithmic noise.
A physical common root would force the integer \(\Delta\) to
vanish. The balanced Euclidean theorem then rules it out.
[Exact linked noise and resolvent estimates](MULTIVARIATE_PRIVATE_UNIT_HIGHER_LABEL_LINKED_RESOLVENT.md)
[Balanced Euclidean theorem](MULTIVARIATE_BOTH_PRIVATE_EUCLIDEAN_DESCENT.md)

## Modulus estimate with the weak row cap

The geometric case \(V=0\) is immediate, so assume \(V>0\).
Write \(N=n+1\). If the physical pair has a common nonreal root,
after conjugation its angle \(\theta\) satisfies
\(\pi/n<\theta<\pi/2\), by the both-private angle boxes.
Put \(b=n\sin\theta\), \(\alpha=1/\sin\theta\), and let
\(x_n,x_N\) denote the log radii of the associated
endpoint-normalized roots. The linked noise is
\[
    Z=\log(G_N/G_n)+\sum_i d_i\log Y_i,\qquad
                       Y_i\sim{\rm Beta}(na_i+1,a_i),
\]
independently, where \(G_m\) is the endpoint scale.
In each orientation \(\log T_n\) has the law of
\(\log T_N+Z\) with independence. Use \(M(t)=\mathbb E e^{tZ}\)
and \(\Lambda=\log M\).

The exact mean bound \(0\le\mathbb E Z\le V/(4nN)\) remains
valid. As \(|d_i|/a_i\le V\) by integrality of the nonzero
gaps, the proof of the linked trigamma bound now gives
\[
 \Lambda''(t)\le \frac{V}{nN(1-V|t|/n)^2},\qquad
 |\Lambda'(t)|\le \frac{V}{4nN}
           +\frac{V|t|}{nN(1-V|t|/n)},\qquad |t|<n/V.
                                                               \tag{1}
\]
The indicated open domain gives finite exponential moments
and their derivatives for both orientations. Exactly the
same estimates hold on the conjugate angle. The packet bound
at any fixed log radius has relative error at most
\[
       \varepsilon=\frac{\mathbb E Z}{\sin\theta}
              +\frac{M''(\alpha)+M''(-\alpha)}{\sin^2\theta}.
                                                               \tag{2}
\]
For the derivation, \(F(s)=(1-e^{s+i\theta})^{-1}\) satisfies
\(|F'|=\operatorname{Im}F/\sin\theta\),
\(|F''|\le\operatorname{Im}F/[\sin\theta\,\sin(\theta/2)]\)
and \(|(\log\operatorname{Im}F)'|\le1/\sin\theta\).
Apply these bounds to the independent shift \(Z\) in
\(F(x+\log T_N+Z)\) by Taylor's formula.
[Resolvent variation at a horizontal position](MULTIVARIATE_PRIVATE_UNIT_HIGHER_LABEL_LINKED_RESOLVENT.md#resolvent-variation-at-any-horizontal-position)

In the logarithmic modulus equation \(m x_m=K_m(x_m,\theta)\),
the derivative in \(x\) is bounded by \(2/\sin\theta\) for
the difference of the two endpoint packet logarithms.
Their adjacent-order difference at fixed \(x\) is bounded by
\(2[-\log(1-\varepsilon)]\) when \(\varepsilon<1\) bounds
(2) in both orientations. Subtracting the two modulus
equations yields
\[
  N|x_N-x_n|
    \le \frac{b(1+1/n)}{b(1+1/n)-2}
       \left(|x_n|+2[-\log(1-\varepsilon)]\right)
                                                               \tag{3}
\]
when \(b(1+1/n)>2\).
At a physical common root, \(x_N-x_n=\log(G_N/G_n)\).
The endpoint drift estimate differs from
\(\Delta\log(N/n)\) by at most \(V/(nN)\), since the private
unit contributions cancel. Thus if \(\Delta\ne0\), the right
side of (3) plus \(V/n\) is a strict upper bound for
\(|\Delta|\).
[Endpoint drift](MULTIVARIATE_ALL_DEGREE_SIGNED_IMBALANCE.md#an-explicit-linked-radial-condition)

For all \(m\ge2\), Eneström–Kakeya gives the useful individual
radius estimate
\[
                     |x_m|\le V(m-1)/(2m).                    \tag{4}
\]
In fact the endpoint-normalized log coefficient chord satisfies
\(0\le f''\le V/m\), so both end neighbors have coefficient
at least \(\exp[-V(m-1)/(2m)]\) by the quadratic Green bound.
The intervening successive ratios are between the end ratios
by log convexity. This proves (4).

## Uniform bounds for the integer drift

For \(n\ge7\) use \(b\ge b_0=303/100\).
Indeed \(n\sin(\pi/n)\ge \pi-\pi^3/(6n^2)>303/100\),
using \(157/50<\pi<22/7\).
On \(|t|\le\alpha\), (1) gives
\[
 \Lambda''(t)\le (5/4)V/(nN),\qquad
 |\Lambda'(t)|\le (41/100)V/N,\qquad \Lambda(t)\le1/40.
\]
For the last use the mean bound and Taylor's formula.
Consequently \(M\le40/39\) and
\[
 M''(\alpha)+M''(-\alpha)
    \le (80/39)\big[(5/4)+(3/10)(41/100)^2\big]V/(nN)
      <(27/10)V/(nN).
\]
Using \(x=1/n\), bound \(\varepsilon\) in (2) by
\[
 E(x)=\frac{V_0}{1+x}
        \left(\frac{27}{10b_0^2}+\frac{x}{4b_0}\right),
                                      \qquad V_0=3/10 .
\]
This bound decreases with \(x\) and is less than one.
The resulting majorant for the integer drift from (3) is
\[
 \frac{b_0(1+x)}{b_0(1+x)-2}
       \left(\frac{V_0}{2}+2[-\log(1-E(x))]\right)+V_0 x.
                                                               \tag{5}
\]
Both factors of the product are positive, decreasing and
convex for \(0\le x\le1/7\).
For the second this follows from \(E\) decreasing and
convex and from the increasing convex function
\(-\log(1-u)\). The product is convex because the two
first derivatives have the same sign. Adding the affine
term preserves convexity, so the maximum occurs at an endpoint.
At either endpoint use
\(-\log(1-u)\le u+u^2/[2(1-u)]\).
Both endpoint bounds are strictly less than one.

For \(4\le n\le6\), use \(b\ge14/5\) because
\(n\sin(\pi/n)\ge2\sqrt2\), and
\(N/n\ge7/6\). The contraction factor in (3) is at most
\(49/19\). On \(|t|\le\alpha\),
\[
 \Lambda''(t)\le(13/10)V/(nN),\qquad
 |\Lambda'(t)|\le(47/100)V/N,\qquad \Lambda(t)\le1/35 .
\]
This gives
\(M''(\alpha)+M''(-\alpha)\le (29/10)V/(nN)\)
and \(\varepsilon\le101/1000\), using \(N\ge5\) and
\(n/N\le6/7\) in (2). By (4), \(|x_n|\le5V/12\).
The drift majorant is bounded by
\[
 \frac{49}{19}\left(\frac{5V_0}{12}
                  +\frac{202}{899}\right)+\frac{V_0}{4}<1.
\]
The rational term bounds \(2[-\log(1-\varepsilon)]\)
by \(2\varepsilon/(1-\varepsilon)\).

For \(n=2,3\), combine (4) directly with the physical endpoint
drift. One obtains
\[
 |\Delta|<V\left[\frac{N(n-1)}{2n}+\frac{N-1}{2}+\frac1n\right]
                                \le 19V/6\le19/20
\]
whenever \(\Delta\ne0\). The first and second moments have
disjoint roots because the first has a real root and the
second is strictly positive on the real axis. In every
remaining case the strict integer bounds give \(\Delta=0\)
and the balanced Euclidean theorem applies.

## Exact scalar checks

The following rational calculations verify the constants
used in the two log-MGF estimates and in the integer bounds.

~~~python
from fractions import Fraction as F

v = F(3, 10)
b = F(303, 100)
assert F(157, 50) - F(22, 7)**3/(6*7**2) > b
assert (1-v/b)**-2 <= F(5, 4)
assert F(1, 28) + 1/(b-v) <= F(41, 100)
assert v/(4*b*8) + F(5, 8)*v/b**2 <= F(1, 40)
assert F(80, 39)*(F(5, 4)+v*F(41,100)**2) < F(27, 10)
assert F(27,10)/b**2 > 1/(4*b)
for x in [F(0), F(1, 7)]:
    e = v*(F(27,10)/b**2+x/(4*b))/(1+x)
    assert e < 1
    assert b*(1+x)/(b*(1+x)-2)*(v/2+2*e+e**2/(1-e))+v*x < 1

b = F(14, 5)
assert b*b < 8
assert (1-v/b)**-2 <= F(13, 10)
assert F(1,16)+1/(b-v) <= F(47,100)
assert v/(4*b*5)+F(13,20)*v*F(6,7)/b**2 <= F(1,35)
assert F(35,17)*(F(13,10)+v*F(47,100)**2) < F(29,10)
assert v/(4*b*5)+F(29,10)*v*F(6,7)/b**2 <= F(101,1000)
assert F(49,19)*(5*v/12+F(202,899))+v/4 < 1
assert 19*v/6 == F(19,20)
~~~
