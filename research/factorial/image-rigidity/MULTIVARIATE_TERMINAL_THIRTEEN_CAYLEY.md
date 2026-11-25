# Terminal-thirteen Cayley separation for light shared rows

Let two distinct ordinary exponent columns have only shared unequal rows.
Write \(u\) for the order-twelve reflected outer mean after positive
endpoint normalization and removal of the binomial coefficients. Then
\(0<u<1\). The following conditions on \(u\) and the largest order-twelve
Gauss reciprocal \(x=d/(12\ell+h)\) suffice for the order-twelve and
order-thirteen factorial moment polynomials to have disjoint root
arguments.

| \(u\) | Reciprocal cap | Projection cells |
|---|---:|---:|
| \([19/20,1)\) | Unrestricted | 196,768 |
| \([23/25,19/20]\) | \(1/90\) | 143,877 |
| \([22/25,23/25]\) | \(1/70\) | 195,968 |
| \([171/200,22/25]\) | \(1/60\) | 104,274 |
| \([19/25,171/200]\) | \(1/55\) | 275,217 |

The assertion holds under independent positive rescalings. The original
polynomials are coprime over every characteristic-zero field by rational
Bézout.

## Scales and grouped coefficient bounds

Write \(w\) for the order-thirteen outer mean and put
\(p=1-u,\ q=1-w,\ r=q/p,\ R=144/143\). If \(c\) bounds the
lower-order reciprocals, the logarithmic ratio lies in
\[
 R\le \frac{\log w}{\log u}
 \le L=1-\frac{\log(1-c/13)}{\log(1+11c)}.
\]
The lower bound is the
[matched outer-ratio inequality](MULTIVARIATE_MATCHED_OUTER_RATIO_MONOTONICITY.md).
For the upper bound use \(x'=12x/(13-(h/d)x)\) and
\(1+12x'\le(1+11x)/(1-x/13)\). The ratio
\(-\log(1-x/13)/\log(1+11x)\) increases with \(x\), since its
numerator is convex and its denominator is concave, with both zero at
\(x=0\). The limiting value of \(L\) at \(c=0\) is \(R\).
The reciprocal bound can always include \(c\le(u^{-2}-1)/11\)
by the product formula for \(u^{-2}\).

For a subband \(p\in[p_-,p_+]\) take an outward rational upper
bound on \(L\) and use the ratio enclosure
\[
 \frac q p \in
 \left[
     R\left(1-\frac{(R-1)p_+}{2(1-p_+)}\right),
     L\left(1-\frac{(L-1)p_-}{2}\right)
 \right].
\]
For \(1<a<2\) and \(1-p\le v\le1\), the derivative of \(v^{a-1}\)
lies in \([a-1,(a-1)/(1-p)]\). Bound \(v^{a-1}\) by integrating
from \(v\) to \(1\), and then integrate those bounds to get
\[
 a\left(1-\frac{(a-1)p}{2(1-p)}\right)
 \le \frac{1-(1-p)^a}{p}
 \le a\left(1-\frac{(a-1)p}{2}\right).
\]
The checker verifies \(R\le L<2\). The auxiliary reference ratio
\(r_*(p)=R(1-(R-1)p/2)\) belongs to the same enclosure because
\(a(1-(a-1)p_-/2)\) increases for \(1\le a\le2\) and
\(p_-\le .24\).

The first ray uses subbands of width \(1/1000\). Its reciprocal
cap also uses the integer row bound. For a unit-gap row, keep the
largest \(1/(12\ell+1)\) with \(\ell\ge1\) at most \(c\).
For a gap \(d\ge2\) with largest reciprocal \(x_1\), each of the
other row reciprocals is at least \(x_1/(1+x_1)\). In particular
\(u^{-2}\ge(1+11x_1/(1+x_1))^2\), so \(x_1\le(1-u)/(12u-1)\).
The maximum of the two cases, each capped by \(c\), gives the
refined bound. These bounds are evaluated at the smallest \(u\) of
the subband.

At order \(j\), let \(h=p\) or \(q\) and \(\lambda=\sqrt h\).
Use the [weight-nine grouped Cayley model](MULTIVARIATE_TERMINAL_TWELVE_CAYLEY.md#weighted-cayley-model)
at \(j=12,13\). An outward upper bound \(q_0\) on the two deficits
gives \(C=(1-q_0/2)/(1-q_0)\). The normalized tanh factors
\(a_i=T_i/h\) satisfy \(\sum_i a_i\le C\) and \(a_i\le A\), with
\(A\) sharpened using the reciprocal cap and the lower deficit
endpoint. The checker keeps all log-moment products of total weight
at most nine and uses the absolute weight-ten remainder bound.
All initial cancellations in a grouped direction are checked before
division by \(\lambda^j\).

Where needed, the projection subdivides the signed weight-two
moment \(P=\sum_i\epsilon_i a_i^2/2\) into eight subintervals of
\([-CA/2,CA/2]\). The
[conditional moment bounds](MULTIVARIATE_TERMINAL_TWELVE_CAYLEY_LOWER.md#condition-the-logarithmic-moments-on-the-leading-skew)
then sharpen each retained odd moment from below and each retained
higher even moment on both sides. Subdivision of \(P\) can be
performed within a failed projection cell, since all eight
subintervals are checked there. The moments of the actual
logarithmic corrections can all be multiplied by
\(s\in[0,1]\) at fixed \(p,r\). The conditional bounds remain
valid at \(P_s=sP\), by \(s\ge s^{d-1}\) for each odd weight
\(d\ge3\) and \(s(CA\pm2P)\le CA\pm2sP\) for the even
weights. The weight-ten error bounds hold for the entire path.

## Root counts and replay

The scaled Cayley disk has radius \(13/5\). The checker proves
leading coefficient dominance and excludes zeros on eleven
rays in the physical left half-plane. The rays have arguments
\(\pi-2\arctan(\sqrt p\,b)\), with distinct positive piecewise
affine levels \(b\). The inverse ray equation is the
[grouped Cayley height relation](MULTIVARIATE_TERMINAL_ELEVEN_MIDDLE_CAYLEY.md#ray-certification-and-continuation).
On the band where the disk can extend outside the physical left
half-plane the checker clips to \(0\le\Re Z\le1/\lambda\),
as required by a point on such a ray. Its reflection covers
the other half of the inverse image. The coefficient homotopy
has positive real coefficients in the original variable.
Its Cayley roots therefore avoid the points mapping to zero
or infinity, and conjugation gives the corresponding lower rays.

First take \(s\) from one to zero at the actual \(p,r\). Then
take \(r\) to \(r_*(p)\) at fixed \(p\) and follow the reference
curve to the smallest \(p\) of the band. Dominance on the disk
boundary and strict projections preserve the sector counts.
The anchor checker evaluates signs between successive
imaginary-axis heights, dividing by the structural zero factor
at odd order. At \(p=0\) it uses the limiting scaled Cayley
polynomials and horizontal inverse ray images. The alternating
signs account for all six positive heights of each reference,
in orders thirteen and twelve starting from the largest height.
Order thirteen has the additional structural root at zero in
the scaled coordinate.

At positive deficit the sectors adjacent to the positive real ray
contain only order-thirteen roots. The sectors adjacent to the
negative real ray, joined across that ray, contain one order-thirteen
root and two order-twelve roots. Conjugation fixes the odd-order
root on the negative real ray. The Gamma integral of the even
power of a nonzero binomial makes the order-twelve moment strictly
positive on the real axis, so its two roots in this joined sector
are nonreal. Thus the sector counts imply disjoint root arguments.

Replay the [exact checker](verify_multivariate_terminal_thirteen_cayley.py)
from the repository root. This verifies 916,104 exact projection
cells across 1,056 order, ray and subband checks.

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_thirteen_cayley.py
python3 research/factorial/image-rigidity/verify_multivariate_terminal_thirteen_cayley.py flat 0 -1 0 15
~~~
