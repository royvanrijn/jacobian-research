# Shared terminal-nine angular separation in verified upper bands

Let two distinct ordinary exponent columns have only shared unequal
rows. Let \(u\) be the order-eight outer geometric mean after positive
endpoint normalization and removal of the binomial factors. If
\[
 u\in[.880,.895]\ \cup[.900,.905]\ \cup[.910,.915],
\]
the order-eight and order-nine polynomials have disjoint root
arguments under independent positive rescalings. The decimals
are exact rationals.
This holds in arbitrary dimension and for arbitrary positive
exponents. Rational Bézout gives coprimality in characteristic zero.

The [marked-row theorem](MULTIVARIATE_TERMINAL_NINE_UPPER_MARKED_ROWS.md)
covers profiles containing either \((3,4)\) or \((4,5)\), including
opposite orientations. The rest of this proof excludes those rows.

## Bounds on ordinary reciprocals

For a row with smaller exponent \(\ell\) and gap \(d\ge1\), use
\(x_h=d/(8\ell+h)\), \(1\le h\le d\). The squared outer mean
of this row bounds \(u^2\) from above:
\[
               u^2\le\prod_{h=1}^d \frac{8\ell+h}{8\ell+7d+h}.
\]
Let \(u_-\) be the lower endpoint of a tested interval and let
\(\ell_-\) be the least positive integer satisfying
\((8\ell_-+1)/(8\ell_-+8)\ge u_-^2\).
If \(d=1\) this gives \(x_h\le1/(8\ell_-+1)\). At \(u\ge22/25\)
the rows with \(\ell=1,2\) are impossible for \(d=1\). Since
\(\ell=3,4\) are already covered, every remaining unit row has
\(x_h\le1/41\).

For \(d\ge2\) put \(y=d/(8\ell)\). The displayed product gives
\(u\le ((1+y)/(1+8y))^{d/2}\le(1+y)/(1+8y)\), so
\[
             x_h\le \frac{1-u_-}{8u_--1}\le\frac3{151}<\frac1{41}.
\]
Also \(u^{-2}\) is the full product of \(1+7x\) over all the
shared factors, so \(x_h\le(u_-^{-2}-1)/7\). The checker uses
the minimum of this last bound, \(1/41\), and the maximum of
the separate unit and nonunit bounds. Write it as \(c\).
The matched order-nine reciprocals are at most \(8c/(9-c)\).

## Interval coefficient model and ray exclusion

Apply the [grouped logarithmic coefficient
model](MULTIVARIATE_TERMINAL_NINE_PRIVATE_LINEAR.md#grouping-the-shared-logarithmic-corrections)
to the fully shared profile. At the midpoint \(u_0\), use the
symmetric reference multipliers
\(u_0^{e_{8,k}}\) and \(u_0^{(64/63)e_{9,k}}\), where
\(e_{j,k}=k(j-k)/(j-1)\). The model preserves each of the six
linked linear directions and every quadratic product separately.
It bounds the omitted logarithmic terms and the exponential
Taylor remainder per coefficient.

Put \(p=1-u=\lambda^2\), \(t=(\lambda X-1)/(\lambda X+1)\).
The transformed polynomial is the sum of the normalized
physical coefficients against \((\lambda X-1)^k(\lambda X+1)^{j-k}\).
For each direction the exact Cayley matrix first combines reflected
weights, then multiplies by the common reference coefficients.
Each resulting coefficient of degree \(h\) is scaled by \(\lambda^h\).
The error array is transformed as independent coefficient intervals.
Near the upper end the checker brackets the reference powers on
the dyadic grid \(2^{-52}\) by comparing exact integer powers.

The [exact checker](verify_multivariate_terminal_nine_upper_bridge.py)
proves leading dominance on \(|X|=2\) and excludes seven rays
throughout the disk. Their positive levels \(b\), in decreasing
order, are piecewise affine functions of \(u\) with rational knots.
Their physical arguments are \(\pi-2\arctan(\lambda b)\). At
\(X=x+iy\) the ray equation is
\[
                    (y-b)(1+pby)=-pbx^2,\qquad 0<y\le b .
\]
Start with \(y\in[0,b]\) and intersect three times with the interval
for \(b-pbx^2/(1+pby)\). It suffices to check \(0\le x\le2\).
Reflecting \(X\) to \(-\overline X\) reverses the normalized
real coefficients under the Cayley map. That reverses the skew
signs and preserves all bounds in the model.

For each tested rectangle, first enclose the real and imaginary
parts of the base polynomial and of each transformed direction
by dyadic intervals. Test a rational linear functional on these
complex rectangles before combining their scalar intervals.
The candidate tests use the real and imaginary coordinates and
normals to the centers of the direction rectangles. The upper band
also uses normals to the monomial rectangles to account for the
coefficientwise remainder. Additional tests use the projections
\(\Im(P(X)\overline X^{\,h})\), \(0\le h\le j\).
A rectangle is accepted only when one of the exact interval
projections has a strict sign. All other rectangles subdivide
in \(u\) or \(x\). The root disk, sign exclusions, knot ordering,
and anchor are checked independently.

## Sector counts

At each fixed \(u\) continue along the straight coefficient segment
from the actual arrays to the symmetric reference with ratio
\(64/63\). The model is convex as an affine enclosure using the
independent scalar, quadratic-product and remainder intervals.
It contains both arrays. Continue the two reference arrays within each connected certified
band to its left endpoint. The checker evaluates the same Gaussian
anchor signs at all of those endpoints. The leading estimate keeps the scaled
roots in \(|X|<2\), which maps into the physical left half-plane
since \(4(1-u)<1\). The seven rays are disjoint inside this region
and continue without a root crossing.

At the anchor set \(Z=t+t^{-1}\). Divide the octic by \(t^4\),
and divide the nonic by \((1+t)t^4\). In each case the quotient
is a polynomial of degree four in \(Z\), evaluated by
\[
 t^{h+1}+t^{-(h+1)}=Z(t^h+t^{-h})-(t^{h-1}+t^{-(h-1)}).
\]
The rational abscissae are \(2,-2\) and the seven values
\(2(pb^2-1)/(1+pb^2)\). Alternating sign changes place the eight
upper anchor roots in successive sectors of orders
\(9,8,9,8,9,8,9,8\). The extra nonic root is \(-1\).
These sector counts persist, with the last upper sector and its
conjugate joined across the negative real axis. In that joined
sector the sole nonic root stays real by conjugacy. The actual
octic has no real root, since the eighth power of the difference
of two distinct real monomials has strictly positive factorial
moment. This completes the angular separation.

The replay certifies 741,493 projection cells for the displayed bands.
Run it from the repository root.

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_nine_upper_bridge.py
~~~
