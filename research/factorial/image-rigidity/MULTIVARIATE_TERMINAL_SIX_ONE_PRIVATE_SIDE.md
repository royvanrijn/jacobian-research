# Terminal-six angular separation for one private side

Suppose two nonnegative integer exponent columns have private
variables on exactly one side. Their ordinary factorial moment
polynomials of orders five and six have disjoint root arguments,
even after separate positive scalings of the argument. In particular
they are coprime over every characteristic-zero field.

The cases of private degree at least three were proved by the
[private-degree criterion](MULTIVARIATE_TERMINAL_SIX_PRIVATE_ANGLES.md).
The new cases have a single private gap \(1\), two private gaps \(1,1\),
or a single private gap \(2\). We link the two coefficient arrays
through a single shared outer mean. The ordinary row structure
strengthens the coefficient bounds exactly where necessary for
the private gap-one case.

## Shared Gauss factors

Orient the private rows to increase and discard the equal rows.
Let \(P\) be the sum of the private gaps. For each shared row let
\(a\ge1\) be the lower exponent and \(d\ge1\) its absolute gap.
Write \(S_0=\sum d^2/(a+d)\) for the shared curvature.
By the [weighted curvature criterion](MULTIVARIATE_TERMINAL_SIX_PRIVATE_ANGLES.md),
we may assume
\[
 S_0<24/5\quad(P=1),\qquad S_0<14/5\quad(P=2).
 \tag{1}
\]
The \(d\) factors of a shared row obtained from Gauss multiplication
have reciprocals
\[
 x_5=\frac{d}{5a+h},\qquad x_6=\frac{d}{6a+h},\qquad 1\le h\le d .
\]
Let \(u,w\) be the outer geometric means of the shared normalized
arrays at orders five and six, respectively. The local contributions
to \(-2\log u\) and \(-2\log w\) are \(\log(1+4x_5)\) and
\(\log(1+5x_6)\). If \(d/a\le R\) in all shared rows, then
\[
 \frac{25}{24}\le\frac{5x_6}{4x_5}
       \le\frac54\,\frac{5+R}{6+R}=:\Lambda,
 \qquad u^\Lambda\le w\le u.
 \tag{2}
\]
Indeed \(x_6/x_5=(5a+h)/(6a+h)\), and
\(1+\lambda(q-1)\le q^\lambda\) for \(q\ge1,\lambda\ge1\).
Multiplying gives (2).

For \(P=2\), a row with \(a=1\) has \(d\le3\), since
\(4^2/(1+4)=16/5>14/5\). When \(a\ge2\) and \(d>(5/2)a\),
its curvature exceeds \(25a/14\), contradicting (1).
Consequently \(x_5\le1/2\), \(x_6\le3/7\), and \(d/a\le3\).

For \(P=1\), a row with \(a=1\) has \(d\le5\). For \(a=2\),
\(d\ge7\) would contribute at least \(49/9>24/5\). For \(a\ge3\)
and \(d>3a\), the curvature exceeds \(9a/4\). Thus \(d\le3a\)
for \(a\ge2\), and altogether \(x_5\le5/6\), \(x_6\le5/7\),
and \(d/a\le5\).

For \(P=1\) and \(u\ge3/8\) there is a sharper estimate.
For \(d\le2\) we already have \(x_5\le1/3\). If \(d\ge3\) and
one Gauss factor has \(x_5>1/3\), then \(5a+1<3d\).
All \(d\) Gauss factors of that row then have \(x_5>1/4\).
The row alone contributes \(u_{\rm row}<2^{-d/2}\le2^{-3/2}<3/8\).
Other rows have outer means at most one, so this is impossible.
Hence \(x_5,x_6\le1/3\) in this range. The largest reciprocal has
\(h=1\), so \(d/a\le(5+1/a)/3\le2\). These bounds give
\[
\begin{array}{c|ccc}
 &x_5\text{ cap}&x_6\text{ cap}&\Lambda\\ \hline
 P=2&1/2&3/7&10/9\\
 P=1,\ u\le3/8&5/6&5/7&25/22\\
 P=1,\ u\ge3/8&1/3&1/3&35/32 .
\end{array}                                                        \tag{3}
\]

## Coefficient envelopes

Write the shared arrays after dividing out binomial coefficients as
\[
 (1,u e^\alpha,v e^\chi,v e^{-\chi},u e^{-\alpha},1),\qquad
 (1,w e^\beta,y e^\eta,z,y e^{-\eta},w e^{-\beta},1).
\]
Put \(U=-\log u\), \(W=-\log w\), and let
\[
 V=\log(v/u^{3/2}),\qquad
 Y=\log(y/w^{8/5}),\qquad Z=\log(z/(yw^{1/5})).
\]
These three defects are nonnegative. For a single factor the
defect is \(G_a(s)=as-\operatorname{artanh}(a\tanh s)\) with
\(a=1/2,3/5,1/5\), respectively, and \(s\) its contribution to
\(U\) or \(W\). This is the
[geometric mean formula](MULTIVARIATE_TERMINAL_SIX_COEFFICIENT_LINKS.md#cross-order-transfers).
Since \(G_a\) is convex and vanishes at zero, \(G_a(s)/s\)
is nondecreasing.

For clarity, the local unsigned skew magnitudes are as follows.
Write \(g_\ell=\log(1+\ell x_j)\), using the order of the particular
array. The quintic bounds are
\[
 |\alpha|_{\rm loc}=\frac{2(g_1+g_2+g_3)-3g_4}{10},\qquad
 |\chi|_{\rm loc}=\frac{4g_2-g_1-g_3-g_4}{10}.
\]
The sextic bounds used here are
\[
 |\beta|_{\rm loc}=\frac{g_1+g_2+g_3+g_4-2g_5}{6},\qquad
 |4\beta-5\eta|_{\rm loc}
   =\frac{3g_1-2g_2-2g_3+3g_4-g_5}{2}.
\]
Changing a row orientation changes signed skews by a sign.
Each displayed magnitude is convex in \(\log(1+(j-1)x_j)\)
and vanishes at zero. For the quintic this is the
[local skew calculation](MULTIVARIATE_TERMINAL_FIVE_ANGULAR_SEPARATION.md)
in equation (2). For the sextic outer skew this follows because
\(g_1,\ldots,g_4\) are convex as functions of that variable and
\(g_5\) is linear. For the linked sextic skew it is the positive
second derivative in the
[within-order bounds](MULTIVARIATE_TERMINAL_SIX_COEFFICIENT_LINKS.md#within-order-bounds).

We bound each magnitude divided by its local \(s\) at the
corresponding upper endpoint in (3), then sum. The triangle
inequality controls the signed sums. In the following table
the first three bounds are multiples of \(U\), and the last
five multiples of \(W\). All decimals denote exact rationals.
\[
\begin{array}{c|ccc|ccccc}
 &|\alpha|&|\chi|&V&|\beta|&|4\beta-5\eta|&|\eta|&Y&Z\\ \hline
 P=2&.14&.07&.04&.15&.03&.13&.04&.02\\
 P=1,\ u\le3/8&.18&.085&.06&.20&.055&.18&.07&.04\\
 P=1,\ u\ge3/8&.11&.055&.025&.13&.017&.11&.03&.015
\end{array}                                                        \tag{4}
\]
The \(\eta\)-column uses the two skew columns to its left.
Every other entry follows from an inequality between rational
powers, checked by raising both sides to a common positive
integer power in the
[certificate](verify_multivariate_terminal_six_one_private_side.py).

In any row of (4), denote the constants for
\(|\alpha|,|\chi|,V,|\beta|,|\eta|,Y,Z\) by
\(A,C,D,B,E,H,J\). By (2), the shared coefficients lie in the
following intervals. Each displayed exponent interval
\([a,b]\) means the multiplier belongs to \([u^b,u^a]\).
\[
\begin{array}{c|c|c}
 \text{order}&\text{indices}&\text{exponent interval}\\ \hline
 5&1,4&[1-A,1+A]\\
 5&2,3&[3/2-D-C,\,3/2+C]\\
 6&1,5&[1-B,\,\Lambda(1+B)]\\
 6&2,4&[8/5-H-E,\,\Lambda(8/5+E)]\\
 6&3&[9/5-H-J,\,9\Lambda/5].
\end{array}                                                        \tag{5}
\]
Multiplying the corresponding bounds by the exact private
coefficients gives all the coefficient enclosures used below.
The same \(u\) governs both orders.

## Moving rays and sector counts

The [certificate](verify_multivariate_terminal_six_one_private_side.py)
uses four continuous rays in each private partition. Their
cosines for \(P=2\) are affine functions of \(u\):
\[
\begin{array}{c|cccc}
 (1,1)&.84-.26u&.60-.50u&-.15-.26u&-.73-.07u\\
 (2)&.84-.19u&.60-.49u&-.15-.21u&-.73-.063u .
\end{array}
\]
For \(P=1\) the first cosine is piecewise affine on the
intervals with endpoints \(0,1/5,3/8,1\), taking the values
\(.84,.735,.605,.145\) there. The second cosine is
\(.55-.82u\). The third is \(-.15-.64u\) up to \(3/8\), and
\(-.39-.352(u-3/8)\) afterwards. The fourth is \(-.70-.16u^2\).
They are ordered from largest to smallest throughout.

We allow each coefficient to vary independently inside its
interval from (5), times its private coefficient. On radii
\(1/5\) and \(5\), the constant and leading terms respectively
dominate for every such polynomial, by the private-only upper
caps. Thus all roots lie between these circles. The exact
projections (5) in the
[fixed-ray proof](MULTIVARIATE_TERMINAL_SIX_PRIVATE_ANGLES.md#exact-ray-inequalities)
apply on each moving ray. The checker bounds the powers in
(5) by rational intervals with denominator \(10000\); an
integer power comparison certifies each endpoint rounding.
It splits the two-dimensional parameter and radius rectangle
into rational boxes and verifies one uniform strictly signed
projection in every box. The full replays use \(69748\), \(1324\)
and \(820\) boxes for the partitions \((1),(1,1),(2)\).

At \(u=1\) the boxes collapse to the private-only polynomials.
The checker counts their roots in the sectors by exact
rational Rouché disks, as in the
[small shared-curvature certificate](MULTIVARIATE_TERMINAL_SIX_SMALL_SHARED.md#exact-angle-certificates).
The first, third and fifth upper sectors contain one sextic
upper root each, and the second and fourth each contain one
quintic upper root. To carry the counts to a given actual
profile, first use a canonical path from \(u=1\) down to its
value of \(u\). In (5) choose shared coefficients with
exponents \(1,3/2,3/2,1\) for the quintic and
\(1,8/5,9/5,8/5,1\) for the sextic. These exponents lie in
every relevant interval, on both sides of the \(3/8\) split.
At fixed \(u\), interpolate from this canonical array to
the actual array within the coefficient boxes. No root
crosses a ray along either path.

The sectors about the real axes are counted as full
conjugate-symmetric sectors. In the final sector the only
quintic root is real, whereas the actual sextic is strictly
positive on the real axis by the factorial integral. All
shared angles are therefore excluded. Ordinary factorial
moments have rational coefficients, so Bézout over
\(\mathbb Q\) proves coprimality over every characteristic-zero
field.
