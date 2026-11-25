# Escaping fields at the exact real cusp of a fully shared support

Let \(a,b\) be distinct fixed ordinary fully shared exponent
profiles. Write \(\Delta=\sum_i(b_i-a_i)\) and let \(\rho(q)\)
be the curvature of the
[conditional saddle equation](MULTIVARIATE_SHARED_SUPPORT_CONDITIONAL_SADDLE.md#stationary-equation-and-its-critical-points).
Suppose \(\max\rho=1\). Denote the real cusp by \(q_c\),
its conditional coordinate by \(s_c\), and the logarithm
of its scaled external parameter by \(\beta_c\). Put
\(\beta=\log(n^\Delta t)-\beta_c\), using the local
branch at the positive cusp scale.

**Theorem.** A fixed complex neighborhood of the real cusp
scale excludes consecutive common roots for all sufficiently
large orders. This includes both equal and unequal total degrees.

The [bounded quartic window](MULTIVARIATE_SHARED_SUPPORT_CONDITIONAL_SADDLE.md#bounded-scaled-windows-at-the-cusp)
reduces the new assertion to \(\beta\to0\) with
\(n^{3/4}|\beta|\to\infty\). The following contour comparison
holds for every such sequence, including fields of complex
argument close to the two-sector transition.

## A scaled phase and its exact contour

Use the analytic scalar phase \(\phi(s,\beta_c+\beta)\)
from the conditional integral. The stationary jet there gives
\[
 \phi(s_c+h,\beta_c+\beta)-\phi(s_c,\beta_c+\beta)
 =g h^4+K\beta h+O(|h|^5+|\beta||h|^2+|\beta|^2|h|),
 \qquad g<0,\quad K=q_c(1-q_c)>0.
 \tag{1}
\]
The error is analytic on a fixed bidisk, with the indicated
factors. Let \(\epsilon=|\beta|^{1/3}\), choose
\(\ell=(K/(-4g))^{1/3}\), put \(c=K\ell>0\), and take
\(h=\ell\epsilon z\), \(M=n c\epsilon^4\),
\(\theta=\arg\beta\). After subtracting its value at \(s_c\)
and dividing by \(c\epsilon^4\), the phase is
\[
 \Psi(z)=e^{i\theta}z-z^4/4
          +O\bigl(\epsilon(|z|^5+|z|^2)+\epsilon^3|z|\bigr).
 \tag{2}
\]
The limit and its derivatives are uniform on every fixed
\(z\)-disk. Here \(M\to\infty\) on the escaping scale.

For this fixed profile the complete real conditional contour
can first be restricted to \(|s-s_c|\le r_0\), with a fixed
strict phase gap on the omitted integral for all small
complex \(\beta\). Indeed the real cusp has a unique global
maximum at \(s_c\), compact real upper level sets, and the
same absolute integral estimate as the
[real saddle proof](MULTIVARIATE_SHARED_SUPPORT_CONDITIONAL_SADDLE.md#a-local-exclusion-around-the-positive-real-axis).
On the retained contour, after a common nonzero order
factor has been removed, the analytic amplitude is
\(A_n(s)=A_0(s)+O(n^{-1})\), uniform on a fixed complex disk
about \(s_c\), with \(A_0(s_c)>0\). This is the transverse
Gaussian expansion in the conditional integral.
The real interval can be shifted by \(O(\epsilon)\) in
height, keeping its ends connected to the original endpoints.
The end connectors retain a strict physical phase gap.
One can also deform the two halves separately after choosing
a split point of strict loss on such a shifted line.

For clarity, all the shifted horizontal contours below
have bounded \(z\)-heights. Uniformly for bounded
\(|\Im z|\) and \(|\Re z|\) sufficiently large with
\(|\ell\epsilon z|\le r_0\), (2) gives
\(\Re\Psi(z)\le-c_1|\Re z|^4\), by taking \(r_0\)
sufficiently small. The term \(\epsilon |z|^5\) is
bounded by a small multiple of \(|z|^4\) there.
Consequently the scaled tails have a uniform quartic
majorant. The local pieces and connectors stay in the
analytic conditional strip.

## One-saddle and two-saddle contours

By conjugation take \(0\le\theta\le\pi\). For the limiting
phase \(G_\theta(z)=e^{i\theta}z-z^4/4\) and
\(0\le\theta<\pi/2\), the saddle \(z_+=e^{i\theta/3}=u+iv\)
is the unique maximum of the real phase on its horizontal
line. In fact on that line, up to an additive constant,
\[
 \Re G_\theta(x+iv)=-x^4/4+3v^2x^2/2+x\cos\theta.
\]
A negative \(x\) gives strictly smaller value than \(-x>0\).
The derivative has exactly one positive real zero, \(x=u\),
by its cubic with coefficients \(-1,3v^2,\cos\theta\)
in degrees \(3,1,0\). It has negative second derivative
at \(u\) since \(0\le\arg z_+<\pi/6\).
For \(\pi/2<\theta\le\pi\) take the negative saddle
\(-e^{-i(\pi-\theta)/3}\) and reflect the same argument.
On any closed angular region separated from \(\pi/2\),
these maxima and all gaps off fixed saddle neighborhoods
are uniform. For (2), a simple saddle persists near the
model one and the horizontal line through it has the same
local negative curvature and strict global gap.
The Gaussian estimate on this actual contour is uniform
and has nonzero leading term. In particular the lower-order
polynomial is eventually nonzero in those regions.

At \(\theta=\pi/2\), the relevant limiting saddles are
\(z_+=\sqrt3/2+i/2\) and \(z_-=-\sqrt3/2+i/2\).
On the line of height \(1/2\), the two horizontal half-lines
have unique respective maxima at those points, with negative
curvatures. Both have real phase \(-3/8\).
At the midpoint \(i/2\) the real phase is \(-33/64<-3/8\).
Choose a fixed sufficiently small angular neighborhood of
\(\pi/2\), let both simple saddles of (2) continue, and
shift the entire local contour first to height \(1/2\)
in the \(z\)-coordinate. Split there at \(i/2\), and
shift each half to the height of its own saddle,
using short connectors near the midpoint.
The midpoint connectors have a uniform strict gap from
both saddles, by the displayed strict inequality.
On each shifted half-line strict loss away from its
own saddle persists, and the tails obey the quartic
bound following (2). In a fixed small horizontal
neighborhood of either saddle, the real phase drops
by at least \(c_2|z-z_j|^2\). Gaussian comparison on
each half therefore gives the actual two packets
\(I_{n,L}\) and \(I_{n,R}\), with nonzero Gaussian
amplitudes. More precisely, after the common normalization
and the Jacobian, each integral is comparable in modulus
to \(\epsilon M^{-1/2}e^{n\Re\phi(s_j,\beta_c+\beta)}\).
Here and below one can divide out the common
\(e^{n\phi(s_c,\beta_c+\beta)}\) from the comparison.
On each half the integral of the absolute value of
its Gaussian majorant has the same size. This follows
by substituting \(\sqrt M(z-z_j)\) in the local integral;
all other parts are smaller by a uniform strict
scaled or physical gap.

## Consecutive packet comparison

For a fixed original \(t\), the next logarithmic field is
\(\beta'=\beta+\Delta\log(1+1/n)\). At order \(n+1\)
retain exactly the same two local contours as at order
\(n\). Their identities as a deformation of the real
conditional integral are independent of the order.
The exterior chain stays below a fixed physical gap
at \(\beta'\) as at \(\beta\).
On a common analytic neighborhood the change in the
integrand, apart from the amplitude ratio, is precisely
multiplication by \(e^{D_n(s)}\), where
\[
 D_n(s)=(n+1)\phi(s,\beta_c+\beta')-n\phi(s,\beta_c+\beta).
\]
The amplitude ratio differs from one by \(O(n^{-1})\)
uniformly, since \(A_0\) is nonzero on the local disk.
For \(\Delta\ne0\) Taylor expansion gives
\[
 D_n(s)=\phi(s,\beta_c+\beta)
    +\Delta\,\phi_\beta(s,\beta_c+\beta)+O(n^{-1})
\]
with the same error bound on \(s\)-derivatives.
It follows that \(|D_n(s)-D_n(s_j)|\le C\epsilon|z-z_j|\)
on the fixed local disk. For \(\Delta=0\) the increment
is exactly \(\phi(s,\beta_c+\beta)\), so the difference
is \(O(\epsilon^4 |z-z_j|^2)\) near \(z_j\),
and \(O(\epsilon^4(1+|z|^4))\) over the retained
contours by (2).

These bounds permit comparison of each packet and its
successor by the same contour. After dividing by
\(e^{D_n(s_j)}\), insert
\(e^{D_n(s)-D_n(s_j)}-1\) into the lower integral
and use \(|e^u-1|\le |u|e^{|u|}\).
The local Gaussian width \(M^{-1/2}\) in \(z\)
gives an error \(O(\epsilon M^{-1/2})\) relative to
the packet when \(\Delta\ne0\), and
\(O(\epsilon^4 M^{-1})\) when \(\Delta=0\).
On bounded pieces away from the saddle the same
polynomial weights are absorbed by the strict scaled
loss, and on unbounded scaled pieces by the quartic
majorant. The connectors at the fixed physical
endpoints have a fixed phase gap. Including the
amplitude change, the result is
\[
 \frac{I_{n+1,j}(\beta')}{I_{n,j}(\beta)}
   =e^{D_n(s_j)}
      \begin{cases}
       1+o(\epsilon), & \Delta\ne0,\\
       1+o(\epsilon^4),& \Delta=0 .
      \end{cases}                                          \tag{3}
\]
The common scalar prefactors have been set consistently
in the two packet definitions. This comparison uses
\(n^{-1}=o(\epsilon^4)\), which follows from \(M\to\infty\).

For unequal total degrees,
\[
 D_n(s_R)-D_n(s_L)
    =\Delta K(s_R-s_L)+O(\epsilon^2+n^{-1})
\]
by stationarity and the mixed derivative \(K\).
The leading displacement has magnitude bounded below by
a positive multiple of \(\epsilon\), uniformly in the small
angular neighborhood. For equal degrees the difference
is exactly the stationary phase difference. In the scaled
model at \(\theta=\pi/2\) this is \(3i\sqrt3/4\),
as \(G_{\pi/2}(z_j)=3z_j^4/4\). Therefore the
unscaled difference is bounded below in modulus by
a positive multiple of \(\epsilon^4\), on a sufficiently
small angular neighborhood and for \(\epsilon\) small.
In both cases the two exponential multipliers in (3)
differ by more than their comparison errors.

At a common root the two packets cancel at each order
up to the exponentially small exterior remainder.
Every packet has size at least an inverse polynomial
factor times the common centered exponential and
\(e^{-C M}\), by the uniform Gaussian estimate.
Since \(M=o(n)\), the exterior remainder is negligible
even after the packet comparisons and division by
\(\epsilon^4\). Cancellation at both orders would make
their two adjacent ratios agree to that precision.
This contradicts (3). All escaping angular sequences
have now been excluded. The bounded quartic windows
supply the rest of the local neighborhood, proving
the theorem.

## Eventual coprimality for every fixed ordinary support

**Theorem.** For every fixed pair of distinct nonnegative
integer exponent vectors \(a,b\), there is a terminal
threshold \(N_0(a,b)\) such that the consecutive ordinary
factorial moment polynomials are coprime for all
\(N=n+1\ge N_0(a,b)\). The same threshold works over
every characteristic-zero field.

Remove any equal rows. If there is a private variable
on each side, the
[two-private-side theorem](MULTIVARIATE_TWO_MONOMIAL_EVENTUAL_BORDER.md)
gives the assertion. If there is only one private side,
use the [stable endpoint theorem](MULTIVARIATE_ONE_PRIVATE_SIDE_RICE_EVENTUAL.md),
the [subcritical endpoint theorem](MULTIVARIATE_HALF_GAUSSIAN_ENDPOINT_ESCAPING_FIELDS.md),
and the [critical linear-private endpoint theorem](MULTIVARIATE_CUBIC_ENDPOINT_ESCAPING_FIELDS.md).
These exhaust the private degree and linear transition cases.

For fully shared supports with \(\max\rho>1\), use the
[supercritical theorem](MULTIVARIATE_SHARED_SUPERCRITICAL_EVENTUAL_COPRIMALITY.md).
For \(\max\rho<1\), use the
[unequal-degree theorem](MULTIVARIATE_SHARED_SIMPLE_FOLD_AND_REGULAR_CONTOURS.md#consequences-away-from-the-real-cusp)
or the [subcritical equal-degree theorem](MULTIVARIATE_SHARED_SUBCRITICAL_EQUAL_DEGREE.md).
For \(\max\rho=1\), the roots lie in the fixed scaled
annulus at large orders. Every noncusp real limit is excluded
by the [conditional local comparison](MULTIVARIATE_SHARED_SUPPORT_CONDITIONAL_SADDLE.md#a-local-exclusion-around-the-positive-real-axis).
The real cusp limit is excluded by the local theorem above.
A nonreal limit has only regular upper stationary points
(up to conjugation). The
[actual regular contour](MULTIVARIATE_SHARED_SIMPLE_FOLD_AND_REGULAR_CONTOURS.md#regular-values-and-nonzero-coefficients)
has at least one nonzero saddle coefficient, and either gives a
dominant term or requires equality of the consecutive
phase multipliers at limiting common roots.
The [phase separation results](MULTIVARIATE_REGULAR_SADDLE_PHASE_SEPARATION.md)
exclude that equality in both total-degree cases when
\(\max\rho=1\). Compactness of the fixed annulus proves
eventual coprimality.

The polynomials have integer coefficients, so rational
Bézout gives the field assertion. This theorem is
support-dependent. A bound adequate for profiles that
vary with the order, and coprimality below the eventual
threshold for each profile, remain separate obligations.
