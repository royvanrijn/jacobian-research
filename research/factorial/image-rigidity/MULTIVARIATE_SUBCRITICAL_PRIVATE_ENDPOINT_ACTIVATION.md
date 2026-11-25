# Rice activation for a subcritical private linear endpoint

Let \(a,b\) be fixed ordinary factorial exponent columns
with exactly one private side, oriented so
\[
 \sum_{a_i=0} b_i=1,\qquad b_i>0
\]
on all varying rows. Delete equal rows, set
\(d_i=b_i-a_i,\ \Delta=\sum_i d_i\), and put
\[
 K_a=\prod_{a_i>0}a_i^{d_i},\qquad
 \sigma_a=\sum_{a_i>0}d_i^2/a_i .
\]
Assume \(\sigma_a<1\). The case where there is no varying
shared row is a constant and a linear variable after
deleting equal rows, with coprimality from truncated
exponentials. In the rest of this note there is a
varying shared row.

**Theorem.** On every fixed compact subset of
\(\Im\zeta>0\), with \(\zeta=n^\Delta t\), the actual
consecutive moment pair has no common zero for all
sufficiently large orders. The conjugate assertion holds
in the lower half-plane. For such profiles any unexcluded
sequence of common zeros must converge to the positive
real scale \(1/K_a\) with unbounded half-Gaussian
scaled field.

The private left contour can include the entire right
saddle contribution. The inclusion is determined by the
sign of the imaginary difference between the saddle and
boundary phases. This avoids a vertical-line bound on the
pinched left contour when \(\arg\zeta\) is small.

## The endpoint arm and the saddle barrier

Use the [exact Rice contour](MULTIVARIATE_RICE_CONTOUR_AND_REAL_BOUNDARY.md)
with its upper phase \(\phi\) and lower phase \(\phi_-\).
Write \(q_*(\zeta)\) for the unique upper saddle and
\(S(\zeta)=\phi(q_*;\zeta)-\phi(0;\zeta)\).
The boundary phase \(\phi(0;\zeta)\) is real and is
independent of \(\zeta\). The
[right contour proof](MULTIVARIATE_ONE_PRIVATE_SIDE_RICE_EVENTUAL.md#the-actual-right-contour)
provides its actual descending curve, with one end
upper and one lower, and the relative expansion
\(n^{1/2}J(\zeta)e^{nS(\zeta)}(1+O(n^{-1}))\)
after division by the pure \(a\) moment. The proof of
that right period does not assume the left vertical
estimate.

At \(q=0\) the private linear term cancels the logarithm
of zero in \(X\). Both phase branches extend analytically
to zero, with derivatives
\[
 k_+=\log(\zeta K_a),\qquad k_-=k_+-2\pi i .
\]
The branch for \(k_+\) has imaginary part
\(\theta=\arg\zeta\in(0,\pi)\). Run descending arms
from zero using the upper and lower phases, respectively.
They initially point into their respective half-planes
because the vertical component of
\(\dot q=-\overline{\phi'(q)}\) equals \(\Im\phi'(q)\).
The lower arm necessarily goes to lower infinity.

For the upper arm, the real boundary, nonaccumulation,
and outward-region arguments for the right contour apply
also here. It can reach lower infinity by crossing the
real right gap, reach upper infinity directly, or
approach the upper saddle. The last possibility requires
\(\Im S=0\) and \(\Re S<0\), by the constant imaginary
part and strict real phase loss. These alternatives are
exhaustive by the same escape estimates: finite real
branch points and non-gap limits are excluded, and a
bounded trajectory away from the real cuts must have a
stationary limit point.

Let \(R^+\) be the part of the actual right saddle
contour in the upper half-plane, from its real-gap
crossing to upper infinity. It is simple: its two
saddle arms cannot intersect at a noncritical point
by uniqueness of the descent flow, and they join
smoothly in the local quadratic coordinate.
It separates the upper half-plane into left and right
sides. Near zero the boundary arm is on the left.
If \(\Im S\ne0\), the two paths have distinct imaginary
phase levels and cannot intersect.

On the right real gap the imaginary height of the
upper phase, relative to its value at zero, is
\(\pi-(\pi-\theta)q\). Thus \(R^+\) starts at
\[
 x_R=\frac{\pi-\Im S}{\pi-\theta},
\]
and a boundary exit would have to be at
\(x_B=\pi/(\pi-\theta)\).
When \(\Im S>0\) one has \(x_B>x_R\). The boundary
arm cannot cross the gap on the right side of \(R^+\),
so it goes to upper infinity. When \(\Im S<0\), an
upper escape would put both paths at high positive
levels \(y\). At these heights the horizontal segment
between their positions is in the preferred escaping
half-plane if \(\Delta\ne0\). The large-height
inequality \(\Im\phi'>c>0\) holds on this segment,
by the outward-region estimate. It holds on the
whole segment also when \(\Delta=0\), by the
constant derivative limit. Therefore
\(\Im\phi(x+iy)\) increases strictly along the segment
from left to right. The heights zero and \(\Im S<0\)
would place the boundary arm to the right of the
saddle arm, contradicting separation. The boundary
arm must instead cross the right gap and escape
downwards. Here left and right at sufficiently large
height agree with the horizontal order, since the
outward vertical velocity gives just one intersection
of \(R^+\) with each such line.

If \(\Im S=0\), a nonlanding boundary arm could neither
cross the gap nor escape upwards, since these
possibilities would give, respectively, the same
crossing point or the same large-height point as
\(R^+\). Uniqueness excludes such intersections
away from the saddle. Thus the arm lands at the
saddle and \(\Re S<0\). In particular a point of
equal real phase has a nonzero imaginary phase
difference. This also follows from \(S'(w)=-q_*(w)\)
in the strip variable \(w=-\operatorname{Log}\zeta\):
\(iS'\) has positive real part, so \(S\) is injective.
It has boundary limit zero at the real scale \(1/K_a\),
where \(X'(0)=\sigma_a-1<0\) gives the upper inverse
approaching \(q=0\).

## Exact contour classes and the local endpoint integral

Use the pinching left vertical line \(q=-c/n+i y\)
with \(0<c<1\). The passage at \(-c/n\) is between the
coefficient pole zero and the first private Gamma pole
\(-1/n\). For fixed \(\zeta\) with \(\Im S\ne0\)
form a contour \(B\) along the two boundary arms, oriented
from lower infinity towards zero and then out on the upper
arm, with a small replacement through that left passage.
If the upper arm escapes upwards, \(B\) is in the same
contour class as the left vertical line. If it crosses the
right gap downwards, add the upward right saddle contour.
Then the right-gap crossings cancel, and the sum of the
two contours is in the class of the left vertical line.
These claims can be made on domains truncated at heights
\(\pm Y\): join ends at the same height, slide the crossings
within their pole-free gaps, and deform the remaining
segments in the open half-planes. There are no poles in
either half-plane, and the crossing signs and order give
the stated contour identities. The horizontal joining
integrals vanish as \(Y\to\infty\) by the same
reflection-Stirling and preferred-direction estimates as
for the actual right contour. Both the saddle and
boundary escaping arms are descending curves.

In the notation in which the total moment is the right
vertical minus the left vertical, this yields the pure
boundary contribution plus the right saddle for
\(\Im S>0\), and the pure boundary contribution alone for
\(\Im S<0\). The sign of the boundary contribution is
the one in (4) of the
[one-private contour note](MULTIVARIATE_ONE_PRIVATE_SIDE_RICE_EVENTUAL.md#the-pinching-left-contour).
The contour classes are constant on small parameter
neighborhoods away from \(\Im S=0\), by finite-time
continuous dependence and the outward regions.

It remains to show that \(B\) gives the actual boundary
amplitude even at parameters where the full vertical
Barnes bound fails. With \(z=nq\), the private factor in
the reflected integrand is
\(\Gamma(-z)\Gamma(1+z)=-\pi/\sin(\pi z)\).
Retain the central segment of the original vertical line
\(\Re z=-c,\ |\Im z|\le C\log n\). Use horizontal
connectors at those two heights to the boundary arms.
For \(q\) in a fixed small disk about zero, all active
Gamma arguments and \(1+n(1-q)\) have positive real
part comparable to \(n\). On the upper or lower side
with \(|\Im z|\ge C\log n\), the dominant sine exponential
gives the analytic phase \(\phi\) or \(\phi_-\) with its
local bounded amplitude and an overall factor of order
\(n\) in the \(dq\) variable. On the connectors the real
phase relative to its value at zero is bounded above
by \(-c_1|\Im q|\). Indeed at the endpoint on the
original line its linear part is
\(-|\Im k_\pm|\,|\Im q|+O(1/n)\), and on the boundary
arm it is strictly negative with the same order of
magnitude, by the initial descent direction.
On the intervening horizontal segment the linear term
has its maximum at an endpoint, and the quadratic
remainder is \(O((\log n)^2/n^2)\).
The same strict linear bound holds on the arm out to
a fixed small height. Beyond that height there is a
fixed real phase loss, by descent and the escape
estimates. Choosing \(C\) sufficiently large, the two
connectors and the remaining arms give \(O(n^{-2})\)
relative to the boundary exponential. All bounds hold
locally uniformly in the parameter when the end types
are fixed.

On the central segment, the fixed-shift Gamma expansion
of \(T_n(z)\) from the one-private contour note applies
uniformly for \(|z|=O(\log n)\).
The *model* integrand
\(\Gamma(-z)\Gamma(1+z)(-\zeta K_a)^z\) has
uniform exponential decay on the full line
\(\Re z=-c\). Extending the truncated comparison
integral to that line gives
\[
 B_{\rm phys}(\zeta)
     =\frac1{1-\zeta K_a}+\frac{B_1(\zeta)}n+O(n^{-2}),
 \tag{1}
\]
where \(B_{\rm phys}\) includes the negative sign from
subtracting the left vertical integral and the common
Rice prefactor divided by the pure \(a\) moment.
This uses only the local fixed-shift expansion of
\(T_n\), so no uniform bound on \(T_n\) along the
untruncated left vertical line is required.

At a parameter with \(\Im S=0\), the upper boundary
arm approaches the saddle and the real saddle value
is strictly below the boundary value. To form a
reference contour, stop that arm inside a small disk
of the saddle and join to the outgoing upper saddle
arm inside the disk. The real phase remains strictly
below the boundary value on this detour and beyond
it. This contour represents the left vertical integral
and has expansion (1). Its bounded pieces keep the
same gap under a small change of parameter, and its
outward tail estimates are uniform. The right saddle
is exponentially smaller throughout that neighborhood.
Thus the full moment is nonzero there for all
sufficiently large orders.

## Order comparison on the active arc

Any limiting complex common zero must now satisfy
\(\Re S=0\) and \(\Im S>0\), because away from
equal moduli there is a single dominant term, and for
\(\Im S<0\) the saddle cancels exactly.
With \(\Delta=0\), the \(n^{1/2}\) on the saddle
and the constant power on the boundary exclude
consecutive cancellation as in the one-private note.
For \(\Delta\ne0\), the order increment would
force \(\Re q_*=0\) at the limiting point.

This last possibility is on the inactive arc.
Indeed let \(q_*=iy\), \(y>0\), and put
\(\ell_i=a_i/|d_i|>1\) for active rows; the strict
inequality follows from \(\sigma_a<1\).
The equality of leading moduli reads
\[
 \sum_{a_i>0}a_i\log(1+y^2/\ell_i^2)=\log(1+y^2).
\]
Set \(H(t)=t-\arctan t\). The ratio
\(H(t)/\log(1+t^2)\) is strictly increasing for
\(t>0\), since it is the mean of \(s/2\) against
\(2s/(1+s^2)\) on \([0,t]\).
The equality of moduli thus implies
\(\sum_{a_i>0}a_i H(y/\ell_i)<H(y)\).
But the imaginary phase difference is
\[
 \Im S=
   \sum_{d_i<0}a_i H(y/\ell_i)
   -\sum_{\substack{a_i>0\\d_i>0}}a_i H(y/\ell_i)
   -H(y)<0 .
\]
Here the final \(-H(y)\) uses the single private
linear exponent in \(\Delta\).
This contradicts activity and completes the complex
compact exclusion.

The fixed scaled annulus and the real noncritical
neighborhoods were proved in the
[one-private contour note](MULTIVARIATE_ONE_PRIVATE_SIDE_RICE_EVENTUAL.md#eventual-exclusion)
and the
[private boundary analysis](MULTIVARIATE_PRIVATE_SUPPORT_SADDLE_GEOMETRY.md#continuous-endpoint-windows).
The remaining sequences for this profile class must
therefore approach the positive transition \(1/K_a\)
with \(n^{1/2}|\log(\zeta K_a)|\to\infty\).
