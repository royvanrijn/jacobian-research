# Effective simple-fold comparison for fixed shared subcritical bands

Use the chord-centered field \(w\), curvature \(W\), and real
stationary curvature \(\rho\) of the
[regular contour cover](MULTIVARIATE_EFFECTIVE_SHARED_REGULAR_COVER.md).
Let \(K\) be a closed bounded field rectangle strictly inside
\(-\pi<\Im w<0\), with computably specified bounds.

**Theorem.** For every positive rational \(W_-\le W_+\) and rational
\(0<\kappa<1\), a terminating finite search gives a uniform
cutoff for consecutive common roots of fully shared ordinary
profiles with \(W_-\le W\le W_+\) and \(\max\rho\le1-\kappa\),
provided both fields belong to \(K\). The integer total-degree
difference is arbitrary. The search includes every simple-fold
scale. Its constants are unevaluated.

## Testing a fold chart

For each limiting profile in the retained and diffuse compact
family, the subcritical stationary map has one interior
simple critical point \(q_c\), with \(w_c=X(q_c)\).
The [uniform cubic transport proof](MULTIVARIATE_UNIFORM_SUBCRITICAL_EQUAL_DEGREE.md#uniformity-at-the-interior-fold)
gives a finite contour isolating this fold, with a strict
phase loss on its exterior chains. The
[active-arc and fold proof](MULTIVARIATE_UNIFORM_SUBCRITICAL_BANDS.md#the-elementary-airy-class-and-its-transverse-band)
shows the oriented local class is one elementary Airy contour
with nonzero integer multiplier. For that orientation the
coordinate \(\lambda\) satisfies
\[
 \lambda(w_c)=0,\qquad |\Im\lambda'(w_c)|>0 .
\]
At the exceptional row the horizontal ray is inactive, so the
rotation selected by the actual contour still has nonzero
imaginary derivative.

The following finite tests certify such a neighborhood. On a
rational disk about a proposed \(q_c\), isolate the zero of \(X'\)
by a linear Rouché test as in the regular cover. The simple
critical point makes this succeed for suitable disks and phase
enclosures. Denote the centered upper phase by \(\phi(q;w)\),
so that \(\phi_q=X-w\). At the critical field the map
\[
 v=(q-q_c)\left(\frac{3(\phi(q;w_c)-\phi(q_c;w_c))}
                         {(q-q_c)^3}\right)^{1/3}
\]
has nonzero derivative. A strict disk enclosure of the quotient
about its nonzero central value selects one of its three roots.
An inverse \(q=T(v)\) on a smaller disk can be certified by
linear Rouché on the preimage disk, uniformly for \(v\) in the
range disk. The inverse is analytic and its Taylor coefficients
have interval enclosures by Cauchy's formula.

To quantify dependence on \(s=w-w_c\), remove \(\phi(q_c;w_c)\)
and solve
\[
 v^3/3-sT(v)=a(s)+(v+h(v,s))^3/3-\lambda(s)(v+h(v,s)).
\]
On the Banach algebra of power series with norm
\(\sum_{k\ge0}|f_k|r^k\), the linearized map in
\((a,\lambda,h)\) is \(a-\lambda v+v^2h\). It is invertible
by removing the constant and linear terms and shifting by
two places. Use the product norm
\(|a|+r|\lambda|+r^2\|h\|_r\), in which this inverse has
norm one. The nonlinear remainder \(vh^2+h^3/3-\lambda h\)
is Lipschitz with a computable constant tending to zero as
the product ball radius decreases at fixed \(r>0\).
The analytic norm of \(T\) is bounded by Cauchy's estimate
on a larger disk. A direct contraction, with Lipschitz
constant at most \(1/2\), therefore solves the equation on a
computable small disk in \(s\). On smaller \(v,s\) disks
Cauchy's estimates give derivative bounds and iteration gives
interval enclosures. In particular \(\lambda'(w_c)=T'(0)\);
select the rotation specified by the oriented contour.
Shrinking the \(s\) disk keeps the derivative in \(v\) of
\(u=v+h\) close to one and makes \(u\) analytically invertible
on a smaller disk. This gives the chart
\(\phi=C(w)+u^3/3-\lambda(w)u\) and computable derivative
bounds for \(q(u,w)\).

For the contour witness, use two fixed cubic decay sectors
in this disk corresponding to the elementary Airy class.
The exits on a fixed small circle in \(u\) have a strict
negative real height relative to \(C(w)\) when \(|\lambda|\)
is small enough. Join them to the finite exterior chains in
those sectors. Approximate the remaining exterior pieces by
rational polygons and check their strict loss and winding as
in the regular cover. If an exterior chain visits the same
local segment with opposite orientation, those integrals
cancel exactly. A witness can record oriented integer
multiplicities and the net two exits. All these tests are
open in the phase and the field. The qualitative fold contour
supplies a witness at each limiting fold by first truncating
its exiting rays inside the decay sectors and then approximating
the strict exterior paths. No computation of a contour
homotopy is required.

The chart may be shrunk so its field disk is convex,
\(|\Im\lambda'|\) is bounded below by a verified \(c>0\),
and the exterior loss is uniformly positive. The rectangle
\(K\) need not contain the full disk, since the witness
only needs to stay inside the open strip.

## Bounds on each certified chart

On the local contour the Stirling error and its derivatives
have computable constants by the
[shared regular estimates](MULTIVARIATE_EFFECTIVE_SHARED_REGULAR_COVER.md#uniform-expansion-and-the-common-root-certificate).
The limiting amplitude after the change to \(u\) is nonzero.
Let \(U=\mathrm{Ai}\) in the chosen elementary orientation
and use the joint norm for \(J=j^{-1/3}U(j^{2/3}\lambda)\)
and \(L=-j^{-2/3}U'(j^{2/3}\lambda)\) from the
[actual cubic proof](MULTIVARIATE_SHARED_SIMPLE_FOLD_AND_REGULAR_CONTOURS.md#a-uniform-joint-norm-for-the-actual-cubic-class).
It bounds the absolute contour mass \(T_j\) by
\[
 T_j\le C\,(|J|+|L|)
     \begin{cases} j^{1/3},&|j^{2/3}\lambda|\le R_0,\\
                    |\lambda|^{-1/2},&|j^{2/3}\lambda|>R_0.
     \end{cases}
\]
All these constants are computable. On the compact scaled disk,
the minimum of \(|U|+|U'|\) is certified positive by interval
subdivision, since Airy's equation makes their simultaneous
vanishing impossible. On the other range use the
[sectorial Airy contour bounds](EXPLICIT_AIRY_CONTOUR_ENVELOPE.md#3-explicit-primary-airy-error-bounds-when-xige4).
The two permissible rotations in the joint norm give leading
derivative signs that are opposite after division by a square
root. Choose \(R_0\) so their explicit relative errors are less
than \(1/4\). The strict decay sectors also give an exponential
bound for the truncated tails. Shrinking the physical chart makes
the exterior loss larger than the joint norm's exponential rate.

Interpolation of the analytic amplitude at \(u=\pm\sqrt\lambda\)
and division by \(u^2-\lambda\) have computable Cauchy bounds
on smaller disks. For the divided remainder use one integration
by parts, gaining \(1/j\), and the joint norm above. The exterior
contribution is exponentially small relative to the joint norm.
The constant coefficient of the interpolant is bounded away
from zero and the linear coefficient is uniformly bounded.
Absorbing the \(O(j^{-2/3})\) relative error therefore gives,
at a root, \(|U(\xi)|\le C_0j^{-1/3}|U'(\xi)|\) for
\(\xi=j^{2/3}\lambda\), with \(C_0\) computable.

For clarity, set \(\gamma=-U(\xi)/U'(\xi)\).
The denominator is nonzero at such a root by Airy's equation.
Choose the chart and \(j\) so \(|\gamma^2\xi|\le1/4\) and
\(|\gamma|^3\le1/16\). Along the backward homotopy
\(U(z)+t\gamma U'(z)=0\) from \(z(1)=\xi\), simple
zeros satisfy \(z'(t)=-\gamma/(1-t^2\gamma^2 z(t))\).
While \(|z-\xi|\le2|\gamma|\), the denominator has
modulus at least \(5/8\), so this bound continues to \(t=0\).
Since \(U\) has only negative real zeros,
\[
                  |\Im\lambda(w)|\le 2C_0/j .
\]
At two roots in a common convex witness disk,
\[
 c\,|\Delta|\log(1+1/n)
   \le |\Im(\lambda(w_N)-\lambda(w_n))|\le 4C_0/n .
\]
Thus one need only test integers \(|\Delta|\le8C_0/c\).

For every nonzero integer in that finite range use the
[same-contour two-order comparison](MULTIVARIATE_SHARED_SIMPLE_FOLD_AND_REGULAR_CONTOURS.md#comparing-both-orders-on-the-same-cubic-contour).
At the fold its interpolation matrix tends to a matrix with
determinant
\[
 (A(q_c)q_u(0,w_c))^2
     e^{\phi(q_c;w_c)+\Delta q_c}\,\Delta q_u(0,w_c).
\]
This has a computable positive modulus bound from the witness.
The finite phases, the Gaussian sectors, and the analytic
amplitude estimates give a computable perturbation bound, by
shrinking the field disk and increasing the cutoff. The
\(O(T_n/n)\) remainder divided by \(|J|+|L|\) is at most
\(C n^{-2/3}\) on both scaled ranges. This excludes all nonzero
integers in the range. The upper order is evaluated on the lower
order's contour by the same winding identity, and changing the
parameter and exponent contributes precisely the analytic
multiplier in that comparison.

For \(\Delta=0\) and \(|n^{2/3}\lambda|\) large, that determinant
has main term \(-\frac23 e^{C(w)} A_0(v)A_0(-v)\lambda\)
with \(v=\sqrt\lambda\), and absolute error
\(O(|\lambda|^4+n^{-1})\) with computable constants.
The inverse matrix has norm at most \(C/|\lambda|\) there.
The resulting joint-norm relative error is bounded by
\(C_1/(n|\lambda|^{3/2})\) plus an exponentially decreasing
term. Increasing the fixed scaled cutoff excludes this range.

On the remaining compact scaled disk use the
[real-order interpolation on the truncated rays](MULTIVARIATE_SHARED_SIMPLE_FOLD_AND_REGULAR_CONTOURS.md#bounded-scaled-fields-for-equal-total-degrees).
Its absolute integral and parameter derivative bounds reduce to
integrals of a fixed polynomial times \(e^{-c r^3+C r}\) on
\(r\ge0\), all computably bounded by splitting at a rational
large enough radius. The real-order Stirling derivative is
bounded via a complex-order Cauchy disk of fixed relative
radius, since all local Gamma arguments remain in strict
sectors there. Isolate the finitely many Airy zeros on a
slightly larger disk with a zero-free boundary. For example,
subdivide using the Taylor expansion from Airy's equation and
certify simple isolating disks by Rouché. The zeros are simple
and none equals zero, so subdivision gives positive derivative
and separation bounds and a positive function bound outside
their disks. The local implicit branches then satisfy
\(\xi(\epsilon)=\xi_0+O(\epsilon)\),
\(\xi'(\epsilon)=O(1)\), with computable constants for
\(\epsilon=n^{-1/3}\). In physical coordinates the derivative is
\(2\epsilon\xi_0+O(\epsilon^2)\). Projection onto \(\xi_0\)
bounds the displacement below by a computable multiple of
\(n^{-5/3}\). The exterior perturbations are exponentially
small in the scaled coordinate, and distinct branches stay
separated after the adjacent scaling. This excludes the
remaining equal-degree roots.

## Closing the compact subcritical interior

Use the paired field cells of the regular cover with equal
imaginary part. Off the closed actual active arc there is a
dominant-term witness. At any pair on the closed arc the
fields agree by its strictly monotone imaginary height.
For a regular pair use a common convex two-saddle chart as
before. At the fold use a common convex fold chart with the
positive losses and the three comparisons above. The phase
and interval tests persist in neighborhoods of the limiting
profiles and fields. Enumerating regular, dominant, and fold
witnesses in parallel with the profile and field nets therefore
terminates by compactness. Taking the maximum of the finite
local thresholds gives the claimed cutoff.
