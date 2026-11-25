# Uniform regular comparison for one private side

Orient an ordinary factorial profile with private variables on
exactly one side so the private rows are \((0,e_i)\), \(e_i\ge1\).
Delete equal rows and let the remaining shared rows be \(a_i,b_i>0\),
with \(d_i=b_i-a_i\). Write
\[
 P=\sum_{\rm private}e_i,\quad p=\#\{\text{private rows}\},\quad
 S_{\mathrm{curv}}=P+\sum_{\rm shared}\frac{d_i^2}{\max(a_i,b_i)},
 \quad
 \sigma=\sum_{\rm shared}\frac{d_i^2}{a_i},\quad
 \Delta=P+\sum_{\rm shared}d_i .
\]
**Theorem.** On every fixed upper band \(S_{\mathrm{curv}}\le K\),
the profiles with \(P\ge2\) have a uniform eventual terminal
cutoff for consecutive coprimality. The same holds for \(P=1\)
and \(\sigma\ge1+\kappa\), with any fixed \(\kappa>0\).
If \(P=1\) and \(\sigma\to1\) along a sequence of common roots
with unbounded terminals, both centered fields
\[
 \xi_j=\log(j^\Delta t)+
        \sum_{\rm private}e_i\log e_i+\sum_{\rm shared}d_i\log a_i,
                 \qquad j=n,n+1,                         \tag{1}
\]
must tend to zero. The total-degree difference is arbitrary.

For \(P=1\) and \(\sigma\to\sigma_0<1\), the proof also
excludes any common-root limit with \(\Im\xi_j=\theta_j\to
\theta>\theta_{\sigma_0}\), where
\[
 \theta_s=\sup_{y\ge0}\{\arctan y-\arctan(sy)\}
          \quad(0\le s<1).
\]
Any common-root limit at height zero then must have both
fields tending to zero.

## Compact phases and the left integral

Suppose there were a sequence of counterexamples in one of the
stated scopes. We use subsequences throughout. The private rows
have a fixed allocation after extraction, since \(P\le K\).
The
[bounded-private-degree root annulus](MULTIVARIATE_BOUNDED_PRIVATE_ROOT_ANNULUS.md)
bounds the normalized root moduli above and away from zero.
The real parts of (1) are therefore bounded. In fact the
endpoint log center is
\[
 \log G_j=\Delta\log j+\sum_i(F(b_i)-F(a_i))
               +O((V_{\rm sh}+P\log j)/j),\qquad
 F(x)=x\log x-x,\quad F(0)=0,
\]
by the
[endpoint-center estimate](MULTIVARIATE_GROWING_CURVATURE_ENDPOINT_EXPANSION.md#comparison-of-the-orders).
Here \(V_{\rm sh}=\sum_{\rm shared} d_i^2/\min(a_i,b_i)\)
is bounded. Subtracting the row constant in (1) changes
the leading center by at most \(P+O(V_{\rm sh})\), by
integrating \(\log(a_i+d_iq)-\log a_i\) on shared rows.
Conjugation permits \(0<\theta_j=\Im\xi_j<\pi\), since
the even-order moment is strictly positive on the real
physical axis. The two-order field shift is exact and real,
\[
                         \xi_{n+1}-\xi_n
                                   =\Delta\log(1+1/n).    \tag{2}
\]

Write
\[
 B(q)=\sum_{\rm shared}d_i\log(1+d_iq/a_i),\qquad
 X(q)=(P-1)\log q+\log(1-q)+B(q).
\]
The phase normalized to value zero at the pure \(a\) boundary
is
\[
 \Phi(q;\xi)=\xi q+(P-1)(q\log q-q)
                +\int_0^q\{\log(1-z)+B(z)\}\,dz .          \tag{3}
\]
Use the upper logarithms and the conjugate lower logarithms
in the two respective half-planes. The lower Rice phase
includes the extra term \(-2\pi i q\).
Equation (3) follows from the
[exact Rice formula](MULTIVARIATE_RICE_CONTOUR_AND_REAL_BOUNDARY.md#exact-contour-and-its-phases)
after subtraction of its pure \(a\) value. In particular
\(\partial_q\Phi=X+\xi\) on the upper side.

By
[retained and diffuse compactness](MULTIVARIATE_UNIFORM_SADDLE_MAP_COMPACTNESS.md#retained-and-diffuse-limits),
\(B\) and its derivatives converge on compact complex sets
off the cuts, and across bounded real gap charts, to the
retained shared sum plus \(\tau q\) for \(\tau\ge0\).
The same estimates hold on \([0,1]\) for \(B\). One can see
this boundary-centered form by integrating the curvature
transform from zero. The summable retained tail is
analytic across any bounded real interval after discarding
finitely many poles. Every fixed real bounded interval
has no accumulations of retained poles. The limit of
\(\sigma=B'(0)\) is the retained endpoint curvature plus
\(\tau\). We write \(X_*,\Phi_*\) for limiting phases.
The same compactness theorem gives exactly one simple
upper saddle \(q_*(\xi)\) for \(0<\Im\xi<\pi\), with
analytic convergence of nearby-profile saddles on compact
neighborhoods. The lower phase has no stationary point,
since its vertical descent velocity is strictly negative.

For \(P\ge2\) the pinched left vertical integral, relative to
the pure \(a\) moment, is
\[
 {\cal L}_j(\xi)=1+O(j^{-1}) .                              \tag{4}
\]
This is the actual integral on \(\Re q=-c/j\), \(0<c<
\min_{\rm private}1/e_i\), with its sign in the residue formula.
Use the
[Mellin–Barnes formula](MULTIVARIATE_ONE_PRIVATE_SIDE_RICE_EVENTUAL.md#the-pinching-left-contour)
and shift \(z=jq\) from real part \(-c\) to \(5/2\).
Its shared ratio product at any fixed real \(s=\Re z\)
is bounded uniformly in the profiles, after the factors
\((ja_i)^{d_i s}\) have been removed. Indeed the logarithmic
error is \(O((1+\lvert s\rvert)^2 V_{\rm sh}/j)\), by
Taylor expansion of the log Gamma functions, since
\(\sum |d_i|/a_i\le V_{\rm sh}\) and
\(\sum d_i^2/a_i\le V_{\rm sh}\). At imaginary height \(v\)
the absolute shared numerators are bounded by their real
values. The reciprocal falling Gamma factor costs at most
\(e^{\pi|v|/2}\) by the Euler product comparison in the
cited proof. The remaining fixed private and Barnes Gamma
factors give an exponential rate at most
\[
        -\pi P|v|/2+(\pi-\theta_j)v.
\]
This gives integrable exponential tails on compact height
intervals of \(0<\theta_j<2\pi\). The residues at zero,
one and two and the shifted integral give (4), uniformly
locally in \(\xi\). In particular the portion of the left
integral with \(|\Im q|\ge Y\) is bounded by
\(e^{-c_0 jY}\) times at most a power of \(j\), for
large fixed \(Y\) on such a field disk.

For \(P=1\) the same pinched left integral satisfies
\[
        {\cal L}_j(\xi)=(1-e^\xi)^{-1}+O(j^{-1})            \tag{5}
\]
at any fixed limiting height \(\theta>\theta_{\sigma_0}\)
when \(\sigma\to\sigma_0<1\), and at every positive interior
height when \(\sigma_0\ge1\). It also holds at the negative
real parameter boundary \(\theta=\pi\) in all cases.
Here is the uniformity check on the Euler-product estimate.
Relative to the unshifted real arguments \(j\) and \(ja_i\),
replacing the lower arguments of the Euler integrals by
the real Gamma arguments on \(\Re z=-c\) and comparing a
sum with its integral costs at most
\[
 C\log(1+C(v/j)^2)
   +C\sum_{\rm shared}(1+|d_i|)
       \log\left(1+C\bigl(|d_i v|/(ja_i)\bigr)^2\right)
       =O((1+V_{\rm sh})|v|/j),
\]
using \(\log(1+x^2)\le2|x|\). The real Gamma ratios cost
\(O((1+V_{\rm sh})/j)\) as in (4). For the main Euler
integrals the derivative of their difference, in \(y=|v|/j\),
is the negative of
\[
       \sum_{\rm shared}|d_i|\arctan(|d_i|y/a_i)-\arctan y.
\]
By subadditivity of \(\arctan\), this is bounded below by
\(\arctan(\sigma y)-\arctan y\) before taking the negative.
Thus the growth of the ratio is at most
\(\exp((\theta_\sigma+O(1/j))|v|)\), with
\(\theta_s=0\) for \(s\ge1\).
For \(\sigma_0>0\), \(\theta_\sigma\to\theta_{\sigma_0}\)
by the uniform Lipschitz bound for \(\arctan(\sigma y)\)
near \(\sigma_0\); at zero the universal bound \(\pi/2\)
suffices. The original Barnes factors have exponential
rates \(-\theta_j v\) for \(v>0\) and
\(-(2\pi-\theta_j)|v|\) for \(v<0\).
This proves fixed exponential decay on the specified
disks and the same tail conclusion as for (4).
For \(|v|\le C_1\log j\), the exact Stirling ratios differ
from one by \(O((1+|v|)^2/j)\) uniformly, the estimate
following directly from Taylor expansion at the same real
arguments. The integrable exponential majorant and a
sufficiently large \(C_1\) prove (5) by the
[Mellin inversion](MULTIVARIATE_ONE_PRIVATE_SIDE_RICE_EVENTUAL.md#one-private-linear-variable).

## The actual right saddle on compact field disks

In what follows the interior height is one covered by (4)
or (5). The domain for the right contour joins the two
half-planes through the limiting pole-free right gap
\((1,p_R)\), with \(p_R=+\infty\) if there are no retained
negative rows. The exact two phase charts agree there in
real part and in derivative. Every fixed compact
neighborhood of a point in the gap avoids finite-profile
poles for sufficiently late profiles. On such compact sets
and in the open half-planes the exact Rice integrand, after
division by the pure \(a\) moment, is
\[
 j^{(p+1)/2}A(q)e^{j\Phi(q;\xi)}(1+O(j^{-1})),             \tag{6}
\]
with the appropriate lower version of \(\Phi\). Here \(A\)
is order-independent, analytic and bounded away from zero
and infinity on each selected compact chart. To check the
power and amplitude, Stirling for \(j!\) and the two Gamma
factors of the rational Rice ratio leaves \(j^{1/2}\),
including the \(z=jq\) Jacobian. Each private Gamma factor
contributes \(j^{1/2}\sqrt{2\pi e_i q}\). The shared
amplitude is \(\exp(\frac12\sum_{\rm shared}\log(x_i(q)/a_i))\),
and the remaining denominator gives
\(\sqrt q\sqrt{q-1}\), with the branches of the Rice formula.
The normalization cancels all powers involving \(\Delta\);
the remaining exponential is (3).
The sum of shared Stirling errors is uniformly \(O(j^{-1})\)
on these charts, as in the
[uniform Rice estimates](MULTIVARIATE_RICE_CONTOUR_AND_REAL_BOUNDARY.md#exact-contour-and-its-phases).
Across the right gap the two rational Gamma arguments
are positive on the real diameter. The two phase arguments
differ there by an integer multiple of \(2\pi i\), which
gives the same integer-order kernel.

Choose a fixed \(c_R\) in the limiting right gap. One may
start with the finite Rice rectangle with sides at
\(-c/j,c_R\) and sufficiently large fixed heights \(\pm Y\).
Its horizontal panels have strictly smaller phase than the
saddle and the boundary, uniformly after passage to the
limit. Here is the upper height check, uniformly in bounded
\(\Re q\). Each term in \(\Im B_*(q)\) is nonnegative.
The private and rational terms give the lower bound
\[
        (P-1)\arg q+\arg(1-q)+\theta
                   \ \longrightarrow\ (P-2)\pi/2+\theta.
\]
This is positive for \(P\ge2\). For \(P=1\), a retained
shared row contributes at least \(\pi/2-o(1)\);
a positive diffuse mass instead gives \(\tau\Im q\).
In the flat case \(P=1,\sigma_0=0\), the height under
consideration exceeds \(\pi/2\), so the bound is again
positive. The lower chart has uniformly negative vertical
velocity. Integration of these bounds from fixed large
height gives arbitrarily low horizontal panels on the
bounded real range. By (6) the loss on a fixed such panel
persists for the finite profiles. The omitted tails of
the left vertical integral, from (4) or (5), can be
made smaller than either desired exponential scale by
increasing \(Y\) further. Thus it remains to identify
the saddle coefficient of the finite right vertical arc.

The limiting descent flow and compact deformation of this
arc are justified as follows. At smooth upper real cut
points outside the right gap the imaginary derivative of
\(\Phi_*\) is strictly positive. A backward trajectory
from an approaching point then crosses that real interval,
which excludes forward accumulation there. At \(1\) and
at a finite \(p_R\), the
[real branch-point signs](MULTIVARIATE_RICE_CONTOUR_AND_REAL_BOUNDARY.md#where-the-upper-descending-contours-can-meet-the-real-axis)
exclude a constant-imaginary-phase trajectory approaching
the point. At zero with \(P>1\), and at other retained real
poles, the logarithmic expansion is of the form
\[
   \Im(\Phi_*(b+x+iy)-\Phi_*(b))
       =c\,y\log|x+iy|+(B_0+c\arg(x+iy))x+O(y+|x+iy|^2),
\]
where \(c\ne0\) and both real boundary imaginary derivatives
are positive. On the level that could approach \(b\), the
formula forces the real side of sign \(c\) and
\(|x|\ge |x+iy|/2\) sufficiently close to \(b\).
The real descent velocity then has sign \(c\), pointing
away from \(b\). Entry through a fixed small semicircle
cannot lead to further approach, since \(|x|\) increases
as long as the trajectory stays in that neighborhood.
At \(P=1,q=0\) the phase is regular and its upper vertical
velocity is \(\theta>0\). The limiting diffuse term and
analytic retained tail are real on each local real chart
and leave the sign checks unchanged. The lower flow
cannot return upward after crossing the gap.

A trajectory above a fixed stopping level cannot escape
to infinity. Below the axis the vertical speed has a
positive absolute lower bound. If both signs of retained
gaps occur, one of each contributes at least
\(\pi-o(1)\) to the sum of upper row arguments at large
\(|q|\), by their two supplementary limiting angles.
All other row contributions are nonnegative, so the
vertical speed is bounded positively away from zero.
If all retained gaps are positive, the real part of the
trajectory is bounded above, since for large positive
\(\Re q\) the horizontal descent points left. Indeed
all centered row real logarithms there are nonnegative,
the private term grows, the diffuse term has the same
sign, and the real rational logarithm is bounded.
When the real part is bounded above and \(|q|\) is large,
each fixed positive row has argument at least
\(\pi/2-o(1)\). Two units of total retained gap therefore
suffice for strictly positive upper vertical speed.
If there is only the private linear row and \(\tau>0\),
at large height the diffuse term gives the speed, and at
bounded height and far left the private row has argument
tending to \(\pi\). In the completely flat \(P=1\) case,
\(\theta>\pi/2\) already suffices.
These bounds at infinity rule out escape above a fixed
level by the finite descent energy identity. On the
portion outside a large disk, speed bounded away from
zero makes the total time finite, and Cauchy–Schwarz
makes the total length finite. Bounded trajectories
staying above the level must have a stationary accumulation
point, since the finite cut limits have been excluded.

Apply this to the finite right arc with endpoints below
a strict sublevel of the single saddle. In a local analytic
quadratic coordinate take a rectangle whose exit faces
lie below that level. This rectangle together with the
strict sublevel is forward invariant. Every point of the
arc enters this set in finite time, as otherwise it would
have a stationary accumulation point outside the set.
Compactness gives one finite flow time for the arc.
The deformed arc is therefore an integer multiple of a
descending saddle segment plus paths in the strict
sublevel. Indeed in the rectangle the sublevel components
are the two descending sides, and path pieces with
endpoints on them can be joined to the ends of the
descending segment within these components. A transverse
perturbation gives a finite crossing decomposition, with
all other joins at strictly lower level. This is the
[finite-contour argument](MULTIVARIATE_UNIFORM_SUPERCRITICAL_SUPPORTS.md#finite-contour-deformation-for-the-limiting-profile)
applied to an arc with both endpoints in the sublevel.
The homotopy is compact in the gap domain. Thus it
transports to the exact finite-profile contour. Its
saddle segment can be moved slightly to the actual
nearby saddle. The finitely many sublevel pieces retain
a strict exponential loss.

The integer coefficient is locally constant in the field.
One may use the same compact homotopy on a nearby field
disk. Different sufficiently high truncations of the
initial right vertical arc differ by strict-sublevel
segments, since the vertical loss on that fixed line holds
uniformly beyond large fixed height. Different gap lines
can be joined across the gap with high horizontal panels.
The expansions agree up to strict exponential loss, so
the integer coefficient is independent of these choices
by the nonzero oriented Gaussian. At the negative real
parameter the gap saddle is real and simple. The real part
of \(X_*\) runs strictly from \(-\infty\) to \(+\infty\)
on the right gap and its derivative there is positive.
At the outer endpoint this follows from the nearest
retained negative pole if finite, and otherwise from
the private positive logarithm or the diffuse term.
The vertical line through the gap saddle has strictly
decreasing phase upward and downward away from it:
at \(\theta=\pi\),
\(\Im \Phi'_{*,+}>0\) above and \(\Im\Phi'_{*,-}<0\)
below. The local gap Gaussian therefore has coefficient
one up to orientation. The same contour works on a
small complex disk through the negative real parameter.
Continuation from this anchor through the connected
height region of (4) or (5) shows that the right coefficient
is nonzero throughout that region.

We have proved the local exact two-term formula
\[
 \frac{M_j(t)}{M_j(0)}=
     {\cal L}_j(\xi)+j^{p/2}J(\xi)e^{j {\cal S}(\xi)}
                                (1+O(j^{-1})),            \tag{7}
\]
where \({\cal S}(\xi)=\Phi(q(\xi);\xi)\) at the actual
finite-profile saddle, and \(J\) is the locally nonzero
order-independent Gaussian amplitude. The exponentially
small outer errors have been absorbed into the right
relative error by the choice of truncation and sublevel.
The field and profile dependence in (7) is analytic
locally with bounded derivatives. The same formula holds
on a small disk about each negative real parameter.

## The positive real boundary and the adjacent comparison

On \(0<q<1\) the limiting real phase from (3) has
\[
 (1-q)\Phi_*''(q;\xi)=R(q)-1,\qquad
 R(q)=(1-q)\bigl((P-1)/q+B_*'(q)\bigr).
\]
Here \(R'(q)=-(P-1)/q^2-\sum d_i^2 b_i/x_i(q)^2-\tau\),
with summation over retained shared rows. For \(P\ge2\),
\(R\) decreases from infinity to zero. The derivative of
the phase starts and ends negative infinite, and there
can be just one interior maximum. At its first appearance
it is below the boundary. The interior saddle value then
strictly increases with the field, with derivative \(q>0\),
and eventually exceeds zero. Thus the boundary and
interior maximum coexist at exactly one real field.
For \(P=1,\sigma_0>1\), the same argument starts at
\(R(0)=\sigma_0>1\). The coexistence is at \(\Re\xi<0\),
since at \(\xi=0\) the phase increases from zero.
For \(P=1,\sigma_0\le1\) the phase is strictly concave
inside \((0,1)\), with the sole boundary-to-interior
transition at \(\xi=0\).

The real analysis gives uniform local expansions on small
complex disks at the real stable points. For a stable
boundary, sum the terms on \(k\le\eta j\) with fixed small
\(\eta>0\). After dividing by the pure \(a\) moment,
the shared normalized factorial ratios differ in their
logarithms from their leading linear term by
\(O((k+k^2)/j)\), by the same Taylor bounds as on the
left contour. For \(P\ge2\), use
\((e_i k)!/e_i^{e_i k}\le(k!)^{e_i}\) and
\(\binom jk\le j^k/k!\) to bound the terms by
\[
                       C^k(k!)^{P-1}/j^{(P-1)k}.
\]
Successive ratios of this envelope are below \(1/2\)
for \(k\le\eta j\) after choosing \(\eta\) small.
The sum is \(1+O(j^{-1})\). For \(P=1\) the stable
boundary requires a fixed \(\Re\xi<0\). The private
factor cancels \(k!\) in the binomial coefficient,
so the terms have geometric envelope with a fixed
ratio smaller than one for \(\eta\) small. The coefficient
at \(k\) is \(e^{k\xi}(1+O((k+k^2)/j))\) for
\(k=O(\log j)\), by the falling factorial estimate
and the shared ratios. The sum is
\((1-e^\xi)^{-1}+O(j^{-1})\).
In either case the endpoint coefficient is uniformly
nonzero.

For an interior stable contribution take a small interval
around its limiting nondegenerate maximum. Use the
[cotangent sum deformation](MULTIVARIATE_VARYING_PROFILE_SUBCRITICAL_SADDLE.md#the-lattice-sum-and-its-contour)
with half-grid endpoints. The summand interpolation from
Stirling carries power \(j^{(p-1)/2}\) relative to the
pure \(a\) term, so the sum and local Gaussian give total
power \(j^{p/2}\). The interval endpoints and connecting
contours retain a strict phase loss on a small complex
disk. This gives an actual saddle expansion with nonzero
analytic leading amplitude and uniform \(O(j^{-1})\)
relative error. Outside the retained interval and endpoint
neighborhood the limiting real phase has a fixed strict
gap. Uniform convergence and an upper bound for the
summand by \(Cj^{p/2}e^{j\Re\Phi(q;\xi)}\) make the
discarded sum exponentially small. This upper bound
follows from \(\binom jk\le e^{j(-q\log q-(1-q)\log(1-q))}\),
the shared Stirling ratios, and the private factorial
upper bound \(C\sqrt{1+e_i k}\,\exp(F(e_i k))\).
Consequently at regular real fields of strict dominance
there are no zeros in a neighborhood for large orders,
and at positive real coexistence the two-term formula
(7) holds with the just proved endpoint expansion.
For the continuous transition, only \(\xi=0\) remains
possible as a limiting real field.

It remains to exclude common roots in the indicated
regular two-term neighborhoods. At each positive interior
height where (4) or (5) holds, the upper saddle stays in
\(\Re q>0\). On \(q=iy\), \(y>0\),
\[
 \Im X_*(iy)=(P-1)\pi/2+
          \sum_{\rm retained}|d_i|\arctan(|d_i|y/a_i)
                  +\tau y-\arctan y.
\]
For \(P\ge2\) it is positive. For \(P=1\) the sum and
diffuse term together are at least \(\arctan(\sigma_0 y)\)
by subadditivity and the slope bound on \(\arctan\).
Thus there is no root of \(X_*=-\xi\) on the axis
at the interior heights of (5). The branch starts
to the right at the negative real gap and cannot cross.
The stationary phase satisfies
\(\partial_\xi {\cal S}_*=q_*\), so its real part
is strictly increasing with the real field on each
such height. At height \(\pi\) the same holds with
\(q_*>1\). There is therefore at most one limiting
equal-modulus point per height, and a zero in (7)
must limit to such a point. At height zero the
real analysis similarly gives at most one possible
limit. At a common root the two limiting fields
have the same height and hence agree. By (2),
\(\Delta/n\to0\).

At a regular two-term point, cancelling the terms in (7)
gives \(j\Re{\cal S}(\xi_j)=-(p/2)\log j+O(1)\).
Dividing the two zero equations gives a ratio whose
right side has modulus tending to one, because the
amplitudes are bounded, nonzero and their fields tend
to the same point. The real exponent on the left is
\[
 \Re\bigl((n+1){\cal S}(\xi_{n+1})-n{\cal S}(\xi_n)\bigr)
  =\Re {\cal S}(\xi_n)+\Delta\Re q(\xi_n)
          +O(|\Delta|/n+\Delta^2/n).
\]
Since \(\Re q\) has a positive limit and \(\Delta/n\to0\),
a nonzero integer \(\Delta\) gives a contradiction.
If \(\Delta=0\), the fields and leading amplitudes
coincide exactly. Division then gives
\(\Re {\cal S}(\xi_n)=O(n^{-1})\), incompatible with the
\(-(p/2)\log n\) term in the single-order equation.
This excludes the regular cases claimed above and completes
the uniform theorem. Coprimality transfers over every
characteristic-zero field by rational Bézout.
