# Eventual coprimality for fixed fully shared supercritical supports

Let \(a,b\in\mathbb Z_{>0}^r\) be distinct exponent columns,
with equal rows removed. Write \(d_i=b_i-a_i\),
\(\Delta=\sum_i d_i\), \(x_i(q)=a_i+d_iq\), and
\[
 \rho(x)=x(1-x)\sum_i\frac{d_i^2}{x_i(x)}\qquad(0<x<1).
\]
The term supercritical means \(\max\rho>1\).

**Theorem.** For every fixed fully shared supercritical
ordinary exponent profile, the factorial moment polynomials
\(M_n(t),M_{n+1}(t)\) are coprime for all sufficiently large
\(n\). There is no restriction on the two total degrees.
The threshold depends on the profile and is unevaluated.
The result transfers from \(\mathbb Q\) to every
characteristic-zero field.

The main additional ingredient is a contour deformation
in each complex region of unequal saddle moduli. It
isolates the larger saddle with a fixed nonzero Gaussian
coefficient, including when its descending trajectory
passes through the smaller saddle. This gives zero-freeness
there. Where the moduli agree, the two actual contributions
have different consecutive multipliers.

## The two regions of strict dominance

Use the upper Rice phase \(\phi_+\), the lower phase
\(\phi_-\), and \(\zeta=n^\Delta t\) from the
[exact contour](MULTIVARIATE_RICE_CONTOUR_AND_REAL_BOUNDARY.md#exact-contour-and-its-phases).
For \(\Im\zeta>0\), set \(w=-\operatorname{Log}\zeta\),
so that \(-\pi<\Im w<0\). The stationary equation is
\(X(q)=w\), a proper two-sheet map with no branch point
in this case. Choose \(x_0\in(0,1)\) with \(\rho(x_0)>1\).
The [vertical divider argument](MULTIVARIATE_REGULAR_SADDLE_PHASE_SEPARATION.md#equal-total-degrees-at-and-above-the-real-cusp)
uses only this inequality and positive integer exponents,
independently of \(\Delta\). It gives \(\Im X(x_0+iy)>0\)
for \(y>0\). The two real inverse branches at the boundary
\(\Im w=-\pi\) lift into opposite sides of the divider.
Indeed \(X'>0\) on the two real gaps \((p_L,0),(1,p_R)\),
where \(p_L,p_R\) denote the nearest exterior row poles,
or the infinite endpoints when no such poles exist.
There is one inverse in each gap for each target real
part on that boundary. As its real coordinate follows the
gap left to right its image real part increases from
\(-\infty\) to \(+\infty\). At an infinite endpoint this
follows from the logarithmic term \(\Delta\log|q|\)
and the fact that all gaps have the same sign. Thus label
the two interior inverse branches \(q_L,q_R\) with
\(\Re q_L<x_0<\Re q_R\).

Let \(D(w)=\phi_+(q_R;\zeta)-\phi_+(q_L;\zeta)\), where
the common additive constant in the stationary phase
cancels. Differentiating the stationary values gives
\(D'=q_L-q_R\), because \(\phi_+\) depends on \(w\)
by the term \(-wq\). In particular \(\Re D'<0\)
throughout the strip.

For a fixed height \(-\pi<\eta<0\), as \(\Re w\to+\infty\)
the two solutions tend to \(0,p_R\) when \(p_R\) is finite.
These solutions follow directly from the terms \(-\log q\)
and \(-D_R\log(p_R-q)-i\pi\), respectively, where
\(D_R=\sum_{i:p_i=p_R}|d_i|\). The inverse near each
endpoint is upper for the stated height, and the two
solutions exhaust the strip degree. If \(p_R=+\infty\),
all gaps are positive, and the right branch instead
tends to infinity with argument \((\eta+\pi)/\Delta+o(1)\)
by the leading term \(\Delta\log q-i\pi+\mathrm{const}\).
Here \(\Delta\ge2\), since a single increasing unit-gap
row has \(\rho(x)\le1/4\). Thus its real part diverges
positively. The stationary phase up to a common constant is
\(F_a(q)-\Delta q\) as in the contour formula. It remains
bounded at \(0\), tends to negative real infinity at a
finite \(p_R\) since some \(x_i\) vanishes, and also tends to
negative real infinity along the infinite right branch since
\(-\Delta\Re q\) dominates the logarithmic terms. Therefore
\(\Re D\to-\infty\) at the right end of each horizontal line.
At the left end \(\Re D\to+\infty\) by reversal of the columns
and \(q\mapsto1-\overline q\). Under this reversal the roles
of the two saddles interchange and the real phase difference
changes sign. In fact \(F_a^{\,\mathrm{rev}}(1-\overline q)
=\overline{F_b(q)}\) and \(F_b=F_a+X\), so all extra terms
are common to the two branches.

The equal-modulus locus \(E\) is consequently a single
smooth graph over \(\eta\). Its two complements in the
strip are connected, with a fixed larger saddle on each.
On \(E\) the Cauchy–Riemann equations give
\[
 \frac{d}{d\eta}\Im D=\frac{|D'|^2}{\Re D'}<0 .
\]
At the lower boundary the unique equal-modulus point has
\(D=i\pi\): on the two real gaps the stationary phases
have imaginary parts zero and \(\pi\) up to a common
constant. The real phase difference is strictly decreasing
from positive to negative infinity there by the same
endpoint limits. The graph extends to this point by
local inversion and the nonzero real derivative.
At the upper boundary the two real stable maxima at the
[positive coexistence scale](MULTIVARIATE_SHARED_SUPPORT_CONDITIONAL_SADDLE.md#a-local-exclusion-around-the-positive-real-axis)
have equal stationary value and straddle the divider.
Their derivatives \(X'\) are negative, so they are the
limits of the upper branches. The graph extends to their
equal-modulus point with \(D=0\), again by local inversion
and the nonzero real derivative. In particular \(0<\Im D<\pi\)
on the open graph.

## Descending paths and Stirling bounds at infinity

By column reversal it is enough to treat \(\Delta\ge0\).
Work on the domain consisting of both open half-planes joined
through the two real gaps. The real phase \(h=\Re\phi\) and
its gradient agree through the gaps as in the
[equal-degree contour proof](MULTIVARIATE_RICE_CONTOUR_AND_REAL_BOUNDARY.md#actual-two-saddle-expansion-on-the-equal-degree-graph).
The deformation uses \(\dot q=-\overline{\phi'(q)}\).
For fixed interior \(w\), descending trajectories cannot
accumulate at a finite real cut. On a non-gap interval the
upper vector points strictly into the upper half-plane,
and the lower vector strictly into the lower half-plane.
At the real branch points use the
[boundary exclusion](MULTIVARIATE_RICE_CONTOUR_AND_REAL_BOUNDARY.md#where-the-upper-descending-contours-can-meet-the-real-axis).
For a farther pole, its possible upper level curves approach
tangentially from the permitted side. Upon entering a small
half-disk a trajectory on that level has horizontal distance
at least half the radius from the pole, and the horizontal
component then points away. This excludes accumulation even
along a subsequence. A crossing through a gap enters the
lower half-plane and never returns, by the strictly negative
vertical velocity of the lower phase.

At infinity there are constants \(C_\pm\) with
\[
 \phi_\pm'(q)=\Delta\log q+C_\pm+O(|q|^{-1})
\]
in the respective open half-planes, with
\[
 \Im C_+=\theta+\pi(D_--1),\qquad
 \Im C_-=\theta-\pi(D_-+1),
 \qquad \theta=\arg\zeta,\qquad
 D_\pm=\sum_{\pm d_i>0}|d_i|.
\]
The expansion follows from factoring out \(q\) in \(X\)
with the principal logarithms. Its error is uniform up to
the real directions in either half-plane. After integration
\(\phi_\pm(q)=\Delta(q\log q-q)+C_\pm q+O(\log(1+|q|))\).
On a fixed vertical line in each half-plane the real phase
thus has strict linear loss at infinity. The rates per
unit height on the upper and lower sides are
\(\theta+\pi((D_++D_-)/2-1)\) and
\(\pi((D_++D_-)/2+1)-\theta\), respectively. They are
positive, including at \(\theta=\pi\), since the profile
is supercritical and \(D_++D_-\ge2\).

When \(\Delta=0\), the upper and lower vertical velocities at
infinity point outward and are bounded away from zero.
The horizontal velocity is bounded. Thus any unbounded
trajectory must eventually move to arbitrarily large
positive or negative height. For \(\Delta>0\), the horizontal
velocity is strictly negative outside a large circle, with
magnitude comparable to \(\log|q|\), and the vertical velocity
is bounded there. The upper vertical velocity is bounded away
from zero positively when \(\Re q\le0\) and \(|q|\) is large:
its leading lower bound is
\(\theta+\pi(D_--1)+\Delta\pi/2>0\).
If \(D_-=0\) this uses \(\Delta\ge2\). The lower vertical
velocity has a negative bound when \(\Re q\le0\) at infinity.
After a trajectory reaches a sufficiently negative real
coordinate it stays in this region, in the same half-plane,
and its height diverges. An unbounded forward trajectory has
its real coordinate bounded above, since the horizontal
velocity outside the large circle is negative. If it reaches
arbitrarily large height before entering the far left region,
bounded vertical speed and the strictly negative horizontal
speed take it into that region before its height can return
to the large circle. Hence every unbounded trajectory
eventually escapes via one of these two stable regions.
The endpoint type, upper or lower infinity, is open in the
parameter whenever the trajectory has reached such a region.

For \(\Delta>0\), along an escaping trajectory one has
\(|\Re q|=O(|\Im q|\log(2+|\Im q|))\), uniformly on compact
parameter paths for which the stable regions are reached
in bounded time. Indeed outside the left threshold the
horizontal speed is comparable to \(\log|q|\) and the
outward vertical speed is bounded above and below.
First \(|\Im q|=O(1+|\Re q|)\) by these speed bounds.
Integrating \(d|\Im q|/d|\Re q|\ge c/\log(2+|\Re q|)\)
then gives the claimed growth. For \(\Delta=0\) the
linear bound of the equal-degree proof applies.

For quantitative bounds on a Rice contour at large heights,
the relative Stirling expansion remains valid uniformly on
the required horizontal connectors, even when the arguments
approach a negative real direction. For a row with
\(\Re z\le0\), where \(z=n x_i(q)\), apply
\(\Gamma(1+z)=-\pi/(\Gamma(-z)\sin\pi z)\) and use the
right-half-plane Stirling formula for \(-z\). The sine on
either open half-plane has a single leading exponential
with relative error \(O(e^{-2\pi|\Im z|})\) at large
\(|\Im z|\). Combining the two expansions gives the same
phase and amplitude as in the fixed-sector formula, with
the branches identified on an overlapping sector.
The \(q\)-ratio is reflected as a pair when \(\Re q\le0\),
using the exact gap reflection identity also off the real
line. For \(0\le\Re q\le1\), at large height both factors
are already in fixed sectors. For \(\Re q\ge1\) they are
in the right half-plane. Thus on a large-height connector
the error and all amplitudes are polynomially bounded in
\(|q|,n\) after the common scalar normalization of the
Rice formula. This suffices for exponential tail bounds.
[Gamma reflection and sectorial Stirling](https://dlmf.nist.gov/5.5.E3)

When \(\Delta>0\), for a connector to a vertical line of
bounded real part one has \(\partial_{\Re q}h>0\) at
large height. The rightmost endpoint at the vertical line
therefore gives a uniform upper bound on \(h\), with
strict linear loss. When \(\Delta=0\), the phase is linear
up to the logarithmic error. Its loss along the connector
is bounded by the larger of the two endpoint values of
that linear phase. At the descending endpoint the velocity
tends to \(-\overline{C_\pm}\), which gives strict linear
loss in the height, and the vertical line has the same
property. Together with the growth bounds above, this
proves that connector integrals vanish as the height tends
to infinity. It also proves \(h\to-\infty\) on each
escaping descent end.

## The two actual saddle periods at equal modulus

At an interior point of \(E\), both descending curves from
each simple saddle escape to infinity. To verify this,
\(h\) is strictly decreasing as soon as the curve leaves
the saddle, so neither candidate stationary point can be
a limit. A bounded trajectory away from the cuts would
have a stationary limit point, since boundedness of \(h\)
from below makes \(\int|\phi'|^2\,dt\) finite. The finite
boundary exclusions and the infinity alternatives above
leave only escape. No curve can leave the domain at a finite
time by these same boundary exclusions; at infinity the
velocity is at most logarithmic in \(|q|\).

At the negative-axis endpoint each saddle is real in its
respective gap. The Hessian \(X'\) is positive there and
one descending direction goes up, the other down. At
\(\theta=\pi\) the upper velocity has strictly positive
vertical component on the entire upper half-plane, since
\(\Im X>-\pi\). The lower velocity is strictly negative
vertically there. Neither trajectory can accumulate at
another stationary point because there is no stationary
point in the open half-planes at this parameter. The
infinity alternatives give one upper and one lower end.
The two paths, parametrized by their heights, are homologous
to the original two vertical lines by horizontal interpolation
that crosses the real line in the respective gaps.
The Rice formula is valid at this real negative scale with
\(\operatorname{Log}(-\zeta)\) real.

The two end types and the orientation of each path persist
locally along \(E\) by finite-time continuity of the flow and
the outward regions at infinity. The graph is connected.
Along each compact parameter subpath from the negative end,
truncate at large positive and negative heights and retain
the horizontal connectors to its original vertical line.
Continuous dependence and eventual monotonicity in height
make these truncated paths vary continuously, without
crossing any exact integrand pole, since all poles lie
outside the domain joining the half-planes through the
two gaps. The connector bounds show that the limiting
integrals are the actual Rice periods, with their original
relative orientations.

Near each interior point of \(E\), the two saddle phases
have distinct imaginary parts. Therefore no descending
upper trajectory can connect them at a nearby parameter:
the imaginary phase is constant along such a trajectory.
A trajectory that has crossed into the lower half-plane
cannot return. The same contour homotopies and end types
hold on a sufficiently small complex parameter disk.
The local quadratic coordinate at each simple saddle
gives a nonzero analytic Gaussian amplitude. Off a fixed
saddle neighborhood the descending contour has strict
real-phase loss, and at infinity it has the tail bounds
already proved. These provide an actual two-saddle expansion
with uniform relative error tending to zero separately
for each contribution.

When the two orders use the same \(t\), one has
\(w'=w-\Delta\log(1+1/n)\) at order \(n+1\).
The ratio of the phase multipliers between the orders
tends to \(\exp(D-\Delta D')\), uniformly on any smaller
regular parameter disk, by Taylor expansion in \(w\).
If \(\Delta\ne0\), its log modulus on \(E\) is
\(-\Delta\Re D'\ne0\). If \(\Delta=0\), the imaginary
part of \(D\) lies strictly between zero and \(\pi\)
there. Thus the multiplier ratio differs from one
throughout \(E\). At a zero on a small disk the two
saddle contributions must cancel. Dividing these relations
for consecutive zeros, the nonvanishing analytic amplitudes
at \(w,w'\) contribute a ratio tending to one. This
contradicts the fixed nonunit multiplier ratio.
Hence each interior point of \(E\) has an eventual local
common-root exclusion.

## Isolation of the larger saddle through transitions

Fix a parameter in either component of the complement of \(E\).
Denote the larger saddle by \(p\), with height \(h_p\),
and the smaller saddle by \(r\), with height \(h_r<h_p\).
Use a local analytic coordinate \(u\) with
\(\phi_+(q)=\phi_+(p)-u^2\).
Take a small rectangle \(R=\{|\Re u|<2b,\ |\Im u|<b\}\)
contained in this coordinate disk, with \(h_r<h_p-2b^2\).
Put \(c=h_p-2b^2\) and \(S_c=\{q:h(q)<c\}\).
The flow in the coordinate is
\(\dot u=2|u'(q)|^2\overline u\). It enters through the
horizontal sides of the rectangle and can exit only
through the vertical sides, on which \(h\le h_p-3b^2<c\).
Since \(S_c\) itself is invariant under descent,
\(S_c\cup R\) is forward invariant.

Start with the compact exact Rice loop from (1), for example
a finite rectangle extending in both half-planes through
the two real gaps. Every point of the loop enters \(S_c\cup R\)
in finite time. Otherwise its height stays at least \(c\),
so it cannot escape to infinity by the escape and tail
estimates. It cannot accumulate at the real cut by the
boundary and gap analysis. The remaining compactness
argument for the descent then gives a stationary
accumulation point. It is not \(r\), whose height is
less than \(c\), and it is not \(p\), which is inside \(R\).
This is a contradiction. The entry set \(S_c\cup R\)
is open and invariant, so continuity and compactness of
the original loop give a common finite flow time at which
the entire loop lies in this set. Flowing the loop in the
domain gives an exact contour homotopy for the Gamma
integrand for each \(n\). This is the only use of the
global descent in this step. It does not require a
descending trajectory from \(p\) to avoid \(r\).

Inside the \(u\) rectangle the sublevel intersection consists
of two connected pieces, given by \(x^2-y^2>2b^2\)
with \(x>0\) and \(x<0\), respectively, where \(u=x+iy\).
Each is convex, from \(x>\sqrt{2b^2+y^2}\) or its reflection.
Subdivide the deformed loop into finitely many pieces
lying in \(S_c\) or in \(R\), by compactness and the open
two-set cover. A piece through \(R\) that joins endpoints
in the same sublevel component can be replaced by a path
inside that component. A piece joining opposite sides is
homologous, with endpoints kept fixed, to sublevel paths
and one copy of a fixed real segment \([-b_1,b_1]\)
in the \(u\) coordinate, where \(\sqrt2 b<b_1<2b\).
The segment is oriented by the crossing direction.
These homotopies take place inside the rectangle, which
is simply connected. If the whole loop is in the rectangle,
its integral vanishes and no segment is needed.
After collecting orientations the exact contour consists
of an integer number \(\nu\) of copies of this saddle
segment and a finite chain in \(S_c\).

At a nearby parameter the same chain and slightly
perturbed saddle segment give a representative with a
strict phase gap on the chain outside the segment.
The analytic saddle and its quadratic coordinate vary
analytically on a smaller parameter disk. Move the
segment and join its endpoints by sufficiently short
paths in the original sublevel parts of the rectangle.
Every other part lies in a fixed compact subset of a
strict sublevel at the original parameter, so the
gap persists uniformly after shrinking the disk.
The exact integral, after the common scalar
normalization, therefore has the expansion
\[
 e^{n\phi_+(p;\zeta)} n^{-1/2}
       \bigl(\nu\,A_p(\zeta)+o(1)\bigr)
        +O\!\bigl(e^{n(h_p-\epsilon)}\bigr)
\]
locally uniformly, where \(A_p\) is analytic and nonzero
and \(\epsilon>0\). In this formula the heights are
evaluated at the current parameter. The amplitude is the
usual Gaussian factor times the analytic nonzero amplitude
of the Rice integrand at \(p\).

For a fixed orientation of the saddle coordinate the integer
\(\nu\) is unique: two different integers would give different
leading asymptotics for the same exact integral.
In particular the property \(\nu\ne0\) is locally constant
on the region of strict dominance, even across points where a
descending path from the larger saddle connects to the smaller
one. In each of the two connected regions near the graph the
already established two-saddle expansion has a nonzero
coefficient at the larger saddle. The larger exponential
dominates there at any fixed point off \(E\).
Thus \(\nu\ne0\) in the entire region on each side.
This proves eventual local zero-freeness of the individual
polynomials away from \(E\) in the open upper half-plane
of the scaled parameter.

## From the local conclusions to a fixed-profile terminal

For fixed fully shared profiles all scaled roots lie in a
[fixed annulus](MULTIVARIATE_SHARED_SUPPORT_CONDITIONAL_SADDLE.md#a-fixed-annulus-for-possible-common-roots)
at large orders. Every nonzero real scale has an eventual
local consecutive exclusion by the real-scale analysis for
fixed fully shared supports, since the strict inequality
\(\max\rho>1\) omits the real cusp. The positive real
coexistence point is included there. Every scale in the open
upper half-plane has the exclusion just proved, either by
a single actual dominant saddle or near the graph by the
consecutive phase comparison. Complex conjugation gives
the same statement below. If common roots existed at
arbitrarily large terminals, some sequence of their scaled
values would converge in the fixed annulus, contradicting the
local exclusion at its limit. This proves eventual
coprimality. The polynomials have integer coefficients, so
the resulting rational Bézout identities extend to every
characteristic-zero field.
