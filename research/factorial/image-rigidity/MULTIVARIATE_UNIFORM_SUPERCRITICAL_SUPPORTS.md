# Uniform eventual coprimality on bounded supercritical support bands

Let \(a,b\) be distinct ordinary nonnegative exponent columns.
Remove equal rows and assume both entries of every remaining
row are positive. Put
\[
 d_i=b_i-a_i,\quad m_i=(a_i+b_i)/2,\quad
 r_i=d_i/(2m_i),\quad W=\sum_i d_i^2/m_i,\quad
 \rho(q)=q(1-q)\sum_i\frac{d_i^2}{a_i+d_iq}.
\]
The dimension and the exponent columns may vary with the
terminal. Within each consecutive pair the profile is fixed.

**Theorem.** Fix \(W_0<\infty\) and \(\kappa>0\).
There is a terminal \(N_0(W_0,\kappa)\) such that for every
fully shared ordinary profile with
\[
 W\le W_0,\qquad \max_{0<q<1}\rho(q)\ge1+\kappa,
\]
the moments at orders \(N-1,N\) are coprime whenever
\(N\ge N_0(W_0,\kappa)\). This holds over every
characteristic-zero field. In particular it covers the
entire band \(4+4\kappa\le W\le W_0\), since \(\rho(1/2)=W/4\).
The threshold is unevaluated.

The proof allows the full retained and diffuse closure of the
bounded-curvature profiles. In that closure exact finite Rice
contours can still be deformed using the limiting phase.
One then uses the actual finite-profile phases on these
contours. Convergence of the phases requires no rate in the
terminal order.

## Centering, compact scales, and uniform local expansions

For each profile let \(x_i(q)=a_i+d_iq\) and put
\[
 C_0=\sum_i(b_i\log b_i-a_i\log a_i-d_i),\qquad
 G(q)=\sum_i[x_i(q)\log x_i(q)]_{\rm c},
\]
where the subscript subtracts the chord through zero and one.
For order \(n\), use the external variable
\[
 w=-\operatorname{Log}(n^\Delta t)-C_0,\qquad
 \Delta=\sum_i d_i .
\]
For nonreal \(t\) in the upper half-plane take
\(-\pi<\Im w<0\). At a common root the same \(t\)
gives \(w'=w-\Delta\log(1+1/n)\) for order \(n+1\).

The endpoint center \(\mu_n\) of the
[uniform root annulus](MULTIVARIATE_UNIFORM_ROOT_ANNULUS.md#shared-rows-and-large-degree-imbalance)
satisfies
\[
 \mu_n=\Delta\log n+C_0+O(W/n)
\]
on the band. Indeed, use
\(\mu_n=\sum_i\int_{a_i}^{b_i}\psi(1+nx)\,dx\)
and \(0\le\psi(1+y)-\log y\le1/y\) for \(y>0\).
The latter follows by integrating the positive-real trigamma
bounds from infinity. The error is at most
\(n^{-1}\sum_i|\log(b_i/a_i)|=O(W/n)\).
The annulus bounds apply because
\(\sum d_i^2/\min(a_i,b_i)\le W/(1-R)\)
with \(R=\max_i|r_i|\) uniformly below one on the band.
Thus \(w,w'\) for possible consecutive common roots
have bounded real parts, and \(\Delta/n\) is bounded.

After a scalar normalization the upper
[Rice phase](MULTIVARIATE_RICE_CONTOUR_AND_REAL_BOUNDARY.md#exact-contour-and-its-phases)
is
\[
 f_+(q,w)=G(q)+(q-1)\log(1-q)-q\log q-wq .
 \tag{1}
\]
The lower phase is \(f_-(q,w)=f_+(q,w)-2\pi i q\),
with lower principal determinations in (1).
The real phases and their derivatives join across the
pole-free gaps to the left of zero and right of one.

The exact integrands, divided by their respective
exponentials in (1), have analytic expansions of the
form \(A(q)(1+B(q)/n+O(n^{-2}))\) on the compact
contours used below. Their constants, and all derivatives
on smaller compact neighborhoods, are uniform in the
profiles. The coefficients \(A,B\) depend on the profile
and its contour determination, with \(A\) bounded and
bounded away from zero. To see this, normalize the
individual Gamma factors by their values at \(a_i\)
and use sectorial Stirling through the \(1/(12n x_i)\)
term. On a fixed compact neighborhood avoiding the
limiting row cuts along the given contour, all small
\(|r_i|\) satisfy \(x_i(q)/m_i=1+O(|r_i|)\).
Row differences of the logarithmic amplitudes are
\(O(|d_i|/m_i)\). The Stirling remainders after
subtracting the normalization point have sum
\(O(W/n^2)\), by Cauchy's estimate in the row argument
and \(|d_i|\ge1\). The finitely many retained row
types away from \(r=0\) have uniform sector margins.
On real gaps the ratio of the binomial Gamma factors
uses the exact [reflection across the left gap](MULTIVARIATE_RICE_CONTOUR_AND_REAL_BOUNDARY.md#actual-two-saddle-expansion-on-the-equal-degree-graph),
or the direct positive-argument ratio on the right.
The resulting fixed-contour bounds persist for all
sufficiently near finite profiles in the closure.
[Stirling and reflection formulas](https://dlmf.nist.gov/5.11.E1)

## The limiting phase and its real-boundary geometry

Suppose a counterexample sequence has \(n\to\infty\),
and pass to a retained and diffuse subsequence as in the
[varying-profile phase estimates](MULTIVARIATE_VARYING_PROFILE_PHASE.md#retained-measures-and-the-diffuse-term).
The limiting phase \(G_*\) is the sum of the retained
row chord functions and \(\tau q(q-1)/2\), where \(\tau\)
is the diffuse mass. The retained row list can be
infinite, but \(\sum_i w_i<\infty\) and its row poles
have no finite accumulation. Indeed each fixed
\(|r_i|\ge\delta>0\) admits only finitely many retained
rows, including their multiplicities. Uniform convergence
of \(G\) and its derivatives holds on compact sets of
either open half-plane, across \([0,1]\) for \(G\) itself,
and across the limiting real pole-free gaps. This follows
either by integrating the limiting curvature transform
or directly by expanding the small-\(|r_i|\) chord
functions on each bounded \(q\)-set. The same expansion
extends the infinite tail analytically across any
bounded real interval.

The limit stationary map is
\[
 X_*(q)=G_*'(q)+\log(1-q)-\log q .
\]
The difference from midpoint centering in the
[stationary compactness theorem](MULTIVARIATE_UNIFORM_SADDLE_MAP_COMPACTNESS.md)
is just a real constant, namely the negative of the
integral of the centered row transform from zero to one.
This constant converges, by convergence on \([0,1]\).
Thus the limit is a proper two-sheet map over
\(\mathcal T=\{-\pi<\Im w<0\}\).
The real limiting curvature \(q(1-q)G_*''(q)\) is
strictly concave by the positive row and diffuse formulas,
and its maximum is at least \(1+\kappa\).
Choose \(q_0\in(0,1)\) where the curvature exceeds one.
The vertical-divider inequality of the
[fixed-profile proof](MULTIVARIATE_SHARED_SUPERCRITICAL_EVENTUAL_COPRIMALITY.md#the-two-regions-of-strict-dominance)
passes to the limit at this fixed \(q_0\), since nearby
positive profiles also have curvature exceeding one there.
In particular \(\Im X_*(q_0+iy)\ge0\) for every \(y>0\).
An equality at some height is harmless because targets
have strictly negative imaginary part.

The real pole-free gaps are \((p_L,0)\) and \((1,p_R)\),
where \(p_L,p_R\) are the nearest retained row poles,
with \(-\infty,+\infty\) when there is no pole on that side.
On each gap \(X_*'>0\), and the upper imaginary height
is \(-\pi\). The real value runs from negative infinity
to positive infinity in both gaps. At a finite outer
endpoint the retained logarithm gives the divergence.
At an infinite endpoint either the diffuse affine term
diverges in the required direction or all retained gaps
have the one sign giving the endpoint logarithmic
divergence. Each gap therefore lifts a real inverse into
the upper half-plane near the lower boundary of
\(\mathcal T\). By the divider, the two proper branches
never meet. Label them \(q_L,q_R\) on the left and right
of \(q_0\), respectively.

The upper limit phase \(f_*\) is given by (1). The phase
difference \(D=f_*(q_R,w)-f_*(q_L,w)\) satisfies
\(\Re D'=\Re(q_L-q_R)<0\). For fixed height,
as \(\Re w\to+\infty\) the left inverse tends to zero
by the term \(-\log q\). Then \(\Re D'\le-q_0/2\)
eventually. As \(\Re w\to-\infty\) the right inverse
tends to one and \(\Re D'\le-(1-q_0)/2\) eventually.
This gives one equal-modulus point at each height and
two connected complementary regions. On the lower
boundary the equal point is unique by the same real
gap monotonicity and these derivative bounds. Its phase
difference is \(i\pi\). On the upper boundary the real
function \(G_*(q)-q\log q-(1-q)\log(1-q)-w q\)
has two nondegenerate maxima at its unique coexistence
value. This follows from strict concavity of the
curvature by the [real stationary comparison](MULTIVARIATE_SHARED_SUPPORT_CONDITIONAL_SADDLE.md#a-local-exclusion-around-the-positive-real-axis).
Their phase difference is zero. At both boundary points
the nonzero derivative \(\Re D'\) permits local extension
of the graph. On the open equal-modulus graph \(E\),
\[
 0<\Im D<\pi,\qquad
 \frac{d}{d\Im w}\Im D=\frac{|D'|^2}{\Re D'}<0.
 \tag{2}
\]

Here the real phase joins through the two gaps as for
finite profiles. The real-cut exclusions apply to the
limit as well. At zero, one, and the nearest finite row
poles the local logarithmic coefficients and one-sided
imaginary derivatives have exactly the signs in the
[Rice boundary argument](MULTIVARIATE_RICE_CONTOUR_AND_REAL_BOUNDARY.md#where-the-upper-descending-contours-can-meet-the-real-axis).
Farther retained poles are isolated, with both boundary
phase derivatives of positive imaginary part above
the cuts. The infinite analytic tail and diffuse term
are real on the real axis locally, so do not change those
sign checks. A descending trajectory cannot land at
any finite real branch point. It can cross from
upper to lower only through a gap, and below it
has strictly negative vertical velocity. These are
statements for each fixed limiting profile and interior
\(w\); uniform estimates in the selected finite
contour neighborhoods then follow from compactness.

## Finite contour deformation for the limiting profile

On the union of the open half-planes and the two real
gaps use descent for the real part of the limiting phase.
For every fixed interior \(w\), a trajectory that stays
above a finite level cannot escape to infinity. To prove
this it suffices to bound its speed away from zero
outside a large disk along the trajectory. Below the real
axis this follows from \(\Im f_{-,q}<\arg(t)-\pi<0\).
Above the real axis the following alternatives apply.

* If retained gaps of both signs occur, take one of each.
  At large \(|q|\) their two angular contributions to
  \(\Im G_*'\) sum to at least \(\pi-o(1)\), uniformly
  up to the real rays. This is the angular law for
  \(\log(1+r(2q-1))\); the two limiting arguments are
  \(\arg(2q-1)\) and \(\pi-\arg(2q-1)\), multiplied by
  positive integers. All other terms, including the
  diffuse contribution, have nonnegative imaginary part.
  Thus the upper vertical velocity \(\Im f_{+,q}\)
  is bounded positively away from zero.
* If all the retained gaps have one sign, assume it is
  positive; reversal gives the other case. The horizontal
  component of descent points left for sufficiently large
  positive \(\Re q\), since each real row contribution is
  nonnegative there, one grows with \(\Re q\), and the
  diffuse contribution has the same sign. The real
  logarithm of \((1-q)/q\) is bounded in that range.
  Hence a forward trajectory has its real coordinate
  bounded above. If \(\tau>0\), at large height the
  diffuse term gives positive upper velocity. At bounded
  height and far to the left one retained row gives
  angular contribution tending to at least \(\pi\).
  If \(\tau=0\), a finite set of retained rows has total
  absolute gap at least three. Otherwise all gaps would
  sum to at most two and \(\rho<1\), by the elementary
  bound \(q(1-q)d_i^2/x_i(q)\le d_i^2/4\), strict when
  the sum of squares is four. At infinity with real
  coordinate bounded above this finite set contributes
  at least \(3\pi/2-o(1)\). Thus in all these cases
  the upper velocity is again bounded positively away
  from zero along the far part of the trajectory.
* With no retained rows the diffuse term is
  \(\tau(q-1/2)\), \(\tau>0\). The horizontal component
  bounds the real coordinate in both directions, and
  the vertical component is positive at large height.

Indeed if the height were bounded below, descent gives
a finite integral of the square of the speed. On the
part outside the disk, the speed lower bound then makes
the total time finite, and Cauchy–Schwarz makes the total
length finite there. This contradicts unbounded escape,
also if a trajectory were to reach infinity in finite
time. Finite cut accumulation is excluded as above. A
bounded trajectory staying above the relevant stopping
level therefore has a stationary accumulation point,
by the finite energy identity.

We can now use the
[finite sublevel contour argument](MULTIVARIATE_SHARED_SUPERCRITICAL_EVENTUAL_COPRIMALITY.md#isolation-of-the-larger-saddle-through-transitions)
with the limiting phase, without requiring convergence
of any finite-profile descent paths. Start with a
compact Rice loop, crossing the real gaps at its two
vertical sides and enclosing \([0,1]\). At strict
unequal heights flow it to a rectangle about the larger
limiting saddle or to a strict sublevel above the smaller
one. At equal heights use disjoint small rectangles about
both saddles and a common strict sublevel. In either
case the union of the chosen rectangles and the sublevel
is forward invariant, because each exiting face of a
rectangle is below the stopping level. Every starting
point of the compact loop enters that union in finite
time: if not, the preceding argument gives a stationary
accumulation point, and each such point is either inside
a rectangle or below the stopping level. Openness and
compactness give one finite flow time for the loop.
The paths in each rectangle are then replaced by oriented
copies of its local saddle segment and paths in the
sublevel components. This is the finite integer
crossing decomposition of the cited argument.

The complete finite homotopy lies in a compact subset of
the domain formed by the two half-planes and the two
gaps. For all sufficiently late finite profiles of the
sequence it avoids every exact pole: the retained cuts
on any bounded real range are stable, and the remaining
row poles are outside that range. The two limiting
saddles are simple and the corresponding finite-profile
saddles are uniformly close with all local parameter
derivatives. At the parameter of the homotopy replace
the two chosen saddle segments by their finite-profile
quadratic segments with short connectors at their ends.
By local convergence of the phases, the rest of the
homotopy still has a strict phase loss from the active
saddles used in the decomposition. On a sufficiently
small parameter disk the same construction and strict
loss apply, uniformly in late profiles. Therefore the
actual statements are Gaussian expansions for the
finite-profile stationary values, with their finite
integer crossing coefficients. Uniform analytic
Stirling and Gaussian integration give a relative
\(O(n^{-2})\) remainder after the first relative \(1/n\)
correction. The expansion holds at both orders with the
same coefficient functions for each chosen profile.
No approximation of a stationary exponential by the
limiting one has been made.

The nonzero integer coefficients are anchored at the
negative real boundary of the equal-modulus graph.
Near its lower endpoint the two limiting saddles lie
in the two real gaps, with positive real Hessians.
Use finite vertical segments through these points
at \(\arg t=\pi\). On the upper half of each segment,
\(\partial_y h=-(\pi+\Im X_*)<0\); below it the
phase decreases with downward motion. The horizontal
closing sides can be made strictly lower than both
saddles by choosing a large fixed height. In fact at
bounded real part a nonempty finite selection of
retained rows contributes a positive lower bound to
\(\Im G_*'\) for large positive height, uniformly over
the interval of the closing side. In the purely diffuse
case the diffuse term does so. The lower phase similarly
decreases uniformly at large negative height. These
finite contours give both nonzero contributions at
the lower endpoint and on an interior neighborhood.

The larger coefficient at strict dominance, and both
coefficients near an equal-height point, are locally
unique up to orientation. To see uniqueness, evaluate
any two candidate contour expansions with the same
orientation along the approximating sequence at a
fixed parameter where just the compared saddle has
larger limiting height. Its exponential dominates all
other contributions, so distinct integer coefficients
would contradict equality of the original periods.
Near equal height use this argument on the two sides
of the graph to compare both coefficients. The
anchored coefficients thus stay nonzero along the
connected equal-modulus graph. The larger coefficient
stays nonzero throughout each of the two connected
strict-dominance regions. In particular any limiting
nonreal zero lies on \(E\).

## Local expansions at the two real parameter boundaries

At every finite point of the lower boundary
\(\Im w=-\pi\), the same fixed vertical-loop argument
used at the equal point gives two actual periods, one
in each real gap. There is a strict outer phase loss
along each segment; the horizontal closing panels can be
taken low enough relative to both real saddles. Analytic
variation on a small disk, passing through the gap
charts, gives the same uniform Gaussian expansions as
above, including the first relative correction.
Hence a zero tending to the lower boundary must tend
to its unique real-scale equal-modulus point.

For \(\Im w=0\), use the coefficient sum and the real
entropy phase rather than moving the Rice loop to the
top boundary. The coefficient divided by the binomial
coefficient has the analytic interpolation with logarithm
\(n(G(q)-wq)+H_0(q)+H_1(q)/n+O(n^{-2})\).
The normalizing scalar here is the constant endpoint
coefficient; \(H_0\) and \(H_1\) are the Stirling log
amplitude and inverse-argument terms relative to that
endpoint. Multiplication by the interpolated binomial
gives the phase
\[
 G(q)-q\log q-(1-q)\log(1-q)-wq .
\]
The compact limiting real phase has one nondegenerate
global maximum except at the unique coexistence
parameter, where it has two separated nondegenerate
maxima. All are strictly inside \((0,1)\).
For a fixed boundary point, keep small disjoint real
intervals around its global maxima, with endpoints
at neighboring half-grid points. The entropy upper
bound on the binomial coefficient gives a uniform
strict exponential loss for the actual sum outside
those intervals, by compact convergence and the
bounded log amplitudes. This loss persists on
a small complex parameter disk.

On each interval apply the
[cotangent summation proof](MULTIVARIATE_VARYING_PROFILE_SUBCRITICAL_SADDLE.md#the-lattice-sum-and-its-contour).
Here it needs only the nondegenerate maximum and the
strict endpoint gap within the retained interval.
The vertical walls have cotangent modulus at most
one, and on lines of small fixed height \(h\) the
cotangent differs from its constant boundary value
by \(O(e^{-2\pi n h})\). The short horizontal contours
can be moved through the actual saddle, with endpoint
connectors below the strict phase gap. Uniformity
follows from the convergence of the analytic phases
near the limiting maximum and the strict negative
real curvature there. Expanding the Stirling amplitude
and the exact local quadratic integral by one more term
gives the uniform relative \(O(n^{-2})\) remainder
after the first relative \(1/n\) correction. The
intervals both give positive Gaussian amplitudes on
the real axis, since each represents a real maximum
with the original positive summand. This gives a
single nonzero leading term off coexistence, and two
actual terms with positive real amplitudes near
coexistence. Consequently a limit of zeros at the
upper boundary can occur only at the upper endpoint
of \(E\).

## Uniform comparison of consecutive orders

Consider a proposed sequence of common roots at
orders \(n,n+1\), with the profiles and \(n\) as
above. By positivity on the even order no common
root is real; take the upper half-plane by conjugation.
Pass to limits of the two centered parameters \(w,w'\)
in the closed strip. They have the same imaginary
part. Each must lie on the corresponding closed
equal-modulus graph of the same limiting profile,
by the actual expansions just proved. There is one
such point at each height, so their limits agree
and \(\Delta/n\to0\).

Use a common parameter disk about that limiting point
and the two actual saddle expansions there. Denote the
finite-profile stationary phase difference by \(D_n^*(w)\),
which is fixed by the profile and parameter, independently
of the polynomial order. Orient it as the right minus
left stationary value of the limiting branches.
The Gaussian coefficients are analytic, uniformly bounded
away from zero and infinity relative to the common scalar,
with bounded derivatives on a smaller disk. The stationary
difference and its derivatives converge to \(D,D'\).
At a zero, cancellation gives an equation of the form
\[
 e^{nD_n^*(w)}=L_n(w)\bigl(1+B_n(w)/n+O(n^{-2})\bigr),
 \tag{3}
\]
where \(L_n\) is the negative ratio of the two actual
leading amplitudes and the subscript denotes the profile.
The same functions \(D_n^*,L_n,B_n\) apply at order
\(n+1\) with parameter \(w'\). In particular
\(\Re D_n^*(w)=O(n^{-1})\).

If \(\Delta\ne0\) along a subsequence, divide the two
zero equations. The real part of the exponent on the
left is
\[
 \Re\bigl((n+1)D_n^*(w')-nD_n^*(w)\bigr)
  =-\Delta\,\Re(D_n^*)'(w)
                 +O\!\left(n^{-1}+|\Delta|/n+\Delta^2/n\right).
\]
The limiting derivative is strictly negative, including
the endpoints of the graph, by the branch divider.
The leading term is therefore bounded away from zero
when \(|\Delta|\ge1\) and dominates the error if
\(|\Delta|\to\infty\). But the modulus ratio on the right
of (3) tends to one because \(w'-w\to0\), a contradiction.

If \(\Delta=0\), then \(w'=w\). The adjacent quotient in
(3) yields \(e^{D_n^*(w)}=1+O(n^{-2})\), since the
common first relative term cancels. This is impossible
in the interior of the graph or at its lower endpoint
by (2) and the lower boundary value \(i\pi\).
At its upper endpoint the difference tends to zero,
so \(D_n^*(w)=O(n^{-2})\) and \(e^{nD_n^*(w)}\to1\).
The ratio \(L_n(w)\) instead approaches the negative
real half-line away from zero, by positivity and the
uniform local analytic bounds of the two positive-scale
Gaussian amplitudes. This contradicts (3).
Thus no counterexample sequence has unbounded terminals.
For each pair the resulting rational Bézout identity
transfers to every characteristic-zero field.
