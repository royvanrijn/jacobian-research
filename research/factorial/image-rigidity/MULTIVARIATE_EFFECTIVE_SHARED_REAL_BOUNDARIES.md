# Effective all-field comparison in fixed shared curvature bands

Write \(W,\rho,w\) for the shared curvature, real stationary
curvature, and chord-centered field of the
[shared saddle comparison](MULTIVARIATE_UNIFORM_SUPERCRITICAL_SUPPORTS.md).
In particular \(w_N-w_n=-\Delta\log(N/n)\) at a common root.

**Theorem.** For positive rational \(W_-\le W_+\) and
\(0<\kappa<1\), a terminating finite search computes a uniform
consecutive coprimality cutoff for every fully shared ordinary
profile on \(W_-\le W\le W_+\) satisfying either
\(\max\rho\le1-\kappa\) or \(\max\rho\ge1+\kappa\).
It covers every physical field and every integer degree
difference. The cutoff is unevaluated.

To bound the real coordinate, let
\(R=\max(1/2,W_+/(2+W_+))<1\). This bounds every row
\(|r_i|=|d_i|/(2m_i)\) because \(2|d_i r_i|\le W_+\) and
\(|r_i|\le |d_i|/(2+|d_i|)\). The
[endpoint annulus](MULTIVARIATE_UNIFORM_ROOT_ANNULUS.md#shared-rows-and-large-degree-imbalance)
then bounds the ordinary endpoint-normalized log radius
by \(W_+/[2(1-R)]\). The chord-centered real coordinate
differs from its negative by at most \(W_+/[j(1-R)]\)
by the explicit endpoint center estimate in the shared
comparison. Thus any possible root at \(j\ge2\) has
\(|\Re w|\le M=2W_+/(1-R)\).
Use retained and diffuse profile nets and pairs of fields
with equal heights in the closed rectangle
\(\{|\Re w|\le M,\ -\pi\le\Im w\le0\}\).
For interior heights use the
[regular](MULTIVARIATE_EFFECTIVE_SHARED_REGULAR_COVER.md)
and [fold](MULTIVARIATE_EFFECTIVE_SHARED_SUBCRITICAL_FOLD.md)
witnesses. The remaining local witnesses are given below.

## Lower boundary

At \(\Im w=-\pi\), the [real-gap comparison](MULTIVARIATE_UNIFORM_SUPERCRITICAL_SUPPORTS.md#local-expansions-at-the-two-real-parameter-boundaries)
gives a saddle in each real gap, with real parts
\(q_L<0<1<q_R\) and positive real Hessians. This geometry
applies to both stated bands. The two stationary heights
are equal at a unique real scale because the difference
has real derivative \(q_L-q_R<-1\) in the orientation
right minus left, and at equality their phase multiplier
is \(-1\). The two-gap vertical-rectangle construction
gives actual nonzero coefficients and strict exterior
losses. The [subcritical active-arc proof](MULTIVARIATE_UNIFORM_SUBCRITICAL_EQUAL_DEGREE.md#actual-coefficients-and-the-active-arc)
uses the same gap construction.

These are finite rational contour-witness tests. Isolate
each simple saddle, choose a vertical direction of
quadratic decay, and enclose the finite exterior paths
with the strict loss relative to the leading height.
When both terms are needed the paths are below both
heights. The contour identity uses winding one at \(1/2\)
and zero at the external Gamma poles, as in the regular
cover. Across the left gap the rational factorial ratio is
exactly
\[
 \frac{\Gamma(j(q-1))}{\Gamma(j q+1)}
   =(-1)^{j+1}\frac{\Gamma(-j q)}{\Gamma(j(1-q)+1)} .
\]
Across the right gap both original arguments are positive.
The exact reflection, continued on the left neighborhood,
and the right ratio give uniform Stirling errors on these
crossings. The logarithm of the physical parameter with
its sign reversed is analytic through the lower boundary.
Thus a small disk through it has the ordinary analytic
two-term expansion with its constant phase factor included
in the appropriate exponential. Stirling error and the
local Gaussian bounds give a computable relative
\(O(1/j)\) error. In particular the coefficients used at
both consecutive orders are the same analytic functions.
The proof of these estimates on the two half-planes and
across the gaps is the
[exact Rice expansion](MULTIVARIATE_RICE_CONTOUR_AND_REAL_BOUNDARY.md#exact-contour-and-its-phases).

If one field has unequal boundary heights it has a single
dominant term on a small neighborhood. If both fields
have equal heights they coincide by the uniqueness of the
scale. The common convex chart has positive lower bounds
on \(|\Re D'|\) and \(|e^D-1|\), since at the center the
former is at least one and the latter is two. These
strict tests suffice for the two-order comparison in
the regular cover.

## Upper boundary

At every limiting field on \(\Im w=0\), the real entropy
phase on \([0,1]\) has one nondegenerate global maximum,
apart from the unique supercritical coexistence scale
where it has two nondegenerate maxima of the same height.
In the subcritical band the strict curvature inequality
makes the phase strictly concave. In the supercritical band
the [real stationary geometry](MULTIVARIATE_UNIFORM_SUPERCRITICAL_SUPPORTS.md#the-limiting-phase-and-its-real-boundary-geometry)
gives precisely the described coexistence. Keep small
disjoint intervals in \((0,1)\) around the one or two
global maxima. Off them the phase has a strict loss by
compactness. This also holds for the exact coefficient
sum at nearby complex fields by the entropy upper bound
on the binomial and the order-independent log amplitude
bounds.

A certificate can specify these intervals, isolate their
saddles, and check losses outside by interval subdivision
of \([0,1]\), using continuity of the entropy term at
both endpoints. Deform the cotangent summation panels of
the [lattice proof](MULTIVARIATE_VARYING_PROFILE_SUBCRITICAL_SADDLE.md#the-lattice-sum-and-its-contour)
inside arbitrarily small analytic neighborhoods of the
retained intervals. Their endpoints may move to the
nearest half-grid points at large order. The vertical
cotangent modulus is then at most one, the horizontal
cotangent is its limiting constant with computably
exponential error, and the short connectors lie in
strict sublevels. In a one-maximum chart this provides
a computable zero-free threshold. In a two-maximum
chart, each saddle has the same finite expansion for all
orders,
\[
 e^{jP_i(w)}A_i(w)
       \bigl(1+B_i(w)/j+O(j^{-2})\bigr),
\]
after removal of the common nonzero normalization.
At the real center both \(A_i>0\) by the Gaussian
expansion of the original positive coefficient sum.

For the \(O(j^{-2})\) bound, use Stirling through the
inverse-argument term on the compact saddle disks.
The real parts of the normalized row arguments stay
uniformly positive there, and the reciprocal powers sum
to a computable bound from \(\sum_i1/m_i\le W_+\).
On the local straight Gaussian segment put
\(q-q_s=j^{-1/2}v\) in the selected descent direction.
On \(|v|\le j^{1/12}\), Taylor expansion of the analytic
phase and amplitude through order \(j^{-3/2}\) gives the
first relative correction, with a remainder bounded by
\(j^{-2}\) times a fixed polynomial times the absolute
Gaussian. Its coefficients come from Cauchy bounds on
the strictly larger saddle disk. Odd terms at orders
\(j^{-1/2},j^{-3/2}\) integrate to zero on the symmetric
segment. The remaining segment has strict quadratic
loss and is \(O(j^{-2})\) relative to the leading local
scale with computable constant, since its tail is
bounded by a polynomial times \(e^{-c j^{1/6}}\).
The exterior loss is fixed. These steps certify the
error bound using only finite interval inequalities and
elementary Gaussian integrals.

At a pair of upper boundary fields, a one-term chart
excludes a root unless both fields lie at coexistence.
There they coincide. Choose a convex neighborhood with
\(|\Re D'|\ge c>0\), since \(D'\) is the difference
of the two distinct real saddle coordinates up to sign.
For \(\Delta\ne0\) use the same real modulus comparison
as for the regular cover. For \(\Delta=0\) both fields
are identical. The two-saddle zero equation in that
neighborhood has the analytic form
\[
 e^{jD(w)}=L(w)\bigl(1+B(w)/j+O(j^{-2})\bigr),
 \qquad L=-A_1/A_2 .
\]
Shrink the disk so \(|D|\le1/4\) and \(\Re L\le-\eta<0\).
Dividing the two zero equations cancels their common
first correction up to \(O(n^{-2})\), yielding
\(e^{D(w)}=1+O(n^{-2})\).
On \(|D|\le1/4\) this forces \(|D(w)|=O(n^{-2})\) with
a computable constant, hence \(\Re e^{nD(w)}>1/2\)
above a computable order. The zero equation instead
forces \(\Re e^{nD(w)}<0\) by the negative amplitude
ratio and the bounded errors. This certifies the
upper endpoint with equal total degrees.

## Finite search on the closed field strip

For pairs with an interior height, the regular and fold
arguments give either a one-field exclusion or a common
chart whenever both fields are on the respective active
locus. At a boundary height both fields are on that
same boundary, and the local alternatives above apply.
All witnesses have strict open tests and persist for a
small profile cell and field neighborhood. The
admissible profile and paired-field family is compact;
phase derivatives on all fixed contour neighborhoods
are enclosed by the finite retained and diffuse nets.
Dovetail the tolerances, rational witnesses, and
refinements, testing the closed band conditions on
\(\max\rho\) by its bounded real derivative.
A finite subcover ensures termination. The maximum of
the certified thresholds then applies to every common
root candidate in the required fixed bands.
