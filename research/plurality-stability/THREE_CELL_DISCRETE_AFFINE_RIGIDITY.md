# Discrete affine symmetries of planar three-cell maximizers

Fix a balanced planar three-cell Gaussian maximizer \(A\)
at \(0<r<1\). A nonsimplex \(A\) cannot be invariant under
an invertible affine map whose two eigenvalues have modulus
less than one. No such maximizer can be invariant under
a fixed-point affine map with one contracting eigenvalue
and one neutral eigenvalue of modulus one, or under a
nontrivial unipotent affine map with a fixed point. These
statements also allow the affine map to permute the labels.
Take a finite power to remove that permutation. The argument
also applies to the inverse map whenever it has the stated
linear type.

The proof uses the analytic adjusted pair contrasts
\[
 F_{ij}=T_r(1_{A_i}-1_{A_j})-(\beta_i-\beta_j)
\]
of the actual maximizer. Each contrast is nonconstant by
injectivity of Gaussian smoothing. Its zero set contains
all actual \(i,j\) jump arcs. The essential cell supports
are invariant under any label-preserving affine symmetry,
so every iterate of a point at a pair wall still gives
pair equality, even if the iterate is a singular junction.
We also use the reciprocal-gradient restriction (G20),
full active rank, and the exclusion of an individual
cylindrical cell.
[Maximizer structure](THREE_CELL_GLOBAL_OBSTRUCTIONS.md)

## Strict contractions about a fixed point

Translate coordinates so that the fixed point of the map
is the coordinate origin. This is a change of coordinates,
so the Gaussian conditional sample has the extra constant
drift associated with the new center.

First suppose a power of the linear part is scalar \(aI\),
with \(0<a<1\). At any actual interface point with vector
\(z\ne0\) from the fixed point, the associated contrast
vanishes at all \(a^n z\). Analytic continuation on the
line makes it vanish along that line. For each pair, the
possible directions are among the zeros on the circle of
the first nonzero homogeneous Taylor term of its contrast
at the fixed point. This is a finite set, since the
nonzero homogeneous polynomial cannot vanish on the
whole circle. Thus the walls lie on a finite union of
rays. In each complementary open sector the indicators
have zero distributional gradient, so the cells are
conical about that point. Translated-cone rigidity makes
the balanced maximizer the centered simplex.

If the eigenvalues are nonreal, a real linear change of
coordinates makes the map a rotation followed by scalar
contraction. A rational rotation angle reduces to the
scalar case. For an irrational angle take an actual
pair-interface point different from the fixed point.
On dividing the stationary zero identities for its
iterates by \(a^{nm}\), with \(m\) the first nonzero Taylor
degree of that pair contrast, the remainders tend to
zero. Since the iterated directions are dense on the
circle in these coordinates, the leading homogeneous
term would vanish on that circle. This is impossible.

For distinct real eigenvalues, take an even power and
write the contraction in its eigenbasis as
\((x,y)\mapsto(a x,a^q y)\) with \(0<a<1\), \(q>1\).
For each contrast let \(P(x,y)\) be the sum of its
Taylor monomials of least weighted degree, assigning
degrees \(1,q\) to \(x,y\). At an actual pair-wall point
\((x,y)\), evaluate the contrast on its iterates and
divide by \(a^{n\deg P}\). Truncating the Taylor
series at any ordinary degree exceeding \(\deg P\)
leaves a remainder tending to zero after division.
It follows that \(P(x,y)=0\).

On each sign half-plane of \(x\) the zero set of each
weighted leading polynomial is a finite union of graphs
\(y=c|x|^q\). Indeed \(P(\pm1,c)\) is a nonzero
polynomial, since two distinct monomials of the same
weighted degree have different \(y\)-degrees. The cells
therefore have constant labels on the connected strips
between these graphs. A transition between different
labels has stationary pair equality all along its graph.
If there is a positive transition, take the largest
one in that half-plane. At observations
\((x,c|x|^q)\) with \(|x|\to\infty\), the conditional
sample lies on the same sign half-plane with probability
tending to one and its rescaled height tends almost surely
to \(r^{1-q}c>c\). This follows from conditional mean
\(r(x,y)\) plus the constant drift, and bounded Gaussian
displacements for each fixed noise realization. The
argument works in this eigenbasis even if it is oblique.
The source limit has one deterministic label, an
impossibility at an active pair wall. Indeed in the
presence of a deterministic limiting label all other
adjusted scores are strictly smaller, since all pair
price differences are strictly below one in magnitude.
A smallest negative transition is excluded identically.

The only remaining labels are constant on the coordinate
quadrants. All three occur, so a label is repeated.
If it occupies adjacent quadrants, it is a half-plane,
contrary to the individual-cell restriction. Otherwise
the same label occurs on opposite quadrants. Across
each straight axis ray the corresponding tangential
derivatives of the adjusted scores agree at the fixed
point. Matching the two coordinate derivatives therefore
forces all three gradients there to coincide, since the
opposite quadrants already share both derivatives.
This contradicts (G20) along any of the incident pair
arcs.

The remaining strict-contraction case, after taking a
power, has a real Jordan block
\((x,y)\mapsto a(x+\tau y,y)\) with \(0<a<1\),
\(\tau\ne0\). At a pair-wall point use the first nonzero
ordinary homogeneous Taylor term \(P_m\) of its contrast.
The iterates give
\[
 P_m(x+n\tau y,y)=O(a^n(1+n)^{m+1}).
\]
The left side is a polynomial in \(n\), hence identically
zero. If \(y\ne0\) this would make \(P_m\) vanish on a
whole line at height \(y\), which forces the homogeneous
polynomial to vanish identically. Thus all actual
interfaces lie in the eigenline, which cannot separate
three positive labels.

## One contracting direction with a fixed line

After taking a power, the eigenvalues are \(a\in(0,1)\)
and \(1\), with linearly independent eigenvectors.
Let \(H\) be the neutral fixed line, and use orthonormal
coordinates \((x,y)\) with \(x\) parallel to the contracting
eigenvector. In these coordinates \(H=\{x=ky+b\}\) for
some real \(k,b\). The map acts on each horizontal line
by contraction about its intersection with \(H\).

For almost every fixed \(y\), no pair contrast vanishes
identically on its horizontal line. To obtain an exceptional
null set, for each nonzero contrast choose a fixed
horizontal coordinate where its restriction in \(y\)
is not identically zero. Its discrete zero set contains
all exceptional heights for that contrast. A genuine
label change in one of the two rays of a nonexceptional
line would give stationary zeros accumulating at its
intersection with \(H\). One-dimensional analyticity
would force the contrast identically zero on the line.
Thus on either side of \(H\) the label depends only on
\(y\). Each of these profiles is piecewise constant with
locally discrete transitions, by restriction of the
analytic winner rule to a generic line parallel to \(H\).
At least one profile is nonconstant, since all three
labels have positive mass. A transition of that profile,
say at height \(c\) for a pair \(i,j\), gives a regular
horizontal pair arc away from \(H\), and stationary
equality continues across its whole supporting line.

If \(k=0\) the ray profiles are invariant under every
positive orthogonal dilation of the coordinate normal
to \(H\), leaving its tangent coordinate fixed. This
is the rank-one affine flow already excluded for
noncylindrical maximizers. If \(k\ne0\), let
\(F_+(t),F_-(t)\) be the \(i,j\) indicator contrasts in
the two profiles. At height \(c\), deconvolution along
the horizontal direction in the stationary zero
identity says that
\[
 \int_{\mathbb R}
 \bigl[F_+(t)1_{\{x>kt+b\}}+F_-(t)1_{\{x<kt+b\}}\bigr]
        p_\sigma(t-r c)\,dt
\]
is constant in \(x\), where \(\sigma=\sqrt{1-r^2}\)
and \(p_\sigma\) is the normal density of variance
\(\sigma^2\). The horizontal conditional Gaussian
convolution is injective on bounded functions by its
nonvanishing Fourier multiplier. Differentiation of
this constant integral gives \(F_+=F_-\) almost
everywhere because the density is positive and \(k\ne0\).
The three possible values \(1,-1,0\) distinguish all
three labels. Hence the profiles agree and the partition
is cylindrical, a contradiction.
[Rank-one orthogonal affine-flow exclusion](THREE_CELL_FIRST_CONTACT_RIGIDITY_TESTS.md#semidefinite-affine-flows)

## Shears about a fixed line

In orthonormal coordinates about a fixed point the
shear is \((x,y)\mapsto(x+c y,y)\) with \(c\ne0\).
By invariance, almost every slice of each indicator in
the horizontal coordinate is periodic with period
\(|cy|\). Except for the null set of heights at which
some pair contrast vanishes identically, strict winners
coincide with the labels on the open intervals between
the discrete zeros of the pair contrasts.

The slices at sufficiently small \(|y|>0\) are pure,
meaning that exactly one label occurs, almost
everywhere in height. Otherwise choose a sequence
of generic heights \(y_n\to0\) with mixed labels.
Some pair \(i,j\) then has a genuine horizontal-slice
transition for infinitely many of those heights.
By periodicity each transition repeats with spacing
\(|c y_n|\). Factor the analytic contrast by \(y^m\)
where \(m\) is the first order whose normal Taylor
coefficient on the fixed line is not identically
zero. For every \(x_0\) there are periodic roots
\((x_n,y_n)\to(x_0,0)\), hence that first coefficient
vanishes at \(x_0\). This contradiction proves purity
near the fixed line.

Let \(p_i(t)\) be the fraction of a horizontal period
occupied by label \(i\) at height \(t\). For every
fixed observation height \(y\), the score \(U_i(x,y)\)
as \(x\to+\infty\) tends to
\[
             \int p_i(t)\,p_\sigma(t-r y-\delta)\,dt,        \tag{AS1}
\]
where \(\delta\) is the constant normal drift in the
shifted coordinates. To prove this, convolve the
periodic slice with the horizontal normal density.
On \(\eta\le |t|\le M\) its Fourier series is uniformly
absolutely summable, since its \(n\)-th coefficient
has Gaussian damping \(\exp(-2\pi^2\sigma^2 n^2/(c^2t^2))\).
For each nonzero mode the integration in \(t\) has
phase equal to a nonzero multiple of the horizontal
conditional mean divided by \(t\). After the change
of variable \(u=1/t\), the Riemann–Lebesgue lemma
makes that term vanish in the limit. Letting
\(\eta\downarrow0\) and \(M\to\infty\) controls the
discarded part by the normal probability in height.
The remaining constant Fourier coefficient is
exactly (AS1).

At any generic mixed height both occurring labels
win at arbitrarily large \(x\), by periodicity.
They must tie in the limiting adjusted scores.
If labels \(i,j\) mix on a positive-measure set of
heights, the two smoothed profiles in (AS1) differ
by their price difference everywhere, by real
analyticity in \(y\). One-dimensional Gaussian
injectivity now gives \(p_i-p_j\) constant almost
everywhere. Pure slices near the fixed line force
the constant to be \(-1,0\), or \(1\); either
nonzero choice contradicts positive Gaussian mass
of both labels. Thus \(p_i=p_j\). If two different
pairs mix on positive-measure sets, this gives all
three profiles identically \(1/3\), contradicting
purity near the fixed line. If only one pair mixes
on a positive-measure set, the remaining label
occurs only on pure slices up to null sets and is
cylindrical. If none mixes, every label is
cylindrical. All possibilities are excluded.

## Translations and rational correlations

For \(r\) irrational or \(r=p/q\) in lowest terms
with \(p\ge2\), a planar balanced maximizer has no
nonzero translation period, even up to permutation
of its labels. The argument concerns an exact
translation, without a fixed point.

Remove a possible permutation by taking a power,
and write \(L>0\) for the resulting period in its
direction. Gaussian smoothing makes each adjusted
score periodic with period \(L/r\) in that
direction. At almost every point the strict
winning label persists at all its translates by
integer multiples of \(L\). It is therefore still
a winner after averaging the scores on that orbit.
In the irrational case this average is the
continuous mean in the phase of the score, so
the averaged scores are independent of the
coordinate along the period. The inequalities
are strict for almost every transverse coordinate
and almost every starting phase: for such a phase
the original strict inequalities hold at the
starting point and nonstrict inequalities hold
on its dense orbit closure. Continuity makes the
integrals of both gaps strictly positive. Thus
the partition is cylindrical.

In the rational case the average is over the
\(q\) distinct phases. The strict inequalities
persist on this finite orbit, and the averaged
scores now have period \(L/p\). Their almost
everywhere unique winner inherits that period.
Repeating yields periods \(L/p^n\) tending to zero.
Continuity of translations in local \(L^1\)
makes the indicators invariant under every
translation in that direction, again contradicting
full active rank. The averaging leaves the
reciprocal integers \(r=1/q\) open.

## Uniform separation from affine images

Fix a compact correlation band in \((0,1)\) and a
compact family of invertible affine maps each covered
by the stated symmetry exclusions. Nonsimplex ties
have label disagreement with their pullbacks by those
maps bounded below by a positive constant, even
allowing a label permutation. Otherwise choose a
convergent subsequence of maps and permutations and
use strong maximizer compactness in the active plane.
Uniform local isolation keeps the limit nonsimplex.
Affine pullbacks pass to the limit in Gaussian label
measure: on compact input sets the pulled-back laws
have uniformly bounded densities on a fixed compact
image set, and continuous approximations to the limit
indicators handle the varying maps there. Tightness
controls the discarded input tails. Vanishing
disagreement would give an exact excluded symmetry.
[Local comparison and maximizer continuation](POSITIVE_CORRELATION_CONTINUATION.md)
