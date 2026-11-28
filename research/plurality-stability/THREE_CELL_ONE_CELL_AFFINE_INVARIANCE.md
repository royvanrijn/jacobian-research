# Affine expansion invariance of one maximizing Gaussian cell

Let \(A=(A_1,A_2,A_3)\) be a balanced planar Gaussian
noise-stability maximizer at \(0<\rho<1\). Suppose one
cell \(A_i\) is invariant up to null sets under an invertible
affine map \(F(x)=p+M(x-p)\). If all eigenvalues of \(M\)
have modulus greater than one, or all have modulus less
than one, then \(A\) is the centered regular simplex and
\(p=0\).

The proof includes similarities, nonscalar Jordan blocks
and expanding maps whose eigenvalues admit a power-law
resonance. The same affine map need not preserve the two
other cells.
[Exact affine symmetries of the entire partition](THREE_CELL_CONTRACTING_AFFINE_SYMMETRY.md)

## The score bound used on curved arcs

For actual three-cell maximizers with positive masses,
complete adjacency implies \(1/3+\beta_i>0\) for every
centered price. Indeed at a regular interface of the two
other labels, both adjusted scores are at least the adjusted
\(i\)-score, so \(U_i-\beta_i\le1/3\). The strict positivity
of \(U_i\) yields the claim. At every winning point of \(i\)
the reverse inequality \(U_i-\beta_i\ge1/3\) holds. See the
separate pressure propagation argument for uniformity and
the consequences for general isolated ends.
[Pressure threshold and phase recurrence](THREE_CELL_PRESSURE_CONTRACTION_PROPAGATION.md)

## The analytic boundary of a single invariant cell

Use coordinates centered at \(p\), and replace the map
by its inverse if needed so all eigenvalues of \(M\)
have modulus greater than one. Write \(\Gamma_i\) for the support of the
distributional gradient of \(1_{A_i}\). It is preserved by
\(M\). Locally it is the closure of the analytic jump arcs
of cell \(i\). The actual-maximizer reciprocal-gradient
bound (G20) guarantees that every pair score on a regular
arc has nonzero gradient at an endpoint of that arc,
including the fixed point. Thus each analytic half-arc reaching
the fixed point continues as an analytic graph past its endpoint.
There are finitely many such half-arc germs by the local
analytic stratification. The support description follows
from local finite perimeter and the same stratification:
an open region with no actual jump arcs has zero
distributional derivative of its indicator.
[Regular arcs, local stratification and endpoint gradients, G20](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#regular-interface-arcs-and-isolated-triple-junctions)

The fixed point is in \(\Gamma_i\). Indeed invariance under the
contracting inverse sends positive-volume parts of both
the cell and its complement into every sufficiently small
neighborhood of the fixed point. Pass to a positive power of \(M\) fixing all
the finitely many half-arc germs. The power still expands and
fixes the oriented tangent ray of every germ. In real
eigenvalue cases we may take it even, so all eigenvalues of
that power are positive.

- If the original eigenvalues are a nonreal conjugate pair,
  the power fixing a real ray has both eigenvalues equal
  and positive. The original linear map is diagonalizable
  over the complex numbers, so that power is a positive
  scalar matrix. The same scalar alternative holds when
  the real eigenvalues have equal modulus and the map is
  diagonalizable.
- For a positive scalar matrix, any boundary graph
  \(v=h(u)\) over its tangent has
  \(h(\lambda u)=\lambda h(u)\), with \(h(0)=h'(0)=0\).
  Comparison of analytic coefficients makes it zero.
- For a nonscalar Jordan block with repeated positive
  eigenvalue \(\lambda>1\), the tangent of each fixed germ
  must be its eigenline. Use coordinates with that line
  horizontal and write the map as
  \((u,v)\mapsto(\lambda u+c v,\lambda v)\) with \(c\ne0\).
  The graph equation is
  \(h(\lambda u+c h(u))=\lambda h(u)\).
  A first nonzero term of order \(m\ge2\) would require
  \(\lambda^m=\lambda\), impossible. Hence all germs
  are straight also in this case.
- For two distinct positive eigenvalues \(1<\lambda<\mu\),
  use their eigenvectors as coordinates \((u,v)\). The
  basis need not be orthogonal. A vertical-tangent graph
  \(u=h(v)\) has \(h(\mu v)=\lambda h(v)\) and is
  necessarily zero. A horizontal-tangent graph
  \(v=h(u)\) has
  \[
                           h(\lambda u)=\mu h(u).            \tag{LA1}
  \]
  Comparison of analytic coefficients shows it is zero
  unless \(\mu=\lambda^m\) for an integer \(m\ge2\). In
  that case every graph is either zero or a monomial
  \(v=c u^m\). If an earlier positive power was taken to
  fix the germs, the same resonance relation holds for the
  original eigenvalue moduli.

The graph equalities on one-sided germ intervals suffice,
since both sides are analytic at zero. Iterating the inverse
of the powered map sends every point of \(\Gamma_i\setminus
\{0\}\) into the local germs. Conversely the images of every
local curve segment stay in \(\Gamma_i\). Thus each
straight germ extends to its full half-ray. In the resonant
case a nonzero monomial germ extends to the whole corresponding
half-curve. In particular the global support is contained
in the finite union of the resulting axes and power-law
half-curves, even though the other two cells can have
additional interfaces. These conclusions are about the
boundary of the one invariant cell.

## An extreme monomial forces a large pressure gap

Suppose there is a nonzero monomial boundary half-curve.
On its chosen sign half-plane of \(u\), write thresholds
as \(v=c|u|^m\). Pick a nonzero coefficient \(c\) of largest
absolute value for that sign. Along its curve the conditional
Gaussian source has mean equal in those coordinates to
\((\rho u,\rho c|u|^m)\) plus the fixed drift coming from
\((\rho-1)p\). Its coefficient relative to the mean horizontal
coordinate tends to \(\rho^{1-m}c\). This lies strictly outside the finite
set of coefficients on that half-plane. In fact the
conditional mean is at distance at least \(c_0|u|\) from
the entire boundary support for all sufficiently large
\(|u|\). To verify this, for \(u'\) within \(\eta|u|\) of
\(\rho u\), the values \(c'|u'/u|^m\) of every boundary
coefficient \(c'\) stay a fixed positive distance from
\(\rho c\) if \(\eta>0\) is small enough. The vertical
separation is then of order \(|u|^m\). Any point with
horizontal coordinate further away is already separated
in the first coordinate. The opposite sign half-plane
and the vertical axis are also a distance comparable to
\(|u|\) away. All distances in this fixed eigenbasis
are comparable to the physical Euclidean distances.

Consequently the ball about the conditional mean of radius
\(c_0|u|/2\) belongs to the invariant cell or to its
complement up to null sets. This follows from the zero
distributional gradient in that ball. The conditional
Gaussian tail outside it tends to zero, so along the
monomial curve the unadjusted score \(U_i=T_\rho1_{A_i}\)
tends to either zero or one along subsequences. At every
point of this boundary the adjusted \(i\)-score is
maximal and ties some other adjusted score, by continuity
from both sides and the winner rule. If \(U_i\) tends to
one, each other unadjusted score tends to zero. At least
one fixed other label ties \(i\) on a subsequence, forcing
a price contrast of one. This is impossible since both
cells have positive measure, so their unadjusted contrast
is everywhere strictly between \(-1\) and \(1\).

If \(U_i\) tends to zero, maximality of the adjusted
\(i\)-score at those observations gives
\(U_j-\beta_j\le U_i-\beta_i\) for every \(j\). Summing, and
using \(\sum_j U_j=1\) and \(\sum_j\beta_j=0\), gives
\(U_i-\beta_i\ge1/3\). Taking the limit, and using each
\(U_j\ge0\) in the pointwise inequality, yields
\[
              \beta_i\le-1/3,\qquad \beta_i\le\beta_j
                     \quad\text{for every }j.               \tag{LA2}
\]
In particular the total price gap to the other labels is
at least one, so some other label is priced at least
\(1/2\) higher than \(i\).
This contradicts the score bound (PP1), so there is no nonzero monomial
boundary. No pressure assumption beyond actual three-cell
maximality is needed.

## Identification and uniform separation

If all the boundary arcs of \(A_i\) are straight rays,
its indicator is constant on each open sector cut out by
them, since its distributional derivative vanishes there.
Hence it is homogeneous about \(p\). Consider its
first-order tangent label there. If only two strict
winner sets approach \(p\), their score difference has
nonzero gradient by (G20), so their positive-angle
tangent regions are half-planes. If all three approach,
their score gradients are distinct by (G20), and the
positive-angle tangent regions are the normal cones of
the extreme points of the convex hull of the three
gradients. Each of these is a connected sector of
aperture at most \(\pi\). A collinear middle label can
have a cusp of zero angle. To see the tangent convergence,
in each direction with a unique strict linear winner the
Taylor expansion gives that label at all sufficiently
small positive radii. The exceptional directions form
a finite set. When only one label approaches, the
tangent is constant.

Homogeneity and positive masses rule out a full or empty cell.
A half-plane is excluded since each maximizing cell uses
both of its active directions. Thus \(A_i\) is an exact
translated pointed wedge with vertex \(p\). The
translated-wedge rigidity theorem puts \(p=0\) and
identifies the entire partition as the centered regular
simplex.
[Each maximizing cell uses both active directions](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#each-cell-uses-both-active-directions)
[Translated wedge rigidity](THREE_CELL_TRANSLATED_WEDGE_RIGIDITY.md)

There is a uniform form on compact positive-correlation bands,
bounded centers \(p\), and compact sets of invertible matrices
whose eigenvalues are all uniformly outside or all uniformly
inside the unit circle. The Gaussian symmetric difference
between \(A_i\) and \(F A_i\) has a positive lower bound over
nonsimplex balanced maximizers and all those affine maps.
Otherwise pass to a subsequence with converging correlations,
centers and matrices, fixed label, and strong Gaussian
convergence to a maximizing partition, using planar compactness.
The limit remains nonsimplex by uniform simplex isolation.
Continuity of Gaussian pullback by these converging affine
maps follows by bounded continuous approximation and uniform
bounds on the transformed densities on compact sets, with
uniformly small tails. The limiting cell would be invariant
under the limiting map, contradicting the theorem.
[Simplex isolation and strong compactness](POSITIVE_CORRELATION_CONTINUATION.md)
