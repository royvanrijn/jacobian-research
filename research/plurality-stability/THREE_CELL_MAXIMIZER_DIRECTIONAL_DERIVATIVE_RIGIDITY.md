# Directional derivatives distinguish the active planes of three-cell maximizers

Fix \(0<\rho<1\). An individual cell of an ordinary
three-cell Gaussian noise-stability maximizer with positive
prescribed masses is never annihilated in its active plane
by a finite product of nonzero directional derivatives.
Thus in a finite pointwise linear identity between label
vectors of actual three-cell maximizers, with nonzero
coefficients, every active plane that occurs must occur
at least twice.

In particular, let \(A,B\) be maximizers with the same
prescribed masses and distinct active planes. If \(C,D\)
are also maximizers at the same masses and correlation and
\[
                    1_A+1_B=1_C+1_D,                           \tag{DD1}
\]
then \(\{A,B\}=\{C,D\}\) as labeled partitions up to null
sets. The conclusion does not assume that a partition is
conical or simplicial. It gives strict cost for every
nontrivial balanced complementary patch of two maximizers
in distinct planes.

## An individual maximizing cell has no directional product equation

Dimension reduction gives each maximizer an active plane
and, for three labels, each individual cell in that plane
depends on both coordinates. The reduced interfaces are
locally finite analytic arcs. A regular pair arc continues
at its finite accumulations until it meets a genuine triple
junction. At a noncollinear genuine triple, a given incident
cell has two nonparallel boundary rays. At a collinear
genuine triple, all incident arcs have the same tangent
line. These properties hold for positive prescribed
masses, including unequal masses.
[Heilman–Tarter, Theorem 1.9](https://doi.org/10.1017/fms.2021.56)
[Full active direction for each cell and local jump geometry, G20](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#regular-interface-arcs-and-isolated-triple-junctions)

Suppose \(f=1_{A_i}\) has \(P(D)f=0\) in distributions,
where \(P(D)=\prod_{j=1}^q D_{\tau_j}^{m_j}\), the
\(\tau_j\) are unit vectors in distinct unoriented
directions, and \(m_j\ge1\). Any repeated parallel factors,
including opposite orientations, reduce to this form up
to a nonzero scalar. At a regular arc straighten the
boundary locally to \(t=0\) with a smooth change of
coordinates. The indicator is a step in \(t\). If
\(M=\sum_jm_j\), the coefficient of its highest normal
delta derivative \(\delta^{(M-1)}(t)\) in \(P(D)f\)
is, up to a nonzero smooth factor, the restriction
of \(\prod_j(\tau_j\cdot\nabla t)^{m_j}\) to the arc.
All terms with derivatives on coefficients or fewer
normal derivatives have lower normal order. Thus at each
regular point some \(\tau_j\) is tangent to the interface.
Continuity of the tangent into this finite set of
directions forces each connected arc to be straight.

At a finite point \(p\) in the closure of these boundary
arcs, the local boundary of the cell consists of finitely
many straight segments issuing from \(p\). Write \(n_j\)
for a chosen normal to direction \(\tau_j\). On a sufficiently
small neighborhood, let \(\mu_{j,+}\) and \(\mu_{j,-}\)
be arclength on the rays in directions \(\tau_j\) and
\(-\tau_j\). With jump signs \(\alpha_{j,\pm}\in\{-1,0,1\}\)
measured using the same \(n_j\) on both rays,
\[
 D_w f=\sum_j (w\cdot n_j)
     \bigl(\alpha_{j,+}\mu_{j,+}+\alpha_{j,-}\mu_{j,-}\bigr)
\]
locally for any constant \(w\). A zero sign encodes an absent
arc. This first derivative has no point contribution,
since \(f\) has locally finite perimeter and the jump
measure is carried by its reduced boundary. Locally
\(D_{\tau_j}\mu_{j,+}=\delta_p\) and
\(D_{\tau_j}\mu_{j,-}=-\delta_p\), by integration along
a ray with compactly supported tests. Commuting \(D_w\)
through \(P(D)\) therefore gives
\[
 0=\sum_j(w\cdot n_j)(\alpha_{j,+}-\alpha_{j,-})
                 (P/D_{\tau_j})(D)\delta_p.                    \tag{DD2}
\]
Here \((P/D_{\tau_j})(D)\) denotes the product with one
copy of \(D_{\tau_j}\) removed.

The corresponding polynomials
\(\prod_l(\tau_l\cdot\xi)^{m_l}/(\tau_j\cdot\xi)\)
are independent. Factor out the common
\(\prod_l(\tau_l\cdot\xi)^{m_l-1}\); at a nonzero point of
\(\tau_j^\perp\), only the \(j\)-th remaining polynomial
is nonzero. Polynomial independence is equivalent to
independence of the differential operators applied to
\(\delta_p\), by testing arbitrary Taylor jets at \(p\).
Choose \(w\) with all \(w\cdot n_j\ne0\). Equation (DD2)
makes \(\alpha_{j,+}=\alpha_{j,-}\) for every \(j\).
Every jump thus continues through a finite endpoint
along the same full line, with its jump sign.

The continuation is global. For a line supporting an
initial arc, the local continuation just proved makes
its jump points relatively open in the line, including
finite accumulations. They are relatively closed by
the local finite analytic graph description. Thus the
line carries the jump throughout. Two nonparallel
lines would cross, giving four incident rays in two
directions for the same cell. A regular pair point has
one boundary curve. At a noncollinear genuine triple
an individual cell has only two incident rays, and
the collinear alternative has only one tangent line.
No point allows the proposed crossing. Thus all
boundary lines are parallel. Distributional derivatives
of \(f\) along their common tangent vanish, making the
cell cylindrical and contradicting its full active
dimension.
[Each maximizing cell uses both active directions](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#each-cell-uses-both-active-directions)

## Every active plane in an identity must repeat

Suppose a finite linear combination of maximizing label
vectors is constant and some active plane \(V\) occurs
exactly once. For each other distinct plane \(W\), choose a
vector in \(W^\perp\) whose projection onto \(V\) is
nonzero. Such a vector exists because the planes have
equal dimension and are distinct. Apply the product of
the corresponding distributional derivatives to the
identity. Every term except the one in \(V\) vanishes.
On the term in \(V\) the derivatives act through their
nonzero projections and contradict the obstruction above.
The case of a family of only one term is impossible
because a positive-mass cell indicator is nonconstant.
This proves the active-plane occurrence assertion. The
argument permits a constant on the other side and any
nonzero real coefficients.

## A difference supported on one common coordinate is zero

Let \(A,B\) be actual planar three-cell maximizers at
the same prescribed masses and correlation. Suppose
\(d=1_A-1_B\) depends only on a coordinate \(s\), with
orthogonal coordinate \(t\). Almost everywhere, either
\(d(s)=0\) or \(d(s)\) is a nonzero label difference
\(e_i-e_j\). In the second case the fiber is purely
label \(i\) in \(A\) and purely label \(j\) in \(B\),
apart from null sets.

Since \(A\) uses both active directions, some regular
pair arc \(i/j\) has a tangent with a nonzero \(s\)
component. Otherwise all regular walls would be tangent
to the \(t\)-direction and all label indicators would
be invariant in that direction. On a small neighborhood
of the chosen arc its projection onto the \(s\)-axis
is an interval and the arc makes those fibers contain
both labels on positive-measure parts. Hence \(d(s)=0\)
for almost every \(s\) in that interval. In particular
the two partitions agree on an open slab up to null
sets, and the regular arc in the slab is a common
\(i/j\) interface.

On that arc subtract the two stationary pair equations.
Writing \(\beta^A,\beta^B\) for the label prices yields
\[
 T_\rho(d_i-d_j)(s)
      =(\beta_i^A-\beta_j^A)-(\beta_i^B-\beta_j^B).
                                                               \tag{DD3}
\]
The left side is a one-variable analytic function, so
it is constant on the whole line. By Gaussian
injectivity and equality of the prescribed masses,
\(d_i=d_j\) almost everywhere. A nonzero three-label
difference has all three coordinates distinct, so
this forces \(d=0\).

## Complementary patches of two distinct planes

Apply active-plane occurrence to (DD1). The planes
of \(C,D\) must be \(V_A,V_B\), respectively, after
interchanging \(C,D\) if necessary. Their differences
\(1_A-1_C=1_D-1_B\) are invariant along both
\(V_A^\perp\) and \(V_B^\perp\). Hence the common
difference depends only on \(V_A\cap V_B\).
Indeed invariance under the sum of the complementary
subspaces follows by smoothing and passing to local
\(L^1\). If the intersection is zero-dimensional,
the difference is constant and must be zero by balance.
If the intersection is a line, (DD3) forces it to
be zero. Thus \(A=C\) and \(B=D\).

Let \(E\) be a measurable cut with
\(\int_E(1_B-1_A)\,d\gamma=0\), and write \(e=1_B-1_A\).
The indicator vectors \(1_A+1_Ee\) and
\(1_A+1_{E^c}e\) are ordinary balanced competitors
with the given masses. The sum of their deficits
relative to the common maximum is
\[
              2\langle 1_E e,T_\rho(1_{E^c}e)\rangle_\gamma .
                                                               \tag{DD4}
\]
A zero deficit would make both patches actual maximizers
and give (DD1). Thus (DD4) is strictly positive
whenever both parts of the cut contain positive-measure
disagreement. The statement leaves two maximizers in
the same active plane outside its scope.
[Complementary patch identity and the simplex antecedent](THREE_CELL_COMPLEMENTARY_PATCH_RIGIDITY.md#strict-cost-of-a-nontrivial-cut)


## Uniform loss with separated planes and positive patch masses

Let the common prescribed mass vector vary over a compact
subset of the positive mass simplex and \(\rho\) over a compact
subset of \((0,1)\). With the notation of (DD4), suppose that
\[
 \|P_{V_A}-P_{V_B}\|_{\mathrm{op}}\ge\delta,\qquad
 \min\bigl(\gamma(E\cap\{A\ne B\}),
           \gamma(E^c\cap\{A\ne B\})\bigr)\ge\delta>0.            \tag{DD5}
\]
Then the cost in (DD4) has a positive lower bound depending
only on these compact sets and \(\delta\), uniformly in the
ambient dimension.

Otherwise take a sequence with vanishing cost. Its prescribed
masses and correlations converge along a subsequence. Each
complementary patch is a near-maximizer, since the two
individual deficits are nonnegative. Apply dimension-free
near-maximizer compactness to each of the two outputs and
apply fixed-plane compactness to each input in its own active
coordinates. The approximating cylindrical outputs add at
most four dimensions to the four dimensions of the input
planes. After padding and rotating a common span into
\(\mathbb R^8\), pass to convergent frames for the four planes.
A fixed planar indicator tuple pulled back by convergent
orthonormal frames converges in Gaussian label distance,
by approximation in \(L^1(\gamma_2)\) by bounded continuous
functions and the fixed Gaussian marginal in each frame.
The original two inputs and two output limits are therefore
actual maximizers at the limiting masses and correlation,
and the indicator identity passes to them. Each input
still has full intrinsic dimension two, so the projector
bound preserves distinct planes. The lower bounds on
changed masses make either patched limit have positive
distance from both input limits. This contradicts (DD1).
[Dimension-free near-maximizer compactness with moving parameters](GAUSSIAN_NEAR_MAXIMIZER_COMPACTNESS.md)
[Finite-dimensional compression for the simplex case](THREE_CELL_COMPLEMENTARY_PATCH_RIGIDITY.md#strict-cost-of-a-nontrivial-cut)


## The difference from a simplex cannot satisfy such an equation

Let \(A\) be a balanced actual maximizer in the plane and let
\(T\) be any centered regular three-sector partition in that
plane. The reference \(T\) need not be globally optimal at
the chosen correlation. If a finite product \(P(D)\) of
nonzero directional derivatives satisfies
\[
                         P(D)(1_A-1_T)=0,                   \tag{DD6}
\]
then \(A=T\) or \(A=-T\) up to null sets, where \(-T\)
denotes the spatial reflection through zero with the same
label names. Both possibilities can occur in (DD6).

First every regular arc of \(A\) is straight. Away from the
reference rays this follows from the leading normal delta test
used for (DD2). An analytic arc with an open part on a
reference line lies on that line throughout its connected
regular continuation. Any other arc has an open part disjoint
from the reference rays and is therefore straight by the
delta test and analytic continuation. Thus at each finite
point the jump sets of \(A,T\) are locally finite unions of
straight segments.

For a direction not among the factors of \(P(D)\), the
jump vector of \(1_A-1_T\) must vanish at any regular
point of a line in that direction, again by the leading
normal delta test. Near a finite vertex, discard such
segments, since their vector jumps cancel. The remaining
directions are characteristic for \(P(D)\). Apply (DD2)
componentwise to the difference, treating the jump entries
as signed real coefficients. It shows equality of the two
jump vectors on the opposite rays in each characteristic
direction, with the same chosen normal on both. For all
noncharacteristic directions the difference has zero jump
on each ray separately. Thus the jump of \(1_A-1_T\)
cannot have any unpaired endpoint at a finite point.

At any point other than zero, \(T\) is locally constant
or has a single full straight interface. Its own opposite
jump vectors therefore agree. The endpoint condition
for the difference shows that the same is true for \(A\)
at such points. There can be no genuine triple there.
A noncollinear triple has unpaired rays, and a collinear
triple cannot be made of distinct straight rays all
tangent to one line.

If \(A\) has no genuine triple at zero, the same local
regularity makes every interface a complete line.
Nonparallel lines would cross, so the partition would
be cylindrical, contradicting its active dimension two.
Consequently \(A\) has a noncollinear genuine triple at
zero. Its three straight rays continue to infinity
because no other point can end them. Their sector angles
are all below \(\pi\), so every affine line not through
zero crosses at least one of these rays at a positive
distance from zero. An additional interface would have
to be such a complete line and would create a prohibited
crossing. The three cells of \(A\) are therefore just
three sectors at the origin. Balance makes each angle
\(120^\circ\).

Finally apply the endpoint condition at zero to each
label. For a \(120^\circ\) wedge each of its two
unoriented boundary lines carries a jump on only one
half-ray. Its opposite-ray jump difference is nonzero
on exactly those two directions. Thus \(A_i\) and
\(T_i\) have the same two unoriented boundary lines,
and the opening angle forces \(A_i=T_i\) or
\(A_i=-T_i\). Choices with different signs for two
labels give overlapping sectors of positive angle,
because \(T_i\) and \(-T_j\) have center directions
separated by \(60^\circ\) when \(i\ne j\).
All signs must be the same. The antipodal example
does satisfy (DD6) for a suitable product: if \(T_i\)
is defined by \(L_1,L_2>0\) for linear forms
\(L_1,L_2\), then up to null sets
\[
 1_{T_i}-1_{-T_i}=1_{\{L_1>0\}}+1_{\{L_2>0\}}-1 .
\]
A product with one tangent derivative for each
reference wall annihilates the entire vector difference.
