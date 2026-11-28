# One translated wedge determines a balanced maximizing partition

Fix \(0<\rho<1\). Suppose one cell of an actual balanced
planar three-cell Gaussian noise-stability maximizer is
a translated pointed wedge, up to null sets. Then the
partition is the centered regular simplex. The wedge
need not have a prescribed aperture or a vertex at the
Gaussian origin. On every compact positive-correlation
band, each cell of any nonsimplex balanced maximizer
stays a positive Gaussian symmetric-difference distance
from every translated wedge, including degenerate
intersections of at most two half-planes.

## The vertex and the two pressure inequalities

Write the wedge of label \(i\) as \(K=p+W\), where
\[
                    W=\{x:n_1\cdot x\ge0,\ n_2\cdot x\ge0\}
\]
for nonparallel inward unit normals \(n_1,n_2\). Its
supporting offsets are \(c_\ell=n_\ell\cdot p\).
The two outward boundary rays issue from \(p\).
Write \(U_a=T_\rho1_{A_a}\) for the score of label
\(a\), with centered prices \(\sum_a\beta_a=0\).

Each boundary ray has an opposing label with an adjusted
score equal to that of label \(i\) along the full supporting
line. Indeed the winning-score rule supplies at least
one of the other two labels tied at each boundary point,
and a finite union of discrete zero sets of analytic
restrictions could not cover a ray. The supporting
ray has a unique such label: otherwise all scores would
have constant differences on the line, making \(U_i\)
constant there, contrary to its limits below.
[Stationary score limits on a convex supporting ray](THREE_CELL_CONVEX_SUPPORTING_RAYS.md#the-two-score-limits-on-a-supporting-line)

Let the opposing label on line \(\ell\) be \(j_\ell\).
Along the forward ray the other wedge inequality recedes
to infinity in the conditional Gaussian. Thus \(U_i\)
tends to
\[
         q_\ell=\Phi\!\left(-\frac{(1-\rho)c_\ell}
                                      {\sqrt{1-\rho^2}}\right).
\]
In the reverse direction \(U_i\to0\). Since
\(U_{j_\ell}-U_i=\beta_{j_\ell}-\beta_i\) along that
full line, the reverse limit and \(0\le U_{j_\ell}\le1-U_i\)
give
\[
              0\le\beta_{j_\ell}-\beta_i\le1-2q_\ell .
                                                               \tag{TW1}
\]
In particular \(q_\ell\le1/2\), so \(c_\ell\ge0\) and
\(p\in W\).

The two opposing labels are distinct. If the same label
opposed both rays, its adjusted contrast with \(i\)
would vanish on the two nonparallel supporting lines.
Its gradient would vanish at \(p\), whereas an actual
jump arc of that pair approaches \(p\). This contradicts
the reciprocal-gradient bound (G20). In verifying that
the open arcs are actual jumps near \(p\), use uniqueness
of the opposing label on each full line: the remaining
score ties there only at isolated points, and label \(i\)
occupies the interior side of the wedge. The two line
identities now give a three-score tie at \(p\).
By (TW1), \(\beta_i\) is a minimum price, and hence
\[
                          U_i(p)=\tfrac13+\beta_i\le\tfrac13.
                                                               \tag{TW2}
\]
[Nonvanishing pair gradients at finite endpoints, G20](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#regular-interface-arcs-and-isolated-triple-junctions)

## The conditional wedge score centers the vertex

Put \(t=(1-\rho)/\sqrt{1-\rho^2}=\sqrt{(1-\rho)/(1+\rho)}\).
The conditional source at observation \(p\) has
mean \(\rho p\), and the cone property gives
\[
                       U_i(p)=\gamma_2(W+tp).
\]
Since \(p\in W\) and \(0<t<1\), the set \(W+tp\)
contains \(W+p\), whose Gaussian measure is \(1/3\).
If \(p\ne0\), at least one \(c_\ell\) is strictly
positive because the normals are independent, and
the inclusion is strict on an open set. Thus (TW2)
forces \(p=0\). Balance makes the centered opening
\(120^\circ\), and one-cell simplex rigidity gives
the full simplex.
[One simplex cell determines the full partition, G16](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#one-simplex-cell-determines-the-whole-partition)

There is also a useful version for a potentially
cut-off tip. Suppose a convex cell has two nonparallel
supporting boundary half-rays with distinct opposing
labels as above, and the intersection of their
supporting lines lies in the central maximal-tie
gradient-separation regime. The supporting line
equations make that intersection a triple score tie.
Both contrasts against label \(i\) have nonzero
gradients normal to their respective lines, so those
gradients are independent. There are points
arbitrarily near their intersection where label \(i\)
strictly wins both comparisons. Thus the intersection
belongs to the closed convex representative \(K\).
Each outward ray direction is a recession direction
of \(K\), by taking the limit of line segments from
an arbitrary point in \(K\) to points tending to
infinity along the ray. The two recession directions
span the supporting wedge. Starting at their
intersection fills the whole wedge, so \(K\) is
that wedge, and the preceding proof applies.
[Central maximal-tie gradient separation](THREE_CELL_CENTRAL_INTERFACE_COMPACTNESS.md#uniform-separation-at-every-maximal-tie)

## Uniform distance from affine wedges

Suppose the asserted separation failed on a compact
correlation band. Take actual nonsimplex maximizers
whose one cell approaches a sequence of intersections
of at most two closed affine half-planes in Gaussian
symmetric difference. Fix the label by a subsequence
and pass to an actual maximizer limit by strong
compactness. For each half-plane take a subsequence
of the unit normals and the offsets in the extended
real line. Its indicator converges almost everywhere
to a half-plane indicator or to a constant, and so
does the two-factor intersection. The limiting cell
is thus an intersection of at most two half-planes
up to null sets. Empty or full cells are excluded
by balance. A half-plane or a parallel-sided strip
is excluded by full directional rank of each
maximizing cell. The remaining possibility is the
pointed wedge, already proved to force the simplex.
Uniform local simplex isolation on the band prevents
that limit for nonsimplex maximizers.
[Full directional rank of a maximizing cell](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#each-cell-uses-both-active-directions)
[Strong compactness and uniform local comparison](POSITIVE_CORRELATION_CONTINUATION.md)
