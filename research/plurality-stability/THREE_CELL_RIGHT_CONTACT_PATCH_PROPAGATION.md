# Complementary patches in the simplex plane at a right contact

Let \(r\in(0,1)\) be a correlation at which the balanced regular
three-cell simplex is globally optimal. Suppose distinct correlations
of simplex optimality tend to \(r\) from above. If \(A\) is any balanced
planar maximizer at \(r\), one can orient and label a simplex \(T\)
so that every nontrivial balanced complementary patch of \(A,T\)
has strictly smaller stability. In fact for almost every rotation
of one appropriate cyclic labeling of \(T\) the conclusion holds.

In indicator notation, put \(d=1_A-1_T\), \(d_E=1_Ed\).
For any measurable \(E\) with \(\int d_E\,d\gamma_2=0\) and
\(0<\gamma_2(E\cap\{A\ne T\})<\gamma_2(A\ne T)\), the claim is
\[
     \langle d_E,T_r d_{E^c}\rangle>0 .                       \tag{RP1}
\]
There is no assumption on the connectivity of \(A\) or on its
interface at the origin. The right-accumulation hypothesis enters
only in the final derivative comparison. For arbitrary \(r\)
at which \(T\) is optimal, the geometric argument shows that
an exact complementary recombination into two other maximizers
must give one of them every outward simplex ray as an interface
up to arclength-null sets.

## Tangent partitions and angular interval endpoints

For an actual planar maximizer \(F\), use the open representatives
given by the adjusted analytic scores. At a point \(p\), retain
only labels whose strict winner regions accumulate at \(p\).
The other labels are absent in a neighborhood. If only one label
remains the partition is constant in that neighborhood. With two
labels the gradient of their score contrast is nonzero at \(p\),
by (G20) and the local analytic jump stratification. With three,
the three gradients are pairwise distinct by the same regularity
result. Hence on every closed cone of directions where one of
these gradients has a strict linear advantage, the partition
agrees with that linear winner sufficiently near \(p\).
This follows by the Taylor expansion of the scores. It gives
a tangent partition in local label measure about \(p\).

The planar tangent partition is either constant, binary with
two opposite half-planes, or a three-sector partition whose
positive sector angles are all strictly below \(\pi\). In the
binary case, a third label may occupy a cusp of zero angle.
In fact the positive-angle regions are the normal cones of
the extreme points of the convex hull of the retained gradients.
A tangent consisting of one label entails local constancy of
the original partition.
[Interface regularity and distinct gradients, G20](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#regular-interface-arcs-and-isolated-triple-junctions)

We use a fact about the circle. Suppose the sum of the indicators
of two connected angular intervals is also the sum for two
other connected angular intervals, ignoring endpoints. Suppose
the original endpoints are distinct, their union is a proper
subset of the circle, and both original intervals have positive
length. If they overlap in positive length, both new intervals
must cover the overlap. The union can be cut open to a line
interval. Since the new intervals have overlapping interiors,
their endpoints cannot cancel in the sum. The four original
endpoints, with their signs, must therefore be exactly the four
new endpoints. If the original intervals are disjoint, they
are separated by a gap and each connected component of their
union must be an entire new interval. In either case there
are two nonempty new intervals with no hidden endpoints.

One consequence is that a semicircle \(I\) and a
\(120^\circ\) interval \(J\) with distinct endpoints can be
recombined into two intervals both shorter than a semicircle
only if \(J\) is contained in \(I\). The disjoint case follows
immediately. If they overlap without containment, the endpoint
matching gives either the original two intervals or their
union and intersection. The union is longer than a semicircle.

## An orientation that preserves all three rays at the origin

Suppose \(B,C\) are actual maximizers satisfying
\[
                         1_A+1_T=1_B+1_C.                    \tag{RP2}
\]
The tangent partitions at any point satisfy the same identity.
Choose \(T\) so none of its switching rays agree with a
switching ray of the tangent of \(A\) at the origin, and so
its supporting lines are distinct from the binary tangent
line if \(A\) has one. If the tangent of \(A\) has three labels,
choose the reverse cyclic label order for \(T\). These exclude
only finitely many rotations within the indicated cyclic order.

If the tangent of \(A\) is constant with label \(i\), each of
the other two labels is supported in its \(120^\circ\) sector
of \(T\). Any output tangent using such a label must have
three positive sectors, since a binary tangent would require
a half-plane for that label. If both outputs did so, their
combined label-\(i\) angle would be less than \(2\pi\), whereas
the input angle is \(2\pi+2\pi/3\). Thus one output tangent is
\(T\) and the other is constant \(i\).

If the tangent of \(A\) is binary with labels \(i,j\), the two
source intervals for each of these labels have total angle
\(5\pi/3\). They satisfy the interval fact, so each output
uses both labels. If both outputs also use the third label,
their \(i,j\) intervals are all shorter than a semicircle.
The containment consequence would place \(T_i,T_j\) in the
respective input semicircles. Their common ray would have to
lie on the dividing line of the semicircles, contrary to
the orientation choice. Therefore at least one output is
binary. Its two opposite switching rays have to be visible
source endpoints by the interval fact. The only antipodal
pair of source endpoints is the line of the binary input,
since no reference line is that line. The label ordering of
this output agrees with the input: at least one of \(T_i,T_j\)
overlaps its like-labeled input semicircle, since otherwise
their common ray would again lie on the input line. This
output is therefore the binary tangent of \(A\), and the
other is the reference tripod.

If the tangent of \(A\) has three positive sectors, the
interval fact for each label gives three positive output
sectors in both \(B,C\), with no hidden endpoints. Every
output switching ray is then an input switching ray carrying
its original oriented label transition. To see this, each
input ray changes one label into another, and the four
endpoints for each label are visible with no cancellation.
Both outputs could change at one input ray only by canceling
a third label's endpoints there. In a three-sector output
one finds all three oriented transitions of one cyclic
order. Since \(A,T\) have opposite cyclic orders, an output
cannot mix their transitions. Thus the two output tangents
are again the two input tangents.

In all cases name the output having the tripod tangent \(B\).
The identification of tangents upgrades to \(B=T,\ C=A\)
in a full neighborhood of zero, up to null sets. Indeed the
input switching directions are disjoint. In a conic
neighborhood of any direction, at least one input and its
corresponding output have the same unique tangent label.
They agree there sufficiently close to zero, by the tangent
construction. Equation (RP2) makes the other two agree there
as well. Finitely many such conic neighborhoods cover the
circle, so the agreement holds on a disk. In particular each
adjusted pair contrast of \(B\) vanishes on the full supporting
line of its reference ray, by analyticity along that line.

## Propagation along the rays

Further restrict the rotation so that no adjusted pair score
of \(A\) is constant on a reference supporting line, and
every nonzero reference ray point is either a local constant
point of \(A\) or a regular transverse crossing of an
\(A\) interface. Both requirements hold after excluding
a countable set of angles. For the first, the first nonzero
homogeneous Taylor term of a pair contrast after subtracting
its value at zero has only finitely many real zero directions
on the circle. The contrast is nonconstant by Gaussian
convolution injectivity and the positive distinct cells.
For the second, genuine triple points are locally finite
by (G20). Away from them the regular pair interface has
countably many analytic arc charts. On each chart the points
whose tangent is parallel to the radial direction are
discrete unless the whole chart is in a radial line, since
the signed determinant of the position and tangent is
analytic. A regular radial arc supplies just a single
supporting line to exclude. Where \(A\) has only one
locally accessible label it is locally constant regardless
of scores of absent labels.

The partition \(C\) cannot carry any open reference \(i/j\)
ray segment with that pair as its interface. This would
give a constant \(i/j\) score contrast for \(C\) on the full
supporting line. The contrasts of \(B,T\) are already
constant on the line, so by (RP2) the contrast of \(A\)
would also be constant, contrary to the rotation choice.

At a regular transverse crossing of the two input interfaces,
the input tangents are binary with different dividing lines.
Their output tangents must again be the input tangents in
some order. If the pairs are the same, apply the circular
interval fact to the two overlapping input semicircles for
each label. Both outputs are binary and their respective
boundary lines are the only antipodal pairs of source
endpoints. The overlap of the like-labeled semicircles fixes
their orientation. If the pairs differ, they share a label
whose total input angle is \(2\pi\), with a positive angle
outside its support. Neither output can be constant. Both
cannot have three sectors, since their combined angle for
the common label would be smaller than \(2\pi\). A binary
output must use the common label, since otherwise the other
output would devote its entire circle to it. Its second
label is supported on the single semicircle contributed by
one input, and it fills that semicircle. Thus the binary
output is that input, and the remaining output is the
other input. The same conic-neighborhood argument as at the
origin upgrades this identification to local agreement.
Since \(C\) cannot carry the reference pair interface, at
each such crossing \(B=T\) locally.

Consider a point on the reference \(i/j\) ray where \(A\)
is locally constant. If the third label strictly won for
\(B\) there, \(B\) would locally be that constant label.
It would have to be the constant label of \(A\), by (RP2),
and \(C\) would equal \(T\) locally, which was excluded.
On the ray the third adjusted score of \(B\) can equal the
tied pair only at isolated points. Indeed the contrast on
the line is analytic and there is a strict pair win near
the origin. The normal derivative of the tied contrast is
also nonzero except at isolated points, because it is
analytic on the line and nonzero near the origin. Off these
isolated exceptions the two tied labels are actual winners
on opposite sides. Their ordering has to agree with \(T\)
by (RP2), since giving the opposite order would require
the constant \(A\) to supply two different labels. Thus
\(B=T\) locally at such points of the ray as well.

We conclude that \(B\) carries all three outward reference
rays as their original interfaces, up to arclength-null
sets. This conclusion has used no inequality about the
correlation derivative.

## Strict comparison at right contact

The unobstructed rays give a strict derivative test for
\(B\) unless \(B=T\) almost everywhere. In fact (OA2) is
a positive sum of Gaussian pair interactions of surface
measures, and \(S'_r(B)\ge\mathcal J_r(B)\) by the first
extra-coordinate second variation. Since each reference
surface measure is dominated by that of \(B\),
\(\mathcal J_r(B)\ge\mathcal J_r(T)=S'_r(T)\).
If \(B\ne T\), some surface measure has positive mass off
its reference rays. Otherwise the distributional gradients
of the cells would vanish in the interiors of the reference
sectors, and agreement near zero would give \(B=T\).
Strict positivity of the Mehler kernel therefore makes
\(\mathcal J_r(B)>S'_r(T)\).
[Surface test for the correlation derivative, OA2](THREE_CELL_ORIGIN_ARM_SURFACE_COMPARISON.md#an-unoriented-surface-bound-from-stability)

On the other hand \(s\mapsto S_s(T)-S_s(B)\) vanishes at
\(r\) and is nonnegative at the distinct simplex-optimal
correlations approaching \(r\) from above. Differentiability
at \(r\) gives \(S'_r(B)\le S'_r(T)\). Thus \(B=T\) and
(RP2) gives \(C=A\). For a balanced complementary patch
of \(A,T\), the two outputs satisfy (RP2), and the sum of
their deficits from the common maximum is
\(2\langle d_E,T_r d_{E^c}\rangle\).
Each deficit is nonnegative, and equality would make
both outputs actual maximizers. The strictness statement
(RP1) follows.

For fixed \(r,A,T\) satisfying the hypotheses and this
rotation choice, the strict loss is uniform over cuts
whose two changed masses are each bounded below by
\(\varepsilon>0\). Otherwise the two outputs of a sequence
of cuts would each approach the optimum in the fixed
plane. Strong near-maximizer compactness gives ordinary
maximizer limits, preserving (RP2) and the distances
from both input partitions. This contradicts the proven
reconstruction.
[Strong near-maximizer compactness](POSITIVE_CORRELATION_CONTINUATION.md)

The comparison constrains recombinations of tied maximizers.
Existence of a nonsimplex tied maximizer with this strict
comparison remains open.
