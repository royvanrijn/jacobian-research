# Rigidity of connected three-cell networks with algebraic ends

Let \(A\) be a balanced planar Gaussian maximizer at \(0<r<1\).
Suppose each of the three strict-winner sets is connected and
the closure of its actual jump interfaces is semialgebraic.
If all three cells are unbounded, \(A\) is the regular simplex.
Exactly two unbounded cells are impossible. Thus a potential
nonsimplex maximizer under this reduction must have a bounded
interface network with just one unbounded cell. Its graph is a
compact theta with two trivalent junctions. Its three arcs cannot
all be graphs over the same axis, and each arc is curved.
The exterior algebraic-arc obstruction also excludes this compact
case when the jump network is semialgebraic.
[Connected semialgebraic maximizers](THREE_CELL_CIRCULAR_INTERFACE_OBSTRUCTION.md#connected-semialgebraic-maximizers)

The proof uses the analytic interface structure of actual
maximizers and their exact Gaussian conditional scores on
unbounded ends. The same argument applies if the jump network
compactifies to a finite graph and each end has a convergent
Puiseux parametrization at infinity. Semialgebraicity provides
these hypotheses by triangulation on the sphere and the
one-variable Puiseux expansion after inverting the longitudinal
coordinate.
[Łojasiewicz, Theorem 3](https://www.numdam.org/item/ASNSP_1964_3_18_4_449_0.pdf)
[Mondal–Netzer, Section 4.2](https://arxiv.org/pdf/1305.1215)
[Coste, convergence of the expansions in Section 1.5](https://perso.univ-rennes1.fr/michel.coste/polyens/RASroot.pdf)
[Interface regularity and complete adjacency](THREE_CELL_GLOBAL_OBSTRUCTIONS.md)

## Counting vertices and ends

Take the closure of all actual jump arcs, including their
endpoints, and ignore passive equality lines where the label
does not jump. By local analytic stratification, each face
of this graph has constant label almost everywhere: the
distributional derivatives of the label indicators are
supported on the jump graph. Each face has a strict-winner
point, because the score ties are null. A connected strict
winner set cannot visit two faces. Thus in the connected
case the graph has exactly three faces.

At an accumulating actual arc the pair gradient cannot
vanish by (G20). Hence finite endpoints are absent except
at genuine triples, where the local analytic description in
the maximizer analysis gives valence three or four. Insert
degree-two vertices on any closed regular components.
There is at least one genuine triple. Otherwise each jump
component is a smooth closed curve or a proper embedded
line, because its analytic level continues at each finite
accumulation point and there is no vertex at which it
can stop. Any such \(i,j\) interface separates the plane.
The connected strict winner sets of \(i\) and \(j\) lie on
opposite sides. The connected third label lies on one
side and cannot share an interface with the opposite
label off the separating curve. It cannot have a regular
interface on that curve either, since every point of its
regular arcs has only the two \(i,j\) labels on its sides.
This contradicts complete adjacency.

If at least two cells are unbounded, include their
point at infinity as a graph vertex. Near infinity the
ends are finite in number and each meets a sufficiently
large generic circle once, by semialgebraic monotonicity
or by the Puiseux parametrization. Upon one-point
compactification the graph is a finite graph in the
sphere. If \(C\) is its number of components, Euler's
formula with three faces gives
\[
                 \sum_z(\deg z-2)=4-2C.                    \tag{SN1}
\]
With all three cells unbounded, each large circle must
meet all three strict winner sets, since a connected
unbounded set meets every sufficiently large radius.
Thus \(\deg(\infty)\ge3\). At least one finite triple
gives another positive contribution. Equation (SN1)
leaves exactly one finite triple, of degree three,
with three ends at infinity and \(C=1\). A loop at
the finite triple would trap a bounded winner region,
contrary to connectivity and unboundedness of every
label. It follows that all three edges connect the
triple to infinity. The pair types are distinct at
the trivalent vertex, by the noncollinear or the
one-sided cusp alternative of the local classification.

If only cells 1 and 2 are unbounded, each large circle
sees just those labels, both of them. Its actual
transitions alternate in label, so \(\deg(\infty)\)
is positive and even. Equation (SN1) with a finite
triple now gives \(\deg(\infty)=2\), \(C=1\), and
either two finite trivalent triples or one four-arm
cusp. Both ends have pair type \(1,2\). If their
directions at infinity coincide, the two oriented
transverse switches have opposite signs.

## No sublinear divergence or decaying Puiseux term

Fix a limiting unit direction of one or more ends and
use coordinates in which it is the positive \(x\)-axis.
The Puiseux parametrization makes these ends graphs
\(y=h_j(x)=o(x)\) for large \(x\). Other directions
stay at a positive angle and contribute only an
exponentially small probability when conditioning on
a point of this group. This follows by separating
the direction with a fixed angular sector, and taking
the conditional Gaussian within a fixed small multiple
of \(x\) of its mean.

Suppose some graph is unbounded transversely. Take the
largest Puiseux exponent \(p\in(0,1)\) of such graphs
and write \(h_j(x)=c_j x^p+o(x^p)\) for the whole group.
Observing at \((x,h_j(x))\), the conditional source
horizontal coordinate is \(rx+O(1)\), and its
transverse coordinate in units of \(x^p\) tends to
\(r c_j\). The neighboring source walls in these
units tend to \(r^p c_l\). Thus if \(r^{1-p}c_j\)
is distinct from all coefficients \(c_l\), a single
label has probability tending to one. This is
incompatible with two distinct adjusted winners at
the observed wall. Indeed every price difference is
strictly smaller than one in absolute value, by
positive conditional probability of every cell at
every finite point and the winning equation on
the positive-mass cells. The finite set of coefficients
would therefore be closed under multiplication by
\(r^{1-p}\), which forces all to vanish. This
contradicts the choice of exponent. Every graph in
the group consequently has a finite limit \(b_j\).

Order the graphs from bottom to top, allowing ties
among their limiting heights \(b_j\). Write \(d_j\)
for the vector of the lower label minus the upper label
at switch \(j\). In the three-unbounded-cell case the
pair types are distinct at the single triple. With
two or three of those ends in one direction, the
successive strip labels traverse distinct edges of
the triangle on the labels. Thus \(d_j\cdot d_l=-1\)
for distinct switches in the group. In the
two-unbounded-cell case the possible two switches
in a single direction have \(d_1\cdot d_2=-2\).

For the limiting parallel-strip profile put
\[
 K_{jl}=(d_j\cdot d_l)p_\sigma(b_l-r b_j),
 \qquad a_j=r\sum_l K_{jl},\qquad \sigma=\sqrt{1-r^2},
\]
where \(p_\sigma\) is the centered normal density
of variance \(\sigma^2\). The numbers \(a_j\) are
strictly positive, as follows.

If all offsets are distinct, the limiting adjusted
lower-minus-upper contrast on switch \(j\) must
cross from nonnegative to nonpositive. This follows
by testing the actual winners on both sides at fixed
transverse coordinates and sending \(x\to\infty\).
The derivative of that contrast at transverse
coordinate \(y\), divided by the positive factor
\(r p_\sigma(b_j-r y)\), is
\[
 -2+\sum_{l\ne j}(-d_j\cdot d_l)
                          \frac{p_\sigma(b_l-r y)}
                               {p_\sigma(b_j-r y)} .
\]
The ratios are exponentials of affine functions of
\(y\) with nonzero slope. With at least one other
switch the bracket is strictly convex. If it
vanishes at the switch with nonzero derivative, the
contrast has a strict local extremum, inconsistent
with the sign crossing. If its derivative also
vanishes there, the contrast is locally increasing,
also inconsistent. With one switch its derivative
is immediately negative. The contrast is analytic
and nonconstant, so the sign test concerns genuine
crossing even though winner inequalities hold only
weakly in the limit. Evaluating its negative
derivative at \(b_j\) yields \(a_j>0\).

For two coincident offsets with three distinct strip
labels, the middle label disappears in the limit,
and each of its pair contrasts has negative derivative
\(-r p_\sigma((1-r)b_j)\). If there are only two
alternating labels and the two offsets coincide,
the limiting source is a single label, excluded
by the winner condition as above. For a three-switch
cycle with two coincident offsets, replace the
two collapsed switches by their sum, which is the
incidence of the surviving transition. The remaining
two switches then form a two-label alternating
profile with distinct offsets, for which the
crossing argument just proved strict negative
derivatives. At the collapsed offset the derivative
of each original pair contrast is half the derivative
of this surviving contrast, by the incidence dot
products. At the other offset the third contrast
is unchanged. All three offsets coinciding would
again leave a deterministic source label. These
are all cases.

If any tail differs from its constant height, let
\(x^{-\alpha}\), \(\alpha>0\), be the first nonzero
decaying power among the graphs in this group and
write \(h_j(x)=b_j+v_j x^{-\alpha}+o(x^{-\alpha})\).
The exact conditional formula on the graph \(j\)
has coefficient of this power equal to
\[
                    -a_j v_j+r^{-\alpha}(Kv)_j .           \tag{SN2}
\]
In fact, condition on the horizontal source \(S\)
and integrate the vertical normal. Moving source
switch \(l\) contributes the normal density at its
limiting height times \(d_j\cdot d_l\). The
displacement is \(r^{-\alpha}v_l x^{-\alpha}
+o(x^{-\alpha})\) uniformly on
\(|S-rx|\le x^{1/2}\). Moving the observed switch
changes its conditional mean by \(r v_j
x^{-\alpha}\), which gives the other term.
Outside the specified source interval the
Gaussian tail is smaller than every inverse
power of \(x\). This justifies the expansion in
expectation and applies also to coincident offsets.
Stationarity makes (SN2) vanish.

The matrix \(\operatorname{diag}(a_j)^{-1}K\) has
row sum \(1/r\) and nonpositive off-diagonal entries.
Choose a component of \(v\) of largest absolute
value, and flip the sign of the whole coefficient
vector if needed to make this component positive.
In that row its image is at least \(v_j/r\). But
(SN2) requires its image to be \(r^\alpha v_j\),
a contradiction. Therefore all Puiseux tails in
every group are exactly affine. By analytic continuation
of each regular arm, a tail on a straight line
extends on that line along the entire connected
arm. One may use the vanishing of the linear
coordinate normal to the line on its analytic
parameterizations. At finite limits with only the
associated pair active, (G20) keeps their gradient
nonzero.

## Consequence for the bounded cells

When there are three unbounded cells the three
arms are now rays from the same finite triple.
The partition is invariant under dilations about
that point. Translated cone rigidity identifies it
as the centered regular simplex.

With exactly two unbounded cells, analyticity
extends the stationary \(1,2\) score equality
on the line of either straight exterior arm to
its entire supporting line. If the two directions
are not opposite, the reverse direction of one
arm meets no end. The conditional source there has
a deterministic label at infinity, making the
limit of \(U_1-U_2\) equal to \(1\) or \(-1\).
This contradicts its constant stationary pressure
on that line. When the directions are opposite,
take a common transverse coordinate. Their limiting
upper labels agree, as do their lower labels,
since large semicircles between the opposite
end directions traverse no interface. If their
heights are \(b_1,b_2\), the two pressure limits
are the same strictly monotone function
\(2\Phi((1-r)b/\sigma)-1\) of the height, with
the common choice of sign of the contrast.
Thus \(b_1=b_2\).

The whole \(1,2\) interface lies on this one
line. To see it, in the connected finite graph
its two ends must meet genuine finite junctions:
an entire infinite component with no junction
would be a separate graph component. In the
two-triple case the two ends use the sole \(1,2\)
incident edge of each triple. In the four-arm
cusp case a pair type has at most two
incidences, both used by the two ends. This
also rules out any further bounded component
of this pair. A translation tangent to the
common line has zero normal speed on the
\(1,2\) interface. Orient the other two
interfaces outward from cell 3. Their signed
interaction kernel is strictly positive in
these orientations. The first extra-coordinate
form (G1) is zero on the translation. Taking
absolute values would strictly lower the
form if its speed changed sign on positive
measure sets. Hence the Gaussian distributional
derivative of \(1_{A_3}\) along this translation
has one sign. A positive-mass indicator with
one-sided monotone slices has unbounded support
in that direction. This contradicts the assumed
bounded third cell.
[Translated cone rigidity and translation sign test](THREE_CELL_GLOBAL_OBSTRUCTIONS.md)

## The remaining compact network

If only one cell is unbounded, every active interface
lies in the bounded closures of the two bounded cells.
There is no vertex at infinity, so (SN1) again gives
\(C=1\) and a total degree surplus of two, because
there is a finite triple and thus positive surplus.
A single four-arm cusp would have just two pair
types at the vertex by the local classification of
the thin middle label, and there could be no third
pair type elsewhere in the connected graph. This
contradicts complete adjacency. Therefore there are
two trivalent vertices, each with one incidence of
each pair type. No same-type loop can return to one
vertex, since it would require two incidences there.
Hence the three compact pair arcs join the two
junctions.

No whole pair arc can be straight. Translation tangent to
it has zero speed there, and the other two interfaces,
oriented outward from their common label, have strictly
positive mutual interaction kernel. The absolute-value
argument used above forces that label to be monotone along
the translation. Either the label or its complement is
bounded and has positive mass, which is impossible for
an indicator with nontrivial one-sided monotone slices.

Suppose these arcs are graphs over a common axis.
Their projections have the same two endpoints, and
the graphs are ordered between them. Outward along
a transverse line the labels then occur in a
three-switch cycle: the single unbounded cell, the
two bounded cells in order, and the unbounded cell
again. Orient each normal in the positive transverse
direction. The associated incidences have mutual
dot product \(-1\). Translation in that same
direction has nonnegative speed on every arc and
positive speed on a set of positive surface
measure on each. Its first extra-coordinate form
is zero by the translation identity. Flipping its
sign on one whole arc preserves the weighted local
term and the three self-interactions. It strictly
increases the interaction with the other two arcs
because their cross kernels have negative sign.
This violates (G1). Thus at least one arm
backtracks in any projection that separates the
junctions.

A compact theta with backtracking and nonalgebraic exterior arcs
remains open under the finite-graph hypotheses. The cases without
connected winner sets or the stated finite graph and end expansions
also remain open.
