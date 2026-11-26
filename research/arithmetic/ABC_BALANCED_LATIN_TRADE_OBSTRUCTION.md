# Balanced label cancellations in dense endpoint graphs

In the complete-label reduction of the balanced \(abc\) endpoint,
each selected edge is an actual additive triple and each pair of
labels determines at most one edge. Both properties matter to the
small configurations that this reduction supplies. Pair uniqueness
and near-quadratic edge counts permit a graph with no balanced
six-cycle. Exact addition supplies a balanced eight-edge pattern
with simultaneous arithmetic parallelograms in the three letters.
[Complete-label reduction](ABC_DIFFUSE_FACTOR_GRAPH.md#the-complete-label-reduction-close-to-minimum-mass)

## Balanced cycles and matching supports

Let \(H\) be a simple linear tripartite \(3\)-uniform hypergraph
on parts \(A,B,C\). Linear means that two different edges share
at most one vertex. View \(H\) as a simple bipartite graph on
\(A,B\), with each edge \((a,b)\) carrying its \(C\)-label \(c\).
This projection has a proper edge coloring, since two equal
colors at \(a\) or at \(b\) would violate linearity. Call a simple
even cycle in the projection balanced if its two alternating
edge lists have the same multiset of \(C\)-labels.

A balanced four-cycle is impossible. Each edge of one alternating
list is adjacent to both edges of the other and therefore has no
partner of its own color.

A balanced six-cycle is exactly the projection of a \(3\times3\)
grid in \(H\) with the inherited parts. For each edge of the
first alternating list, its two adjacent edges have different
colors, so its matching color on the other list must occur on
the unique nonadjacent edge of that list. The three colors must
be distinct. If two first-list edges shared a color, each would
have a different nonadjacent partner of that color. One of
these partners is adjacent to the other first-list edge,
contradicting properness. Hence each alternating list lifts to
three pairwise disjoint hyperedges, with the same three vertices
in each part. Between the two lists every two edges intersect:
the adjacent pairs intersect on \(A\) or \(B\) and the
nonadjacent pairs on \(C\). This is the grid.

Conversely a grid consists of two sets of three pairwise disjoint
edges with the same nine-vertex support. Both sets project to
matchings on the same three labels in \(A\) and on the same
three labels in \(B\). Their edges in this projection are distinct
by linearity. The union is a simple bipartite graph of degree two
on six vertices. It is a single cycle, since every component would
need at least four vertices. Its alternating color multisets agree
because both lifted matchings cover the same \(C\)-vertices.
The same projection argument shows more generally that two
disjoint hypergraph matchings on the same vertex support give a union
of alternating cycles in each pair projection. The joint color
multisets agree across the union. Individual cycle color multisets
need not agree if that union is disconnected.

Any such matching equality gives a signed incidence cancellation.
Giving the first matching weight \(+1\) and the second weight
\(-1\) makes the sum of the incident weights zero at every vertex.
For any assignment of nonzero numerical factors to labels, the
products of those factors across the two edge lists agree, on each
of the three parts. This concerns factors determined by labels.
A cofactor that varies between edges at a given outside label
does not cancel by this argument.

## An obstruction from near-quadratic grid-free graphs

Kong, Shangguan and Wang, Theorem 1.1, construct for every
sufficiently large \(N\) a linear tripartite \(3\)-uniform
hypergraph on \(N\) vertices with at least
\[
                    N^2\exp\!\bigl(-C(\log N)^{2/3}\bigr)
\]
edges and with neither a loose triangle nor a \(3\times3\) grid.
They also prove for the same construction that \(e\) distinct
edges cannot span at most \(e+3\) vertices for every
\(3\le e\le11\), in Remark 1.2 and Appendix A.
[Kong, Shangguan and Wang, arXiv:2611.05261v2, Theorem 1.1 and Appendix A](https://arxiv.org/html/2611.05261v2#S1.SS1)

Putting \(N=\lfloor n^c\rfloor\) gives three part caps \(n^c\)
and edge exponent \(2c-o(1)\), with the same pair uniqueness
as the outside factor projection. All three parts also have
size \(n^{c-o(1)}\): the edge count is at most the product of
any two part sizes. This gives a combinatorial obstruction
to forcing balanced six-cycles from those exponents. The
absence of loose triangles also excludes triangles of the pair
projections assembled from three different edges.

The full range of forbidden configurations in the source
has two further consequences.

- There is no nonempty signed cancellation, with arbitrary
  nonzero weights on its support, supported on at most seven
  distinct edges. Every vertex of its support would have
  degree at least two, since an isolated incident edge has
  nonzero weight. With \(e\) support edges there are at most
  \(\lfloor3e/2\rfloor\) vertices, which is at most \(e+3\)
  for \(3\le e\le7\). Supports of size one or two are
  impossible by linearity.
- A signed cancellation supported on eight edges, if present,
  must use exactly twelve vertices, all of degree two, and
  no other edge can be supported on those twelve vertices.
  The vertex bound follows from \((11,8)\)-freeness and the
  absence of degree one. Any additional edge on the same
  vertices violates \((12,9)\)-freeness. In particular a
  balanced eight-cycle with four distinct colors, if present,
  has no additional hyperedge on its twelve labels.

The KSW construction uses finite-field labels and does not supply
primitive positive integer values satisfying the same edges with small
radicals. Its obstruction applies to use of the part exponents and
linear incidence alone.

## Short signed cancellation from common-label pairs

There is a complementary forcing statement. Let \(N\) bound
each of \(|A|,|B|,|C|\), and suppose \(H\) has
\(E=N^{2-o(1)}\) edges. Then there is a nonzero integer
incidence cancellation on at most twelve edges, with weights
of absolute value at most two.

Use an auxiliary simple bipartite graph with classes of
ordered distinct pairs \((a,a')\) from \(A\), and ordered
distinct pairs \((b',b)\) from \(B\). Connect these two
nodes whenever there are edges \((a,b,c)\) and
\((a',b',c)\) of \(H\) with the same \(c\). The center \(c\)
is determined by the node pair, since \((a,b)\) belongs
to at most one original edge. Writing \(\mathbf e_v\) for
a unit vertex vector, the difference of the two original
edge incidence vectors equals
\[
     (\mathbf e_a-\mathbf e_{a'})
        -(\mathbf e_{b'}-\mathbf e_b).
\]
Consequently a path, with the appropriate sign for each
traversed auxiliary edge, telescopes to the difference of
its two endpoint pair vectors.

The auxiliary edge count is
\[
        M=\sum_{c\in C}\deg(c)(\deg(c)-1)\ge E^2/N-E,
\]
and the graph has at most \(2N^2\) vertices. Remove vertices
of degree below \(M/(4N^2)\). At most \(M/2\) edges are
removed, so a nonempty subgraph has minimum degree
\(N^{1-o(1)}\). Starting at one auxiliary node, choose
three-step paths for which all the original labels at their
pair nodes and cherry centers are distinct, counted separately
in their disjoint vertex parts. Each earlier label forbids
at most two choices for the next node. For a prescribed
center the other pair is determined by linearity, and a
prescribed entry of that pair first determines the center
from its incident original edge, then determines the other
entry. Thus only \(O(1)\) neighbors are lost at each step.

There are at least \(N^{3-o(1)}\) such paths and at most
\(2N^2\) endpoints. Two distinct paths have the same endpoints.
Their signed incidence sums cancel by telescoping. The resulting
edge coefficient list is nonzero. Each individual path uses
distinct centers in \(C\), so its signed edge list determines
its pairs of original edges at those centers. The initial pair
node identifies the first cherry uniquely, since the path visits
no original label twice. Its other pair node is determined by
the order at the initial node. Remove the first cherry and
repeat. Thus equal signed lists with the same start would give
identical paths. Each path uses six edges, proving the support
and coefficient bounds.

This is the signed version, for a tripartite graph of the
specified density, of the common-center pair-graph argument
used to force a bounded even configuration.
[Dellamonica et al., *On even-degree subgraphs of linear hypergraphs*, §4](https://homepages.math.uic.edu/~mubayi/papers/tworeg-3g.pdf)

## Eight-edge supersaturation for actual additive values

Suppose now the vertices on each part are distinct integers,
and every retained edge satisfies \(a+b=c\). There is an
absolute constant \(K>0\) such that if \(E\ge K N^{7/4}\),
there are
\[
                             \gg E^{10}/N^{15}
                             \tag{BT1}
\]
ordered pairs of disjoint four-edge matchings supported on
the same twelve vertices. They have the following particular
arithmetic form, for nonzero integer parameters \(h,t,s\):
\[
 \begin{split}
  A_0&=\{a,a+h,a+t,a+t+h\},\\
  B_0&=\{b,b+h,b+s,b+s+h\},\\
  C_0&=\{a+b+h,a+b+h+t,a+b+h+s,a+b+h+t+s\}.
 \end{split}                                                       \tag{BT2}
\]
All four labels in each set are distinct. For each of the
four targets in \(C_0\), the two crossed summand pairs
obtained by putting the \(h\) on either summand are both
original edges. Explicitly, for \(i,j\in\{0,1\}\), they are
\[
 (a+it+h,\ b+js,\ a+b+h+it+js),\qquad
 (a+it,\ b+js+h,\ a+b+h+it+js).
\]
Assign sign \((-1)^{i+j}\) to the first edge in each cell
and the opposite sign to the second. This gives the two
matchings. Thus at \(E=N^{2-o(1)}\) there are at least
\(N^{5-o(1)}\) trades of this form. There are \(O(N^5)\)
choices of the five numerical inputs \((a,a+h,a+t,b,b+s)\),
so the exponent in this conclusion is sharp.

For the proof, the auxiliary edge from the two triples
with a common target satisfies \(a-a'=b'-b\). Each connected
component therefore lies in one fixed nonzero difference
\(h\), and has at most \(N\) nodes on either side. Each
auxiliary node has degree at most \(N\). With \(M\) as
above, prune again at \(\delta=M/(4N^2)\).
At least \(M/2\) edges survive and at least \(M/(2N)\)
left nodes remain incident to an edge. Fix one such left
node. It has at least \(\delta^2\) two-step walks back to
its own side. Exclude terminal nodes whose underlying
two \(A\)-labels intersect the initial pair. There are
at most three such terminal nodes, because the difference
within a pair is fixed and nonzero. Their walks account
for at most \(3N\) choices. Since \(M\gg E^2/N\), the
assumption with sufficiently large \(K\) gives
\(\delta^2\gg N\). At least \(\delta^2/2\) walks remain.

For each allowed terminal node let \(g\) be its number of
intermediary \(B\)-pairs. Given the first intermediary pair,
there are only \(O(1)\) choices for a second intermediary
that reuse any of its \(B\)-labels or either of its two
cherry-center targets. Indeed, the candidate \(B\)-pair is
determined by either entry because the difference \(h\)
is fixed. Prescribing a cherry center and its incident
left node likewise determines the entire intermediary pair
by subtracting an \(A\)-label from the prescribed center.
The two centers along a single walk are distinct, because
the initial and terminal left pairs have different first
labels. It follows that at least \(g^2-O(g)\) ordered
intermediary pairs give a faithful trade at this terminal
node. The number of terminal nodes is at most \(N\).
By Cauchy–Schwarz these counts summed over the terminals
are \(\gg\delta^4/N\), using the lower bound
\(\sum g\ge\delta^2/2\) and enlarging \(K\) to absorb the
linear error. Sum over the \(\gg M/N\) starting nodes.
This gives \(\gg M^5/N^{10}\gg E^{10}/N^{15}\) faithful
ordered auxiliary configurations. Passing to the induced
matching pair loses at most an absolute constant, because
its twelve labels and their two edge lists determine all
the auxiliary nodes up to boundedly many orderings. The
common-target equation gives exactly (BT2), completing
the proof.

The numeric additivity in this argument is the reason the
common-center pair components have at most \(N\) nodes on
each side. Part caps alone only give the quadratic auxiliary
vertex count. At the minimum nonatomic endpoint, the
complete-additive reduction gives precisely the exponent
hypothesis of (BT1), with arbitrarily small extraction loss.
This supplies a fixed-size arithmetic pattern for the
three radical caps. A count that rules out its abundance
under those caps remains separate. The balance in (BT2)
automatically cancels label records on all three letters,
and so cannot by itself produce an off-grid defect on a
different unbalanced output.
