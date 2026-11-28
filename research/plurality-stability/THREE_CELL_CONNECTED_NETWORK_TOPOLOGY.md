# Topology of connected three-cell maximizers

Let \(A\) be a balanced planar Gaussian global maximizer at \(0<r<1\).
Suppose each of its three open strict-winner sets is connected. Let \(G\)
be the closure of the actual jump interfaces, including their genuine
triple junctions. Then \(G\) is connected and has only the following
possibilities.

| Unbounded strict-winner sets | Actual jump graph |
| --- | --- |
| Three | One trivalent junction with one proper ray for each pair of labels |
| Two, with labels 1 and 2 | Two trivalent junctions joined by compact arcs of types \(1,3\) and \(2,3\), with one proper \(1,2\) ray at each junction |
| One | Two trivalent junctions joined by three compact arcs, one for each pair of labels |

Here a proper ray is an embedded copy of \([0,\infty)\) whose tail
eventually leaves every compact subset of the plane. This gives a finite
embedded graph after one-point compactification of the unbounded network.
It supplies the topological part of the connected semialgebraic reduction
without assuming any parametrization at infinity. The asymptotic rigidity
arguments in that reduction impose additional hypotheses on the ends.
[Connected semialgebraic network rigidity](THREE_CELL_SEMIALGEBRAIC_NETWORK_RIGIDITY.md)

## The local graph and its three faces

Use the adjusted analytic maximizing scores \(V_i=T_r1_{A_i}-\beta_i\)
and write \(W_i=\{V_i>\max_{j\ne i}V_j\}\). Each pair difference is a
nonzero analytic function. Indeed, if \(V_i-V_j\) were identically zero,
integrating and using equal volumes would give \(\beta_i=\beta_j\).
Injectivity of \(T_r\) would then give \(1_{A_i}=1_{A_j}\) almost
everywhere, contrary to disjointness and their positive masses. The tie
set has measure zero, and the \(W_i\) give the maximizing cells up to
null sets. The graph \(G\) lies in the closed set where at least two
scores tie for the maximum. In particular, \(G\cap W_i=\varnothing\)
for every \(i\). This disjointness uses the maximizing score equations.
[Maximizer interface and score identities](THREE_CELL_GLOBAL_OBSTRUCTIONS.md)

The local analytic description in (G20) and its subsequent junction
analysis gives the following graph structure. An actual pair arc has
nonvanishing pair-score gradient at any finite accumulation point.
The only finite endpoints are genuine triple junctions. To see the
continuation away from a triple, a sufficiently small neighborhood
contains strict wins from only the two labels of that arc. At the
accumulation point their equality set is a regular analytic curve by
the gradient bound. On each side, the sign of the pair score difference
selects one of those two labels wherever the winner is strict. Thus
the pair jump continues along that curve. At genuine triple junctions
the local graph is a finite star of valence three or four. Its three
incidences in the trivalent case have distinct pair types. In the
four-valent case the cusp description gives exactly two pair types.
[Regular interface arcs and isolated triple junctions](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#regular-interface-arcs-and-isolated-triple-junctions)

More explicitly, the local finite analytic stratification separates
the score ties, the open winners and their closures. A sequence of
actual jumps accumulating at a point in a stratified neighborhood
has a subsequence in a single one-dimensional jump stratum, so
(G20) applies to a connected arc at that accumulation point. The genuine
triple set is closed and locally finite. Every point of \(G\) has a
neighborhood in which \(G\) is a single regular arc or a finite star
of arcs. A finite subcover of a compact subset of \(G\) by these
neighborhoods meets only finitely many graph edges, since every local
branch belongs to a single edge. Here an open edge is a component
after deleting the genuine vertices. The closure of an edge at a
finite end is interpreted with that endpoint adjoined. If an edge end
has finite accumulation, choose a subsequence in one local branch at
the accumulation point. At a regular point the branch continues
through the point, so the subsequence converges to an interior point
in the topology of the edge and cannot approach an end. At a vertex
the half-branch has that vertex as endpoint. In particular, an edge
with an end at infinity is proper at that end. These facts also imply
properness of a nonrepeating infinite edge path, since every compact
set meets only finitely many edges. This local argument allows the
graph to have infinitely many vertices before the global face count
is used.
[Bierstone–Milman, Corollaries 2.8 and 2.11](https://www.numdam.org/item/PMIHES_1988__67__5_0.pdf)

Each connected component of \(\mathbb R^2\setminus G\) has one label
almost everywhere. For a local verification, take a small disc off
\(G\) with a finite analytic stratification by winners and ties. In
the disc minus its finitely many zero-dimensional strata, the label on
the open strata extends locally constantly across each one-dimensional
stratum. Indeed its two sides have the same label, since a jump would
be part of \(G\). The punctured disc is connected, so there is just one
label there. These local assignments agree across overlapping discs
and give the face label. A strict win of a different label in the
face would persist on an open set, violating the almost-everywhere
assertion. Every face contains a strict winner because it is open and
ties have measure zero. Since each \(W_i\) avoids \(G\), is connected
and is nonempty, there are exactly three faces. The strict winners of
the face label are dense in their face, so a face is bounded exactly
when its strict-winner set is bounded.

The same observations give a bound for any closed subset \(H\subset G\).
Each face of \(\mathbb R^2\setminus H\) contains a strict winner, and
no connected \(W_i\) can meet two such faces because it avoids \(H\).
Consequently \(H\) has at most three complementary faces. Moreover, a
bounded face of \(H\) contains a bounded \(W_i\), and distinct bounded
faces use distinct labels. This bound uses strict-winner disjointness
from the jump set and connectedness of the winners.

There is no component of \(G\) made entirely of a regular pair
interface with no junction. Such a component would be a simple
closed curve or a proper embedded line by the local continuation
property. Either separates the plane by the Jordan curve theorem,
applied after one-point compactification in the line case. The two
strict-winner sets for its incident pair lie on opposite sides by
connectedness. The third strict winner lies on one side and cannot
have a regular jump interface with the label on the opposite side.
That interface would have to belong to the separating curve, where
every regular arc has only the original two incident labels. This
contradicts complete adjacency. In particular the nonempty graph
\(G\) has a genuine junction.
[Complete adjacency](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#complete-adjacency)

## Finite subgraphs and proper rays

For a finite embedded graph in the sphere, the number of complementary
faces is \(1+E-V+C\), where \(E,V,C\) are its edge, vertex and
component counts. Write \(b_1=E-V+C\). For a compact finite subgraph
of \(G\) in the plane, there are exactly \(b_1\) bounded faces. The
sphere formula also follows from Alexander duality, which identifies
reduced zeroth homology of the complement with first cohomology of
the graph. A topologically embedded finite graph is compact and
locally contractible as required by this theorem.
[Hatcher, *Algebraic Topology*, Theorem 3.44](https://pi.math.cornell.edu/~hatcher/AT/AT.pdf)

If a compact connected finite subgraph has \(k\ge 1\) disjoint
proper rays attached, with interiors avoiding the original subgraph,
its compactification has first Betti number increased by \(k-1\).
To see this, append \(k\) abstract intervals with separate terminal
points, and identify those terminal points at infinity. Properness
makes the map to the sphere continuous; the resulting finite graph
embeds because it is compact and the map is injective. Euler's
formula therefore applies even if an end meets large circles more
than once. The faces on the sphere are exactly those of the planar
complement. More generally, with rays attached to several disjoint
compact connected cores, the increment is the sum of \(k-1\) over
the ray-bearing cores. Identifying the terminal points also merges
those cores into one component, giving the stated increment.

In a tree component, each edge leaving a finite subtree extends to
a proper ray disjoint from the subtree except at its start. At every
finite vertex the path can continue without backtracking, since there
are no degree-one vertices. It cannot revisit a vertex or meet a
second ray constructed from another outgoing edge, as either event
forms a finite cycle. The same extension works outside a connected
skeleton containing the unique cycle in a component. A path leaving
the skeleton cannot meet it again, since the new path and a finite
path inside the skeleton would give a second cycle. The local
finiteness and properness of edges give properness in the plane.
Connectivity of any component gives finite paths between its vertices:
it is locally path connected, and a compact path image meets finitely
many edges.

## Three unbounded strict winners

A finite cycle in \(G\) would bound a face of that cycle containing a
bounded strict-winner set. Thus \(G\) is a forest. A genuine junction
and its outgoing rays give a skeleton with at least three rays. Four
rays would give at least four faces by Euler's formula. If two genuine
junctions were in the same component, the finite path between them
would have at least two remaining outgoing edges at each endpoint.
Extending them gives four rays, again impossible. If there were a
different component, it would have a genuine junction and hence at
least three disjoint proper rays. Combining this skeleton with the
first one gives more than three faces. It follows that the graph is
connected with a single trivalent vertex and three proper rays. The
incident pair types at the vertex are distinct.

## Two unbounded strict winners

There is exactly one finite simple cycle. In fact any compact finite
subgraph of \(G\) with first Betti number at least two would enclose
two bounded faces and hence two bounded strict winners. Conversely,
the bounded face of the full graph has compact boundary in \(G\).
By local finiteness and properness of edges, finitely many compact
edge segments and their endpoints cover this boundary. Their union
is a compact finite subgraph of \(G\). Its complement has a bounded
component because the covered boundary separates the original bounded
face from the exterior. Thus this subgraph contains a cycle. Two
distinct simple cycles would give a finite subgraph of first Betti
number at least two, by independence of the two distinct cycle
chains over \(\mathbb F_2\).

A component disjoint from the cycle would be a tree with a genuine
junction, and therefore contain a skeleton of three proper rays.
The union of that skeleton and the disjoint compact cycle would have
at least four faces. It follows that \(G\) is connected. It is also
unbounded because there are two unbounded faces, whereas the exterior
of any disc containing a bounded graph is connected.

Starting from the cycle, extend an outgoing edge to a proper ray.
An outgoing edge exists by connectedness and unboundedness. A finite
connected core containing the cycle can have at most two disjoint
proper rays attached, since with \(k\) rays its compactification has
\(k+1\) faces. Choose such a skeleton attaining the largest feasible
\(k\). If some edge were missing, connectedness gives an outgoing
edge from a finite point of the skeleton. Incorporate the compact
initial segment of an existing ray into the core if needed to put
that point in the core. The outgoing edge extends to another proper
ray disjoint from the skeleton except at its start, by uniqueness of
the cycle. This contradicts maximality. The whole graph is therefore
this skeleton, and its three faces force \(k=2\).

The compactified graph has a degree-two vertex at infinity and a
total finite degree surplus \(\sum_v(\deg v-2)=2\), by Euler's
formula. It has either two trivalent junctions or one four-valent
junction. In the latter case every edge has an incidence at that
junction, since there are no other junctions anywhere in the
connected planar graph. A four-valent junction has only two pair
types, contrary to complete adjacency. Thus there are two trivalent
junctions. Outside a large disc the bounded third strict winner is
absent, so both rays have type \(1,2\). Each junction has one
incidence of this type and one incidence of each of the other two
types. Hence each has one ray, and the two compact arcs connect the
junctions because a loop would contribute two incidences of the same
type at one junction.

## One unbounded strict winner

Every genuine jump arc has at least one bounded strict winner on its
sides. Thus \(G\) lies in the union of the closures of the two bounded
strict winners. It is compact and therefore a finite graph. Euler's
formula with three faces gives \(\sum_v(\deg v-2)=4-2C\), for \(C\)
components. Since there is a genuine junction and all vertices have
degree at least two, this forces \(C=1\) and total degree surplus two.
The single four-valent case is impossible by the same pair-type
adjacency argument. The two trivalent vertices are joined by one
compact arc of each pair type, by distinctness of their incidences.

The three compact arcs cannot all be graphs over a common axis.
Their projected intervals would share their two endpoints and
their order would be fixed on the open common interval. The labels
encountered along a transverse line would then follow the same
three-switch cycle from the unbounded label back to itself. Orient
each interface toward the positive transverse direction. The three
label incidence vectors have pairwise dot product \(-1\).
Translation in that direction has nonnegative normal speed on each
arc and positive speed on a set of positive interface measure on
each. Its first extra-coordinate form is zero by the translation
identity. Reversing its speed on one whole arc preserves the local
term and the three self-interactions, and strictly increases the
interaction with the other two arcs by positivity of the Mehler
kernel and the negative incidence products. This contradicts
nonnegativity of the first extra-coordinate form.
[Translation identity and nonnegativity](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#translation-equations-on-an-unconstrained-space-of-speeds)

## Straight arcs give affine separation

Suppose an open actual \(i,j\) arc lies on a line \(L\). The adjusted
score difference \(V_i-V_j\) vanishes along the arc and therefore on
all of \(L\) by real analyticity. Both \(W_i\) and \(W_j\) avoid \(L\).
At a regular point of the arc their strict winners lie locally on
opposite sides, so connectedness puts each entire strict-winner set
in its corresponding open half-plane. Thus \(L\) affinely separates
the two labels up to null sets.

At any positive correlation the affine-separation theorem makes
every strict-winner component unbounded. Consequently a graph with
any bounded strict winner has no open straight actual arc, including
on its unbounded ends. For \(r\ge2/3\), affine separation also makes
\(L\) pass through the origin and gives a finite continuous height
in the exclusive-fiber representation. This range includes every
possible first transition by the explicit comparison through
\(339/500>2/3\). With a second straight arc of a different pair type,
the straight-arc theorem in the exclusive-fiber representation
gives the centered simplex.
[Affine-separation theorem](THREE_CELL_AFFINE_SEPARATION_RIGIDITY.md)
[Comparison through 339/500](THREE_CELL_CENTRAL_SHIFT_BOOTSTRAP.md)
[Straight-arc rigidity](THREE_CELL_EXCLUSIVE_STRAIGHT_ARC_RIGIDITY.md#a-straight-arc-propagates-across-the-base-label)

## Finitely many strict-winner components suffice for a finite graph

Without assuming the winners are connected, suppose their total
number of connected components is a finite number \(N\). Then the
jump graph is still finite after including the point at infinity
when necessary. The local analytic graph description is unchanged.
For every closed subset of the jump graph, the complement has at
most \(N\) faces. Each face contains a strict winner, and each
connected strict-winner component is disjoint from the graph.

Every graph component contains a cycle or contains two disjoint
proper rays from a common point. In the cycle-free case, continue
in both directions along an edge, using the absence of finite
degree-one ends and the local finiteness argument above. Select a
cycle or the union of two such rays in each of any finite collection
of graph components. In the sphere the union has first Betti number
equal to the number selected, including when multiple pairs of rays
are joined at infinity. The face bound gives at most \(N-1\) graph
components.

Choose a compact finite subgraph with largest first Betti number,
which is at most \(N-1\) by the same face bound. It contains every
simple cycle. Otherwise adding a missed cycle would add an
independent cycle, since it uses an edge segment outside the chosen
subgraph, after subdivision. Within each graph component, join the
chosen pieces by finitely many finite paths to make a connected
compact core. Use a single point for a component with no chosen
piece. Any outgoing edge from a core continues to a proper ray.
It cannot return to the core or meet another such ray, because
this would create a cycle outside the core.

The number of attached rays in any finite union of skeletons is
bounded. For \(k\) rays attached to one core, their first Betti
number increment in the sphere is \(k-1\) if \(k>0\).
There are finitely many cores and the total first Betti number
is at most \(N-1\). Choose a skeleton attaining the largest
possible total number of rays. If it missed an edge, a finite
path from the core would leave the skeleton. Include the initial
compact ray segment up to the first exit in the core if needed.
Then the exiting edge continues to an additional proper ray
disjoint from the skeleton except at its start, since a further
meeting makes a cycle outside the core. This contradicts maximality.
The whole graph is therefore the chosen finite collection of
compact cores and proper rays.

Conversely, every face of the full jump graph contains exactly one
strict-winner component. Indeed, the face has some label \(i\) almost
everywhere by the local label assignment above. Every point of the
face at which \(V_i=V_j\), \(j\ne i\), is a dormant interior tie. By
(FC4) it is an isolated zero of that pair contrast within the face.
Each pair's zero set is relatively closed, hence locally finite there.
Removing the union of these sets from the face leaves a connected
set. Two points of the complement can be joined by a compact
polygonal path in the open face, detouring in small disks around
the finitely many zero points that path meets. The complement is
exactly the strict-winner set in that face. Thus the number of
strict-winner components equals the number of faces. In particular,
finiteness of the strict-winner components and finiteness of the
compactified jump graph are equivalent.
[Dormant interior ties, (FC4)](THREE_CELL_FIRST_CONTACT_VARIATIONS.md#point-transfers-at-a-dormant-interior-score-tie)

## A dormant score along an entire interface component

Suppose all three adjusted scores agree on an open actual jump arc
between labels \(i,k\), and the third label \(j\) has no strict winner
near the arc. Then the common equality persists on a whole component
of \(G\), a simple closed regular curve or a proper regular line,
and label \(j\) remains dormant along it.

The connected regular \(i,k\) edge containing the arc has analytic
pair level and nonvanishing pair gradient. The score of \(j\) agrees
with the other two along that edge by analytic continuation. At any
point of its interior, \(j\) is absent locally by the definition of
a genuine triple. In a transverse coordinate pointing toward the
\(k\)-side, the derivative of \(V_j\) is between those of \(V_i\)
and \(V_k\). Indeed these last two labels win on their respective
sides and the scores all agree on the wall. If a bound were an
equality, the corresponding dormant pair contrast would have zero
gradient. The half-ellipse insertion inequality (CC1) makes its
Hessian positive definite. But the contrast vanishes identically
along the wall, so its second derivative in the tangent direction
at a zero-gradient point is zero. Both bounds on the middle
gradient are therefore strict.
[Half-ellipse insertion test, (CC1)](THREE_CELL_CENTRAL_INTERFACE_COMPACTNESS.md#inserting-a-missing-label-next-to-a-regular-wall)

Any finite endpoint of the edge would have to be a genuine triple.
There the \(i,k\) pair gradient is still nonzero, so the same equality
with \(j\) extends through the endpoint on the analytic \(i,k\) level
curve. The transverse derivative of \(V_j\) lies between those of the
other two by continuity, and it lies strictly between them by
distinctness of the three gradients at a genuine triple. The order
of the derivatives persists in a neighborhood. Moving transversely
from the common equality curve then gives either \(i\) or \(k\) a
strict win on each side. Hence \(j\) has no local strict wins, contrary
to the endpoint being genuine. The edge has no finite endpoints
and makes a whole junction-free graph component.

Under connectedness of the three strict winners, a junction-free
component is impossible by complete adjacency and separation. This
also excludes the open common-score wall that is necessary for
unbounded junction counts in a fixed central disk satisfying the
radius hypothesis of the central interface compactness theorem,
provided the limiting maximizer has connected strict winners.
Without connectedness, an unbounded dormant component still forces
all labels to be unbounded. Normalize the prices to have sum zero.
The common adjusted score on the component is \(1/3\), so each
unadjusted score is a positive constant along it. The score of a
bounded label tends to zero at infinity, since its conditional source
Gaussian has mean \(\rho x\) and fixed variance. Any dormant component
in the presence of a bounded label must therefore be a closed curve.
[Central junction accumulation alternative](THREE_CELL_CENTRAL_INTERFACE_COMPACTNESS.md#collapse-of-components-and-accumulation-of-junctions)

## Inserting the dormant label along a wall

Let \(\Gamma\) be a dormant regular component of type \(i,k\) as above.
Orient its unit normal \(n\) from \(i\) toward \(k\) and set
\[
 b_i=-\partial_n(V_i-V_j),\qquad
 b_k=\partial_n(V_k-V_j),\qquad a=b_i+b_k.
\]
The two \(b\)'s are positive by the strict middle-gradient result.
Use the Gaussian surface measure \(d\nu=\phi_2\,ds\). For a smooth
nonnegative speed \(f\) of compact support on \(\Gamma\), write
\[
 \alpha=b_k/a,\qquad q=e_j-\alpha e_i-(1-\alpha)e_k,
 \qquad \lambda=\int_\Gamma q f\,d\nu.
\]
For small \(t\ge0\), insert label \(j\) in a strip of width
\(t\alpha f\) on the \(i\)-side and \(t(1-\alpha)f\) on the \(k\)-side.
Compact support permits a regular tube with only those two original
labels. Translate the whole partition and the two strips by
\(t v+O(t^2)\) to preserve the three masses exactly. The velocity \(v\)
is the unique solution of \(B s_v=-\lambda\), where \(s_v=v\cdot N\)
is the translation speed on the full graph and \(B\) is its mass-flux
map. Full moment rank gives the inverse and the mass equations follow
by the implicit function theorem. Indeed the changed masses are smooth
normal-coordinate integrals in the strip widths and the translation,
with a smooth signed-integral extension through zero width.
[Translation mass map and normal speed](THREE_CELL_GLOBAL_INTERFACE_BOUNDS.md#a-global-matrix-bound-on-volume-flux)

For the objective priced by the original mass multipliers, the
translation contributes \((\rho^{-1}-1)\|s_v\|_a^2t^2+o(t^2)\) by
(G2). The signed strip changes, divided by \(t\), converge as measures
on the fixed compact tube to \(q f\,d\nu\). Their quadratic interaction
therefore has coefficient \(\iint M_\rho(x,y)(q(x)\cdot q(y))f(x)f(y)
d\nu(x)d\nu(y)\), because the kernel is continuous on the tube.
The two linear score costs, divided by \(t^2\), are
\(\int[b_i\alpha^2+b_k(1-\alpha)^2]f^2\,d\nu\). For the
translated source the score on a translated point is evaluated at
the unshifted point plus \((1-\rho^{-1})(t v+O(t^2))\). Its resulting
additional linear cost at order \(t^2\) vanishes, because the two
transverse terms are proportional to
\(-b_i\alpha+b_k(1-\alpha)=0\). All these expansions are on compact
supports with smooth score kernels and regular normal coordinates.
Maximality of the exactly balanced competitor gives
\[
 \int_\Gamma \frac{b_i b_k}{a}f^2\,d\nu
 \ \ge\ \iint_{\Gamma^2} M_\rho(x,y)(q(x)\cdot q(y))
                             f(x)f(y)\,d\nu(x)d\nu(y)
             +(\rho^{-1}-1)\|s_v\|_a^2.                  \tag{NT1}
\]
The priced translation expansion is the same as for the interior
point transfer (FC4).
[Priced translation in (FC4)](THREE_CELL_FIRST_CONTACT_VARIATIONS.md#point-transfers-at-a-dormant-interior-score-tie)

Write \(q=q_0+(1/2-\alpha)(e_i-e_k)\), with
\(q_0=e_j-(e_i+e_k)/2\) orthogonal to \(e_i-e_k\) and of squared
norm \(3/2\). Positive semidefiniteness of the Mehler kernel gives a
nonnegative contribution from the second summand to the integrated
quadratic form. The local factor satisfies \(b_i b_k/a\le a/4\). For
nonzero nonnegative \(f\), the vector \(\lambda\) is nonzero because
its \(j\) coordinate is \(\int f\,d\nu>0\), so the translation term is
strictly positive. Hence
\[
       \int_\Gamma a f^2\,d\nu>6\iint_{\Gamma^2}M_\rho(x,y)
                                       f(x)f(y)\,d\nu(x)d\nu(y).
                                                        \tag{NT2}
\]
By absolute values and density, the same-pair interaction operator
\(2M_\rho\) restricted to this component has norm at most \(1/3\) in
the weighted space with local term \(\int a f^2d\nu\). For example,
if \(s_v^\Gamma=1_\Gamma s_v\) and \(s_v^{\rm out}=s_v-s_v^\Gamma\),
the translation identity forces
\[
 \rho\langle s_v^\Gamma,Ks_v^{\rm out}\rangle
                  \ge (1-\rho/3)\|s_v^\Gamma\|_a^2 .     \tag{NT3}
\]
[Positive Mehler factorization](THREE_CELL_GLOBAL_INTERFACE_BOUNDS.md#a-global-matrix-bound-on-volume-flux)
[Translation identity (G2)](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#translation-equations-on-an-unconstrained-space-of-speeds)

## A radius bound for compact dormant walls

If \(\Gamma\) is closed, let \(E\) denote its Jordan interior and
\(L\) the length of \(\Gamma\). Define the arclength mean squared
radius \(m_2=L^{-1}\int_\Gamma |x|^2ds\) and put
\[
 \sigma=\sqrt{1-\rho^2},\qquad
 \kappa=\frac{1-\rho}{2(1+\rho)},\qquad
 P(z)=2\Phi(z)-1.
\]
Then every closed dormant component satisfies
\[
 m_2>\frac{1}{\kappa}\log\!\left(
       \frac{3P(e\sqrt{\pi/(2\rho)})}{\rho^{3/2}}\right).
                                                        \tag{NT4}
\]
Here \(e\) is the base of the natural logarithm. In particular, no
such component lies entirely in the closed disk whose squared radius
is the right-hand side.

Use \(f=\phi_2^{-1/2}\) on \(\Gamma\) in (NT2). It gives
\(\int_\Gamma a\,ds>6\iint_{\Gamma^2} k_\rho(x,y)\,ds_x ds_y\), where
\[
 k_\rho(x,y)=\frac{1}{2\pi\sigma^2}
       \exp\!\left[-\frac{\kappa}{2}(|x|^2+|y|^2)
                         -\frac{\rho}{2\sigma^2}|x-y|^2\right]
\]
is the joint Gaussian density divided by the square root of its
two marginal densities. The convolution derivative bound gives
\(a\le 2\phi_1(0)\rho/\sigma\). Also choose the sign of \(F=U_i-U_k\)
to make the inside label its positive side at the wall. The outer
normal derivative of this contrast is \(-a\), even if other labels
occur inside \(E\). Euclidean divergence and the conditional Gaussian
Laplacian bound therefore give
\[
 \int_\Gamma a\,ds=-\int_E\Delta F\,dx
              \le 4\rho^2 |E|/(e\sigma^2),\qquad
 |E|\le L^2/(4\pi).
\]
The derivative bounds follow from
\(\mathbb E|Z_1|=2\phi_1(0)\) and
\(\mathbb E||Z|^2-2|=4/e\) for standard planar Gaussian \(Z\).
[Conditional Gaussian Laplacian bound](THREE_CELL_ISOLATED_SUBCLUSTERS.md#minimum-diameter-in-a-central-disk)

For a periodic arclength parametrization \(x(s)\), the chord from
\(x(s)\) to \(x(s+u)\) has length at most \(|u|\) for
\(-L/2\le u\le L/2\). For each such shift, Jensen's inequality
bounds the average of the first exponential factor below by
\(e^{-\kappa m_2}\). Put \(\ell=L\sqrt\rho/(2\sigma)\). Integrating the
remaining Gaussian in \(u\) gives
\[
 \iint_{\Gamma^2} k_\rho\,ds_x ds_y
 \ge \frac{L e^{-\kappa m_2}P(\ell)}{\sqrt{2\pi}\sigma\sqrt\rho}.
\]
Combining this bound with the two separate upper bounds for
\(\int a\,ds\) yields
\[
 e^{-\kappa m_2}
   <\min\left\{\frac{\rho^{3/2}}{3P(\ell)},\
              \frac{\rho^2\sqrt{2/\pi}\,\ell}{3e P(\ell)}\right\}.
\]
The first function decreases, while the second increases because
\(P(\ell)/\ell\) is the average of a decreasing positive density on
\([0,\ell]\). They cross at \(\ell=e\sqrt{\pi/(2\rho)}\), giving (NT4).

For a compact component of small spatial variance an alternate bound
can be stronger. Let \(p=L^{-1}\int_\Gamma x\,ds\) and
\(v=L^{-1}\int_\Gamma|x-p|^2\,ds\). Jensen's inequality on the full
product of the normalized arclength measures gives
\[
 \iint_{\Gamma^2} k_\rho\,ds_xds_y
 \ge \frac{L^2}{2\pi\sigma^2}
    \exp\!\left[-\kappa|p|^2-(\kappa+\rho/\sigma^2)v\right],
\]
since two independent points of that measure have mean squared
distance \(2v\). Use the divergence and isoperimetric bounds again to
obtain the necessary condition
\[
           \kappa|p|^2+(\kappa+\rho/\sigma^2)v
                        >\log(3e/\rho^2).                  \tag{NT5}
\]
