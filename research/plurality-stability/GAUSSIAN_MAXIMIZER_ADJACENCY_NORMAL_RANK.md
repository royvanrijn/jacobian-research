# Adjacency and normal rank at Gaussian noise maximizers

A set of pair interfaces whose deletion leaves a bipartite graph must
carry normals spanning all but at most one active direction of a
Gaussian noise-stability maximizer. If the deletion isolates any label,
those normals span the full active space. This extends the known odd
cycle obstruction. For a
full-dimensional four-cell maximizer, if a single pair is absent and the
other five pairs are present, the interface joining the two labels
outside the absent pair cannot have parallel normals. These statements
apply at a single positive correlation, including at a tie with the
simplex.

## The translation kernel restricted to a bipartite subgraph

Fix \(0<r<1\) and a maximizing ordinary partition \(A\) with prescribed
positive masses in an admissible dimension \(d\ge k-1\). Let \(G\) have
vertices the cell labels and one edge for each pair interface of positive
Gaussian surface measure. The distributional gradients of the cell
indicators are the oriented pair-interface measures. If \(G\) were
disconnected, the union of one component's cells would have zero
distributional gradient and measure strictly between zero and one.
Thus \(G\) is connected.

For each edge \(e=(i,j)\), let \(d_e=e_i-e_j\) be its label vector,
\(N_e\) the unit normal out of \(A_i\), and \(\nu_e\) its Gaussian
surface measure. Put \(a_e=|\nabla T_r(1_{A_i}-1_{A_j})|\) on the
interface, and let \(K\) be the tagged scalar interaction with
kernel \((d_e\cdot d_f)M_r(x,y)\), where \(M_r\) is the relative Mehler
kernel. The translation equations and stability in one unused
coordinate give
\[
      L=\frac1r a-K\succeq0,\qquad L(v\cdot N)=0
                          \quad(v\in\mathbb R^d).            \tag{AR1}
\]
The inequality means nonnegativity of the quadratic form on
\(L^2(a\,d\nu)\). It first follows on compact regular speeds by
multiplying by the linear Gaussian mode in the unused coordinate,
with a compact odd cutoff and its limit. First-order mass fluxes
then vanish automatically. To extend it, first test indicators
on compact parts of a single tag with approximating regular
speeds. Positivity of its Mehler interaction gives \(a>0\)
almost everywhere on that tag. Exhaustion and the boundedness of
\(a\) and \(\nu\) give \(T_{\sqrt r}\nu_e\in L^2\).
The scalar form extends from bounded speeds by its bound in
\(L^2(a\,d\nu)\). The norm measure is finite, so all translations
belong to this space. Equivalently the second assertion in (AR1)
follows by Gaussian integration by parts applied to the cell
indicators and the stationary normal score relation. Orienting
any edge the other way changes both its incidence vector and its
normal, so both assertions continue to hold.
[Heilman–Tarter, Theorem 1.9 and Lemmas 4.8 and 5.4](https://arxiv.org/pdf/2011.05583v2)
[The extra-coordinate bound for all label counts](GAUSSIAN_ADAPTIVE_SURFACE_HIGH_CORRELATION_LIMIT.md#consequences-for-optimal-values-and-moving-maximizers)

Let \(D\subseteq E(G)\) be a set of deleted edges such that
\(G\setminus D\) is bipartite, possibly disconnected. Let \(V_D\)
be the linear span, in the essential-support sense, of all normal
vectors on the deleted interfaces, and \(V\) the span of the
normals on all interfaces. Then
\[
           \dim V\le\dim V_D+1,\qquad
           V=V_D\ \text{if }G\setminus D\text{ has an isolated label}.
                                                               \tag{AR2}
\]
Here \(\dim V\) is the active dimension of the partition. In fact,
for \(v\perp V\) each directional distributional derivative of
the cell indicators vanishes, so the indicators are constant
along that direction. This can be seen by Fubini on lines or
by convolution and passage to the limit. Conversely dependence
only on a subspace forces all normals into that subspace.

To prove (AR2), restrict the form \(L\) to the retained tags.
Each connected label component with retained edges gives its own
block, since distinct components have disjoint labels. Orient all
edges in such a component from one side of its bipartition to
the other. Their incidence dot products are two for the same edge,
one for distinct edges with a shared endpoint, and zero otherwise.
In particular the scalar kernel in that block is nonnegative.

The null space of that block has dimension at most one. If \(f\)
is a nonzero null vector, taking the absolute value keeps the
diagonal term and can only increase the interaction. Nonnegativity
of the form therefore makes \(|f|\) a null vector too. For
unbounded \(f\) the interaction with \(|f|\) is obtained by
monotone approximation with bounded nonnegative fields and
the form bound in (AR1). The original signed interaction is
absolutely convergent by the same bound for \(|f|\).
The weak equation \(K|f|=(a/r)|f|\) follows by testing against
bounded \(L^2(a\,d\nu)\) functions and using the Cauchy–Schwarz
inequality for the nonnegative form \(L\). The equation and
strict positivity of \(M_r\) imply that \(|f|\) is positive
almost everywhere on any tag where it has a positive-mass
support, and on every tag incident with that tag.
Indeed a positive-measure zero set on a receiving tag could be
tested by its indicator, giving a positive left side and a zero
right side. Connectedness propagates this to every retained tag
in the block. Equality on replacing \(f\) by \(|f|\) now implies
that \(f\) has one sign on each tag and the signs agree across
incident tags. Thus every nonzero null vector has a fixed strict
sign almost everywhere in its block. Two linearly independent
null vectors would permit a nonzero pair orthogonal in
\(L^2(a\,d\nu)\), which is impossible for such vectors.

The equations on the deleted tags tie together the null vectors
of the retained components. Let \(f\) be a global null vector of
\(L\) vanishing on all deleted tags. On a retained component that
has a nonzero block kernel, write \(f=c_C f_C\), where \(f_C\) is
the strictly positive block null vector in the bipartite orientation.
On any component with zero kernel or no edges the coefficient is
fixed at zero.

For a deleted edge \(e\) joining distinct retained components
\(C,C'\), the global equation on \(e\) is \(Kf=0\) there. The only
contributing retained edges touch one of the two endpoints of \(e\).
All retained edges at a given endpoint share its incidence sign in
their bipartite orientation. Thus for either component that has a
coefficient, its contribution has one fixed nonzero sign times that
coefficient. The strict sign follows from positivity of the Mehler
kernel and of \(f_C\) on each incident retained edge. Test this
equation against a nonzero nonnegative bounded field supported on
\(e\). The integrals are finite: each component contribution has a
constant sign and is a bounded form pairing between its null vector
and the test. This gives one scalar relation involving only \(c_C\)
and \(c_{C'}\), with nonzero multiplier for each coefficient present.
In particular a zero coefficient on one side forces zero on the
other, and otherwise one coefficient determines the other.

The graph obtained by contracting retained components is connected
because \(G\) is connected. The scalar relations along a spanning
tree show that the space of global null vectors supported on retained
tags has dimension at most one. If any retained component is an
isolated label with no edges, the zero coefficient there propagates
to every component, and this space is zero. Finally the translations
\(v\cdot N\) for \(v\in V\cap V_D^\perp\) lie in this global null space
by (AR1), and the map from these vectors is injective by the definition
of \(V\). This proves (AR2).

## Geometric exclusions

If no edges need deletion, (AR2) gives
\(\dim V\le1\). All normal directions are then of the form
\(\pm v\) for one unit direction \(v\). The common sign conclusion
in the bipartite orientation makes them all equal after choosing
the sign of \(v\). For each label on either side, all outward
normals have one common sign, opposite on the other side.
Consequently each indicator is monotone in \(v\cdot x\), as
follows from the sign of its directional distributional
derivative. Its positive mass and positive complementary mass
make its cell a half-space in that direction. Two labels on
the same side would have overlapping half-spaces. Hence each
side contains at most one label. A maximizer with \(k\ge3\) has
a nonbipartite adjacency graph and therefore contains an odd cycle.
For \(k=3\) this recovers complete adjacency. The three-cell
argument by absolute values also appears in the adaptive setting.
[Cross-interface cut inequality and the known nonbipartite obstruction](FINITE_REFLECTION_TEST_FOR_MAXIMIZERS.md#cross-interface-normals-require-internal-interactions)
[Complete adjacency from normal saturation in three cells](THREE_CELL_ADAPTIVE_NORMAL_SATURATION.md#complete-adjacency-at-saturation)

Consider a label whose removal leaves a bipartite induced graph.
Delete all its incident edges. The retained graph isolates that
label, so (AR2) says its own interface normals span the full active
space. In particular the cell's indicator cannot be cylindrical
in any active direction. For three cells this applies to each label.

Take \(D\) to be the edges omitted by a maximum cut. If every pair interface has
parallel normals within that pair, (AR2) gives
\[
         |E(G)|-\operatorname{maxcut}(G)\ge\dim V-1,          \tag{AR3}
\]
since each deleted edge supplies at most one normal direction.
For \(k=4\) and \(\dim V=3\), any graph with a missing pair
\(i,j\) has a cut omitting at most one edge: put \(i,j\)
on one side and the other two labels on the other. Thus
a full-active four-cell maximizer with parallel normals on
each pair must have all six pairs present. If exactly one
pair is missing, the crossing subgraph of that cut is the
connected four-cycle. Its only deleted edge is the pair
between the two other labels. Equation (AR2) forces the
normals of that single interface to span at least two
directions, even if the other interfaces have arbitrary
geometry. For a triangle-free five-cell maximizer, the only
possible graph is a chordless five-cycle by the known odd-cycle
obstruction. Deleting any one edge leaves a connected bipartite
graph. Thus each pair interface in this graph must have normals
spanning at least three directions if the active dimension is four.

This obstruction uses the actual partition and the extra-coordinate
stability inequality. A complete adjacency graph is compatible
with both simplicial and curved or disconnected configurations,
so the rank condition leaves general tied maximizers open.
