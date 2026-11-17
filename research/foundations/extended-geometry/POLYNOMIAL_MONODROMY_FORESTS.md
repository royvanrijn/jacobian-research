# Polynomial monodromy forests and resonance source trees

## Result

The source dual graph and all node expansion indices on a simple-branch
resonance stratum are determined by one labelled tree.

Let

\[
 f:\mathbf P^1\longrightarrow\mathbf P^1
\]

be a degree-\(d\) polynomial with simple finite branch points.  Label a
regular fiber by \(\{1,\ldots,d\}\).  The finite branch monodromies are
transpositions

\[
 \tau_1,\ldots,\tau_{d-1}\in S_d,
\]

and their product is the inverse of the \(d\)-cycle at infinity.  Regard
each transposition \((ab)\) as an edge joining the sheet labels \(a,b\).
Those \(d-1\) edges form a labelled tree \(T_f\).

If a subset \(S\) of the finite branch points collides and bubbles on the
stable target, let \(T_f[S]\) be the subforest consisting of the
corresponding edges.  Then:

1. the source components over the new target bubble are the connected
   components of \(T_f[S]\);
2. a component with \(k\) vertices has degree \(k\);
3. its attaching source node has expansion index \(k\);
4. isolated vertices give degree-one identity strands; and
5. inclusions \(S\subset S'\) give refinements of forest components, hence
   the source trees are compatible on nested target strata.

This is the missing combinatorial rule for the residue-resonance part of the
wonderful target graph.  Together with the
[general radial source atlas](GENERAL_RADIAL_SOURCE_ATLAS.md), it determines
the source dual tree, component degrees, and node indices on every stratum
whose generic finite branching is simple.

The algebraic half of `H1-STACK` is supplied by the
[recursive resonance atlas](RECURSIVE_RESONANCE_ATLAS.md): normalized
initial forms place the actual target-flag divisors in framed root-residue
coordinates, make the third-fiber equations automatic, prove contractions,
and extract the simultaneous labelled-node character.  Its selected-factor
log-étale comparison supplies all-degree coverage and effective descent.
The root-side coarse graph is independently identified by
[finite normalization](SOURCE_GRAPH_FINITE_NORMALIZATION.md).

## 1. The Dénes tree

The monodromy at infinity of a polynomial is a \(d\)-cycle.  Simple finite
branch points have transposition monodromy.  Hence

\[
 \tau_{d-1}\cdots\tau_1=\sigma_\infty^{-1}.          \tag{1.1}
\]

The transpositions generate a transitive subgroup because the source is
connected.  Their edge graph is therefore connected.  It has \(d\) vertices
and \(d-1\) edges, so it is a tree.

Conversely, the minimal transposition factorizations of a fixed \(d\)-cycle
are counted by

\[
 d^{d-2},                                           \tag{1.2}
\]

the number of labelled trees.  A direct tree-factorization bijection and
refinements of Dénes's theorem are given by Goulden--Yong,
[Tree-like properties of cycle
factorizations](https://arxiv.org/abs/math/0106039).

The tree \(T_f\) should not be confused with the stable source dual tree.
It is the monodromy tree on the \(d\) sheets of a regular fiber.  Stable
source trees are obtained from its subforests.

## 2. Collision-subforest theorem

### Theorem 2.1

Let \(S\subset\{1,\ldots,d-1\}\), and put

\[
 \sigma_S=\prod_{i\in S}^{\mathrm{path\ order}}\tau_i. \tag{2.1}
\]

The connected components of the source over a disk containing precisely the
branch points in \(S\) are the orbits of

\[
 G_S=\langle\tau_i:i\in S\rangle.                   \tag{2.2}
\]

These orbits are exactly the vertex sets of the connected components of
\(T_f[S]\).  If their sizes are \(k_1,\ldots,k_m\), then

\[
 \operatorname{cycle\ type}(\sigma_S)
 =(k_1,\ldots,k_m).                                 \tag{2.3}
\]

After the disk becomes a target bubble, there is one source component of
degree \(k_j\) for each orbit and one attaching source node of expansion
index \(k_j\).

### Proof

Orbits of the local monodromy group classify connected components of the
cover over the punctured disk.  Since the generators are the edge
transpositions, their orbits are the connected components of the selected
edge graph.

Each nontrivial component is a tree on \(k_j\) vertices and has \(k_j-1\)
edges.  The product of its edge transpositions, in any inherited order, is
a \(k_j\)-cycle.  Thus the boundary monodromy has one cycle of length
\(k_j\) on that component.  Boundary-monodromy cycles are the points over
the attaching node, and cycle lengths are their expansion indices.
\(\square\)

The component is genus zero.  Its selected simple branch points contribute
\(k_j-1\) to ramification and its attaching \(k_j\)-cycle contributes
another \(k_j-1\), so

\[
 (k_j-1)+(k_j-1)=2k_j-2.                            \tag{2.4}
\]

This is Riemann--Hurwitz component by component.

## 3. Nested compatibility

If \(S\subset S'\), every component of \(T_f[S]\) is contained in a unique
component of \(T_f[S']\).  Therefore the source component partition for
\(S\) refines that for \(S'\).  For a nested set of target collisions, apply
this observation along every inclusion.  It produces:

- the source vertices above every target vertex;
- the degree of each source vertex;
- every source edge above a target edge; and
- the expansion index of that edge.

This is exactly the harmonic morphism of trees underlying the tropical
admissible cover.  The compatibility of classical and tropical branch maps
is developed by Cavalieri--Markwig--Ranganathan,
[Tropicalizing the space of admissible
covers](https://arxiv.org/abs/1401.4626).

The forest rule is braid-equivariant.  Changing branch paths applies Hurwitz
moves to the ordered transposition tuple and transports the labelled cover;
the orbit/cycle construction itself is intrinsic.  For an algebraic global
comparison one must still carry out that descent on the root-cluster graph,
rather than fix branch cuts once and for all.

## 4. Maxwell and caustic in one rule

The degree-six source-node partitions now have a uniform explanation.

Two disjoint transposition edges select a forest with component sizes

\[
 (1,1,2,2).
\]

This is the pairwise Maxwell node.  The normalization branch count is

\[
 \frac{1\cdot1\cdot2\cdot2}{\operatorname{lcm}(1,1,2,2)}
 =2.                                                \tag{4.1}
\]

Three disjoint edges give

\[
 (2,2,2),
\qquad
 \frac{2\cdot2\cdot2}{\operatorname{lcm}(2,2,2)}
 =4,                                                \tag{4.2}
\]

which is the triple Maxwell node.

Two adjacent edges give a three-vertex tree and three isolated sheets:

\[
 (1,1,1,3).                                        \tag{4.3}
\]

This is the simple two-critical-point caustic pattern.  Its normalization
branch count is one.  Thus “Maxwell versus caustic” is simply “disjoint
versus adjacent edges” in the same polynomial monodromy tree.

The direct degree-six radial/Maxwell atlases remain valuable algebraic
checks, but their node partitions are no longer isolated formulas.

## 5. Interface with source-vertex rigidity

The forest theorem supplies the combinatorial part of a flag-complete source
enhancement.  For a component with vertex set \(O\):

- its degree is \(|O|\);
- the outer node fiber is one point of multiplicity \(|O|\);
- an inner node fiber is the collection of components of the corresponding
  smaller subforest, with their sizes as multiplicities; and
- when the target zero mark lies on the bubble, its regular preimages carry
  the sheet labels in \(O\).

These are the multiplicities of the target-flag divisors.  Their positions
on the algebraic source component, and one third-fiber equation, determine
the rational map by
[source-vertex rigidity](SOURCE_VERTEX_RIGIDITY.md).

Over \(\mathbf C\), a fixed monodromy tuple and fixed target branch
configuration determine the cover up to isomorphism by the Riemann existence
theorem.  Hence no pointwise discrete choice remains after \(T_f\) is fixed.
The stronger algebraic-family statement is:

> recover this monodromy-tree Hurwitz family functorially from the
> root-cluster graph, express its flag divisors by the recursive residue
> coordinates, and prove compatibility with every wonderful contraction.

The recursive resonance theorem establishes that statement; the component
maps are therefore unique and the
[labelled-node theorem](LABELLED_NODE_SATURATION.md) supplies normalized
branches and actual inertia.

## 6. H1-STACK pipeline

The current comparison can be organized with no repeated case-by-case
Hurwitz search:

1. build the wonderful target tree from branch-value diagonals;
2. read radial source components from cluster multiplicities;
3. decorate every resonance vertex by subforests of \(T_f\);
4. construct the corresponding target-flag divisors in residue coordinates;
5. apply the two-fiber/third-fiber reconstruction;
6. saturate node monoids and compute inertia characters; and
7. take the formal H2 subgroup quotient.

All seven steps are general.  The recursive resonance atlas implements step
4, its selected-factor log-étale comparison supplies global coverage and
descent, and its vertex-centralizer tables feed into step 6.

## Reproduction

Run:

```bash
.venv/bin/python scripts/verify_polynomial_monodromy_forests.py
```

The checker exhausts all 1,441 reduced transposition factorizations of a
fixed cycle in degrees two through six, with counts
\(1,3,16,125,1296=d^{d-2}\).  It checks 43,614 collision subforests,
325,515 nested subset pairs, 151,499 componentwise Riemann--Hurwitz
identities, and recovers the degree-six pairwise Maxwell, triple Maxwell, and
two-edge caustic node partitions.
