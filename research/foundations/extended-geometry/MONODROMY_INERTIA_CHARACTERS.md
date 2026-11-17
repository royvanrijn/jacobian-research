# Monodromy centralizers and resonance inertia

## Result

The generic stack inertia of every radial or simple-branch resonance node is
determined by the same labelled monodromy data that determines its source
tree.

For a connected cover of a punctured target component with monodromy group
\(G\subset S_d\), deck transformations are

\[
 \operatorname{Deck}(f)=C_{S_d}(G),                 \tag{0.1}
\]

restricted to transformations preserving the source labels and node cycles.
Their rotations on the node cycles give the character

\[
 \chi:\operatorname{Aut}_{\mathrm{lab}}(f)
 \longrightarrow\prod_j\mu_{e_j}                   \tag{0.2}
\]

in the [labelled-node inertia theorem](LABELLED_NODE_SATURATION.md).

For polynomial simple-branch resonance components this becomes especially
small:

- a connected \(k\)-sheet transposition tree generates \(S_k\);
- its deck centralizer is trivial for \(k\ge3\);
- for \(k=2\), the centralizer is the quadratic flip \(\mu_2\); and
- one labelled point in a regular fiber kills that flip.

For a radial chain, those endpoint anchors must be matched through every
connector and every node normalization.  Equal cluster multiplicities give
trivial radial inertia, but unequal multiplicities can leave a finite
full-chain subgroup.  On a resonance bubble away from the target zero flag,
every two-vertex forest component is unanchored and may contribute a flip.
On a bubble containing the target zero flag, the fully labelled zero fiber
kills all flips on that vertex; this alone does not kill rotations on
different connector vertices.

This generalizes the degree-six observation that the same node partition can
have trivial radial inertia and nontrivial Maxwell inertia.

## 1. Deck groups as centralizers

Choose a base point in the complement of the target branch flags and label
its \(d\) inverse images.  A deck transformation permutes those sheets and
commutes with continuation around every target loop.  Thus it centralizes the
monodromy image.  Conversely a sheet permutation centralizing every
monodromy permutation descends to a deck transformation.  This proves
(0.1).

If a source point in that regular fiber is labelled, the deck permutation
must fix its sheet.  If a source node is retained individually, the deck
permutation must preserve its monodromy cycle setwise.  On a cycle of length
\(e\), every commuting setwise-preserving permutation is a rotation, giving
one element of \(\mu_e\).  Collecting those rotations is (0.2).

This supplies a direct algorithm for the abstract character in the
labelled-node theorem:

1. generate the vertex monodromy group;
2. compute its centralizer in the symmetric group;
3. intersect with stabilizers of all source labels and node cycles;
4. record the rotation on each node cycle; and
5. take the inverse image of the diagonal \(\mu_L\), where
   \(L=\operatorname{lcm}(e_j)\).

No phase subgroup is declared to be inertia before step 5.

## 2. Simple-branch components

Let a resonance source component correspond to a connected component of the
[polynomial monodromy subforest](POLYNOMIAL_MONODROMY_FORESTS.md) on
\(k\) sheets.  Its edge transpositions generate \(S_k\).  Therefore

\[
 C_{S_k}(S_k)
 =
 \begin{cases}
  S_2\simeq\mu_2,&k=2,\\
  1,&k\ge3.
 \end{cases}                                        \tag{2.1}
\]

The degree-one case is also trivial.  A \(k=2\) component is the map
\(U\mapsto U^2\) after choosing its two branch flags.  Its flip
\(U\mapsto-U\) survives if its regular fibers are unlabelled.  If one sheet
is labelled, the flip exchanges that sheet with the other and is killed.

The target zero flag is the natural anchor.  When it lies on the bubble,
every sheet in its inverse image is a labelled polynomial root, so all
quadratic flips disappear.  When it lies on the other side of the attaching
node, a quadratic resonance component carries no root label and its flip can
survive.

## 3. Closed node formula

Let

\[
 (e_1,\ldots,e_r)
\]

be the forest-component sizes over one target node.  Mark each \(e_j=2\) as
**unanchored** if its quadratic component contains no labelled regular-fiber
point.  All other indices are constrained to have trivial deck phase.  Put

\[
 L=\operatorname{lcm}(e_1,\ldots,e_r),\qquad
 M=\operatorname{lcm}\{e_j:
   e_j\ne2\ \text{or the }e_j\text{-component is anchored}\},       \tag{3.1}
\]

with \(M=1\) if the set is empty.  Then the inertia order on one normalized
lift is

\[
 |I|=\frac{L}{M}\in\{1,2\}.                         \tag{3.2}
\]

Indeed a diagonal phase \(q\in\mathbf Z/L\) must vanish modulo every
constrained index.  The only allowed nonzero coordinates are flips modulo
two.  Formula (3.2) follows.

In particular, if the target zero flag lies on the bubble, every quadratic
component is anchored, \(M=L\), and inertia is trivial.

Away from the zero flag, inertia is \(\mu_2\) exactly when the lcm of all
nonquadratic component sizes is odd and at least one quadratic component is
present.  An even nonquadratic component kills the diagonal flip.  For
example:

\[
\begin{array}{c|c}
\text{node partition} & \text{unanchored inertia}\\ \hline
(1,1,2,2) & \mu_2\\
(2,2,2) & \mu_2\\
(1,1,1,3) & 1\\
(2,3) & \mu_2\\
(2,4) & 1.
\end{array}                                         \tag{3.3}
\]

## 4. Radial connectors

A degree-\(\mu\) power connector has monodromy generated by a
\(\mu\)-cycle, so its unanchored deck group is

\[
 C_{S_\mu}(\langle(1\,2\,\cdots\,\mu)\rangle)
 \simeq\mu_\mu.                                     \tag{4.1}
\]

But a radial connector is not an isolated cover component.  A rotation can
remain an automorphism of the nodal special fiber even when a different
component farther down the chain contains the labelled roots.  Thus labels
do not kill connector rotations componentwise.

The relevant inertia is smaller than the special-fiber automorphism group.
One must match the cycle permutations and rotations from both sides of every
source node and ask whether the resulting action preserves one simultaneous
Kummer normalization branch.  The
[recursive resonance atlas](RECURSIVE_RESONANCE_ATLAS.md) performs this
full-chain calculation.

For an ordered radial scale partition, with \(B_0\) first active,
\(B_0\mid\cdots\mid B_k\), set

\[
 M_j=\operatorname{lcm}_{i\in B_j}\mu_i,\qquad
 L_j=\operatorname{lcm}_{i\in B_j\cup\cdots\cup B_k}\mu_i. \tag{4.2}
\]

The diagonal node phases telescope along each connector chain.  Once the
earlier phases are fixed, endpoint rigidification leaves \(L_j/M_j\)
choices at level \(j\).  Therefore:

> **General radial inertia theorem.**
> \[
>  |I_{\rm radial}|
>    =\prod_{j=0}^{k}\frac{L_j}{M_j}.               \tag{4.3}
> \]

This is a statement about one normalized polynomial lift, not about the raw
automorphism group of its nodal special fiber.  Equal multiplicities make
all factors one.  This explains why the thirteen degree-six
\((2,2,2)\) types have trivial normalized inertia despite special-fiber
automorphism groups of orders up to eight.  Unequal multiplicities need not:
the two strict orders on multiplicities \(2,3\) have inertia orders \(3\)
and \(2\).  The compiler also exhibits an unrigidified nonpolynomial tame
chain with order-six ghost inertia, showing why a componentwise anchor
argument is insufficient.

## 5. Degree-six unification

At a pairwise Maxwell divisor, two disjoint transposition edges and two
isolated sheets give

\[
 (1,1,2,2).
\]

The bubble is away from target zero, so both quadratic components are
unanchored.  Formula (3.2) gives one diagonal \(\mu_2\), while normalization
has two branches.

At the triple Maxwell divisor, three disjoint edges give

\[
 (2,2,2).
\]

All three flips exist, but their intersection with the diagonal phase is
again only one \(\mu_2\); normalization has four branches.

At a radial coordinate point, the same partition \((2,2,2)\) occurs over a
bubble containing target zero.  Its six labelled root sheets kill the three
flips, so inertia is trivial.  This explains the former apparent mismatch
without chart-specific automorphism guesses.

Two adjacent transposition edges give the caustic partition
\((1,1,1,3)\).  The degree-three component has full \(S_3\) monodromy and
trivial deck centralizer, so its inertia is trivial whether or not zero lies
on the bubble.

## 6. Higher-codimension completion

Formula (3.2) gives generic inertia at one target node.  At a stratum with
several nested target nodes, the same global deck transformation may act on
more than one node.  The correct group is obtained by applying the
centralizer-character algorithm simultaneously to all vertex monodromies and
edge cycles; it need not be a naive product of the generic node groups.

The degree-six pair--triple crossings have already passed this simultaneous
test and give \((\mu_2)^2\).  The shared phase engine now implements the
general operation: compose local character tables with multiplicity, split
the resulting phase vector by target node, and require every segment to lie
in its diagonal \(\mu_{\mathrm{lcm}}\).  It distinguishes two independent
node flips (order four) from one connector flip coupled across two nodes
(order two), and rejects off-diagonal flips.

The [recursive resonance atlas](RECURSIVE_RESONANCE_ATLAS.md) now performs
the missing extraction.  It uses full vertex centralizers, retains
permutations of isomorphic disconnected source components, matches their
oriented node-cycle actions across adjacent vertices, and then feeds the
resulting two-sided rotations into the simultaneous diagonal intersection.
Thus the automatic nested-tree character compiler is complete for arbitrary
tame branch-cycle profiles.  Its checker also verifies (4.3) on 76 ordered
equal- and unequal-multiplicity radial charts of total degree at most seven.

## Reproduction

Run:

```bash
.venv/bin/python scripts/verify_monodromy_inertia_characters.py
```

The checker computes deck centralizers for all 1,441 polynomial monodromy
trees through degree six, cyclic connector centralizers through degree
eight, and anchored/unanchored inertia for all 43,614 collision nodes.  It
finds 28 node partitions and recovers the radial, pairwise Maxwell, triple
Maxwell, and caustic conclusions above.  It also verifies independent,
coupled, and off-diagonal simultaneous higher-codimension character tables.
