# Source-vertex rigidity and the flag-atlas reduction

## Result

The source-cover problem does not require a new Hurwitz enumeration at every
boundary vertex.

Let

\[
 f:E\simeq\mathbf P^1\longrightarrow P\simeq\mathbf P^1
\]

be one component of a genus-zero admissible cover.  Choose three distinct
flags \(a,b,c\) on the stable target component \(P\), and let \(t\) be the
unique target coordinate satisfying

\[
 t(a)=0,\qquad t(b)=\infty,\qquad t(c)=1.
\]

Then the two pullback divisors

\[
 D_a=f^*[a],\qquad D_b=f^*[b]                         \tag{0.1}
\]

determine \(t\circ f\) up to a nonzero scalar.  Requiring the prescribed
third fiber \(D_c=f^*[c]\), or merely one prescribed point of that fiber,
fixes the scalar.  Hence the source component and the inverse-image data of
three target flags determine the component map uniquely.

An explicit logarithmic presentation of the source half of `H1-STACK` is therefore
reduced to one exact problem:

> **Flag-extension problem.**  
> Construct, over every stratum of the wonderful stable-target graph, the
> source tree together with the weighted inverse-image divisors of its target
> flags, and prove that this construction is compatible with edge
> contraction.

Once this enhancement is written in root-residue coordinates, the vertex
maps are forced by the elementary
rigidity theorem below.  Node normalization and actual stack inertia are
then supplied by the independent
[labelled-node saturation theorem](LABELLED_NODE_SATURATION.md).

This replaces the vague instruction “select a polynomial Hurwitz class at
every vertex” by a concrete divisor-construction and compatibility problem.
It does **not** by itself provide that explicit atlas.  Abstract existence
and the coarse comparison are supplied independently by
[finite normalization](SOURCE_GRAPH_FINITE_NORMALIZATION.md).

The subsequent
[general radial source theorem](GENERAL_RADIAL_SOURCE_ATLAS.md) solves this
flag-extension problem on every ordered first-scale stratum and for arbitrary
cluster multiplicities.  Accordingly “every stratum” in the remaining
problem can first be replaced by “every nonradial residue-resonance
stratum.”  The subsequent
[monodromy-forest theorem](POLYNOMIAL_MONODROMY_FORESTS.md) also determines
the dual graphs, component degrees, and node indices there.  The subsequent
[recursive resonance atlas](RECURSIVE_RESONANCE_ATLAS.md) constructs the
flag divisors as normalized initial-form equations and proves their
contraction descent, completing the problem isolated here.

## 1. Two-fiber rigidity

### Theorem 1.1

Let \(E\) be a smooth proper geometrically connected genus-zero curve over a
field \(k\).  Let \(D_0,D_\infty\) be disjoint effective divisors of the
same positive degree.  The set of rational functions \(R\in k(E)^\times\)
with

\[
 \operatorname{div}(R)=D_0-D_\infty                 \tag{1.1}
\]

is either empty or a torsor under \(k^\times\).

If \(q\notin\operatorname{Supp}(D_0+D_\infty)\) is a \(k\)-point, there is
at most one such function satisfying \(R(q)=1\).

### Proof

If \(R\) and \(R'\) both satisfy (1.1), then \(R/R'\) has zero divisor.
It is therefore an everywhere regular invertible function on the proper
geometrically connected curve \(E\), hence belongs to \(k^\times\).
Evaluation at \(q\) fixes that constant. \(\square\)

On an affine coordinate \(W\), write

\[
\begin{aligned}
 D_0&=\sum_i m_i[a_i]+m_\infty[\infty],\\
 D_\infty&=\sum_j n_j[b_j]+n_\infty[\infty].
\end{aligned}
\]

Disjointness means that at most one of \(m_\infty,n_\infty\) is nonzero.
Every solution of (1.1) is

\[
 R(W)=
 c\,
 \frac{\prod_i(W-a_i)^{m_i}}
      {\prod_j(W-b_j)^{n_j}},                         \tag{1.2}
\]

and equality of the total degrees gives

\[
 \deg(\text{denominator})-\deg(\text{numerator})
 =m_\infty-n_\infty.                                 \tag{1.3}
\]

Thus (1.2) also records the required order at infinity.

### Relative form

The same proof works over a connected base \(S\).  If
\(\pi:E\to S\) is proper with geometrically connected fibers and
\(\pi_*\mathcal O_E=\mathcal O_S\), two relative rational functions with
the same horizontal zero and pole divisors differ by a unit of
\(\mathcal O_S^\times\).  A section disjoint from those divisors fixes the
unit.  This is the form needed on a boundary chart: the only freedom left
after two fibers have been constructed is one base unit, not a discrete
Hurwitz choice.

## 2. Intrinsic target formulation

Let \(a,b,c\) be distinct target flags.  The normalized coordinate
\(t_{abc}\) on the target is characterized by

\[
 \operatorname{div}(t_{abc})=[a]-[b],
 \qquad t_{abc}(c)=1.                                \tag{2.1}
\]

Given a source component \(E\), two effective degree-\(d\) divisors
\(D_a,D_b\) determine a rational function \(R\) with
\(\operatorname{div}(R)=D_a-D_b\) up to scale.  A proposed third divisor
\(D_c\) is compatible precisely when the scale can be chosen so that

\[
 \operatorname{div}(R-1)=D_c-D_b.                   \tag{2.2}
\]

Equations (1.1) and (2.2) are a choice-free local realization test.  If they
hold, the component map is the unique map satisfying

\[
 t_{abc}\circ f=R.                                   \tag{2.3}
\]

Changing the ordered triple of flags changes the coordinate by the unique
Mobius transformation carrying the old triple to the new one.  Uniqueness
in Theorem 1.1 shows that the reconstructed map is independent of the
chosen triple whenever all flag divisors are compatible.  Retaining every
target-flag divisor is therefore the simplest permutation-equivariant
presentation; any two fibers and one third-fiber equation are the minimal
local input.

## 3. The flag-complete source enhancement

Let \(P\) be a labelled stable target tree.  A **flag-complete source
enhancement** consists of:

1. a nodal genus-zero source curve \(C\) and an assignment of every source
   component \(E_w\) to a target component \(P_v\);
2. a positive local degree \(d_w\);
3. for every flag \(h\) incident to \(P_v\), an effective degree-\(d_w\)
   divisor \(D_{w,h}\) on \(E_w\), with its ramification multiplicities;
4. the identifications of points in node-flag divisors that form the source
   nodes; and
5. the zero-fiber labels inherited from the root-stable curve.

The divisors for distinct target flags are disjoint.  At a target node, the
two source branches being identified have equal multiplicity.  These are
exactly the combinatorial and divisor parts of admissibility, before writing
the component rational functions.

### Proposition 3.1

Fix a flag-complete source enhancement for which the local equations
(1.1) and (2.2) hold at every source component.  There is at most one
admissible-cover map realizing it.  Its automorphism group is exactly the
group of automorphisms of the enhancement preserving the reconstructed
functions.

### Proof

Every stable target component has at least three flags.  Apply Theorem 1.1
and (2.2) to reconstruct the map on each source component.  The node-flag
identifications make these maps agree at source nodes because both branches
map to the same target node.  Thus they glue uniquely.  An automorphism
preserving the enhancement preserves the two determining divisors and the
third-fiber normalization, so it preserves the reconstructed map.  The
converse is immediate. \(\square\)

The proposition is intentionally compatible with nontrivial deck inertia.
For example, \(U\mapsto U^2\) retains the involution \(U\mapsto-U\) when
the source labels do not kill it.  The theorem does not incorrectly identify
the diagonal phase group with inertia; the character test in the
labelled-node theorem still decides which automorphisms survive.

## 4. Consequence for H1-STACK

Let \(B_N^{\mathrm{tgt}}\) be the
[normalized wonderful pullback](BRANCH_GRAPH_WONDERFUL_PULLBACK.md), and
let \(\mathcal E_N\) denote the stack of compatible flag-complete source
enhancements over it.  Forgetting the component maps gives a morphism

\[
 \Phi:
 \Gamma_N^{\mathrm{ACV}}\longrightarrow\mathcal E_N. \tag{4.1}
\]

Proposition 3.1 says that \(\Phi\) has no residual vertexwise choice: on the
realizable locus it is a monomorphism on objects and preserves the actual
automorphism groups.  Its image is cut out by:

- equality of the two divisor degrees on each source component;
- the third-fiber equations (2.2);
- componentwise Riemann--Hurwitz;
- equality of expansion indices across each source node; and
- the edge-gluing equations.

Consequently the source theorem reduced to the following statement.

> **Polynomial flag-extension theorem.**  
> Over the normalized wonderful branch graph, the labelled root-stable
> polynomial family has a unique dominant flag-complete source enhancement.
> Its local divisor equations are compatible with every nested-set edge
> contraction.  Reconstructing its vertex maps and applying labelled-node
> saturation gives the normalized polynomial ACV graph.

This formulation separates three layers that had previously been mixed:

\[
\boxed{
\begin{array}{c}
\text{wonderful target tree and branch residues}
\\ \Downarrow\\
\text{source tree plus weighted inverse-image divisors of target flags}
\\ \Downarrow\quad\text{Theorem 1.1}\\
\text{unique component rational maps}
\\ \Downarrow\quad\text{labelled-node theorem}\\
\text{normalized branches and actual tame inertia}.
\end{array}}
\]

All four lines are now general theorems.  On the complete radial face poset
the second line is the general radial atlas; on nested tame resonance strata
it is the normalized-initial-form construction of the recursive resonance
atlas.

Tropically, the second line is the harmonic source tree with its expansion
factors and target-flag fibers.  Tropical admissible-cover theory supplies
the compatible combinatorial branch map, but algebraic realizability still
requires the residue equations (2.2); see
Cavalieri--Markwig--Ranganathan,
[Tropicalizing the space of admissible
covers](https://arxiv.org/abs/1401.4626).

## 5. Degree-six recovery

All maps in the degree-six `(2,2,2)` atlas are instances of (1.2).

The equal-scale central component has

\[
 D_0=2[0]+2[1]+2[3],
 \qquad D_\infty=6[\infty],
\]

and hence

\[
 f(W)=c\,W^2(W-1)^2(W-3)^2.                          \tag{5.1}
\]

The normalized target scale fixes \(c=1\).  Thus the source-root
cross-ratio \(3\), which separated the two ambient Hurwitz classes in the
direct enumeration, is simply part of the divisor \(D_0\).  Once that
labelled divisor is retained, the competing class cannot occur.

A quadratic cluster tail has

\[
 D_0=[0]+[1],\qquad D_\infty=2[\infty],
\]

so its map is \(cW(W-1)\).  A cyclic square tail has

\[
 D_0=2[0],\qquad D_\infty=2[\infty],
\]

so its map is \(cW^2\).  These are exactly the natural-tail and Maxwell-tail
maps used in the radial and Maxwell atlases.  The former is anchored by its
two labelled zero points; the latter may retain its square deck involution.

The direct degree-six Hurwitz calculation remains useful as an independent
check, but it is no longer the model for the general proof.

## 6. Optimized next experiment

The degree-six experiment should now be generalized by recording less data,
not more cover equations.

For each nested target stratum:

1. compute the harmonic source tree and expansion indices;
2. record the inverse-image divisors of all incident target flags;
3. reconstruct each source-vertex map from two fibers;
4. test one equation of the form (2.2) for every remaining flag;
5. check edge contraction on the divisor data; and
6. apply the node phase/inertia rule only after the coarse source
   enhancement has passed.

This is the smallest modular certificate that can complete `H1-STACK`.  It avoids
re-enumerating ambient Hurwitz classes and keeps the target fan, algebraic
residue realizability, and stack inertia as separate verifiable layers.

## Reproduction

Run:

```bash
.venv/bin/python scripts/verify_source_vertex_rigidity.py
```

The checker exhausts 2,024 disjoint zero/pole divisor profiles in degrees
one through seven, including finite and infinite support, verifies 2,512
permutation invariances, fixes the scalar with a third flag in every case,
and recovers the three degree-six vertex maps above.
