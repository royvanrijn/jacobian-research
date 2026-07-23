# The source graph by finite normalization

## Result

The wonderful stable-target graph is already the complete **coarse**
polynomial admissible-cover graph.  No further source-side blowup can occur.

Let \(X^\circ\) be the fully labelled polynomial root locus and let

\[
 B_N^{\mathrm{tgt}}
 =
 \operatorname{Norm}
 \overline{\Gamma
   (X^\circ\dashrightarrow\overline M_{0,b+2})}
                                                               \tag{0.1}
\]

be the normalized wonderful target graph.  Let

\[
 \mathrm{br}:\overline{\mathcal H}^{\mathrm{adm,lab}}_N
 \longrightarrow\overline M_{0,b+2}                            \tag{0.2}
\]

be a fully marked tame admissible-cover stack with the polynomial discrete
data.  Its branch morphism is proper with finite geometric fibers, and its
coarse branch morphism is finite.  The stack morphism itself need not be
representable at a boundary cover with relative inertia.  The generic
polynomial determines a section of the pullback of (0.2) over
\(X^\circ\).  Take its reduced closure and normalize:

\[
 \mathcal G_N^{\mathrm{poly}}
 =
 \operatorname{Norm}
 \overline{
   X^\circ
   \subset
   B_N^{\mathrm{tgt}}
   \times_{\overline M_{0,b+2}}
   \overline{\mathcal H}^{\mathrm{adm,lab}}_N
 }.                                                            \tag{0.3}
\]

Then:

> **Finite-normalization theorem.**  
> The coarse morphism
> \[
>  G_N^{\mathrm{poly}}\longrightarrow B_N^{\mathrm{tgt}}
> \]
> is an isomorphism.

Thus the source cover contributes stack structure and a universal admissible
cover over the wonderful graph, but contributes no additional coarse
modification.  In particular:

- degree five needs exactly the normalized blowup of \((x^3,y^2)\);
- degree six needs exactly the six-line/four-center wonderful pullback; and
- in every degree, any missing correction after the wonderful target graph
  is stacky, not another coarse fan subdivision.

This closes the corrected coarse H1 comparison abstractly.  The subsequent
[recursive resonance atlas](RECURSIVE_RESONANCE_ATLAS.md) gives the explicit
logarithmic presentation of \(\mathcal G_N^{\mathrm{poly}}\): flag divisors
as normalized initial-form equations, node-root indices, full simultaneous
inertia characters, and contraction charts.  Thus corrected H1 is complete
both abstractly and explicitly.

## 1. Finiteness input

For fixed degree, genera, ramification profiles, and fully marked branch
points, an admissible cover has no continuous moduli over a fixed stable
target.  The branch morphism (0.2) is proper and has finite geometric
fibers.  On coarse spaces it is proper and quasi-finite, hence finite.  This
is the finiteness used below; representability of (0.2) is not required.
Formation of the coarse space commutes with this base change because the
stack is tame in characteristic zero; see Abramovich--Olsson--Vistoli,
[Tame stacks in positive
characteristic](https://arxiv.org/abs/math/0703310), Corollary 3.3.

The Abramovich--Corti--Vistoli twisted-cover stack gives the normalization of
the Harris--Mumford admissible-cover space; see
[Twisted bundles and admissible
covers](https://arxiv.org/abs/math/0106211).  In the tropical comparison of
Cavalieri--Markwig--Ranganathan, the classical branch morphism has the
corresponding Hurwitz number as its degree and is locally toroidal; see
[Tropicalizing the space of admissible
covers](https://arxiv.org/abs/1401.4626), especially the discussion of the
classical branch map and its node-normalization multiplicities.

Full marking matters.  It removes ambiguities caused by permuting branch
points and source sheets before taking the finite quotients used for H2.

## 2. Proof of the theorem

Base change (0.2) along
\(B_N^{\mathrm{tgt}}\to\overline M_{0,b+2}\).  The resulting stack

\[
 \mathcal Z_N
 =
 B_N^{\mathrm{tgt}}
 \times_{\overline M_{0,b+2}}
 \overline{\mathcal H}^{\mathrm{adm,lab}}_N              \tag{2.1}
\]

is proper over \(B_N^{\mathrm{tgt}}\) and has finite geometric fibers.  Its
coarse morphism is finite.

The normalized polynomial determines one point of the generic Hurwitz fiber,
not the whole Hurwitz fiber.  Its graph over \(X^\circ\) is integral and is
generically isomorphic to \(X^\circ\).  Its reduced closure in
\(\mathcal Z_N\) remains proper with finite geometric fibers over
\(B_N^{\mathrm{tgt}}\).  Normalization is finite because these stacks are of
finite type over a characteristic-zero field.  Consequently the coarse
morphism associated to (0.3) is proper, quasi-finite, and generically of
degree one over \(B_N^{\mathrm{tgt}}\).

Pass to coarse spaces.  The induced morphism

\[
 G_N^{\mathrm{poly}}\longrightarrow B_N^{\mathrm{tgt}}  \tag{2.2}
\]

is finite and birational.  The target is normal by definition.  A finite
birational morphism to a normal integral algebraic space is an isomorphism.
This is the standard uniqueness of normalization
([Stacks Project, Tag 0AB1](https://stacks.math.columbia.edu/tag/0AB1)).
\(\square\)

The argument also explains why the order of construction matters.  If one
uses the nonnormal root-stable graph, a finite birational cover can separate
several branches over one boundary point.  Normalizing the wonderful target
graph first absorbs exactly those branches.  The degree-six fourfold Kummer
normalization over the triple Maxwell point is a concrete instance.

## 3. What survives at stack level

The stack in (0.3) need not be isomorphic to its coarse space.  Boundary
covers may have:

- multiple normalized node phases;
- cyclic deck transformations on unanchored connector components;
- ghost automorphisms of twisted nodes; and
- finite stabilizers induced by label quotients.

These do not contradict the finite-normalization theorem.  They are inertia,
not extra coarse points or blowups.

For one target node with source indices \(e_j\), the normalized branches are

\[
 \left(\prod_j\mu_{e_j}\right)/\Delta\mu_L,
 \qquad L=\operatorname{lcm}(e_j),
\]

and actual inertia is the inverse image of \(\Delta\mu_L\) under the
label-preserving automorphism character.  This is precisely the
[labelled-node saturation theorem](LABELLED_NODE_SATURATION.md).

The
[general radial atlas](GENERAL_RADIAL_SOURCE_ATLAS.md) supplies all radial
indices.  The
[polynomial monodromy-forest theorem](POLYNOMIAL_MONODROMY_FORESTS.md)
supplies all indices on nested simple-branch resonance strata.  The
[source-vertex theorem](SOURCE_VERTEX_RIGIDITY.md) proves uniqueness once
the flag divisors are written explicitly.

The [recursive resonance atlas](RECURSIVE_RESONANCE_ATLAS.md) writes those
flag-divisor residue coordinates and contraction maps and computes the
simultaneous full-centralizer character at higher-codimension nested
resonances.  Thus the stack presentation is complete.

## 4. Corrected H2 is now unconditional

The fully labelled stack \(\mathcal G_N^{\mathrm{poly}}\) is canonical:
closure and normalization commute with every root-label permutation.  Put
\(n=N-2\).  The corrected marked and unmarked stacks are

\[
 \mathcal G_N^{\mathrm{mark}}
 =
 [\mathcal G_N^{\mathrm{poly}}/S_{n-1}],
 \qquad
 \mathcal G_N^{\mathrm{unmark}}
 =
 [\mathcal G_N^{\mathrm{poly}}/S_n].                 \tag{4.1}
\]

Subgroup inclusion gives

\[
 \mathcal G_N^{\mathrm{mark}}
 \longrightarrow
 \mathcal G_N^{\mathrm{unmark}}                      \tag{4.2}
\]

finite, representable, etale, and of degree \(n=N-2\).  It is the
normalization in the generic selected-root algebra.  No H1 identification
is now being assumed: (0.3) constructs the repository-specific corrected
graph, while the finite-normalization theorem identifies its coarse base.

Thus corrected H2 is complete.  The old statement over the unmodified
root-stable quotient remains false.

## 5. Corrected coarse H3

The marked coarse graph is normal.  The contracted selected-root incidence
is finite over it, and the closure of the generically distinguished
degree-one point is finite and birational.  It is therefore an isomorphism
by the same normality lemma.  The affine mark has a unique coarse
specialization over every DVR limit.

Thus corrected H3 is also unconditional over (0.3).  This says nothing about
properness of an unrelated relative isomorphism functor; it is only the
finite selected-root statement.

## 6. Simplified completion ledger

The H1/H2 package now separates cleanly:

| Layer | Status |
|---|---|
| coarse target correction | normalized wonderful pullback; complete |
| coarse polynomial source graph | equal to the target graph by finite normalization; complete |
| radial source components and nodes | complete for arbitrary multiplicities |
| simple-branch resonance source trees and nodes | complete by monodromy subforests |
| source-component uniqueness | complete by two-fiber rigidity |
| local node normalization and inertia test | complete |
| generic radial/resonance automorphism characters | complete |
| explicit higher-codimension flag coordinates and simultaneous characters | complete by recursive normalized initial forms and full-centralizer matching |
| corrected H2 subgroup quotient | complete |
| corrected coarse H3 specialization | complete |

This closes the H1/H2 compactification ledger.  Rechecking coarse source
selection or H2 label gluing would duplicate formal consequences of finite
normalization.
