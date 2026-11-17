# External audit: collision ideals and off-diagonal sheets

This note records the repository's audit of *Collision Ideals and
Off-Diagonal Sheets*.  The manuscript names Chloe van der Vlugt as its author.
The public source and Lean development are published by the GitHub account
[`what-social-construct`](https://github.com/what-social-construct).  No public
source checked in this audit establishes that the named manuscript author and
the GitHub account holder are the same person, so this repository records
those two attributions separately.

This is an external-source audit, not a second copy of the paper's proofs and
not an attribution of the repository's earlier off-diagonal constructions to
that paper.

The source reviewed was
[`what-social-construct/collision-ideals`](https://github.com/what-social-construct/collision-ideals)
at commit
[`a409db9922279907493d96f691b5ea9eb71baaf9`](https://github.com/what-social-construct/collision-ideals/tree/a409db9922279907493d96f691b5ea9eb71baaf9).
The 22-page PDF at that commit has SHA-256
`9ffb46c1f4019b76f18363c3459f0c7e77c955a0a9b132338dd17cf26544b6ba`.
The manuscript is CC BY 4.0 and the Lean source is MIT licensed.

## 1. Credit and provenance

Use the following attribution when citing or reusing this audit.

- **Manuscript:** *Collision Ideals and Off-Diagonal Sheets*, attributed on
  its title page to Chloe van der Vlugt, dated 27 July 2026.
- **Lean source and public repository:** published by the GitHub account
  [`what-social-construct`](https://github.com/what-social-construct/collision-ideals)
  and audited here at the pinned commit above.  Commit metadata is evidence
  of repository publication, not evidence connecting that account to a legal
  or real-world identity.
- **AI assistance disclosed by the manuscript:** OpenAI Codex (GPT-5.6 Sol)
  assisted with Lean proof engineering, mathematical and bibliographic
  checking, and manuscript editing.  The manuscript says that its author
  reviewed the work and takes responsibility for all statements, proofs,
  citations, and formalization code.
- **Literature inputs:** the automorphism implication is credited to James
  Ax's Ax--Grothendieck theorem; the finite-etale and purity interfaces are
  sourced by the manuscript to SGA 1, edited by Alexandre Grothendieck and
  Michèle Raynaud, and the
  [Stacks Project, Tag 0BMB](https://stacks.math.columbia.edu/tag/0BMB);
  the mechanization uses Lean 4, credited in the bibliography to Leonardo de
  Moura and Sebastian Ullrich, and Mathlib, credited to the Mathlib Community.
  These are inputs to or infrastructure for the paper, not coauthorship
  claims.
- **Three-dimensional counterexample context:** the manuscript's motivation
  refers to the explicit map announced by Levent Alpöge.  The announcement is
  reported by the Archive of Formal Proofs as crediting Akhil Mathew for
  prompting the question and Claude Fable for work leading to the map.  The
  separate, independently authored
  [Archive of Formal Proofs verification](https://isa-afp.org/entries/Jacobian_Counterexample.html)
  is by Arthur Freitas Ramos, David Barros Hulak, and Ruy Jose Guerra Barretto
  de Queiroz.  That counterexample and its verification are background
  provenance only: the collision-ideals paper does not claim to discover the
  map, and its structural results do not depend on that explicit map.
- **This audit:** the mathematical review, build reproduction, dependency
  audit, and repository integration recorded here were performed for this
  repository with Codex assistance.  They are not peer review, external
  endorsement, or authorship of the external manuscript or Lean project.

The repository had its own off-diagonal and saturation constructions before
this external source was added.  Reusing the paper's verified interface does
not transfer authorship of that prior work, and similarities in terminology
must not be used to make a priority claim.

## 2. What was checked

Let

\[
 S=\mathbb C[x_1,\ldots,x_n,y_1,\ldots,y_n],\qquad
 I_R=(F_i(x)-F_i(y))_i,\qquad
 I_\Delta=(x_i-y_i)_i.
\]

The following algebraic statements in the paper are correct.

1. `I_R` is contained in `I_Delta`, and
   `Obs(F)=I_Delta/I_R` is the kernel of diagonal evaluation on `S/I_R`.
2. With `A=C[x]` and `B=C[F]`, the collision ring is canonically
   `A tensor_B A`; its spectrum is the affine self-fiber product.
3. Over `C`,
   `Obs(F)=0`, `I_R=I_Delta`, injectivity on complex points, and polynomial
   invertibility are equivalent.  The injective-to-invertible implication is
   Ax--Grothendieck.
4. For a planar Keller map, a secant determinant produces an idempotent
   splitting of the collision ring into diagonal and off-diagonal factors.
   In this clopen situation,

   \[
     I_R:I_\Delta
     =I_R+(\delta_F)
     =I_R:I_\Delta^\infty.
   \]

5. For a dominant generically finite map, generic base change identifies the
   collision algebra with `L tensor_K L`.
6. If `L/K` is a separable nonnormal cubic extension and `N` is its normal
   closure, then, as marked `L`-algebras,

   \[
     L\otimes_KL\simeq L\times N,\qquad
     \operatorname{Gal}(N/K)\simeq S_3.
   \]

The cubic conclusion is a conditional structural statement.  It does not
construct the cubic Keller map; when such a Keller map and extension are
supplied, the nonzero generic factor proves that the map is noninvertible.

## 3. Planar status

The paper defines a hidden-inertia locus on a supplied finite-normalization
diagram.  Keller étaleness puts every conjugate center with positive inertia
index on the deleted boundary, while core-freeness ensures that every
ramified divisor moves some conjugate sheet.  Consequently, for the paper's
planar diagram,

\[
 \mathcal H(\mathcal D_F)=\varnothing
 \quad\Longleftrightarrow\quad
 \operatorname{Ram}^{(1)}(Z/Y)=\varnothing.
\]

Purity and finite-etale rigidity then give the conditional implication

\[
 \mathcal H(\mathcal D_F)=\varnothing
 \Longrightarrow I_R=I_\Delta.
\]

The paper does **not** prove that the hidden-inertia locus is empty for every
planar Keller map.  Thus it gives a clean finite-cover reformulation, not a
new proof or unconditional advance on the plane Jacobian problem.  The
repository's plane degree-three exclusion remains the separate
[Orevkov-based result](../plane-jc/JC2_FINITE_NORMALIZATION_FRONTIER.md).

## 4. Lean audit

At the pinned commit:

- `lake build` succeeds with Lean `v4.24.0` and Mathlib commit
  `f897ebcf72cd16f89ab4577d0c826cd14afaafc7`;
- the source contains no `sorry` or `admit`;
- the planar convenience endpoint exposes exactly three named literature
  axioms: branch purity, finite-etale rigidity of the complex affine plane,
  and planar Ax--Grothendieck;
- `PlanarNoHiddenInertia` is an explicit hypothesis, not an axiom hidden
  inside the model;
- the axiom-parameterized theorem `planarVanishing_of` itself has no declared
  project axiom;
- the cubic endpoint `complexThreeCubicS3Collision` has no declared project
  axiom, but requires the residual-field equivalence and the nontrivial
  fixing subgroup as theorem inputs.

The Lean development therefore machine-checks the implication from its
interfaces.  It does not yet machine-check two manuscript front ends:

1. construction of the complete planar Keller collision model together with
   no hidden inertia; and
2. derivation of the residual-field equivalence and nontrivial fixing
   subgroup solely from the phrase “separable nonnormal cubic.”

Those distinctions must remain visible whenever the formalization is cited.

## 5. Repository use

The paper is useful here as a compact, dimension-independent interface among
collision ideals, tensor self-products, secant projectors, and normal-closure
sheets.  We may use its terminology and cite its pinned formalization for
these interfaces.

The repository retains two stricter conventions.

- Outside the clopen Keller situation, the closure of a genuine
  off-diagonal locus is recorded by saturation, as in
  [the decorated-normalization invariant](../extended-geometry/DECORATED_NORMALIZATION_INVARIANT.md).
  The one-step colon `I_R:I_Delta` should not silently replace saturation in
  an arbitrary non-clopen family.
- The generic factor `N` is a function-field statement.  It must not be
  promoted to an affine global component or to a normalization equivalence
  without the corresponding descent, closure, and normality argument.

No status entry in this repository may use the paper to claim unconditional
planar vanishing.  The reusable output is the verified algebraic interface
and its explicitly conditional formal dependency spine.

## 6. Reproduction

```bash
git clone https://github.com/what-social-construct/collision-ideals.git
cd collision-ideals
git checkout a409db9922279907493d96f691b5ea9eb71baaf9
lake exe cache get
lake build
rg -n '\b(sorry|admit|axiom)\b' CollisionIdeals CollisionIdeals.lean
```

For the endpoint dependency audit, import `CollisionIdeals` in a temporary
Lean file and run `#print axioms` on:

```text
CollisionIdeals.collisionDiagonal_ker
CollisionIdeals.obstructionIdeal_eq_bot_iff
CollisionIdeals.Planar.planarVanishing_of
CollisionIdeals.Planar.planarVanishing_assuming_standardGeometry
CollisionIdeals.Planar.planarAutomorphism_assuming_externalLiterature
CollisionIdeals.complexThreeCubicS3Collision
CollisionIdeals.complexThreeCubicS3HiddenInertiaAt
```
