# From one collision to marked-root Keller maps

This repository began with a startling polynomial map circulated publicly in
July 2026 as a three-dimensional counterexample to the Jacobian conjecture.
The first task was modest: check the determinant and the reported collision
exactly.  Those identities held.  The next question—*why does this map
work?*—led from a cubic with one root marked to normalized root covers,
boundary divisors, cancellation constructions, and finally many stably
inequivalent maps in every generic degree.

The repository is both a proof project and a research diary that has now been
cleaned into a small active theorem chain.  Mathematical completion and
external review are deliberately kept separate in [STATUS.md](STATUS.md).

## Where the example came from

The earliest public item located by the repository’s provenance audit is
[a post by Levent Alpöge](https://x.com/__alpoge__/status/2079028340955197566).
Contemporaneous accounts attribute the example to Alpöge and Fable; [Zihan
Zhang’s source audit](https://zzhang-iu.github.io/papers/direct-consequences-jacobian/index.html)
additionally reports that Akhil posed the question and Fable produced the
example.  [Alexis Gallagher’s same-day
explainer](https://jacobianfun.org/jacobian-explained) credits “Levent Alpöge
+ Fable” and develops a weighted-lift interpretation.

The direct announcement has not been recovered in a form that settles the
complete discovery history, original prompt, model conversation, or exact UTC
timestamp.  This repository therefore repeats the contemporaneous attribution
without making a priority claim.  The detailed record is in the
[provenance audit](archive/legacy-notes/PROVENANCE_AUDIT.md).

## The counterexample

Let `u=1+xy` and define `F:A^3 -> A^3` by

\[
F(x,y,z)=\left(
u^3z+y^2u(4+3xy),
y+3xu^2z+3xy^2(4+3xy),
2x-3x^2y-x^3z
\right).
\]

On `x!=0`, put

\[
t=y+1/x,\qquad r=2/x,
\]

and write the target as `(a,b,c)`.  Then

\[
a=t^2+rt/2-ct^3,
\qquad
b=r+4t-3ct^2.
\]

The two coordinate Jacobians multiply to

\[
\det\frac{\partial(a,b,c)}{\partial(t,r,c)}
\det\frac{\partial(t,r,c)}{\partial(x,y,z)}
=\frac r2(-2x)=-2.
\]

Since this is a polynomial identity,

\[
\boxed{\det DF=-2}
\]

everywhere.  Exact substitution also gives

\[
F(0,0,-1/4)=F(1,-3/2,13/2)=F(-1,3/2,13/2)=(-1/4,0,0).
\]

Thus the map is everywhere etale and noninjective.  Over `C` it is a
counterexample to the Jacobian conjecture in dimension three; adjoining
identity coordinates gives examples in every higher dimension.

Dean Cureton independently formalized the determinant, the three evaluations,
a determinant-one rescaling, and the complex specialization in
[Lean 4](verified/LEAN_FOUNDATIONAL_MAP.md).  His formalization is separately authored and
pinned rather than copied into this repository.

## The geometric turn

The counterexample has a shorter coordinate-free construction.  Start from

\[
\mathbb P^1\times\operatorname{Sym}^2(\mathbb P^1)
\longrightarrow\operatorname{Sym}^3(\mathbb P^1),\qquad (x,Q)\mapsto x+Q,
\]

and remove its ramification divisor together with the inverse image of a
tangent non-osculating hyperplane.  The target complement is `A^3`.  The two
boundary lines in each `Sym^2(P^1)` fiber become chords of an auxiliary
conic; residual intersection makes the source an iterated affine-line bundle,
hence `A^3`.  The restricted addition map is étale and generically
three-to-one, so this already proves the counterexample without coordinates.

Coordinates turn that construction into the binary cubic

\[
Q_{a,b,c}(U,V)=cU^3-2U^2V+bUV^2-2aV^3.
\]

The source is the incidence space of this cubic together with one *simple
projective root marked*.  The map `F` forgets the marking.  This immediately
explains why a generic fiber has three points, why a double-root cubic leaves
one affine point, and why a triple-root cubic is omitted.

The marked-projective-root organization was publicly described by
[Andy Jiang](https://x.com/davikrehalt/status/2079175065695035442) and linked
in a same-day MathOverflow update by Qiaochu Yuan.  The repository’s earlier
affine primitive-element calculations reached overlapping reconstruction,
discriminant, and fiber formulas.  It adopts Jiang’s projective organization
with attribution; the repository contribution is the global two-chart proof,
including the root at infinity and the identification of the affine source as
the regular-reconstruction open.

The auxiliary-conic proof comes from a later geometric manuscript prompted by
Semon Rezchikov.  The repository now uses it as the primary existence proof;
the two charts remain the explicit coordinate certificate, including the root
at infinity.  This simplification affects the foundational cubic only and
adds no new cancellation, boundary-invariant, or stable-equivalence result.

This cubic picture became the template:

\[
\text{finite root cover}
\supset
\text{regular-reconstruction open}
\cong
\mathbb A^3.
\]

Roots give generic source points.  Repeated roots produce discriminant
inertia.  When reconstruction develops a pole, the corresponding sheet moves
to the canonical normalization boundary.

## Parallel discoveries and repository contributions

The project developed during a burst of same-day public work.  The following
distinctions matter:

| Topic | Public or independent source already located | What this repository contributes |
|---|---|---|
| Explicit map, determinant, and rational collision | Alpöge/Fable announcement trail | Independent exact implementations, a compact proof, and a maintained reproduction target; no discovery claim |
| Marked projective root | Andy Jiang’s public geometric interpretation | The global two-chart scheme proof, including the root at infinity, plus exact image and nonproperness arguments |
| Auxiliary conic and residual-intersection coordinate | Rezchikov-prompted geometric manuscript | Integration with the marked-root model while retaining explicit coordinate charts and downstream scope boundaries |
| Weighted lift and higher-degree seed family | Alexis Gallagher’s explainer | A normalized-incidence formulation, the regular-reconstruction open, uniform `S_N` monodromy, and later complete boundary calculations |
| Quartic Islands A/B/C | Juntang Zhuang's pinned compilation and exact collision checkers | Independent compact reconstruction showing one canonical and two split weighted seeds, complete boundary signatures, and exclusion from the cancellation normal forms |
| Formal verification of the foundational map | Dean Cureton’s Lean project | Pinned attribution, scope audit, and an optional upstream build target |
| Immediate consequences and source trail | Zihan Zhang’s audit | A separate proof architecture focused on inverse geometry, normalization, and stable equivalence |
| GMC, `(xz)`, and `SU(2)` counterexamples | Christopher D. Long's two arXiv papers | Exact local Gaussian and beta/binomial checks, a complete `SU(2)=S^3` Haar proof, all 18 BCW steps through an independently replayed 79-variable sparse artifact, and a local proof of the fixed-dimensional nonexplicit `GMC(158)` route; no claim that the direct witnesses derive from our map |
| Weighted-seed/Gaussian bridge | Long's Lagrange--Good search architecture and Good's inversion theorem | A repository-derived polynomial determinant correction turning every nonconstant normalized seed into an explicit four-real-Gaussian witness family; locally proved and checked, not externally reviewed |
| Cancellation maps and canonical boundary invariants | No earlier source has been identified in the present audit | The finite cancellation operator, all-parameter reconstruction and collision, boundary-exhaustion theorem, thick intersections, and rigidity results |
| Many stable classes in each degree | No earlier source has been identified in the present audit | The divisor-count theorem giving at least `tau(N-1)` stable classes in generic degree `N>=4` |

“No earlier source identified” is a statement about the current search, not a
claim that no such source exists.  The cancellation and boundary results are
written proofs in this repository but have not received external specialist
review.

## From one map to many

The weighted construction replaces the cubic by

\[
H(W)-BCW+cAC^2
\]

and produces simply ramified `S_N` root covers.  The cancellation construction
uses

\[
C\int_0^T\{1-t(Q-Pt)^m\}^r\,dt-R
\]

and produces a distinguished boundary ramification index `r+1`.

For a fixed generic degree `N`, every proper divisor `r|(N-1)` gives

\[
m=\frac{N-1}{r}-1\ge1.
\]

The resulting cancellation types are pairwise stably inequivalent: their
canonical boundary data contain the ramification index `r+1` and nilpotency
index

\[
m(N-1).
\]

A weighted marked-root map supplies one additional class because its canonical
boundary intersection is reduced.  Therefore:

\[
\boxed{\text{For every }N\ge4,\text{ there are at least }\tau(N-1)
\text{ pairwise stably inequivalent degree-}N\text{ maps.}}
\]

The canonical statement and proof are in the standalone paper
[Marked-Root Keller Maps and Degreewise Stable Multiplicity](papers/marked-root-multiplicity/main.tex).
[Marked-root Keller maps](MARKED_ROOT_KELLER_MAPS.md) supplies the common
framework, while the [five-lemma audit](DEGREEWISE_MULTIPLICITY_AUDIT.md) is a
verification companion.  Neither is a competing canonical theorem source.

Three externally published quartic maps called Islands A, B, and C are
classified in the [quartic-islands audit](extended-geometry/EXTERNAL_QUARTIC_ISLANDS.md).
They are not new cancellation classes: one is the canonical triple-zero
weighted seed and two are split weighted seeds with normalized extra roots
`3` and `-1/2`.  Their boundary data separate Island A from the other two,
and the scheme-theoretic decorated-normalization invariant separates the two
split-root parameters from each other.

## Internal geometry and external consequences

The marked-root, boundary, decorated-normalization, and stable-moduli results
are the repository's internal geometric programme.  Christopher D. Long's
direct Gaussian-moment, `(xz)`, and `SU(2)` counterexamples are independent
external consequences discovered after JC(3).  The distinction is documented
in the [external-consequences and provenance note](extended-geometry/EXTERNAL_CONSEQUENCES_AND_PROVENANCE.md).

Long's JC-to-GMC route does use an exact diagonal normalization of the
foundational map.  The repository now reproduces its 18 stable steps, the
dimensions `3 -> 39 -> 79`, the cubic-homogeneous collision, and the complete
fixed-dimensional implication to `not GMC(158)`.  The proof architecture and
authorship remain credited to BCW, DVEZ, and Zhao.  Long's
much smaller direct witnesses were motivated by the JC(3) announcement but
were not algebraically extracted from that map.  His papers do not review the
repository's Hessian/Fitting-divisor or stable-moduli proofs.

There is now a separate internal bridge from the weighted inverse pencil to
Gaussian moments.  For every nonconstant normalized seed, a polynomial
two-complex-coordinate correction cancels the Lagrange--Good determinant and
produces an explicit `GMC(4)` witness family whose mixed moments read one exact
branch of `H(W)-sW+t=0`.  This construction was prompted by Long's method but
is not one of Long's witnesses and has not received external specialist
review.  See the [weighted Gaussian bridge](extended-geometry/WEIGHTED_GAUSSIAN_BRIDGE.md).

There is also an independent internal closure of the supplied rank-two
Poisson gap.  An exact flux calculation forces the shear `Z -> Z-9Q^2` and
produces compact formulas for `R,T,D,S` with all six canonical brackets, a
complete three-point fiber, and determinant one.  This gives a noninjective
exact symplectic map of `A^4`, a counterexample to the Poisson conjecture for
two canonical pairs, and an inverse-Jacobian `A_4` Weyl consequence.  The
unidentified external abstract remains under provenance audit: the repository
does not claim that its independently derived formulas are the manuscript's
formulas or that the abstract reviews this work.  See the
[rank-two construction and audit](extended-geometry/QUADRATIC_LADDER_AND_POISSON_AUDIT.md).

## What is genuinely new here?

Subject to the provenance limitation above, the repository’s clearest new
mathematics is not the original map or the first marked-root and weighted
ideas.  It is the theory built after them:

- the global normalized-cover formulation with the affine source identified
  as the regular-reconstruction open;
- exact image and nonproperness theorems with denominator-safe projective
  reconstruction;
- functorial reduced, scheme-theoretic, and formal boundary invariants that
  survive polynomial left–right equivalence and stabilization;
- the uniform cancellation construction and its arithmetic parameter family;
- complete normalization-boundary exhaustion for the weighted and
  cancellation maps;
- the uniform weighted-seed/Gaussian bridge, including its polynomial
  determinant correction and exact mixed-moment formula;
- reduced versus nilpotently thick boundary intersections; and
- the divisor-count lower bound `tau(N-1)` for stable-equivalence classes in
  every generic degree `N>=4`.

These claims should be read with [STATUS.md](STATUS.md): “proved in the
repository” and “externally reviewed” are different assertions.

## Reading and reproduction

The foundational map, cubic marked-root model, exact image theorem, and
weighted theorem are the stable core. See [STATUS.md](STATUS.md) for detailed
evidence levels.

The stable-core reading path is:

1. [Foundational Keller map](verified/FOUNDATIONAL_GEOMETRY.md)
2. [Cubic marked-root model](verified/MARKED_ROOT_MODEL.md)
3. [Exact image and nonproperness](verified/IMAGE_AND_NONPROPERNESS.md)
4. [Weighted marked-root theorem](verified/WEIGHTED_SEED_THEOREM.md)

Further families continue in one sequence:
[cancellation construction](cancellation/CONSTRUCTION.md),
[boundary geometry](cancellation/BOUNDARY_GEOMETRY.md), and the
[marked-root framework](MARKED_ROOT_KELLER_MAPS.md), followed by the
[canonical degreewise paper](papers/marked-root-multiplicity/main.tex).

For execution, start with [REPRODUCE.md](REPRODUCE.md).  The core paper is
[papers/core-counterexample/main.tex](papers/core-counterexample/main.tex).
Extended geometry is indexed in
[extended-geometry](extended-geometry/README.md); cancellation
arithmetic, rigidity, and open problems live beside the construction under
`cancellation/`.  Superseded derivations and exploratory tools
remain available through the [archive](archive/README.md).

## Credits and sources

This repository explicitly thanks and credits:

- Levent Alpöge and Fable for the public counterexample trail, following the
  attribution in contemporaneous sources;
- Akhil for posing the question, as reported by Zihan Zhang;
- Andy Jiang for the marked-projective-root organization;
- Alexis Gallagher for the explanatory weighted-lift and seed-family account;
- Juntang Zhuang for the pinned public compilation, Island A/B/C labels,
  expanded maps, and exact quartic collision certificates;
- Dean Cureton for the independently authored Lean formalization;
- Christopher D. Long for the independently authored GMC, `(xz)`, and
  `SU(2)` counterexamples and for the explicit fixed-dimensional GMC route;
- Zihan Zhang for the source and direct-consequence audit; and
- Qiaochu Yuan for connecting the public geometric interpretation to the
  MathOverflow discussion.

Classical and technical inputs—including coincident-root geometry,
Zariski–Nagata purity, normalization, and standard Jacobian-conjecture
reductions—are listed in [SOURCES.md](archive/legacy-notes/SOURCES.md).  If a
missing earlier source or more precise discovery record is found, the
provenance audit and this README should be corrected.
