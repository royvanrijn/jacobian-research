# From one collision to marked-root Keller maps

This repository began with a startling polynomial map circulated publicly in
July 2026 as a three-dimensional counterexample to the Jacobian conjecture.
The first task was modest: check the determinant and the reported collision
exactly.  Those identities held.  The next question—*why does this map
work?*—led from a cubic with one root marked to normalized root covers,
boundary divisors, cancellation constructions, and finally many stably
inequivalent maps in every generic degree.  The strongest degreewise result
is now positive-dimensional: in generic degree `N>=4`, the weighted locus
contains an `(N-3)`-dimensional family of stable classes.

## What is new here: result snapshot

> **Recorded 2026-07-22T19:10:02+02:00 (CEST).**  This timestamp fixes the
> state of the claims summarized here; it is not a priority date.  “New here”
> means that no earlier source has been identified in the current provenance
> audit, not that priority has been established.  “Proved” means proved in the
> repository, not externally reviewed; [STATUS.md](STATUS.md) is the detailed
> evidence ledger.

The per-result dates below are CEST Git commit timestamps for the first commit
containing the stated theorem in its present form.  Where a result was
completed in stages, each material stage is dated.  A working-tree date is
identified as uncommitted and therefore is not represented as part of Git
history.

The original `JC(3)` map, the first marked-projective-root interpretation,
and the first weighted-lift idea are credited to the public and independent
sources below.  The main repository-first results are the theory and new
counterexamples built from them:

| Repository-first result | Added to repository | What was proved here | Why it matters |
|---|---|---|---|
| **Positive-dimensional stable moduli in every degree** | Stable theorem: `2026-07-22T13:40:46+02:00` ([`82ded4b`](https://github.com/royvanrijn/jacobian-research/commit/82ded4b2ecb59308266207aa537234073d0bd1f2)); symplectic transfer: `2026-07-22T13:45:04+02:00` ([`2b52cb9`](https://github.com/royvanrijn/jacobian-research/commit/2b52cb920ae4d714c5f734595077c2a7560ff96b)) | For every generic degree `N>=4`, weighted Keller maps contain an `(N-3)`-dimensional family of stable polynomial left--right classes; the same dimension transfers to exact symplectic maps of `A^6`. | This replaces a list of examples by moduli growing linearly with degree.  See the [decorated-normalization invariant](extended-geometry/DECORATED_NORMALIZATION_INVARIANT.md). |
| **A second, cancellation family** | Degreewise count: `2026-07-21T18:57:58+02:00` ([`07eec55`](https://github.com/royvanrijn/jacobian-research/commit/07eec55222b5e8dce98a63f7dd6b156f3db470fb)); active construction note: `2026-07-21T19:20:10+02:00` ([`38ce4e0`](https://github.com/royvanrijn/jacobian-research/commit/38ce4e04ce2a83b8dc1e0097e276cd09a7d6931b)) | A uniform finite cancellation operator gives noninjective constant-Jacobian maps for all `m,r>=1`, with canonical thick boundary data.  Together with the weighted family it gives at least `tau(N-1)` pairwise stably inequivalent maps in every generic degree `N>=4`. | The bound counts genuinely different construction types, not just parameters inside one family.  See the [construction](cancellation/CONSTRUCTION.md) and [degreewise theorem](papers/marked-root-multiplicity/main.tex). |
| **The first higher factorization slice classified** | Tangent slice: `2026-07-22T15:39:43+02:00` ([`5b2a581`](https://github.com/royvanrijn/jacobian-research/commit/5b2a581541f67dd9027dfd74f3166cc809b07114)); all hyperplanes: `2026-07-22T17:00:04+02:00` ([`a47a60a`](https://github.com/royvanrijn/jacobian-research/commit/a47a60a33214b6e0298888e33debdbd930ba8b84)); integral topology: working tree `2026-07-22T19:05:13+02:00`, not yet committed | No quadratic--cubic binary-form hyperplane complement is `A^5`.  The natural tangent slice is a new nonproper etale fivefold, homotopy equivalent to `S^3`, with vanishing positive-codimension integral Chow groups and `K_0=Z`. | It explains why the next normalized factorization does not automatically yield `JC(5)` while isolating a new cover whose remaining product question is algebraic, not topological.  See the [`(2,3)` slice audit](extended-geometry/QUADRATIC_CUBIC_FACTORIZATION_SLICE.md). |
| **Uniform direct `GMC(4)` counterexamples** | `2026-07-22T13:12:08+02:00` ([`b341fd9`](https://github.com/royvanrijn/jacobian-research/commit/b341fd91e91f9eeea6105d9b3b289ffd22c23f04)) | Every nonconstant normalized weighted seed produces an explicit four-real-Gaussian witness family; its full mixed-moment sequence recovers the seed. | This is a new seed-to-moment construction, separate from Long's smaller independent `GMC(3)` witness.  See the [Gaussian bridge](extended-geometry/WEIGHTED_GAUSSIAN_BRIDGE.md). |
| **Explicit Image and Vanishing counterexamples** | `2026-07-22T18:58:08+02:00` ([`1f06762`](https://github.com/royvanrijn/jacobian-research/commit/1f067624d59b708a661d7284c5613e6f195fa53d)) | A 21-variable cubic-homogeneous Keller collision gives `not SIC(21)` and an explicit 42-variable quartic giving both `not GVC_Delta(42)` and `not VC_(HN,4)(42)`. | The collision supplies an all-orders certificate and a specified inverse coordinate, not a bounded expansion.  See [Image and Vanishing counterexamples](extended-geometry/IMAGE_VANISHING_COUNTEREXAMPLES.md). |
| **A smaller JC-to-GMC route bound** | `2026-07-22T18:58:08+02:00` ([`1f06762`](https://github.com/royvanrijn/jacobian-research/commit/1f067624d59b708a661d7284c5613e6f195fa53d)) | Rank compression and constant-kernel quotienting reduce the locally reproduced BCW route to a 21-variable cubic-homogeneous collision, hence `not GMC(42)`. | This improves the same route from `not GMC(158)`; it is a route-based dimension bound, not a direct Gaussian witness or a minimality claim.  See the [reproduction and optimization](extended-geometry/LONG_SU2_AND_BCW_REPRODUCTIONS.md). |
| **Rank-two symplectic/Poisson descent** | `2026-07-22T13:35:11+02:00` ([`3f1452e`](https://github.com/royvanrijn/jacobian-research/commit/3f1452edc97655560a2cdd1eff36e0761c28169e)) | An explicit noninjective exact symplectic map of `A^4`, with all six brackets and a rational three-point fiber, gives the rank-two Poisson/`DC(4)` counterexample and an inverse-Jacobian `A_4` Weyl consequence. | It descends the foundational phenomenon from three canonical pairs to two; the formulas were derived here independently of an unavailable external abstract.  See the [rank-two audit](extended-geometry/QUADRATIC_LADDER_AND_POISSON_AUDIT.md). |

Two other sharp numbers should be read with different provenance.  The
repository locally reproduces the external Newton-polygon reduction giving
`max(deg P,deg Q)>=125` for any hypothetical plane counterexample, subject to
the published minimal/standard normal form; this is not `JC(2)`, does not say
both degrees are at least 125, and does not assert attainability.  This result
was added at `2026-07-22T16:28:54+02:00`
([`eb08867`](https://github.com/royvanrijn/jacobian-research/commit/eb088670ab2cadb5f94625880085a3f019d86605)).
Conversely, Long's independent direct witness proves the stronger dimension
statement `not GMC(n)` for every `n>=3`; it was first documented here at
`2026-07-22T13:12:08+02:00`
([`b341fd9`](https://github.com/royvanrijn/jacobian-research/commit/b341fd91e91f9eeea6105d9b3b289ffd22c23f04)).
The repository's novelty is the uniform
seed-parametrized `GMC(4)` family and the separate JC-to-GMC route compression,
not Long's result.  See the [plane frontier](plane-jc/DEGREE_FRONTIER_125.md)
and [external-consequences provenance](extended-geometry/EXTERNAL_CONSEQUENCES_AND_PROVENANCE.md).

The repository is both a proof project and a research diary that has now been
cleaned into a small active theorem chain.  Mathematical completion and
external review are deliberately kept separate in [STATUS.md](STATUS.md).

The plane problem is tracked separately in the [JC(2) constraint
program](plane-jc/README.md).  Its current scoped result is an external
Newton-polygon reduction with a local exact reproduction showing that any
hypothetical plane counterexample has larger coordinate degree at least 125,
subject to the published minimal/standard normal form.  This neither proves
JC(2) nor follows from the three-dimensional construction.

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

Between that projective proof and the expanded polynomial lies a short
equivariant coordinate certificate.  Write

\[
 L=aT+bS,\qquad Q=cT^2+dTS+eS^2
\]

and normalize the factor pair by

\[
 X=\{a^2e-abd+b^2c=1,\ ad+bc=1\}\subset\mathbb A^5.
\]

The following formulas are a polynomial isomorphism
`A^3_(a,y,z) -> X`:

\[
\begin{aligned}
b&=1+ay,\\
c&=1-\frac32ay+a^2z,\\
d&=\frac12y-az+\frac32ay^2-a^2yz,\\
e&=-2z+4y^2-4ayz+3ay^3-2a^2y^2z,
\end{aligned}
\]

with polynomial inverse

\[
y=2bd-ae,\qquad
z=2d^2+ce+6bd^2+3bce-\frac92e.
\]

No division by `a` occurs, so this is a global certificate across the
root-at-infinity divisor `a=0`.  It intertwines the residual torus actions

\[
(a,y,z)\mapsto(\lambda a,\lambda^{-1}y,\lambda^{-2}z),\qquad
(a,b,c,d,e)\mapsto(\lambda a,b,c,\lambda^{-1}d,\lambda^{-2}e).
\]

Dropping the fixed product coefficient `ad+bc=1`, multiplication is

\[
G(a,y,z)=(ac,ae+bd,be),\qquad \det DG=-1.
\]

For

\[
A(z_1,z_2,z_3)=(z_1,z_2,-z_3/2),\qquad
B(u_1,u_2,u_3)=(u_3,2u_2,2u_1),
\]

the announced polynomial is exactly

\[
\boxed{F_{\rm original}=B\circ G\circ A.}
\]

Thus the algebraic route is

\[
\boxed{\text{coprime factorization}
\Rightarrow\text{étale normalized multiplication}
\Rightarrow\text{exceptional affine slice}
\Rightarrow\text{announced polynomial}.}
\]

Here “exceptional” means tangent but nonosculating contact type `(2,1)` with
the twisted cubic.  The three mechanisms are distinct:

\[
\begin{array}{rcl}
\text{étaleness}&:&\text{universal coefficient--resultant geometry},\\
\text{noninjectivity}&:&\text{three choices of a linear factor},\\
X\simeq\mathbb A^3&:&\text{the special }(2,1)\text{ hyperplane orbit}.
\end{array}
\]

The complete proof and exact symbolic certificate are in the
[three-proposition factorization model](verified/NORMALIZED_FACTORIZATION_MODEL.md).
The hyperplane-independent argument is isolated as three lemmas in the
[foundational incidence construction](verified/FOUNDATIONAL_INCIDENCE_CONSTRUCTION.md).

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
Semon Rezchikov.  The repository uses it as the projective existence proof;
the normalized factorization model is the intermediate algebraic certificate,
and the two marked-root charts remain the direct reconstruction certificate.
Together they include the root at infinity.  This simplification affects the
foundational cubic only and adds no new cancellation, boundary-invariant, or
stable-equivalence result.

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

For completeness, the étaleness step in that route is intrinsic.  For a
linear form `L` and a quadratic form `Q`, multiplication has only the
infinitesimal relative-scaling kernel `(L,-Q)` on the coprime locus.  The
resultant has bidegree `(2,1)` and is nonconstant in exactly that direction.
Hence `(L,Q) -> (LQ,Res(L,Q))` is étale, and the displayed threefold map is
its normalized base change.

The same factorization mechanism is étale for all unequal factor degrees and
has a unique normalization for consecutive degrees.  The first new case,
quadratic times cubic, is nevertheless not a new `JC(5)` example: its natural
normalized source is a smooth factorial fivefold with trivial Picard and
canonical classes, but its class is `L^5-L^3` rather than `L^5` (equivalently,
its finite-field count is `q^5-q^3`).  The visible determinant-one projection
does not split it: its identity fiber is `x^2+zy^2=1`, of class `L^2+L`, and
the associated Euclidean-addition distribution jumps rank along `a_0=0`.
This rules out the natural product charts, but not an abstract isomorphism
which mixes all coefficient functions.  A generic-plane/nodal-complement
argument now proves that the complex fivefold is simply connected.  The top
boundary homology of the incidence normalization also gives
`H^2=0` and `H^3=Q(-2)` without a full normal-crossings resolution.  A
two-chart integral Gysin calculation finishes the topology: `H^*(X,Z)` is
`Z` in degrees zero and three and vanishes otherwise, so the simply connected
complex fivefold is homotopy equivalent to `S^3`.  Its positive-codimension
integral Chow groups vanish and its algebraic `K_0` is `Z`.  The Euclidean
chart is an explicit affine modification of `A^2 x SL_2`, with a two-component
center and exceptional divisor, but that modification morphism has poles on
the complementary coefficient chart.  Thus the remaining product question
is algebraic rather than topological.  The
[`(2,3)` slice audit](extended-geometry/QUADRATIC_CUBIC_FACTORIZATION_SLICE.md)
records the calculation and the precise remaining modification and
additive-invariant questions.

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
| GMC, `(xz)`, and `SU(2)` counterexamples | Christopher D. Long's two arXiv papers | Exact local Gaussian and beta/binomial checks, a complete `SU(2)=S^3` Haar proof, all 18 steps of Long's conservative 79-variable route, and a repository 17-dimensional quadratic--cubic trace whose rank-compressed 24-variable homogenization has a 21-variable essential quotient and gives `not GMC(42)`; every sparse artifact has an independent replay, with no claim that the direct witnesses derive from our map |
| Weighted-seed/Gaussian bridge | Long's Lagrange--Good search architecture and Good's inversion theorem | A repository-derived polynomial determinant correction turning every nonconstant normalized seed into an explicit four-real-Gaussian witness family; the exact mixed moments recover `1+lambda*H`, making the moment realization injective; supported by a standalone constant-term-safe formal Gaussian--Lagrange proof, locally checked and not externally reviewed |
| Cancellation maps and canonical boundary invariants | No earlier source has been identified in the present audit | The finite cancellation operator, all-parameter reconstruction and collision, boundary-exhaustion theorem, thick intersections, and rigidity results |
| Degreewise weighted stable moduli | No earlier source has been identified in the present audit | A generically finite decorated-normalization invariant with image dimension `N-3`, on an explicitly proved nonempty ordinary boundary-clean open for every generic degree `N>=4` |
| Cross-family stable multiplicity | No earlier source has been identified in the present audit | The complementary divisor-count theorem distinguishing one weighted class and `tau(N-1)-1` cancellation types in generic degree `N>=4` |

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

Every weighted seed has a two-dimensional tangent-map core

\[
(W,\gamma)\longmapsto
\bigl(H'(W)+c\gamma,\,W(H'(W)+c\gamma)-H(W)\bigr),
\qquad J=-c^2\gamma.
\]

On `gamma=0` this is exactly the discriminant normalization, its inverse
equation is `H(W)-sW+t=0`, and its relative Fitting divisor is `(H'')`.  The
three-dimensional Keller map is a weighted suspension of this ramified plane
map: `C=x gamma` supplies the complementary boundary factors and leaves the
constant Jacobian `b_0c`.  The determinant, inverse pencil, reconstruction
poles, discriminant normalization, and Hessian divisor are unified in the
[tangent-map core theorem](verified/TANGENT_MAP_CORE.md).

For cubicization, the
[constant-kernel quotient theorem](verified/CONSTANT_KERNEL_QUOTIENT.md)
identifies every constant translation-invariance direction of the homogeneous
part, conjugates the map to a triangular extension of its essential-input
quotient, and preserves the full fiber scheme.  The theorem is GZ-type in its
use of quotient/section matrices, and its search protocol now scores completed
BCW traces only after rank compression, homogenization, and constant-kernel
quotienting.

After the triangular change `s=H'(W)+c gamma`, the core is simply

\[
(W,s)\longmapsto(s,Ws-H(W)),
\]

the projection of the smooth universal incidence `H(W)-sW+t=0`.  The raw
weighted incidence is its ramified base change by
`(A,B,C) -> (BC,cAC^2)`.

The cancellation maps realize a second suspension mechanism.  For fixed `P`,

\[
(s,Q)\longmapsto
\left(Q,C\int_0^s\{1-t(Q-Pt)^m\}^r\,dt\right)
\]

has Jacobian `-C D^r`, where `D=1-s(Q-Ps)^m`; the birational source chart has
Jacobian `-D^{-r}` and cancels this boundary power.  Thus weighted maps use a
polynomial suspension of simple ramification, while cancellation maps use a
birational suspension of higher boundary ramification.

For the foundational weights this suspension also has an exact
Poisson-square normal form.  In oriented invariant-plane coordinates,
`P=C^2A` and `Q=CB` satisfy `[P/2,Q]=C^2`; the normalized sixteen-monomial
coefficient scheme becomes three univariate weighted-Wronskian layers.  This
is the same band-bracket equation used in the plane `(72,108)` cascade, with
different supports and a different forced layer.  See the
[weighted tangent-suspension bridge](extended-geometry/WEIGHTED_TANGENT_SUSPENSION.md).

The normalized admissible degree-`N` weighted seed space has dimension
`N-3`.  For every `N>=4`, its ordinary boundary-clean locus is nonempty: the
generic discriminant theorem supplies the ordinary nodal-cuspidal open, while

\[
\frac12W^2(1-W)(1+W^{N-3})
\]

explicitly witnesses boundary-cleanness and keeps the boundary mark away from
the cusp and node-branch schemes.  On their nonempty intersection, the
decorated-normalization invariant is generically finite.  Since it is
preserved by stable polynomial left--right equivalence,

\[
\boxed{\text{For every }N\ge4,\text{ weighted degree-}N\text{ maps contain an }
(N-3)\text{-dimensional family of stable classes.}}
\]

This is the main degreewise moduli result.  Its proof is in the
[decorated-normalization theorem](extended-geometry/DECORATED_NORMALIZATION_INVARIANT.md).
The simple `W=1` sheet is intrinsic as an affine component over the second
boundary divisor, but it is disjoint from the discriminant ramification
divisor.  A cross-stratum generator-rigidity lemma would be needed to turn it
into a `r=1` mark on the discriminant normalization and upgrade generic
finiteness to generic injectivity; that lemma is not presently proved.
The exact cotangent lift is polynomially right-equivalent to adjoining three
identity coordinates, so the same invariant and dimension transfer without
loss:

\[
\boxed{\text{For every }N\ge4,\text{ exact symplectic maps of }\mathbb A^6
\text{ contain an }(N-3)\text{-dimensional weighted family of stable classes.}}
\]

See [degreewise exact symplectic moduli](extended-geometry/SYMPLECTIC_STABLE_MODULI.md).

The cancellation calculation gives a different, complementary result.  For a
fixed generic degree `N`, every proper divisor `r|(N-1)` gives

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

Unlike the moduli-dimension theorem, this finite count distinguishes multiple
construction types: the weighted locus from every cancellation type and the
cancellation types from one another.

The finer cancellation prime-intersection diagram depends on the contact
resultant `Res(K_{m,r},L_{m,r})`.  Besides the irreducibility ranges and exact
finite certificates, an endpoint-moment argument now proves uniform
nonvanishing for every `m` in all three columns `r=1,2,3`; the `r=3` proof
uses Schur--Cohn separation of a degree-six endpoint eliminant from the
negative-binomial root disk.  The unrestricted problem for `r>=4` remains
open.  See the
[contact-resultant reduction](cancellation/CONTACT_RESULTANT.md).

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
authorship remain credited to BCW, DVEZ, and Zhao.  A separate repository
optimization reuses exposed factors, reduces the route to dimension 16, and
then rank-compresses the cubic output before homogenizing to 24 variables.
An improved 17-dimensional trace has cubic-output rank six; its 24-variable
homogenization has a three-dimensional constant Jacobian kernel, giving a
21-dimensional cubic-homogeneous Keller collision.  The JC-derived,
nonconstructive-at-the-Gaussian-step route now gives `not GMC(42)`; it is a
certified upper bound, not a minimality claim or a result attributed to Long.  Long's
much smaller direct witnesses were motivated by the JC(3) announcement but
were not algebraically extracted from that map.  His papers do not review the
repository's Hessian/Fitting-divisor or stable-moduli proofs.

The same 21-variable cubic collision is also a direct counterexample
generator.  Its first collision coordinate is `0,1,-1`, which proves that the
first formal inverse coordinate has infinitely many nonzero homogeneous
pieces.  This names an explicit `SIC(21)` witness and, after contraction and a
complex orthogonal change of variables, an explicit 42-variable quartic that
simultaneously defeats the generalized Laplacian Vanishing Conjecture and the
classical Hessian-nilpotent quartic Vanishing Conjecture.  The formulas,
all-order arguments, expanded 876-term quartic, and the remaining
nonhomogeneous 20-variable optimization question are in the
[Image and Vanishing counterexample note](extended-geometry/IMAGE_VANISHING_COUNTEREXAMPLES.md).

There is now a separate internal bridge from the weighted inverse pencil to
Gaussian moments.  For every nonconstant normalized seed, a polynomial
two-complex-coordinate correction cancels the Lagrange--Good determinant and
produces an explicit `GMC(4)` witness family whose mixed moments read one exact
branch of `H(W)-sW+t=0`.  This construction was prompted by Long's method but
is not one of Long's witnesses and has not received external specialist
review.  See the [weighted Gaussian bridge](extended-geometry/WEIGHTED_GAUSSIAN_BRIDGE.md).
The complete mixed-moment sequence determines the branch series, whose
compositional inverse recovers `1+lambda*H` exactly; for a degree-`d`
polynomial, moments through order `d+1` suffice.  On normalized seeds,
`H'(1)=-1` then recovers both `lambda=-h'(1)` and `H`.  This is an injective
exact fingerprint, not a proof of inequivalence under arbitrary
Gaussian-variable transformations.
The load-bearing formal identity, including polynomial maps with nonzero
constant term, is isolated in the
[formal Gaussian--Lagrange lemma](extended-geometry/FORMAL_GAUSSIAN_LAGRANGE_LEMMA.md);
it is the highest-priority theorem here for external specialist review.

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

The same flux method now applies to the full normalized degree-five seed
surface.  In coordinates `(kappa,tau)`, a seed-dependent adapted change fixes
the common quotient coordinate `R=2X-3X^2Q`; cancellation of the complete
negative-`X` principal part uniquely forces an explicit rational quadratic
shear.  The resulting maps

\[
G_{\kappa,\tau}:\mathbb A^4\longrightarrow\mathbb A^4
\]

are exact symplectic and polynomially left--right equivalent to the weighted
maps times `id_(A^1)`.  On the ordinary boundary-clean open this descends the
full two-dimensional stable moduli from three canonical pairs to two.  The
former fixed-third-component family is the `kappa=-9` slice.  See the
[original slice theorem](extended-geometry/DEGREE_FIVE_RANK_TWO_DESCENT.md)
and the
[full surface theorem](extended-geometry/ALL_DEGREE_RANK_TWO_DESCENT_PROGRAM.md).

## Detailed novelty ledger

The snapshot near the top gives the headline results.  The machinery behind
them consists of the global regular-reconstruction open; exact image and
nonproperness theorems; reduced, scheme-theoretic, and formal boundary
invariants stable under left--right equivalence and stabilization; complete
normalization-boundary exhaustion; and the distinction between reduced and
nilpotently thick boundary intersections.  The contribution-by-source table
above and [STATUS.md](STATUS.md) delimit these repository results from
external discoveries and from locally reproduced work.

## Reading and reproduction

The foundational map, cubic marked-root model, exact image theorem, and
weighted theorem are the stable core. See [STATUS.md](STATUS.md) for detailed
evidence levels.

The stable-core reading path is:

1. [Foundational Keller map](verified/FOUNDATIONAL_GEOMETRY.md)
2. [Cubic marked-root model](verified/MARKED_ROOT_MODEL.md)
3. [Exact image and nonproperness](verified/IMAGE_AND_NONPROPERNESS.md)
4. [Weighted marked-root theorem](verified/WEIGHTED_SEED_THEOREM.md)
5. [Constant-kernel quotients and essential cubic input](verified/CONSTANT_KERNEL_QUOTIENT.md)

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
