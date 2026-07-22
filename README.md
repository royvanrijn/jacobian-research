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

## Discoveries so far: an ELI5 timeline

> **Snapshot: 22 July 2026.**  This is a timeline of public announcements,
> external contributions, and results first recorded in this repository—not a
> priority ruling or a claim of peer review.  The announcement appears as 19
> or 20 July depending on timezone.  See the [theorem index](THEOREMS.yml),
> its [rendered status summary](STATUS.md), and the
> [provenance audit](archive/legacy-notes/PROVENANCE_AUDIT.md) for the careful
> version.

- **19–20 July — the counterexample.**  After Akhil Mathew suggested the
  problem, [Levent Alpöge posted](https://x.com/__alpoge__/status/2079028340955197566)
  a formula credited to Fable.  It is a polynomial map of three-dimensional
  space that is locally reversible everywhere, yet sends three different
  points to the same point.  In one stroke this disproved the Jacobian
  conjecture in every dimension at least three; the two-dimensional case is
  still open.

- **20 July — people checked it and explained why it works.**  [Dean Cureton
  formalized the finite certificate in Lean](https://github.com/deancureton/jacobian),
  [Andy Jiang](https://x.com/davikrehalt/status/2079175065695035442) explained
  the map as “choose one root of a cubic, then forget which root was chosen,”
  and Qiaochu Yuan connected that picture to the MathOverflow discussion.
  [Alexis Gallagher](https://jacobianfun.org/jacobian-explained) found the
  weighted-lift viewpoint and a proposed example in every generic fiber
  degree.  A later geometric manuscript prompted by Semon Rezchikov supplied
  the auxiliary-conic construction used here.

- **20 July — the first consequence roundup.**  [Zihan Zhang's audit](https://zzhang-iu.github.io/papers/direct-consequences-jacobian/index.html)
  observed that known implications immediately force failures of the Mathieu
  conjecture for `SU(3)` and of the Gaussian Moments, Vanishing, and Image
  conjectures in some finite dimensions.  These were existence consequences,
  not yet small displayed examples.

- **20 July — `not GMC(158)`, followed by a much smaller direct witness.**
  [Christopher D. Long](https://arxiv.org/abs/2607.18186) tracked the announced
  map through a standard reduction to obtain the route-based bound
  `not GMC(158)`.  In the same paper he gave unrelated explicit polynomials in
  only three real Gaussians, proving the stronger statement `not GMC(n)` for
  every `n>=3`.  In everyday terms: every power of one polynomial can average
  to zero while multiplying by a second polynomial keeps producing a nonzero
  average.

- **20–21 July — the Dixmier conjecture fell too.**  Separate public write-ups
  by [Omniscience Research Agent and Jeff Pickhardt](https://omniscienceproject.com/papers/an-explicit-counterexample-to-the-dixmier-conjecture-in-a-3-jfLENtXF)
  and by [Fable 5 in a session directed and circulated by William
  Mayner](https://github.com/wmayner/dixmier-counterexample) lifted the same
  three-variable map to an injective but non-surjective endomorphism of the
  third Weyl algebra.  Thus the analogous noncommutative Dixmier conjecture is
  false in ranks `n>=3`.

- **21 July — more explicit examples appeared.**  [Long's second
  paper](https://arxiv.org/abs/2607.19012) gave tiny direct counterexamples to
  the `(xz)` conjecture and the Mathieu conjecture for `SU(2)`.  Also,
  [Juntang Zhuang](https://github.com/jzkay12/jacobian_conjecture) published
  three checked degree-four Jacobian counterexamples, now called Islands A,
  B, and C.

- **21–22 July — this repository turned one map into whole families.**  Using
  the public root and weighted-lift ideas above, it proved weighted and
  cancellation constructions in every generic degree `N>=3`.  For `N>=4`
  there are at least `1+(N-1)tau(N-1)-sigma(N-1)` stably different maps
  across the cancellation branches and one weighted class, as well as an
  `(N-3)`-dimensional family of stable classes—not merely a growing list of
  isolated formulas.  See the [degreewise theorem](papers/marked-root-multiplicity/main.tex)
  and [decorated-normalization paper](papers/decorated-discriminant-normalization/main.tex).

- **22 July — the consequences became smaller and more explicit here.**  The
  repository produced a seed-parametrized direct `GMC(4)` family, compressed
  Long's map-derived route from `not GMC(158)` to `not GMC(42)`, constructed a
  four-dimensional exact symplectic/Poisson counterexample, and wrote explicit
  Image and Vanishing counterexamples in dimensions 20, 40, and 42.  These are
  separate from Long's stronger direct `GMC(3)` result.  See the
  [external-consequences note](extended-geometry/EXTERNAL_CONSEQUENCES_AND_PROVENANCE.md),
  [Poisson audit](extended-geometry/QUADRATIC_LADDER_AND_POISSON_AUDIT.md), and
  [Image/Vanishing note](extended-geometry/IMAGE_VANISHING_COUNTEREXAMPLES.md).

- **22 July — the obvious five-dimensional sequel hit an obstruction.**  The
  repository classified the next quadratic-times-cubic root construction.  Its
  natural five-dimensional space is not ordinary affine five-space, even
  though topologically it looks like a three-sphere.  In plain terms, the
  elegant cubic trick does not automatically produce `JC(5)`; see the
  [`(2,3)` slice audit](extended-geometry/QUADRATIC_CUBIC_FACTORIZATION_SLICE.md).

- **22 July — dimension two was narrowed, not solved.**  Building on the
  Guccione–Guccione–Horruitiner–Valqui Newton-polygon program, [Billel
  Helali's exact certificates](https://doi.org/10.5281/zenodo.21479814) exclude
  the last transcribed case below 125.  The local
  reproduction therefore gives `max(deg P,deg Q)>=125` for any hypothetical
  plane counterexample, subject to the published normal-form reduction.  It
  does not prove `JC(2)` or say both degrees are at least 125; see the [plane
  frontier](plane-jc/DEGREE_FRONTIER_125.md).

The repository is both a proof project and a research diary that has now been
cleaned into a small active theorem chain.  Mathematical completion and
external review are deliberately kept separate in the machine-readable
[theorem index](THEOREMS.yml); [STATUS.md](STATUS.md) is its rendered summary.

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
The projective hyperplane-orbit and Grothendieck-class argument is isolated in
the [foundational incidence construction](verified/FOUNDATIONAL_INCIDENCE_CONSTRUCTION.md).

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
| Degreewise weighted stable moduli | No earlier source has been identified in the present audit | A degree-`N-2` decorated-normalization invariant with image dimension `N-3`; full marked-cover faithfulness remains a repair target. The separate Hurwitz/LL paper organizes the admissible-cover closure and compactified rerooting and compares the formal conductor square at simultaneous multicluster collisions |
| Cross-family stable multiplicity | No earlier source has been identified in the present audit | The complementary marked-open theorem distinguishes all cancellation parameter roots, all divisor types, and one weighted class, giving at least `1+(N-1)tau(N-1)-sigma(N-1)` stable classes in generic degree `N>=4` |

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

Finite-order deformation theory gives no additional unfiltered moduli.
The [formal orbit-triviality theorem](extended-geometry/FORMAL_ORBIT_TRIVIALITY.md)
shows that every polynomial deformation of a Keller map over a local Artin
base is uniquely trivialized by an unrestricted polynomial source
automorphism.  Thus dual-number directions in bounded-support ansatz slices
belong to the slice, not to global moduli.  The active continuation is
[complexity-filtered contact](extended-geometry/COMPLEXITY_FILTERED_CONTACT.md),
formal-to-algebraic descent, and global boundary invariants.

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
[decorated-normalization paper](papers/decorated-discriminant-normalization/main.tex).
The coarse decorated-normalization map remains generically etale of degree
`N-2`, reflecting the finite rerooting groupoid.  The proposed refinement by
the full marked incidence cover and its regular-reconstruction open remains a
repair target and is not used for this moduli theorem.
The [Hurwitz--LL compactification](papers/hurwitz-ll-rerooting/main.tex)
realizes the discriminant as a universal critical-value incidence and places
the compactified rerooting groupoid in a marked admissible-cover closure.  Its
normalized-Stein comparison identifies every contracted formal root chart and
the complete multicluster conductor square with the canonical finite
normalization.  On the collision-separating model, rerooting is the degree
`N-2` etale quotient stack obtained by selecting one simple zero-fiber point;
generic transposition ramification appears only after contraction to the
coarse coefficient model.  Projectivizing the remaining critical values gives
the matching marked Hurwitz stratum, and its restricted LL map has exact
degree `(N-2)N^(N-3)`.  On the root compactification, the caustic is
`2 kappa_1`; the caustic and Maxwell closures now both have explicit classes
in the `2N-7`-element invariant boundary basis.
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

For a fixed type `(m,r)`, the inverse polynomial is independent of the
chosen cancellation root, but its maximal reconstruction open is not.  The
open canonically recovers

\[
q=\left.\frac{B}{y^{m+1}}\right|_{A=0},
\]

so the `mr=N-1-r` roots give pairwise distinct stable classes.

A weighted marked-root map supplies one additional class because its canonical
boundary intersection is reduced.  Therefore:

\[
\boxed{\text{For every }N\ge4,\text{ there are at least }
1+(N-1)\tau(N-1)-\sigma(N-1)
\text{ pairwise stably inequivalent degree-}N\text{ maps.}}
\]

Unlike the moduli-dimension theorem, this finite count distinguishes the
selected weighted class from every cancellation branch, the cancellation
types from one another, and all parameter roots within each type.

The finer cancellation prime-intersection diagram depends on the contact
resultant `Res(K_{m,r},L_{m,r})`.  Besides the irreducibility ranges and exact
finite certificates, an endpoint-moment argument now proves uniform
nonvanishing for every `m` in all four columns `r=1,2,3,4`; the `r=3` proof
uses Schur--Cohn separation of a degree-six endpoint eliminant from the
negative-binomial root disk, while the `r=4` proof isolates the exceptional
conjugate pair and excludes its endpoint branch by rational Rouche, angle,
and Bernstein certificates.  The unrestricted problem for `r>=5` remains
open.  See the
[contact-resultant reduction](cancellation/CONTACT_RESULTANT.md).

The canonical statement and proof are in the standalone paper
[Marked-Root Keller Maps and Degreewise Stable Multiplicity](papers/marked-root-multiplicity/main.tex).
[Marked-root Keller maps](MARKED_ROOT_KELLER_MAPS.md) supplies the common
framework, while the [five-lemma audit](DEGREEWISE_MULTIPLICITY_AUDIT.md) is a
verification companion.  Neither is a competing canonical theorem source.

The positive-dimensional theorem is now in
[Decorated Discriminant Normalization and Stable Moduli](papers/decorated-discriminant-normalization/main.tex),
and the compactification and specialized LL/boundary calculations are in
[Hurwitz--LL Compactification of Rerooting](papers/hurwitz-ll-rerooting/main.tex).

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

The same cubic collision is also a direct counterexample generator.  Its
first collision coordinate is `0,1,-1`, which identifies the nonpolynomial
inverse coordinate and gives an exact all-order recurrence for its homogeneous
pieces.  Specializing the identity output descends the witness
to `not SIC(20)` and produces an independently expanded 628-term,
40-variable nonhomogeneous Hessian-nilpotent polynomial that defeats the
generalized Laplacian Vanishing Conjecture.  The homogeneous classical
quartic counterexample remains in 42 variables and has 876 terms.  The
formulas, proofs, recurrences, and exact artifacts are in the
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

The same flux method now applies uniformly to every normalized admissible
degree-`N` seed for `N>=5`.  On `kappa!=-1,-2`, the adapted chart fixes the
common quotient coordinate `R=2X-3X^2Q`; weighted homogeneity forces the
complete negative-`X` principal part onto one universal four-residue line,
which a unique quadratic `Q^2` shear cancels.  The replacement chart at
`kappa=-1` has its own universal residue line and is completed by an `XQ`
shear.  Thus every such seed has an exact symplectic completion

\[
G_H:\mathbb A^4\longrightarrow\mathbb A^4
\]

that is polynomially left--right equivalent to the weighted map times
`id_(A^1)`.  The coarse decorated-normalization map is generically etale of
exact degree `N-2` while retaining the full `(N-3)`-dimensional source image.
It therefore transfers the principal moduli result to pairwise stably
inequivalent exact symplectic maps of `A^4` without using full-cover
faithfulness.  A fixed `kappa=-9` slice has dimension `N-4`.  Degrees five
and six are retained as explicit regression specializations of the
all-degree residue theorem, not as independent theorem sources.  See the
[all-degree rank-two theorem](extended-geometry/RANK_TWO_SYMPLECTIC_DESCENT.md)
and the [degree-five worked example](extended-geometry/DEGREE_FIVE_RANK_TWO_DESCENT.md).

## Detailed novelty ledger

The snapshot near the top gives the headline results.  The machinery behind
them consists of the global regular-reconstruction open; exact image and
nonproperness theorems; reduced, scheme-theoretic, and formal boundary
invariants stable under left--right equivalence and stabilization; complete
normalization-boundary exhaustion; and the distinction between reduced and
nilpotently thick boundary intersections.  The contribution-by-source table
above and the [theorem index](THEOREMS.yml) delimit these repository results
from external discoveries and from locally reproduced work.  The compact
[status page](STATUS.md) renders that index.

## Reading and reproduction

The foundational map, cubic marked-root model, exact image theorem, and
weighted theorem are the stable core. See the [theorem index](THEOREMS.yml)
for evidence levels and [STATUS.md](STATUS.md) for the rendered overview.

The stable-core reading path is:

1. [Foundational Keller map](verified/FOUNDATIONAL_GEOMETRY.md)
2. [Tangent-map core](verified/TANGENT_MAP_CORE.md)
3. [Normalized factorization model](verified/NORMALIZED_FACTORIZATION_MODEL.md)
4. [Foundational incidence construction](verified/FOUNDATIONAL_INCIDENCE_CONSTRUCTION.md)
5. [Cubic marked-root model](verified/MARKED_ROOT_MODEL.md)
6. [Exact image and nonproperness](verified/IMAGE_AND_NONPROPERNESS.md)
7. [Weighted marked-root theorem](verified/WEIGHTED_SEED_THEOREM.md)
8. [Stable normalization functoriality](verified/STABLE_NORMALIZATION_FUNCTORIALITY.md)
9. [Constant-kernel quotients and essential cubic input](verified/CONSTANT_KERNEL_QUOTIENT.md)

Further families continue in one sequence:
[cancellation construction](cancellation/CONSTRUCTION.md),
[boundary geometry](cancellation/BOUNDARY_GEOMETRY.md), and the
[log-geometric suspension bridge](cancellation/LOG_GEOMETRY_OF_SUSPENSIONS.md), then the
[marked-root framework](MARKED_ROOT_KELLER_MAPS.md), followed by the
[finite degreewise paper](papers/marked-root-multiplicity/main.tex), the
[decorated-moduli paper](papers/decorated-discriminant-normalization/main.tex),
and the [Hurwitz--LL paper](papers/hurwitz-ll-rerooting/main.tex).
Explicit reciprocal candidates can be audited with the
[exact reciprocal-link classifier](cancellation/RECIPROCAL_LINK_CLASSIFIER.md).

For execution, start with [REPRODUCE.md](REPRODUCE.md).  The core paper is
[papers/core-counterexample/main.tex](papers/core-counterexample/main.tex).
Extended geometry is indexed in
[extended-geometry](extended-geometry/README.md); cancellation
arithmetic, rigidity, and the research roadmap live beside the construction
under `cancellation/`.  Superseded derivations and exploratory tools
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
