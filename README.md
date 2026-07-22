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

The complete proof and exact symbolic certificate are in the
[three-proposition factorization model](verified/NORMALIZED_FACTORIZATION_MODEL.md).

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
its finite-field count is `q^5-q^3`).  The
[`(2,3)` slice audit](extended-geometry/QUADRATIC_CUBIC_FACTORIZATION_SLICE.md)
records the calculation and the precise remaining cohomology question.

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
| GMC, `(xz)`, and `SU(2)` counterexamples | Christopher D. Long's two arXiv papers | Exact local Gaussian and beta/binomial checks, a complete `SU(2)=S^3` Haar proof, all 18 steps of Long's conservative 79-variable route, a repository shared-factor reduction to dimension 16, and rank-compressed homogenization to 24 variables and `not GMC(48)`; every sparse artifact has an independent replay, with no claim that the direct witnesses derive from our map |
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
nonvanishing for every `m` in both columns `r=1,2`; the unrestricted problem
for `r>=3` remains open.  See the
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
This improves the nonexplicit bound to `not GMC(48)`; it is a certified upper
bound, not a minimality claim or a result attributed to Long.  Long's
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
  determinant correction and exact mixed-moment formula, supported by the
  standalone constant-term Gaussian--Lagrange lemma;
- reduced versus nilpotently thick boundary intersections;
- the `(N-3)`-dimensional weighted stable-class family in every generic degree
  `N>=4`, transferred unchanged to exact symplectic maps of `A^6`; and
- the complementary divisor-count lower bound `tau(N-1)`, which distinguishes
  the weighted locus from multiple cancellation types.

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
