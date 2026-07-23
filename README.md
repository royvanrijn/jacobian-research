# From one collision to marked-root Keller maps

This repository verifies a three-dimensional polynomial Keller map with a
three-point collision, explains it through a tangent-map normal form, and
develops weighted, cancellation, decorated-normalization, Hurwitz, and
rank-two symplectic consequences.  In generic degree `N>=4`, the coarse
decorated normalization already gives an `(N-3)`-dimensional family of stable
classes.  Adding one affine root sheet generically recovers the seed exactly.
The same parameters yield four-real-Gaussian witnesses with an optimal
`(N-3)`-moment algebraic fingerprint; in those moment coordinates, the full
affine nonsurjective partition-lattice geometry and every other functorially
defined seed locus are preserved in their native categories.  In degree six,
the two vertical Ritt loci become explicit sextic moment hypersurfaces, one
coinciding with the all-double nonsurjective component.  The
[Moment--Prony theorem](extended-geometry/MOMENT_PRONY_DETERMINANTAL_GEOMETRY.md)
now derives every omitted partition and every vertical Ritt stratum from
intrinsic log-Hessenberg, saturated Christoffel--Hankel/subresultant, or
Krylov minors in those optimal moments; degree-five, degree-six, and
degree-eight ideals are checked scheme-theoretically.  The required
saturation removes the wrong marked-Fitting collision thickness while
retaining the genuine intersection nilpotents.  The same parameters
also yield inequivalent filtered Weyl endomorphisms.  The nonsurjective seed locus has
dimension `floor(N/2)-1`: quartics are generically exceptional, whereas
generic seeds of every degree `N>=5` are surjective.
The [exact degree-spectrum corollary](verified/GEOMETRIC_DEGREE_SPECTRUM.md)
shows that the geometric degrees of noninvertible Keller maps of
`A^3_C -> A^3_C` are precisely `3,4,5,...`.

For Gaussian moments, Long's explicit five-term polynomial settles GMC
negatively in three real variables and, by adjoining unused coordinates, in
every dimension `n>=3`.  Broad high-dimensional searches are archived.
The [lower-face prime theorem](extended-geometry/TWO_REAL_GMC_LOWER_FACE_THEOREM.md)
proves GMC in two real variables, completing the classification:
`GMC(n)` holds exactly for `n<=2`.  The former
[GMC(2) research program](extended-geometry/GMC2_RESEARCH_PROGRAM.md) now
records the proof route; high-dimensional implication chains and weighted
families remain as examples of logical transport.

## The foundational map in two forms

The same counterexample appears naturally in two polynomial coordinate
systems.  The first is the compact expanded form used in the original
announcement.  Put `u=1+xy` and define

\[
F(x,y,z)=\left(
u^3z+y^2u(4+3xy),
y+3xu^2z+3xy^2(4+3xy),
2x-3x^2y-x^3z
\right).
\]

On `x!=0`, set `t=y+1/x` and `r=2/x`.  For target coordinates `(a,b,c)`,

\[
a=t^2+rt/2-ct^3,\qquad b=r+4t-3ct^2.
\]

The two coordinate Jacobians give

\[
\det\frac{\partial(a,b,c)}{\partial(t,r,c)}
\det\frac{\partial(t,r,c)}{\partial(x,y,z)}
=\frac r2(-2x)=-2,
\]

so `det DF=-2` polynomially everywhere.  Exact substitution gives

\[
F(0,0,-1/4)=F(1,-3/2,13/2)=F(-1,3/2,13/2)=(-1/4,0,0).
\]

The second is the normalized factorization form.  For coordinates `(a,y,z)`,
set

\[
\begin{aligned}
b&=1+ay,\\
c&=1-\frac32ay+a^2z,\\
d&=\frac12y-az+\frac32ay^2-a^2yz,\\
e&=-2z+4y^2-4ayz+3ay^3-2a^2y^2z
\end{aligned}
\]

and define

\[
G(a,y,z)=(ac,\ ae+bd,\ be).
\]

Then `det DG=-1`, and

\[
\begin{aligned}
G(0,0,1/2)
&=G(1,-3/2,-13)\\
&=G(-1,3/2,-13)
=(0,0,-1/4).
\end{aligned}
\]

These are two forms of the same map, not two inequivalent counterexamples.
Indeed, with

\[
A(x,y,z)=(x,y,-z/2),\qquad B(p,q,r)=(r,2q,2p),
\]

one has \(F=B\circ G\circ A\).  The expanded form makes the low-degree
formula and collision immediate; the factorization form exposes the
coprime linear-times-quadratic construction and its resultant.

The independent replays are:

```bash
python3 scripts/verify_counterexample_independent.py
python3 scripts/verify_normalized_factorization_slice.py
```

## Canonical proof path

The common geometric normal form is

\[
(W,s)\longmapsto (s,Ws-H(W)).
\]

The [unifying marked-root thesis](UNIFYING_THESIS.md)
interprets the repository through one local--global mechanism.  A ramified
finite root incidence is suspended so that its Jacobian zero is cancelled by
a complementary vertical or source-chart factor.  The resulting polynomial
map is etale, but its affine source is only the regular-reconstruction open
of the finite normalization.  Noninvertibility is therefore carried by
boundary poles and escape to infinity, not by finite ramification.  The same
formal/local-versus-global-polynomial dichotomy reappears in the stable-orbit
and quantization problems.

The
[boundary-cancelled incidence lemma](cancellation/CONTROLLED_BOUNDARY_SUSPENSIONS.md#1-boundary-cancelled-incidence-lemma)
now supplies the common finite-root, chart-valuation, controlled-divisor,
polynomiality, determinant, boundary-normalization, and collision-fibre
mechanism. Weighted, cancellation, and quadratic-gauge maps are three
explicit instances. Their determinant proofs are centralized there; the
family notes retain only their different polynomiality and boundary data.
The [stable-intersection theorem](verified/QUADRATIC_CANCELLATION_STABLE_INTERSECTION.md)
shows that the cancellation and root-engineered quadratic-gauge families
share exactly one stable polynomial left--right class: the foundational
cubic.

The [tangent-map core](verified/TANGENT_MAP_CORE.md) supplies the weighted
inverse pencil, generic degree, critical divisor, discriminant normalization,
full Hessian Fitting divisor, and incidence form.

The [minimal-boundary gateway and classification program](cancellation/MINIMAL_BOUNDARY_CLASSIFICATION.md)
now defines eight predicates on a finite-normalization diagram: selected
critical boundary, saturated link, boundary monotonicity, ledger completeness,
puncture rank, primitive conormal, noncontraction, and chart
straightenability.  The conditional implications from this package to a
plane core and from a marked plane core to a suspension are theorems.  The
conjectural step is that every numerically boundary-minimal map canonically
supplies a diagram satisfying the package.

Inside the two already identified branches, the exact cubic reduction is
unconditional: the weighted degree-three branch and the `(m,r)=(1,1)`
cancellation branch are both polynomially left--right equivalent to the
foundational map.  The one-place cubic core marking is automatic; an
infinite toric `G_m` core atlas shows why the two-place marking still needs
threefold boundary input.  On the positive branch, primitive quotient and
Stein labels imply the weighted chart, with polynomiality forcing its
`-3/2` coefficient.  Thus the cubic calculation collapses the candidate
branches but does not extract them from the canonical finite-normalization
package.

The [cubic normalization frontend](cancellation/CUBIC_NORMALIZATION_FRONTEND.md)
now gives a second route which does not begin with a suspension.  The
canonical cubic normalization is automatically flat away from a
zero-dimensional Fitting defect.  That defect is now an exact finite test:
flatness at a point is equivalent to scheme-theoretic fiber length three,
and every defect fiber has length at least four.  An exact Koszul
trace-module model has `Fitt_3=(x,y,z)`, and the general nonflat
triple-cover correspondence shows that cubic algebra structure alone does
not exclude it; the Keller boundary must do so.  Every defect has an exact
determinantal presentation, and every reduced defect is automatically the
minimal Koszul case with square-zero length-four fiber `k plus k^3`, so the open
problem becomes exclusion of that boundary/affine-sheet collision.  The
foundational triple-root collision remains allowed because its flat fiber
is the curvilinear length-three algebra `k[epsilon]/(epsilon^3)`.
More generally, intrinsic curvilinearity of every collision fiber already
forces point-flatness: a fiber generator lifts by Nakayama and makes the
local cubic algebra monic and free.  Equivalently, the fiberwise relative
cotangent module is cyclic, or its first Fitting ideal is the unit ideal;
this is already part of the scheme-theoretic intrinsic package.
The generic primitive conormal class extends through all closed collisions
whenever its pure two-dimensional scheme-theoretic ramification support is
`S_2` and its rank-one cotangent module is `S_1`; a finite-length failure is
then impossible.  Equivalently, the remaining closed-point ledger consists
of the two computable modules `Ext_A^2(T,A)` and
`Ext_A^3(Omega_{B/A},A)`.  The canonical `S_2` hull
`C=Ext_A^1(Ext_A^1(T,A),A)` makes this ledger geometric: the two modules
are successively the canonical duals of the finite cokernels
`C/T` and `Omega_{B/A}/T tau`.  Point-flatness is therefore the exact
double-saturation statement that both cokernels vanish.  Local cohomology
couples them more tightly: after `C=T`, the second cokernel is exactly the
closed-point torsion `H_Z^0(Omega_{B/A})`.  For a presentation with image
`N` and `I=Fitt_3(B)`, this torsion is the explicit saturation quotient
`(N:I^infinity)/N`.  Certificate E is now one canonical-bidual test and
one module-saturation test.
The same algebraic mechanism controls the remaining degree-forty-two
Hessian synchronization support.  The shared
[support-saturation principle](verified/SUPPORT_SATURATION_PRINCIPLE.md)
identifies local-cohomology torsion, associated-prime avoidance, positive
grade, and presentation saturation; flatness is only a stronger shortcut.
The active
[cubic closure protocol](cancellation/CUBIC_CLOSURE_ATTACKS.md)
organizes the remaining work into three certificates: those two Ext
modules, the phantom-boundary quotient between the nonproperness and branch
equations, and the graded discriminant obstruction for coefficient gauges.
When the defect vanishes,
Deligne--Faddeev and Quillen--Suslin extract a global binary cubic.  If its
coefficient morphism is affine-linear and the source is the full simple-root
open, the three hyperplane orbits force the tangent-nonosculating slice and
the foundational map.  The remaining cubic problem is thereby reduced to
point-flatness, coefficient linearity, and exclusion of extra simple
boundary.  The accompanying
[gauge-straightening theorem](cancellation/CUBIC_GAUGE_STRAIGHTENING.md)
removes an infinite apparent obstruction: every slice
`C_1=q-3C_0h`, `q!=0`, with `h` invariant under the translation locally
nilpotent derivation is polynomially left--right equivalent to the linear
tangent slice; the opposite shear gives the symmetric `C_2` family.  The
Jacobian formula `1+D(h)` proves these invariant times exhaust all
single-shear polynomial automorphisms.  Unit determinant also makes every
polynomial upper or lower Borel gauge a constant diagonal followed by one
of these classified shears.  For a lower--upper alternating pair, an exact
rank-two Jacobian formula and conjugated locally nilpotent derivation
classify the pair whenever its first factor is itself an allowed
variable-time shear.  Exact coefficient comparison also rules out every
origin-normalized alternating cancellation with two homogeneous-linear
times, and an all-degree support theorem excludes every pair of nonzero
monomial times.  For general multi-monomial times, the graded equation has
no cokernel except for one discriminant line in every fourth degree.
The first apparent bilinear obstruction on that line vanishes identically
by an `sl_2` integration-by-parts identity, so the viable route is to
extract the intrinsic root line and reduce to the proved Borel theorem;
formal first-obstruction chasing is not sufficient.
The saturated-flag reduction makes this precise: a projective
Tschirnhausen flag with scalar mixed coefficient puts the family in
`C_1=q`.  The full simple-root cover is then the cartesian pullback of the
foundational tangent-hyperplane map along
`G=(C_0,C_2,C_3)`.  The remaining coefficient question is the narrower
base-change rigidity assertion that the minimal-boundary pullback with
source `A^3` forces `G` to be a polynomial automorphism.
Coefficient linearity must therefore be tested modulo polynomial
Tschirnhausen gauge, not in a displayed basis.
Separately, the cubic DVR
degree sum shows that the critical divisor has exactly its ramified double
sheet and one affine simple sheet, so any extra simple boundary would have
to create a distinct second nonproperness divisor.  Its exact detector is
the phantom factor `u_F=j_F/delta_F`; proving the nonproperness hypersurface
irreducible makes this factor a unit.

Its [all-degree rational-fiber corollary](verified/ALL_DEGREE_RATIONAL_FIBERS.md)
gives, for every `N>=3`, a Keller map with a complete fiber of exactly `N`
distinct rational points and an open real target neighborhood with `N` real
sheets.  The construction uses explicit integer roots and is uniform in `N`;
the former computation through degree 100 is retained only as a regression.
The [real-sheet spectrum theorem](verified/REAL_FIBER_SPECTRUM.md) sharpens
this to the exact chamber spectrum `N,N-2,...,N mod 2`, proves the minimum is
zero in even degree and one in odd degree, supplies rational targets for every
count, and exhibits the full parity chain by successive fold crossings.

The [standalone universal-monodromy theorem](verified/UNIVERSAL_SYMMETRIC_MONODROMY.md)
proves geometric and arithmetic `S_N` monodromy for every polynomial pencil
`H(W)-sW+t`, including monomial, Chebyshev, decomposable, symmetric, and
critical-value-collision cases; no generic-seed hypothesis is present.  The
[universal weighted-seed theorem](verified/WEIGHTED_SEED_THEOREM.md) applies
it to every admissible weighted map.
Its [effective finite-field Chebotarev corollary](extended-geometry/FINITE_FIELD_CHEBOTAREV.md)
realizes every prescribed cycle/factorization type over all sufficiently
large certified good fields.  It gives the split density `1/N!`, irreducible
density `1/N`, the random-permutation fixed-point law and its effective
`Poisson(1)` limit, plus a deterministic generator for modular Keller-fiber
witnesses and small rational lifts carrying the selected modular fingerprint.
Combining this with the real chambers by constructive weak approximation
gives the
[adelic complete-fiber theorem](verified/ADELIC_FIBER_ENGINEERING.md): every
degree and signature occurs as a complete fiber field of the explicit `F_N`,
and finitely many additional unramified splitting conditions may be imposed
simultaneously at sufficiently large good primes.
The [Hasse-principle fiber theorem](verified/HASSE_PRINCIPLE_KELLER_FIBER.md)
goes in the complementary arithmetic direction: one explicit degree-eight
complete regular fiber has points over `R` and every `Q_p` but no rational
point.  Its integral target is `(12138,-308652,1)`, and the full fiber is the
finite etale scheme cut out by an elementary intersective polynomial.
The [global Sunada construction](extended-geometry/GLOBAL_SUNADA_KELLER_COVERS.md)
uses the point and line actions of `GL_3(F_2)` to give two nonisomorphic
degree-seven inverse covers over one two-dimensional target, with identical
zeta functions at every good fiber.  One Cox-ledger coordinate makes them
determinant-one threefold morphisms to the same target; after a finite
tangent-mark base change they are also pullbacks of stably inequivalent
relative weighted Keller maps.  The sources are boundary complements or
affine bundles, so an absolute polynomial three-space realization remains a
controlled-boundary ledger problem.
The follow-up [Cox-boundary audit](extended-geometry/DAVENPORT_COX_BOUNDARY_OBSTRUCTION.md)
computes the exact pullback `E_3 E_6 (g')^2`, proves stable affine-space
straightening and every coordinate-preserving polynomial suspension
impossible, and finds a one-row unimodular ledger completion whose obvious
`E_3` chart nevertheless creates a new coprime divisor.  The separate
[tangent-mark audit](extended-geometry/DAVENPORT_TANGENT_MARK_CURVE.md)
identifies the marking extension as a punctured rational conic, rules out
affine-line descent of the fixed mark, and gives a unit Gröbner certificate
against every affine-linear moving tangent pair.
Allowing the marks to be linear in a square-root parameter does produce a
new section: the
[proportional tangent audit](extended-geometry/DAVENPORT_PROPORTIONAL_TANGENT_SECTIONS.md)
classifies four exact quadratic sections and gives simultaneous conjugate
marks for both covers over `T=-s^2`.  After dividing a common `s^6`, the
relative Keller seeds extend across `s=0`; three explicit determinant and
endpoint divisors still prevent an absolute affine-space realization.  The
[weighted-glue audit](extended-geometry/DAVENPORT_WEIGHTED_GLUE_OBSTRUCTION.md)
shows that the slope pole alone has an affine-space modification, but the
simultaneously required intercept center is coprime to it.  Their combined
modification makes the boundary coordinate `C` invertible and therefore
deletes, rather than fills, the missing boundary.
The resulting determinant threefold has a reverse affine modification with
Jacobian factor `D^2`, matching the doubled derivative column of the Cox
ledger.  The
[derivative-center audit](extended-geometry/DAVENPORT_DERIVATIVE_CENTER_MISMATCH.md)
shows why this does not yet close the construction: the normalized
Davenport derivative curve is a three-punctured rational curve, whereas the
modification center is `G_m`; their geometric unit ranks remain different
after every affine stabilization.
On the normalized derivative boundary, the unique \(K\)-defined involution
swapping the quadratic puncture pair does have quotient `G_m`.  The
[boundary-involution audit](extended-geometry/DAVENPORT_BOUNDARY_INVOLUTION.md)
locates the remaining failure exactly: that involution does not preserve
the two branches above the node and therefore does not descend to the
singular derivative curve.
The
[node-separation audit](extended-geometry/DAVENPORT_NODE_SEPARATING_AFFINE_MODIFICATION.md)
constructs an explicit affine-plane modification which separates those
branches and satisfies
`pi^*J=(2Y+1)^2 S` with modification Jacobian `(2Y+1)/4`.
It also classifies both one-shear polynomial-coordinate orientations in
every degree by a three-puncture norm calculation.  The first orientation
always deletes a torus point.  The opposite orientation has one unique
torus-clean solution, the affine-plane graph of `T+Y^2`, which separates the
node with the exact doubled pullback but deletes the two branches over the
coordinate origin.  The secant and parabola charts complement one another
on the full normalized boundary, but their explicit reciprocal transition
completes the exceptional fiber to `P^1`, so their natural ambient gluing is
non-affine.  A Newton-polygon argument also closes the remaining alternating
search:
any polynomial coordinate line preserving the full derivative boundary
would force `J(T(t),Y(t))` to be quadratic, whereas successive leading
cancellations uniquely reduce it to a nonzero cubic in `Y(t)`.  Hence no
one-chart polynomial-coordinate separation exists at all; the surviving
routes must modify the explicit exceptional `P^1` or use a
higher-dimensional non-coordinate surgery.
The
[post-coordinate attack program](extended-geometry/DAVENPORT_POST_COORDINATE_ATTACKS.md)
ranks those routes.  Its leading proposal promotes the fixed determinant
constant to the target coordinate
`beta=D^2 U-C R`; the total map is between affine four-spaces, has Jacobian
`D^2`, and replaces the fixed-`beta` `G_m` center by `A^2`, removing the
AS6 unit-rank mismatch at the family level.  The direct product with the
unique parabola chart has Jacobian `D^3/2`, so it is not Keller and cannot
be repaired by polynomial coordinate changes or identity stabilization.
Exposing `D=T+Y^2` in a stabilized target is also impossible: it would make
the degree-seven Davenport generator quadratic over the target.  A shifted
mask `d=s+(T+Y^2)` preserves degree seven but does not cancel the derivative,
and every fixed-`T` constant-Jacobian coupling affine-linear in any number
of mask variables is a polynomial automorphism.  The leading live problem
cannot use the two obvious nonlinear variants either: making
`s+T+Y^2` a target coordinate again forces an explicit triangular inverse,
while correcting `g_T(Y)` by a multiple of `T+Y^2` fails the Keller
equation on that parabola.  Moreover, `beta=D^2U-CR` itself has zero
differential on `D=C=R=0`, so it cannot literally be an output of an
absolute Keller map.  The surviving leading attack is the translated
incidence `beta_hat=S+D^2U-CR`, with `S` replacing an old target coordinate
and `T` mixed directly with `C,R,U`.  The translation is itself a
Jacobian-one source coordinate, so retaining `T,D,C,R` reduces the proposal
exactly to a plane Keller family; parameter-linear dependence on `U` also
splits off and becomes a plane pair.  A genuinely new candidate must
replace `T` and use nonlinear `U`-dependence or essential sheet-dependent
coefficients.  In fact every coefficient curve
`(T+a(Y)U, g_T(Y)+b(Y)U, h(T,Y)+U)` is now excluded in all degrees: its
coefficient equations force a Bézout identity and then factor the remaining
Jacobian into a nonunit depending on `T`.  The complementary class with
`a,b` depending only on `T` leaves the derivative `J` as a factor and also
fails.  Even dependence through the distinguished mixed coordinate
`T+Y^2` factors into a nonunit with nonzero `Y^6` term.  Only a different
genuinely mixed `T,Y` pencil or nonlinear `U`-dependence remains.
For an arbitrary polynomial coefficient pencil `(x,y)`, the translated
Jacobian now factors through the necessary unit
`a(x) G_y-b(x) T_y`.  Every affine-linear pencil fails this gate, so future
searches can reject them before solving any coefficient equations.
Both one-triangular Jung families fail as well: degree dominance handles
every case except `x=T+cY^2+dY+e`, whose nine high-degree coefficient
equations have unit Gröbner ideal.  An affine-in-`U` candidate must now
start at alternating triangular polydegree two.
The [stratified extension](extended-geometry/STRATIFIED_ADELIC_FIBER_ENGINEERING.md)
moves the seed and target together on any rational chart of a selected locally
closed seed stratum. It isolates geometric nonemptiness as the only extra
compatibility condition and gives an exact nonsurjective quintic of omitted
type `2^1 3^1` and trivial Hessian symmetry whose three allowed real
signatures all occur with Frobenius types `(5)` at `7` and `(2,2,1)` at `11`.

The primary dependency chain is:

1. `F1` — foundational Keller collision.
2. `W1` — tangent-map core and weighted suspension.
3. `S1` — stable normalization functoriality.
4. `WB1` — weighted clean-locus intrinsic boundary reconstruction.
5. `C1` — universal cancellation construction.
6. `B1` — complete canonical boundary exhaustion.
7. `P1` — cancellation reconstruction residue and parameter faithfulness.
8. `M1` — finite degreewise stable multiplicity.
9. `D1` — marked Hessian-divisor moduli of dimension `N-3`.
10. `F2` — generic affine-mark faithfulness.
11. `H1-COARSE` — finite-normalization isomorphism with the wonderful
    target graph.
12. `H1-STACK` — explicit recursive logarithmic stack atlas.
13. `R1`, `R2` — all-degree rank-two descent and parameter faithfulness.

The exact scopes, dependencies, checkers, proof types, and assurance states
live only in [`MATH_STATUS.json`](MATH_STATUS.json). [STATUS.md](STATUS.md) is
generated from it. Independent replay, formal verification, and external
review are separate booleans. `artifact_hash` pins the registered checker
content, while `software_lock` names the repository lock manifests used for
its replay; neither field upgrades a claim's review status. False historical
claims are first-class `falsified` entries rather than prose-only corrections.

## Main papers

- [Foundational counterexample](papers/core-counterexample/main.pdf)
- [Discriminant pencils](papers/discriminant-pencils/main.pdf)
- [Decorated discriminant normalization](papers/decorated-discriminant-normalization/main.pdf)
- [Marked-root degreewise multiplicity](papers/marked-root-multiplicity/main.pdf)
- [Hurwitz--LL rerooting](papers/hurwitz-ll-rerooting/main.pdf)

The degree-five and degree-six calculations are worked regressions generated
from the all-degree results, not independent theorem authorities.  Quartic
models and external-map classifications are likewise examples or audits.

## Reproduction

Create the pinned Python environment in [REPRODUCE.md](REPRODUCE.md), then run:

```bash
make check verify-minimal verify-master verify-theorems verify-papers
```

`make check` validates links and the status graph.  `verify-minimal` needs only
the Python standard library.  `verify-master` covers the cancellation chain,
including the unconditional all-`(m,r)` thick-contact formula.  The conditional
contact-resultant refinements are separate from the `M1` proof path.
`verify-papers` discovers and compiles every `papers/*/main.tex`; CI uses the
same discovery rule and archives every resulting PDF.

The full command catalogue, heavier symbolic runs, optional Lean certificate,
and independent Macaulay2 comparison are documented in
[REPRODUCE.md](REPRODUCE.md).

## Research status

The H1--H3 chain has been adversarially audited.  The finite-cover valuative
lemma, quotient-stack map, and local collision algebra are proved, including
simultaneous clusters and relative stabilizers.  The audit exposed that the
original comparison with the unmodified root-stable quotient is false.  The
corrected graph is now constructed at the coarse level: the normalized
wonderful target pullback and finite selected polynomial closure prove
`H1-COARSE` and corrected H2, while finite normalization proves coarse H3.
The `H1-STACK` theorem identifies the selected fs logarithmic factor by the
saturated pushout \((P\oplus_{\mathbf N^E}Q)^{\mathrm{sat}}\); recursive
source screens are charts of that stack.  The standalone
[DVR marking audit](papers/hurwitz-ll-rerooting/dvr-marking-audit.tex)
retains the historical conditional analysis and its post-repair resolution.
The subsequent [branch-scale fan experiment](extended-geometry/BRANCH_SCALE_FAN.md)
proves that the degree-five weighted blowup is the first radial chart of a
general weighted braid subdivision, but also finds a degree-six triple
resonance where the full stable-target graph needs a further nested
refinement.  In the labelled degree-six chart that target refinement is now
constructed explicitly as the normalized pullback of the four-point blowup
of `P^2` giving `Mbar_0,5`; the admissible source expansion and saturated-node
comparison are also exact on the generic equal-scale face and its first
triple-Maxwell refinement.  All thirteen radial scale types now have a
verified degree-six admissible component atlas and saturated node monoids.
The three pairwise Maxwell source-node rings, the triple ring, and all their
radial-boundary intersections are also verified; no further local blowup
center occurs in this chart.  Standard admissible-cover deformation theory
supplies ambient toroidal gluing.  An exact central Hurwitz calculation finds
two ambient cover classes but proves that the labelled source-root
cross-ratio selects the polynomial class; consequently the explicit branch
graph is the labelled coarse admissible-cover graph in this chart.  Stack
inertia is now separated from normalization multiplicity: labelled radial
inertia is trivial, while every pairwise or triple Maxwell target node
retains one diagonal \(\mu_2\).  The canonical algebraic candidate is the
iterated second-root stack along the three pairwise and one triple-Maxwell
boundary divisors.  Its face maps, codimension-two inertia, \(S_3\)-action,
and separate local \((S_2)^3\rtimes S_3\) pair-root quotient are explicit;
smooth tame-stack reconstruction identifies it with the labelled ACV graph
on this chart.  The
[general labelled-node theorem](extended-geometry/LABELLED_NODE_SATURATION.md)
makes the full root-label action canonical on the global graph and
identifies the corrected H2 marked/unmarked quotient formally.  The
[wonderful-pullback theorem](extended-geometry/BRANCH_GRAPH_WONDERFUL_PULLBACK.md)
now constructs the stable-target graph for arbitrary labelled critical-value
collisions.  The
[source-vertex rigidity theorem](extended-geometry/SOURCE_VERTEX_RIGIDITY.md)
then proves that two complete target-flag fibers determine each rational
source-component map up to scale and a third fiber fixes the scale.  Thus
there is no residual vertexwise Hurwitz-class choice.  The
[general radial source theorem](extended-geometry/GENERAL_RADIAL_SOURCE_ATLAS.md)
constructs the source tree, vertex maps, and saturated node partitions on
every ordered first-scale stratum for arbitrary cluster multiplicities.  The
[polynomial monodromy-forest
theorem](extended-geometry/POLYNOMIAL_MONODROMY_FORESTS.md) also determines
the source dual graphs, component degrees, and node indices on every nested
simple-branch resonance stratum, unifying Maxwell and caustic.  The
[recursive resonance atlas
theorem](extended-geometry/RECURSIVE_RESONANCE_ATLAS.md) gives formulas that place
weighted source flag by normalized initial-form equations in framed
root-residue screens, proves contraction compatibility by affine
composition, and automatically matches full vertex centralizers across all
nested nodes.  Its log-étale comparison identifies these screens with the
fs-saturated base change of standard admissible-cover charts.  This supplies
all-degree coverage, algebraization, overlap descent, and the actual
stabilizers independently of the bounded regression checks.  The
[finite-normalization theorem](extended-geometry/SOURCE_GRAPH_FINITE_NORMALIZATION.md)
shows that this does not obstruct the coarse comparison: the wonderful graph
is already the complete coarse polynomial admissible-cover graph.  It also
makes the corrected H2 quotient and coarse H3 specialization unconditional.
The [monodromy-centralizer
theorem](extended-geometry/MONODROMY_INERTIA_CHARACTERS.md) further proves
the full-chain radial formula
\(\prod_j\operatorname{lcm}(\mu_i:i\in B_j\cup\cdots)/
\operatorname{lcm}(\mu_i:i\in B_j)\): equal multiplicities give trivial
inertia, while unequal multiplicities can retain a finite subgroup.  It also
gives the generic anchored/unanchored inertia formula for every simple
resonance node.  Together these results prove `H1-COARSE`, `H1-STACK`, H2,
and coarse H3.  The original comparison with the unmodified root-stable
quotient remains false.

The active continuation queue consists of three LR problems together with
the cancellation, Hessian--Ritt, restricted-minima, and minimal-boundary
frontiers:

- `OP-CR`: cancellation contact resultants in the residual staircase
  `m>=1001`, `7<=r<=X_m/m`; the first six fixed-`r` columns, the 1000
  all-`r` columns `m<=1000`,
  the exact effective Dusart region, the asymptotic region
  `r>=5(mr)^(21/40)`, and all other parameter-irreducibility ranges are
  complete.  The `r=6` column is closed by 7424 rigorous Arb
  branch tubes on `0<=1/m<=1/41` and exact modular gcds for `m<=40`.  An
  exact `r=7` pilot shows that the branch schema survives for the next fixed
  column and proves its eventual tail, but does not make that threshold
  effective or provide a uniform continuation in `r`.
- `OP-LR-REES`: linear target-lift Rees strictness as a finite module/SAGBI
  problem, including a structural cutoff to finitely many torus weights.
- `OP-LR-II`: the minimal pair is now nonzero; classify the remaining
  generator matrices of `II_(F,p,-p)`, prove the finite-weight cutoff, and
  determine whether the minimal class generates the quadratic normal image.
- `OP-LR-NE`: valuative no escape through the finite proper closure of the
  marked normalized-cover `Isom` relation, using `F2`, `H2`, `H3`, and generic
  deck rigidity before reconstructing automorphisms.
- `OP-RITT`: `HR12` now proves the exact all-degree codimension formula and
  completes the degree-eight and degree-twelve diagrams, including the
  all-four Chebyshev intersection in degree twelve.  `HRSYNC` closes the
  degree-thirty global all-six intersection and all fifteen pair schemes,
  including the transported pairs `{2,3}`, `{2,5}`, and `{3,5}`.  `HRCF`
  now proves the all-degree common-right-factor uniqueness step over every
  characteristic-zero base ring.  In degree forty-two, `HR42J` supplies the
  fourth-normal-order defect membership and `HR42K` promotes it to exact
  synchronization on the conormal, regular quadratic/cubic, discriminant,
  generic resultant, and punctured higher-gcd strata.  `HR42R` factors the
  remaining scalar quartic divisors and reduces them to the two explicit
  residual graphs `P_A=0`, `P_B=0`.  Those graphs, their contact-five
  vertex, the odd `w2=0` core, uniform terminal-refinement existence through
  nonreduced thickenings,
  component-completeness, coherence, and boundary gluing remain.
- `OP-RMIN`: the exact current intervals are
  `3<=r_cub<=17`, `3<=nu_cub<=18`,
  `3<=rho_HN,4<=37`, and `6<=n_HN,4<=42`.  Polynomial-gate BCW circuits
  has lowered the rank upper bound from 18 to 17 and the index upper bound
  from 19 to 18, while a separate 44-variable HN lift lowers the exact
  quartic Hessian-rank upper bound from 38 to 37.  Continue circuit-level
  rank/power-rank search and resolve the invertibility-only full-class
  index-three case separately from the known power-linear and
  symmetric-Jacobian theorems.  The
  stronger inverse-degree-nine claim is already false by van den Essen's
  exact dimension-five rank-three example of inverse degree thirteen; see the
  [restricted-minima frontier](extended-geometry/RESTRICTED_MINIMA_FRONTIER.md).
- `OP-SUSP`: minimal-boundary gateway and classification.  Starting from a
  numerically boundary-minimal canonical normalization, construct the
  intrinsic finite-normalization diagram and prove its eight `MBPkg`
  predicates.  In degree three this means extracting the
  coordinate-preserving suspension, primitive conormal and noncontraction
  marks, and positive chart-straightening labels.  The formal package-to-core
  and marked-core-to-suspension implications are already separated in
  `MBP1`; once the witness exists, the two cubic branches collapse to the
  foundational map.
  Alternatively, the cubic-normalization frontend reduces uniqueness to a
  cubic fiber-minimality (equivalently no zero-dimensional flatness defect),
  affine-linearity of the intrinsic binary-cubic coefficient map, and
  absence of extra simple boundary.

In particular, `b_m=34m+1` is an exact source-only profile in one
determinant-normalized target gauge.  It is not target-minimal and is not a
stable LR obstruction.  The better target torus produces a surviving
`v^(6m)S^(4m)` class in the invariant-ring-saturated equivariant target
quotient.  Modulo `gamma` its `v`-degree is `10m`, which already points to an
ordinary source-degree bound around `20m`; the exact torus-gauge slope `24` is
not needed.  The open algebraic step is to identify this as a universal
associated-graded obstruction after arbitrary lower gauges.  Formal
target-lifting now shows that all lower source jets are uniquely forced by the
target jets.  At order two, logarithmic coordinates reduce the new gauge
dependence to the affine-linear bracket class `Xi_2(Y_1)`, followed by one
explicit quadratic coordinate reconstruction.  The entire five-parameter
degree-25 homogeneous weight-zero kernel has now been audited: its order-two
quadratic residue never vanishes.  A filtered-coset BCH lemma now shows that
pairwise opposite-weight expansion is unnecessary at order two if the Rees
degeneration of the target-lift coset is strict.  The next implementation is
to construct the graded target-field, lifted, and normal modules and test
their generators by Gröbner/module membership, not to enlarge the symbolic
jet calculation.  Its quadratic obstruction is the explicit second
fundamental form
`II_F(Y_1,Y_2)=-(DF)^(-1)D^2F[ell_F(Y_1),ell_F(Y_2)]`; opposite weights reduce
to invariant-module maps `II_(F,p,-p)` into the weight-zero normal module.  The
smallest possible pair is already nonzero: the constant fields
`(partial_B,partial_C)` have weights `(1,-1)` and saturated normal symbol
`-146880u^5/7`.  Quadratic strict descent therefore fails, and this class is
the canonical quadratic LR invariant anticipated by `OP-LR-II`.  Valuative no
escape is a separate marked-cover problem, not an automorphism-coefficient
estimate.  None is a bottleneck for stable moduli, which already follow from
decorated normalization and the affine sheet.

Other questions—arithmetic Galois theory, wider quantization,
coefficient-scheme gluing, quadratic--cubic flexibility, and the plane degree
frontier—are retained as parked side programmes in
[STATUS.md](STATUS.md).

## Provenance

The earliest public item located by the repository audit is
[Levent Alpöge's post](https://x.com/__alpoge__/status/2079028340955197566).
Contemporaneous sources attribute the example to Alpöge and Fable; Dean
Cureton supplied a separate
[Lean certificate](https://github.com/deancureton/jacobian), and Alexis
Gallagher developed the weighted-lift viewpoint in a
[same-day explainer](https://jacobianfun.org/jacobian-explained).

The surviving record does not settle the complete discovery history, original
prompt, model conversation, or exact UTC timestamp.  This repository repeats
the contemporaneous attribution without making a priority claim.  See the
[provenance audit](archive/legacy-notes/PROVENANCE_AUDIT.md) for the detailed
record and the canonical sources for precise mathematical claims.

## Repository policy

Proofs belong in canonical sources, statuses belong in `MATH_STATUS.json`, and
open-problem IDs belong in `STATUS.md`.  Papers and notes may cite those IDs
but do not maintain independent continuation queues.  Superseded derivations
remain available under `archive/` and outside primary navigation.
