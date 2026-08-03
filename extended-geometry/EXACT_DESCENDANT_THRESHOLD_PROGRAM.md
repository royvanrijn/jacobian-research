# Exact threshold program for descendant conjectures

## 1. Purpose and current ledger

This note isolates exact-minimum questions from witness-size records.  An
upper endpoint below means that the repository contains an exact witness;
it is not a literature-priority claim.  A lower endpoint is used only when
the corresponding full class is proved safe.

| threshold | rigorous status |
|---|---|
| smallest false Dixmier rank | at most \(3\); `DC_1` and `DC_2` remain open |
| Image/SIC pair dimension | exactly \(2\) |
| unrestricted constant-coefficient GVC dimension | exactly \(3\) |
| ordinary-Laplacian GVC dimension | \(2\le n_{\Delta{\rm GVC}}\le40\) |
| minimum two-pair bidegree-\((4,4)\) coefficient rank | \(2\le r_{\rm SIC}\le5\) |
| homogeneous quartic HN Hessian rank | \(3\le\rho_{\rm HN,4}\le37\) |
| homogeneous quartic HN dimension | \(6\le n_{\rm HN,4}\le38\) |
| genuinely ungraded Keller geometric degree | exactly \(3\) |

The last row is new.  The specialization `A=0`, `gamma=1` of the
[universal cubic testbed](../cancellation/UNIVERSAL_CUBIC_UNGRADED_TESTBED.md#8-an-explicit-genuinely-ungraded-cubic)
has two phantom boundary components.  Their unit lattice, together with the
ramified-normalization Fitting divisor, kills every connected decorated
automorphism.  The pointwise-fixed nonnormal-hypersurface lemma then excludes
every algebraic-torus-equivariant polynomial left--right representative.
The exact geometric-degree spectrum supplies the lower bound three.

The remaining sections record finite theorem-producing attacks.  They are
programs, not proved improvements.

## 2. Dixmier rank: glue primitive normalizer charts

The natural torus reduction already identifies the target primitive chart
with \(A_2\).  The unresolved source calculation should be separated into
two finite questions.

First, in the Ore localization at the primitive semi-invariant
\(\delta_B\), compute the full source normalizer rather than the image
subalgebra.  Filter by differential order and pole order and seek four
generators \(Q_1,P_1,Q_2,P_2\) satisfying the Weyl relations.  A PBW
certificate has a finite stop condition:

1. the associated-graded symbols generate the classical reduced coordinate
   ring;
2. all overlap ambiguities reduce to zero;
3. the resulting normal forms have the Hilbert series of \(A_2\).

Second, repeat on one complementary primitive chart.  The global problem is
then an intersection problem inside the common skew field.  If the two
localized \(A_2\) presentations have an inner or tame Weyl transition and
their pole lattices intersect in a four-generator PBW algebra, the reduced
endomorphism globalizes.  If the intersection has an extra generator or a
missing PBW monomial, the natural-torus route is excluded without another
high-support quantization search.

This chart-gluing test is preferable to asking first whether the classical
source fourfold is abstractly \(\mathbb A^4\): a quantum \(A_2\)
identification is the actual `DC_2` requirement, and the associated-graded
calculation detects precisely where polynomiality is lost.

## 3. Vanishing and Image dimensions: keep three notions separate

The absolute Image pair dimension is already exact:

\[
 \boxed{r_{\rm Image}=2.}
\]

It should not be merged with either GVC dimension.  The two-pair Image
witness has coefficient rank five and is not a Segre point
\(A(\zeta)P(z)\), so it does not give a two-variable constant-coefficient
GVC witness.

Unrestricted constant-coefficient GVC now has the exact failure dimension
three: the
[Hall-envelope theorem](BINARY_GVC_ENVELOPE_CLOSURE.md) proves the binary
case, and the
[homogeneous ternary witness](THREE_VARIABLE_HOMOGENEOUS_GVC_COUNTEREXAMPLE.md)
fails in dimension three.  The following Rees and apolar attacks are retained
as historical route ideas, not as an active dimension frontier.

Before closure, one proposed lower-dimensional attack was a Rees--holonomic
one.  In two variables every fixed leading filtration layer
of a hypothetical nonhomogeneous defect vanishes eventually.  Encode the
entire family

\[
 \Lambda^m(QP^m)
\]

as a diagonal/coefficient extraction of one rational Rees generating
function, with the operator order recorded by a second grading.  Creative
telescoping then produces a recurrence together with a Newton polygon in
the filtration depth.  The exact target is a finite-slope theorem:

> every nonzero holonomic branch has one of finitely many leading Rees
> slopes, and each such slope is already excluded by the split-symbol and
> low-order face theorems.

This would turn “the defect moves to unbounded depth” into a contradiction
for a fixed support, rather than continuing degree-by-degree correction
elimination.  It is not yet a uniform theorem over all supports.

Likewise, the former search for a smaller upper witness in dimensions three
or four proposed an **apolar quotient** rather than identifying a Dvorsky
variable with a linear form:
choose a one-dimensional source subspace and a dual symbol subspace whose
Koszul contraction is acyclic, then take the induced quotient on both
\(\Lambda\) and \(P\).  The required certificate is an all-power identity
in the quotient complex, not a finite moment match.  This strictly contains
the closed symmetry-preserving four-variable slice and gives a finite
rank-stratified elimination problem on the cubic catalecticants.

For ordinary Laplacian GVC, the one-block Schur completion is already
obstructed.  The next finite ansatz should use two auxiliary quadratic
blocks and impose both pure and mixed identities before solving for
coefficients.  Its objective is auxiliary quadratic rank, not raw term
count.

## 4. SIC coefficient rank: a projective good-reduction certificate

The rank-one stratum is excluded and the displayed witness has rank five.
The current rank-two route tries to extract an exact point from the
existential thirteen-moment fiber and then derive a recurrence.  There is a
shorter exclusion route which does not require such a point.

On the projective rank-at-most-two determinantal variety, form the corrected
moment ideal

\[
 I_2=(\mu_1,\ldots,\mu_{12},\mu_{14})
 + I_{\operatorname{rank}\le2}.                      \tag{4.1}
\]

Saturate by the projective nullcone ideal.  If the saturated projective
scheme is empty, every all-moment-zero rank-at-most-two point is unstable,
so no rank-two SIC counterexample exists.  The corrected degrees are
important: their Hilbert numerator passes the necessary nonnegativity test,
whereas degrees \(1,\ldots,13\) provably do not.

The certificate can be obtained at one good prime, provided the computation
uses the complete homogeneous projective model rather than a denominator
chart.  Properness supplies the characteristic-zero implication: a
nonempty generic closed subscheme has a closure whose proper image contains
every good fiber.  Thus an empty complete special fiber excludes a
characteristic-zero point.  The machine certificate must include:

1. the homogeneous determinantal and moment generators over
   \(\mathbb Z_{(p)}\);
2. a finite projective cover or a single homogeneous saturation;
3. exact unit certificates on every chart;
4. a nullcone containment certificate, not merely affine origin support.

The same test can be run successively for ranks three and four.  The first
rank where the corrected projective fiber survives becomes the only stratum
requiring point extraction and an all-order recurrence.

## 5. HN rank and dimension: use the projective kernel bundle

For a homogeneous quartic \(P\), the Hessian is a symmetric matrix of
quadrics

\[
 M:\mathcal O_{\mathbb P^{n-1}}(-2)^n
 \longrightarrow \mathcal O_{\mathbb P^{n-1}}^n.
\]

On a constant-rank stratum its kernel and image are vector bundles.
Nilpotency supplies a flag such as

\[
 \operatorname{im}M^2\subseteq\ker M,
\]

while symmetry identifies the orthogonal relations among the image,
kernel, and their twists.  This produces Chern-class equations before any
quartic coefficients are expanded.

For the first open dimension \(n=6\), stratify by generic Hessian rank and
linear common-kernel dimension.  On each moving-kernel stratum:

1. write the Chern polynomial forced by the nilpotent flag;
2. impose the symmetry duality;
3. compare with the degree-two Pluecker map of the kernel bundle;
4. eliminate impossible Chern data;
5. treat the surviving constant-kernel strata by quotienting the common
   kernel and invoking the known lower-dimensional theorem.

This can raise the HN dimension lower bound without a full coefficient
classification.  The same kernel-bundle degree stratification attacks
Hessian rank three: degree zero is a constant-kernel reduction, and the
first genuine case is Pluecker degree one.

For upper bounds, keep dimension and rank separate.  Cotangent lifts obey

\[
 \operatorname{rank}\operatorname{Hess}(y^TH)
 =2\operatorname{rank}JH+\operatorname{rank}(K^TAK).
\]

Thus the next rank-\(36\) target is exactly either cubic rank \(18\) with
zero kernel excess or cubic rank \(17\) with excess at most two.  New
circuit atoms should be solved from \(K^TAK=0\) as coefficient equations.
The next dimension-\(36\) target is different: it requires cubic source
score \(n+r\le17\), regardless of kernel excess.

## 6. Promotion rule

A threshold changes only after an exact characteristic-zero theorem or
counterexample and an independent finite replay where feasible.  Modular
empty fibers are rigorous only with the complete projective/proper model
specified above.  A bounded moment prefix, sampled matrix rank, affine chart
with an unchecked denominator, or failed search does not change an endpoint.
