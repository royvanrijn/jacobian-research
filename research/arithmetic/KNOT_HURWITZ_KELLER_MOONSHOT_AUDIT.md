# Seven knot, alternating-cover and Keller research targets

Infinite prime augmentation-gap families now exist at every sufficiently
large bridge number. The real paired-cover locus has
a uniform matching-walk classification and complete real-root spectra.
One fixed Keller map carries the degree-ten A10 pencil in dimension three;
arbitrary rational finite-etale families, including the rank-301 Hasse
family, have fixed-map realization in dimension four. The degree-twelve
alternating Hurwitz passport has one absolute component. Known graded
Keller constructions give generic A_n in dimension three for every even
n>=6, with an explicit fixed-field formula below. Rational paired-critical
existence in degree 26, higher complex Hurwitz connectivity and the plane
Keller question remain open.

This document connects the current proofs. Mathematical status and the proof
dependency DAG remain in [MATH_STATUS.json](../MATH_STATUS.json); an exact
computation is not promoted to formal verification, external review, or a
historical priority claim.

| Target | Smallest new case tested | Result and next unresolved case |
| --- | --- | --- |
| Prime augmentation gaps | Unbounded bridge number | Proved for every sufficiently large bridge number, with explicit families and non-effective thresholds; direct prime-family certificates also hold at bridge five and six. |
| Rational square pencils | Degree 26 | Reciprocal models form finite schemes with exact lengths and controlled denominators. Both unbalanced sign splits and every rational-square reciprocal parameter are excluded. General rational existence remains open. |
| Hurwitz classification | Degree twelve | One absolute full-A12 orbit. Chosen-leaf deletions are symmetric or a full two-block wreath group; the entire imprimitive case has one normal form and recursive extraction. The symmetric case, degree fourteen and arithmetic refinements remain open. |
| One fixed Keller map | Degree-ten pencil | Solved uniformly for a reparametrized pencil; infinitely many complete A10 field fibres of one fixed A3 map. |
| Generic alternating Keller group | All even natural actions n>=6; the smallest plane case A4 | Polynomial A3 realizations in all even degrees; regular arithmetic A_(m^2+1) over Q(sqrt(-3)). Every plane A4 action is excluded. General plane realization or exclusion remains open. |
| Hasse algebras | A10 subgroup covers including ramified groups | Infinite rank-301 examples in one fixed A4 map. Rank-256 local solubility has density zero after every solvable arithmetic rational base change; thinner families and nonsolvable changes remain open. |
| Extremal real nodes | All even degrees at least eight | Uniform real-algebraic existence, matching-walk chambers and complete spectra. A rational quadratic-algebra chart has an exact Hermite-matrix test for all-real critical points. Rational models in degrees twelve and 26 remain open. |

## Prime knots and augmentation rank

The [protected-tail theorem](../meridional-rank/representations/UNBOUNDED_PRIME_AUGMENTATION_GAP_AUDIT.md)
and [uniform alternating seeds](../meridional-rank/representations/UNIFORM_PROTECTED_AUGMENTATION_GAPS.md)
give explicit knots in every bridge number with

\[
 \operatorname{ar}_{\mathbb C}(K_{b,k})\le b-1<
       \operatorname{mr}(K_{b,k})=\operatorname{br}(K_{b,k})=b.
\]

Writing d=b-2, their protected Coxeter exponents are
(d+1)(1+d²k), with an explicit full A_(2b) meridional quotient.
Marked two-cusped parents and rigorous interval hyperbolicity certificates
prove infinite prime tails at [bridge five](../meridional-rank/representations/PRIME_FIVE_BRIDGE_AUGMENTATION_GAP_FAMILY.md)
and bridge six. The slopes are respectively (1,-12k) and (1,-20k).
The [threaded-curl theorem](../meridional-rank/representations/UNBOUNDED_PRIME_AUGMENTATION_GAP_FAMILIES.md)
identifies all parents with fillings of one fixed certified hyperbolic
three-cusped manifold. Two Dehn filling applications prove the prime-family
result for every sufficiently large bridge number. The exact port induction
supplies the uniform topology; the bridge and filling thresholds are not
made effective.

There is also an exact obstruction to reusing arbitrary powers of the old
five-strand braid: a full S6 transposition quotient exists precisely when
2 divides the power or 147 divides it. For the corresponding known prime
knots it forces augmentation rank five. The first old power left undecided
by this test is the third power. This limitation of the old power family
does not affect the new uniform protected families.

## Uniform rational existence and Hurwitz spaces

The [uniformity and real-geometry proof](ALTERNATING_HURWITZ_UNIFORMITY_AND_REAL_NODES.md)
establishes several structural facts before any further coefficient search.
The unequal coprime pole orders force indecomposability, ruling out composition
recursions. In admissible degrees m^2+1 with odd m, the split long-cycle
classes are power-rational, so that class-rationality test adds no obstruction.
Ordinary spin lifts do not furnish a tuple sign: a double-transposition
centralizer interchanges its two lifts, and the corresponding reduced spin
kernel is killed.

In degree ten all 17,006,112 transitive fixed-product tuples lie in one
braid orbit. Quotienting the order-nine centralizer gives 1,889,568 absolute
Nielsen points. Thus the explicit rational pencil and the original knot
tuple lie in the same absolute Hurwitz component. Its rational point and
unique geometric component establish the indicated Q-definition, without
proving unirationality or classifying inner components.

The next rational literal-square degree is 26. In the normalized even-H
locus the search is equivalent to one degree-eleven polynomial G, a
squarefree critical polynomial S of degree twelve, and a characteristic
polynomial that is the square of a squarefree degree-six polynomial.
The normalized constant a must itself be a rational square. These are exact
necessary-and-sufficient equations for that symmetric locus. The
[binomial-critical obstruction](PAIRED_CRITICAL_BINOMIAL_OBSTRUCTION.md)
excludes S(u)=c(u^(b-1)-k) uniformly: it forces one double-transposition
branch, all remaining finite branches transpositions, and full symmetric
monodromy. The [cycle-folding theorem](ALTERNATING_HURWITZ_CYCLE_FOLDING.md)
braids every transitive tuple to adjacent doubled-edge blocks with a tree
complement, then coalesces to an odd-degree alternating polynomial tuple.
The [degree-twelve proof](DEGREE_TWELVE_HURWITZ_ORBIT.md) settles the next
complex case. After coalescing the doubled edge, the complete marked
polynomial spaces on eleven sheets have respectively 14,992,384 and
3,748,096 tuples. Each is connected under braids that retain the marked
block. Eight short stabilizer loops connect every attachment, internal
two-block braids connect all splittings, and the literal word (-2,1)
joins the two coalescence types. Thus every transitive degree-twelve tuple
lies in one absolute full-A12 orbit. The unique geometric component is
Q-defined; no rational point is supplied. Degree fourteen is the next
unresolved even case in this sequence.
Four explicit degree-seven stabilizer loops also prove attachment
transitivity in both coalesced types uniformly, by appending paired
leaf edges. The [unmarked polynomial reduction](UNMARKED_POLYNOMIAL_HURWITZ_REDUCTION.md)
then removes the marking condition itself. An explicit full-A_(2r+1)
tuple has a literal cyclic braid stabilizer transitive on its branch
labels, while the other polynomial passport has a unique 3-cycle value.
Thus two unmarked polynomial-cactus connectivity assertions suffice for
the all-degree original passport. Neither remaining assertion is proved
uniformly; degree fourteen remains the next complex case.

The [overlap and leaf-extraction theorem](PAIRED_HURWITZ_LEAF_EXTRACTION.md)
now applies to arbitrary full-alternating polynomial DT input: every such
tuple reaches adjacent factors whose supports intersect in one label.
A binary self-orthogonal-code obstruction produces an odd overlap, and
a five-cycle rotation reduces overlap three to overlap one. Three further
local moves exchange the secondary edges. A separate pure four-move
identity extracts two leaves under stated external-support conditions,
with an explicit block braid repairing a wider class of interactions.
This gives an exact recursion for the cyclic tuples down to degree nine.
It does not yet give arbitrary-input leaf extraction or uniform connectivity.

The [fixed-leaf obstruction](FIXED_LEAF_HURWITZ_OBSTRUCTION.md)
explains one failed route to that induction. Deleting a chosen leaf
commutes with every elementary move that keeps it a leaf. An explicit
full-A_(2r+1) family has imprimitive deletion group C2 wr S_r, which
prevents its chosen factor from gaining a second leaf under that restriction.
A length-r braid does extract two leaves by temporarily moving the chosen
sheet in two factors. Thus this is a uniform strategy obstruction with
a constructive escape, without an unrestricted braid-orbit distinction.

The [imprimitive-deletion normal form](IMPRIMITIVE_LEAF_DELETION_NORMAL_FORM.md)
now classifies every chosen-leaf deletion: its group is either S_(2r)
or C2 wr S_r. In the latter case all block swaps form a projected tree;
a lifted transposition-star braid and a root adjustment normalize the
entire tuple, including all attachments, while retaining the chosen leaf.
For fixed product and leaf, every such input reaches the same literal
endpoint. Its r-2 removable factors give an exact recursion with full
alternating cores down to degree seven. Thus the whole imprimitive case
is closed. Reaching this locus from a tuple with only symmetric deletions
is the remaining gate for the pure-DT polynomial problem. The separate
mixed 3-cycle polynomial problem in the uniform reduction remains open.

The [reciprocal ODE theorem](PAIRED_CRITICAL_RECIPROCAL_OBSTRUCTION.md)
classifies the direct extension of the degree-ten reciprocal construction.
In degree 26 its uniform-sign parameter polynomial has no projective zero
modulo eleven, excluding rational parameters. Remaining reciprocal models
must have mixed-sign factors. A [complete finite-field census and lifting
obstruction](DEGREE_TWENTY_SIX_RECIPROCAL_LOCAL_OBSTRUCTION.md) shows that
normalized reciprocal degree-26 models cannot have good paired reduction
at three. Models with denominators or degenerate reduction remain possible.

The [coupled reciprocal equations](MIXED_RECIPROCAL_COUPLED_ODE.md)
give exact rational systems for every remaining mixed-sign case, including
built-in quadratic descent. The degree-26 square-parameter split (2,10)
reduces to four integer equations in four variables; the other splits
have seven rational variables. The same method completely classifies
normalized rational reciprocal degree ten: the known model is unique.
A projective modulo-eleven obstruction rules out mixed signs even with
nonsquare reciprocal parameter and arbitrary rational denominators.
This classification concerns reciprocal critical pairing, rather than all
rational alternating pencils.

The [uniform finiteness and integrality theorem](RECIPROCAL_PARAMETER_FINITENESS_AND_INTEGRALITY.md)
changes the reciprocal search in every degree 4h+2. For each sign split
p+q=h, the raw normalized equations form a finite complete intersection
of length

\[
                         2h\binom{2h}{p}\binom{2h}{q}.
\]

These lengths include degenerate points and multiplicities. All normalized
complex solutions are algebraic. Their parameters and numerator
coefficients are integral outside primes at most 4h+1, and the reciprocal
parameter is a unit there. The proof rules out every nonzero point at
weighted infinity by reducing to a polynomial whose critical points all
map to zero. It supplies denominator control, without assuming good reduction.

That control yields global exclusions in degree 26. The entire necessary
closure for [split (2,10)](DEGREE_TWENTY_SIX_MIXED_2_10_GLOBAL_OBSTRUCTION.md)
has no Q61 point. For [split (4,8)](DEGREE_TWENTY_SIX_MIXED_4_8_GLOBAL_OBSTRUCTION.md),
a complete modulo-29 census leaves one strong residue point, and a literal
left-null-vector certificate prevents its lift modulo 29 squared.
Both orientations of both unbalanced splits are therefore excluded over Q,
including arbitrary denominators and bad reduction.

For the [balanced split (6,6)](DEGREE_TWENTY_SIX_BALANCED_LOCAL_SQUARE_OBSTRUCTION.md),
two complete censuses show that its labelled strong closure has no Q37
point. Thus every remaining rational reciprocal degree-26 model must have
a nonsquare reciprocal parameter lambda, with (lambda/37)=-1.
Its signed square class is supported on primes at most 25, leaving at most
1024 classes before this residue test. The balanced nonsquare case remains
open. These ansatz exclusions do not exclude nonreciprocal degree-26 models.

The [balanced rational quotient](BALANCED_RECIPROCAL_RATIONAL_QUOTIENT.md)
processes all quadratic descents in one scheme over Q. For every balanced
degree 8k+2, the sign-exchange involution has a unit anti-invariant r;
explicit invariant coordinates give a quadratic finite etale cover of
the quotient. Its raw length is 2k times binomial(4k,k) squared.
In degree 26 this is 290400, with seven raw equations and thirteen
strong equations in seven rational variables. Good rational points
correspond exactly to normalized rational balanced models, including
nonsquare reciprocal parameters. A terminating exact rational-candidate
procedure is proved, but its degree-26 elimination and point list have
not been computed.

The [three-adic construction](DEGREE_TWENTY_SIX_THREE_ADIC_PENCILS.md)
gives a five-parameter family over Z3 with literal square discriminant and
full regular A26 over Q3(s). It uses an involution of the critical etale
algebra. Its tangent-rank proof extends uniformly: any good finite-field
paired model in degree 4h+2, in odd characteristic not dividing 4h+1,
lifts to an (h-1)-parameter local family. The criterion assumes that
finite-field model and supplies no rational descent to Q.

The [quadratic critical-algebra chart](PAIRED_CRITICAL_TRACE_CHART.md)
describes the entire normalized even-H paired locus without reciprocity.
Over its squared-critical-value algebra B, the critical-root algebra has
a canonical presentation

\[
 A=B[u]/\bigl(u^2-(2\rho+w/a^2)u+\rho^2\bigr),
                         \qquad N(u)=a(u-\rho).
\]

Here a and rho are units of B and can vary across its components.
The chart and its inverse use rational operations, with explicit
nondegeneracy conditions and interpolation equations. It gives eighteen
coordinates and thirteen equations in normalized degree 26, and nine
coordinates and six equations in degree twelve. It also gives the uniform
norm condition N(0)=(-1)^(b+1)(n-1)Norm(rho)^2. This is an implicit
coefficient locus; its rationality and rational points are not established.

## A single fixed map carrying the arithmetic pencil

The [fixed-map theorem](../jacobian/monodromy-theory/FIXED_KELLER_MAP_ALTERNATING_FAMILY.md)
works for every polynomial pencil H(T)-sT+a in characteristic zero. A weighted
choice of polynomial profiles in the first output gives an exact identity

\[
 p^8E_{\Gamma(p)}(S)
 =H(\tau+p^3S)-(H'(\tau)-p^5)(\tau+p^3S)+a.
\]

The derivative is a unit in every separable root algebra, and the explicit
inverse reconstructs the entire scheme-theoretic fibre. This proves fullness
for one fixed map, rather than a collection of maps chosen after specialization.

For the rational degree-ten model, the fixed map has determinant one,
coordinate degrees at most 146, generic degree ten, and target curve

\[
 \Gamma(p)=(p,-574772/p^2,-2(38112+p^5)/p^8).
\]

The specialized pencil parameter is 1016544-p^5. Perfectness of A10
preserves its regular group under this solvable radical base change;
Hilbert irreducibility then gives infinitely many rational p with full
A10 field fibres. The ambient fixed map has generic arithmetic and
geometric S10. This resolves the fourth requested target in dimension three.

## Generic alternating Keller realization

The [generic architecture audit](../jacobian/monodromy-theory/ALTERNATING_GENERIC_KELLER_AUDIT.md)
proves that arbitrary polynomial coefficient profiles still yield S_N when
the independent quadratic target tilt remains. The old Laurent extension
also retains its normalization-unit obstruction after stabilization. A
different polynomial extension with the same abstract group is not excluded.

The audit additionally gives a uniform odd-prime counterexample to a
field-theoretic primitive-monomial criterion stated in an inspected
prime-degree Keller preprint. The examples have nonconstant Jacobian. They
invalidate that stated auxiliary criterion, without contradicting the
preprint's Keller conclusion. Thus that criterion cannot supply a checked
universal exclusion here.

The graded constructions of Kistner–Shaska in the preprint
[arXiv:2608.02863v1](https://arxiv.org/pdf/2608.02863v1), supply the different
polynomial extensions sought in dimension three. The
[explicit cyclic-cone theorem](../jacobian/monodromy-theory/UNIFORM_EVEN_ALTERNATING_CONE_KELLER_MAPS.md)
gives a compact two-term seed in every even degree n>=6, a polynomial
Jacobian determinant -1, and full ambient geometric A_n. For every odd m>=3,
n=m^2+1 has coefficients over the same field Q(sqrt(-3)) and regular
arithmetic A_n, with component degrees (3n-3,6n-9,3).

The [uniform odd-shear theorem](../jacobian/monodromy-theory/UNIFORM_ALTERNATING_SHEAR_KELLER_MAPS.md)
also gives geometric A_(m(m+1)) over Q(i) for every odd m>=3, with exact
arithmetic constant field. Its degree-twelve member is the preprint's model.
Both constructions have individual three-cycle finite inertia. They do not
give the knot passport of paired transpositions or rational paired-critical
degree-26 pencils. The generic alternating existence question is solved in
dimension three over number fields; a plane realization or universal plane
obstruction remains open. No first-existence or novelty claim is intended.

For each admissible degree, the same cyclic-cone map also has infinitely
many complete A_n field fibres over Q(sqrt(-3)) on the literal affine
line (Y_0,Z,1), for a fixed good nonzero rational Y_0. Exact unit-based
reconstruction proves that these are entire fibre algebras. Its ambient
group is A_n as well. This uses the three-cycle pencil over the quadratic
field, rather than the original rational double-transposition pencil.

At the smallest plane group,
[every faithful transitive A4 action is excluded](../jacobian/monodromy-theory/PLANE_A4_MONODROMY_EXCLUSION.md).
Its possible degrees are 4,6,12. Point-stabilizer normal generation excludes
6 and 12; Domrina's unconditional published four-sheet theorem excludes 4.
The degree-five primary theorem reported in the literature has not yet
been inspected here, so it is tracked as a source-verification obligation,
not claimed to be a mathematically open result.

## Hasse counterexamples with ramification control

The [subgroup-cover and arithmetic proof](ALTERNATING_SUBGROUP_COVERS_AND_HASSE.md)
first separates element coverage from decomposition-group coverage. Three
small subgroup classes cover every element of A10 but miss a possible
double-transposition-inertia local group C2 x C6. Adding the two-subset
stabilizer repairs every group in that local domain and gives total degree

\[
                  10+45+120+126=301.
\]

A Hensel neighbourhood at the finitely many exceptional primes and an
integral Hilbert specialization argument give infinitely many pairwise
linearly disjoint A10 splitting fields. The product of their four
stabilizer fields is locally soluble at every place and has no rational
point. A perfect-matching construction gives a uniform local cover in
every even degree, conditional on having the required rational regular cover.

The pinned maximal-subgroup calculation proves rank 301 optimal for covering
the **whole abstract local-group domain**. Three classes can cover that domain
with minimum rank 1075. These results leave room for a rank-256 arithmetic
family with more restricted actual ramification. The
[zero-branch theorem](ALTERNATING_PENCIL_RAMIFICATION_DENSITY.md) identifies
a density-7/16 set of primes at which odd parameter valuation obstructs
rank-256 local solubility. Its locally soluble integer parameters have density
zero in every arithmetic progression. This rules out an all-parameter
progression construction while allowing thinner infinite families.

The [nonzero-branch residual theorem](NONZERO_BRANCH_RESIDUAL_S6_AND_SOLVABLE_PULLBACKS.md)
strengthens the obstruction uniformly under parametrization. At a nonzero
branch, the residual sextic over the quartic branch field has full S6,
certified by exact Frobenius types (6), (5,1), and (2,1,1,1,1).
Its perfect A6 subgroup survives every solvable residue-field extension.
Riemann--Hurwitz leaves odd contact over the four-element branch orbit,
so every rational base change with solvable arithmetic monodromy still
has rank-256 local solubility of density zero in every integer progression.
This includes all affine power substitutions, even s=t^2, which kills
the original zero-branch inertia. The proof supplies infinitely many
linearly disjoint full-A10 local failures in every progression.

There is also a conditional transfer in every degree n=2 modulo four,
n>=10: a regular rational A_n pencil with a branch orbit of size at least
four and residual S_(n-4) has the same obstruction for its odd-subset
and middle-partition element cover, of rank 2^(n-2). This assumes the
stated rational cover and residual group; it constructs no degree-26 model.

The [fixed finite-etale family theorem](../jacobian/monodromy-theory/FIXED_KELLER_MAP_FINITE_ETALE_FAMILIES.md)
places all these rank-301 algebras in complete fibres of one fixed
Jacobian-one A4 map, with generic S301. Its rational target curve records
the parameter. A fixed A3 map for the product family remains open.

## Extremal real geometry

The [rational degree-eight theorem](ALTERNATING_HURWITZ_UNIFORMITY_AND_REAL_NODES.md)
gives an infinite rational degree-eight family with geometric A8, eight real
critical points, four real paired branch values and four real collinear
ordinary nodes. Its arithmetic group is S8, with constant quadratic
resolvent Q(sqrt(-7)); a literal rational discriminant square in degree eight
is excluded by the known leading-coefficient obstruction.

The [matching-walk classification](REAL_PAIRED_LAURENT_MATCHING_WALKS.md)
now proves real-algebraic maximal-node existence in every even degree n>=8,
including even H and arbitrary negation-symmetric real branch sets.
For each ordered branch set it classifies the all-real chambers by a path
matching plus one auxiliary leaf, toggled twice. Every spectrum is an initial
segment {2,6,...,M}; all allowable maxima are constructed. In the odd-cover
locus the maximum is at most n-2 for n divisible by four and n-4 otherwise.
Without odd symmetry, n=2 modulo four also permits totally real fibres.

There are 63,576 all-real covers over a fixed ordered degree-twelve branch
set, and 269,668,672 odd real covers over a fixed symmetric degree-twenty-six
branch set. These are real chamber counts, not rational-point counts.

The next degree without a rational real paired-node coefficient model is
twelve. For literal squares over Q the next admissible degree remains 26.
These are distinct rational existence questions.

In the quadratic chart, all critical points are real exactly when the
two Hermite trace matrices for w and w+4rho a^2 are positive definite.
For rational data this is an exact test by 2h rational principal-minor
inequalities. It describes the full real locus in that chart without
choosing algebraic roots, while leaving the rational-point question open.

The [constant-secant obstruction](SCALAR_SECANT_PROPORTIONAL_SCALE_OBSTRUCTION.md)
rules out two simple coefficient slices in every degree 4h. With scalar a,
the scale rho=cM is impossible over R, even when critical points are
nonreal; a double-zero identity forces a geometric-sum equation whose
real and unit-circle cases are impossible. If both a and rho are scalar,
the model is impossible even over C, by Euler coefficients and the
paired-product norm. A genuinely affine scale rho=r0+cM with both
coefficients nonzero remains outside these exclusions.

The [degree-twelve reciprocal obstruction](RECIPROCAL_DEGREE_TWELVE_OBSTRUCTION.md)
now excludes that ansatz globally over Q, using a projective modulo-31
eliminant certificate, and excludes all-real reciprocal degree-twelve
models even over R. A uniform weighted-Rolle lemma controls the possible
sign splits. These are structural restrictions on reciprocal critical
pairing, not an exclusion of general degree-twelve real pencils.

The [uniform endpoint-sign theorem](ODD_K_REAL_RECIPROCAL_OBSTRUCTION.md)
strengthens the real conclusion: every all-real good reciprocal model
of degree 4k requires even k. Thus the entire degree class 4 modulo 8
is excluded within reciprocal pairing. Its proof combines a negative
first critical value with the forced alternating order of the two sign
types; it uses no rationality or elimination hypothesis. The known general
real-algebraic paired models in those degrees must have nonreciprocal
critical pairing.
