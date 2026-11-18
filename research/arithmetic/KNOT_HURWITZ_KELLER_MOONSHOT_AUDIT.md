# Seven knot, alternating-cover and Keller research targets

Infinite prime augmentation-gap families now exist at every sufficiently
large bridge number. The real paired-cover locus has
a uniform matching-walk classification and complete real-root spectra.
One fixed Keller map carries the degree-ten A10 pencil in dimension three;
arbitrary rational finite-etale families have fixed-map realization in
dimension four. A different, known alternating cover now supplies uniform
Hasse families of rank 2^(n-2), including the minimum rank 256 in degree ten.
The degree-fourteen alternating Hurwitz passport has one absolute component. Known graded
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
| Rational square pencils | Degree 26 | The reciprocal ansatz is excluded over Q. The complete affine-trace slice lies on an explicit curve of genus at least 1051, so its normalized rational models form a finite set; their existence remains undecided. A good number-field point of degree at most twelve is certified. Normalized rational even-numerator degree-ten models also form a finite set. General rational degree-26 existence remains open. |
| Hurwitz classification | Degree fourteen and an infinite D5-lobe family | One absolute full-A14 orbit. Ordered one-vertex D5 lobe trees reach the canonical wreath-leaf form in at most quartic cost. The distance grows quadratically for the common-root bouquets. Imprimitive leaf deletion and specified coalesced attachment holonomy are classified. Degree sixteen and arbitrary-tuple uniform connectivity remain open. |
| One fixed Keller map | Every admissible degree | One fixed rational A3 map for each degree m^2+1, odd m>=3, has infinitely many complete A_n field fibres. The paired degree-ten pencil also has its own fixed map. |
| Generic alternating Keller group | All even natural actions n>=6; the smallest plane case A4 | Polynomial A3 realizations in all even degrees; regular arithmetic A_(m^2+1) over Q(sqrt(-3)). Every plane A4 action is excluded. General plane realization or exclusion remains open. |
| Hasse algebras | Uniform admissible degrees | A different concentrated-inertia cover gives minimum-rank-2^(n-2) Hasse families and one fixed A4 map per degree. Numerical optimality follows from Maroti's prior covering theorem; the optimal subgroup family is classified uniformly. The original double-transposition pencil retains its rank-256 density obstruction after every solvable rational base change. |
| Extremal real nodes | All even degrees at least eight | Uniform real-algebraic existence, matching-walk chambers and spectra. Affine trace is excluded in degree 4h>=12. The quadratic spectral slice has an exact pairing criterion, real-solubility gate and an unavoidable accessory divisor. General rational models in degrees twelve and 26 remain open. |

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

Arithmetic is more restrictive in the even-numerator degree-ten locus.
The [degree-ten finiteness theorem](DEGREE_TEN_AFFINE_TRACE_CURVE_ARITHMETIC.md)
proves that its normalized rational models form a finite set. When the
canonical scale is nonscalar, the rank-two value algebra makes its trace
automatically affine in that scale. The resulting explicit quartic cover
of the slope line has geometric S4 and at least thirty simple branch
points, so its normalization has genus at least twelve. Faltings gives
finitely many rational points; the scalar-scale boundary is covered by
the earlier reciprocal classification. This gives neither an effective
point list nor finiteness of unnormalized rational twists. It is a
statement about this symmetric coefficient locus, not the dimension or
unirationality of the full degree-ten Hurwitz component.
For each model its pencil parameter remains free and can yield infinitely
many Hilbert specializations.

The known [Arala alternating lines](ARALA_LAURENT_PENCILS_AND_KELLER_FIBRES.md)
already give rational literal-square pencils with regular A_n in every
admissible degree. In Laurent coordinates they are
H(z)=P(z+c)-P(c), a=P(c), where P(X)=X^n-nX^(n-1) and c=(n-2)^2/n.
Their two finite inertia types are an (n-1)-cycle and a three-cycle.
The high critical multiplicities distinguish them from the required
paired-simple-critical pencil. An index-additivity proof shows that
coalescing finite values of the paired passport, with degree, pole type
and smooth genus-zero source fixed, cannot produce its pure three-cycle.
Thus the known construction does not close the degree-26 paired problem.

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
Q-defined; no rational point is supplied.
The [degree-fourteen theorem](DEGREE_FOURTEEN_HURWITZ_ORBIT.md) now extends
this finite classification by complete polynomial censuses on thirteen
sheets. Their absolute orbit sizes are 446,265,625 for six DTs and
89,253,125 for a marked three-cycle followed by five DTs. Exact character
inner products over all 101 cycle types give independent total bounds,
and collision-free breadth-first traversal reaches those bounds. The
uniform folding, marking and attachment theorems then give one absolute
full-A14 orbit. This supplies no rational point or unirationality result.
Degree sixteen is the next unclassified even case in this sequence.
Four explicit degree-seven stabilizer loops also prove attachment
transitivity in both coalesced types uniformly, by appending paired
leaf edges. The [unmarked polynomial reduction](UNMARKED_POLYNOMIAL_HURWITZ_REDUCTION.md)
then removes the marking condition itself. An explicit full-A_(2r+1)
tuple has a literal cyclic braid stabilizer transitive on its branch
labels, while the other polynomial passport has a unique 3-cycle value.
For each remaining original degree n=2r+2≥16, it suffices to prove both
unmarked polynomial connectivity assertions on 2r+1 sheets. Marked
and unmarked connectivity are equivalent for each polynomial space.
This gives a sufficient criterion for the original passport.

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
alternating cores down to degree seven. For the residual pure-DT problem
with r≥7, transitivity is equivalent to every orbit having a representative
with an imprimitive chosen-leaf deletion. Each such representative reaches
the same normal form. The residual mixed 3-cycle problem also remains open.

The [full coalesced holonomy theorem](FULL_COALESCED_ATTACHMENT_HOLONOMY.md)
classifies the attachment action of the specified polynomial bases.
On 2r+1 sheets its outside group is S_(2r-3) in the DT case, and
A_(2r-2) in the mixed three-cycle case for r>=4. At r=3 the mixed
group is exactly the regular Klein four group, certified on the complete
196-state marked move graph. All these groups remain attainable by pure
literal loops preserving each ordered splitting. Their ordered multiple-choice
orbits are therefore explicit. These are actions at the same coalesced
interface; they do not identify the separate leaf-extension interface or
prove arbitrary-tuple connectivity.

The [normalization-distance theorem](UNBOUNDED_HURWITZ_LEAF_NORMALIZATION_DISTANCE.md)
shows why a bounded local repair cannot close that gate. For every r>=4,
explicit full-A_(2r+1) tuples have symmetric deletion at every initial leaf,
and every word reaching a wreath leaf has length at least floor(r/2)-1.
Independent five-sheet dihedral subgroups obstruct an invariant pairing;
one elementary move can expose at most one such subgroup. A fresh buffer
does not remove the growing bound. This also rules out repairs confined
to a bounded set of original factor colors, even with unlimited internal
braids. It permits a rule propagating through linearly many factors and
does not distinguish unrestricted Hurwitz components.

The [root-preserving absorption theorem](D5_LOBE_BOUNDED_CORE_ABSORPTION.md)
now closes this reachability question for the entire unbuffered D5 bouquet
family. Two sweeps isolate one fixed thirteen-sheet core; a finite
connectivity theorem and a final slide restore the canonical wreath star
while retaining the common root. The construction can therefore absorb
every successive lobe. The finite controller words are proved to exist
with explicit length bounds; their literal words have not been recovered.

The [quadratic-distance theorem](QUADRATIC_HURWITZ_BOUQUET_DISTANCE.md)
gives matching asymptotic bounds. For a bouquet with t lobes and 2t factors,
let d_t be the least ordinary signed adjacent Hurwitz word length reaching
any wreath leaf, allowing common sheet relabelings for free. Then

\[
 \left\lceil\frac{(2t-1)\max(0,2\lceil t/2\rceil-3)}8\right\rceil
 \le d_t\le
 6656+(t-2)\,446265984+14t^2-52t+48.
\]

In particular d_t=Theta(t^2). Each move enlarges the total historical
support span of all sheets by at most eight. Every endpoint block has
two jointly occurring sheets whose original lobe positions are far apart;
these disjoint pairs force quadratic total growth. Thus the whole family
lies in one known orbit, but a linear adjacent-word normalization algorithm
is impossible. Sharp distance for the buffered variant, arbitrary polynomial
tuples and degree-sixteen connectivity remain open.

The [ordered one-vertex lobe-tree theorem](NATURAL_D5_LOBE_TREE_HURWITZ_REACHABILITY.md)
allows each adjacent natural D5 pair to attach at any one old sheet and
to share any vertex of its own path. A full twist on the old prefix
moves that sheet to the canonical root, and at most two moves root the
new pair before bounded-core absorption. For t >= 2, all these tuples
lie in the absolute orbit of W_(2t) with full group A_(4t+1). The
general upper bound is O(t^4). The bouquet subfamily gives a
quadratic worst-case lower bound for the larger class.

A [two-leaf absorption theorem](UNIFORM_TWO_LEAF_DT_ABSORPTION.md)
allows one final double transposition with two fresh leaves at any pair
of distinct old sheets. This gives reachability for the buffered bouquets
and buffered ordered trees. The theorem applies to every canonical
component with r ≥ 3. At r = 2, the five adjacent-anchor extensions
stay in the other full-A7 orbit. For r ≥ 7 the pure-DT connectivity
problem is now equivalent to reaching a representative with a removable
two-leaf factor. The proven degree-thirteen case starts the induction.

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
two complete censuses exclude every labelled Q37 point, supplying the
square-parameter case. The
[complete reciprocal obstruction](DEGREE_TWENTY_SIX_RECIPROCAL_GLOBAL_OBSTRUCTION.md)
now excludes the nonsquare case too. For every one of the eighteen
nonsquare residues of eta modulo 37, the thirteen rational-quotient
equations generate the unit ideal. Independent direct Euler/Chebyshev
reconstruction binds every equation, and ordinary exact Singular checks
confirm all eighteen eliminations. Uniform integrality then excludes
the entire strong quotient over Q37, including every denominator pattern
and degenerate residue. It also excludes points over every finite local
extension with residue field F37, regardless of ramification.

All sign classes are therefore excluded: no good rational reciprocal
odd-Laurent model exists in degree 26. Rational normalization of its unique
zero-level double root loses no model. This is an obstruction to extending
the degree-ten construction through constant-product critical pairing;
general nonreciprocal degree-26 existence remains open.

The [odd-numerator affine-trace theorem](ODD_AFFINE_CRITICAL_TRACE_CURVES.md)
provides a uniform replacement for reciprocal coefficient searches.
For N of degree 2h+1 with its zero-level critical root normalized to one,
theta=s+t*rho forces a trace polynomial to equal
2(4h+1+k_h(t)R)q(R)^2. A universal Euler identity makes this equivalent
to a linear polynomial operator on q. Its triangular monic kernel gives
one accessory equation Phi_h(s,t)=N(1)=0, of degree at most 2h in s,
with an exact converse on explicitly listed nonvanishing conditions.
The secant may vary. No good model has slope minus two in characteristic
zero; in degrees ten and twenty-six that is the only possible rational
exception to the recurrence, including every rational denominator.

In degree 26 the retained curve has degree twelve in s. A good simple
point modulo 29 lifts to Q29 and to a number field of degree at most
twelve, giving full regular A26 and nonconstant paired critical product.
This proves the good open is nonempty. No rational Q-point has yet been
found or excluded. The theorem classifies the affine-trace slice of the
odd Laurent locus, rather than every degree-26 polynomial pencil.

The [degree-26 arithmetic theorem](Q26_HIGH_GENUS_AFFINE_TRACE_CURVE.md)
now proves finiteness in that whole normalized rational slice. Its primitive
accessory equation has degree twelve in s and 144 in t, with no discarded
vertical component. Exact Frobenius witnesses give arithmetic S12 over t.
A modular discriminant identity, independently bound by 3,169 extension-field
evaluations and a degree bound, certifies 2,124 simple characteristic-zero
branch points. Geometric S12 and Riemann--Hurwitz give genus at least 1051.
Faltings therefore gives a finite rational point set, without an effective
list. It also excludes a nonconstant rationally parametrized family of
normalized affine-trace models. Whether any of those rational points lies
in the good open remains undecided.

Here S12 acts on accessory roots over the slope t; reconstructed pencils
have A26 monodromy over their independent parameter w. The finiteness
theorem allows infinitely many w-specializations of any rational cover
that is eventually found. Non-affine canonical trace remains outside
this arithmetic restriction.

The [small-prime branch theorem](Q26_SMALL_PRIME_LOCAL_BRANCHES.md)
rules out local obstruction at two, three and the real place.
Its integral infinity charts are G(X,2V) and 3^18G(X/9,9V), both reducing
to X^11(X+1) at their respective primes. Simple roots give infinitely
many good local models, including nonintegral three-adic parameters.
Two simple real infinity roots also give good real coefficient models.
This supplies neither a rational point nor an all-real critical-point
model. Solubility at every prime remains undetermined.

The [local companion through prime 29](Q26_LOCAL_SOLUBILITY_THROUGH_PRIME_TWENTY_NINE.md)
extends these necessary-condition tests to every prime at most 29.
The direct identity
\[
 5^{-40}F(1056+1200X,238)\equiv3X+4X^2\pmod5
\]
has two simple roots. Replacing 238 by 238+5^41V preserves the entire
reduction for every V in Z5, giving two families with infinitely many
good models. Six nonsingular affine residues handle the other new primes.
Thus none of the tested places obstructs the rational affine-trace
problem. Rational nonemptiness and local solubility at larger primes
remain undecided.

The [uniform denominator theorem](Q26_UNIFORM_TRACE_PARAMETER_DENOMINATORS.md)
adds an arithmetic restriction on every local point. Its general
valuation criterion accounts for all relative denominator orders of
the canonical trace parameters s,t. At prime 43, both parameters are
integral at every affine Q43 point; a joint Frobenius certificate gives
the same conclusion at a positive-natural-density set of primes.
There are nevertheless infinitely many good Q43 models.
The smooth projective normalization has no Q43 or rational point at a
pole of either coordinate. This concerns s and t, rather than all
reconstructed coefficients, and it leaves rational existence undecided.

The [balanced rational quotient](BALANCED_RECIPROCAL_RATIONAL_QUOTIENT.md)
processes all quadratic descents in one scheme over Q. For every balanced
degree 8k+2, the sign-exchange involution has a unit anti-invariant r;
explicit invariant coordinates give a quadratic finite etale cover of
the quotient. Its raw length is 2k times binomial(4k,k) squared.
In degree 26 this is 290400, with seven raw equations and thirteen
strong equations in seven rational variables. Good rational points
correspond exactly to normalized rational balanced models, including
nonsquare reciprocal parameters. A terminating exact rational-candidate
procedure is proved, but its degree-26 elimination and raw rational point
list have not been computed. Its good rational locus is now known to be
empty by the separate 37-adic obstruction; degenerate raw points are a
different question.

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

The [Arala transfer](ARALA_LAURENT_PENCILS_AND_KELLER_FIBRES.md)
makes this uniform over Q in every admissible degree n=m^2+1, odd m>=3.
One fixed A3 map per degree carries infinitely many complete degree-n
field fibres with splitting group A_n. Its coordinate degrees are at most
20n-54, or 466 in degree 26. These are root-field fibres of the concentrated
cover, with ambient generic group S_n. They do not realize its Laurent
extension as the generic extension of that polynomial map.

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
family with more restricted actual ramification. For the original
double-transposition pencil, the
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

There is now an infinite rank-256 family from a different cover. The
[concentrated-inertia theorem](ARALA_LAURENT_PENCILS_AND_KELLER_FIBRES.md)
uses Arala's known cover and works in every admissible degree n=m^2+1.
The odd-subset fixed fields together with the unordered middle-partition
field have total rank 2^(n-2). Powers of the long-cycle inertia have a
unique fixed letter; nontrivial three-cycle inertia has invariant support
of size three. Their normalizers therefore lie in the selected stabilizers.
A denominator-safe Hensel progression handles every exceptional prime,
and integral Hilbert irreducibility gives infinitely many pairwise linearly
disjoint A_n splitting fields. The product algebras have points over every
completion and no rational point. These ranks are globally minimum in
every admissible degree, including composite half-degrees such as fifty.
At n=10 the minimum is 256. This different family attains it without
contradicting the obstruction for the original cover.

The [priority and equality addendum](ALTERNATING_WEIGHTED_COVERS_PRIORITY_AND_OPTIMAL_HASSE.md)
identifies the exact prior result: Maroti's Theorem 4.1 already implies
the minimum total coset degree 2^(n-2) for every n=2 modulo four, n>=6,
and its proof gives the same upper-cover family. A separate written
argument classifies equality for every such n>=14, extending the
[prime-half proof](MINIMUM_TOTAL_DEGREES_ALTERNATING_PRIME_HALVES.md).
Every other primitive or imprimitive action exceeds the total degree
budget; two-cycle witnesses force each constituent. This concerns total
degree, not the number of conjugacy classes. Numerical optimality is prior
work; priority of the uniform equality statement and the precise arithmetic
transfer remains unresolved.

There is also a conditional transfer in every degree n=2 modulo four,
n>=10: a regular rational A_n pencil with a branch orbit of size at least
four and residual S_(n-4) has the same obstruction for its odd-subset
and middle-partition element cover, of rank 2^(n-2). This assumes the
stated rational cover and residual group; it constructs no degree-26 model.

The [fixed finite-etale family theorem](../jacobian/monodromy-theory/FIXED_KELLER_MAP_FINITE_ETALE_FAMILIES.md)
places all these rank-301 algebras in complete fibres of one fixed
Jacobian-one A4 map, with generic S301. Its rational target curve records
the parameter. The same compiler places each new rank-2^(n-2) family in
one fixed A4 map, with generic S_(2^(n-2)). A fixed A3 map for these product
families remains open.

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
paired-product norm. These endpoint proofs leave a genuinely affine scale
rho=r0+cM with both coefficients nonzero untreated.

In degree twelve the
[full rational scalar-affine theorem](DEGREE_TWELVE_SCALAR_AFFINE_OBSTRUCTION.md)
now excludes that remaining slice over Q as well, including nonreal critical
points. A rational normalization leads to a uniform square/Hermite
condition. Five retained arithmetic circuits prove unit identities for
the generic parameter and every rational exception. Exact factor-product
identities and projective finite-field no-root certificates account for
all parameter denominators; the replay needs no Groebner or factorization
oracle. Thus a rational even-H degree-twelve model must have a nonconstant
secant element or a genuinely quadratic scale in its cubic value algebra.
The unrestricted rational real-node problem remains open.

The [uniform real affine-trace theorem](UNIFORM_REAL_AFFINE_TRACE_OBSTRUCTION.md)
now excludes a larger slice in every degree 4h>=12, without requiring real
critical points. The canonical pair trace cannot satisfy theta=s+t*rho
for real scalars s,t, even if the secant varies. The critical jets force
a polynomial square; a positive real polynomial forces slope t=-2.
A reflection identity then gives a difference operator whose leading
coefficient is j(j-2)(j-4)/3. It cannot annihilate the required degree
j=2h>=6 polynomial. The degree-eight model gives the sharp j=4 boundary.
Scalar secant and affine scale imply affine trace, so that entire slice
is excluded for real coefficients in these degrees. The earlier degree-twelve
arithmetic certificates provide a separate exact check; unrestricted
nonaffine-trace rational models remain open.

The [quadratic spectral criterion](QUADRATIC_SPECTRAL_PAIRING_AND_ELLIPTIC_DIVISORS.md)
treats a next slice with nonconstant secant: theta=t0+t1*R+t2*R^2 and
squared value M=kappa(R-eta). In degree twelve it gives five exact
equations in the six q/theta coefficients when eta is fixed and kappa
is recovered. The trace difference is q^2 times an accessory polynomial
of exact degree five. On the smooth elliptic spectral curve the resulting
divisor relation is 6B+5B-11B=0, an identity that gives no torsion order
for B. The natural norm also has simple forced zeros and a nonconstant
right side, so it supplies no constant-norm Pell shortcut.

The same note gives the complete real-solubility test for the square-root
spectral cover. With alpha=t0*t2 and b=t1^2/4-1-alpha, positive t0 or t2
ensures real points. If both are negative, the smooth cover has real
points exactly when b>0 and b^2>4alpha. This is necessary for an all-real
critical pencil; rational solubility of the cover is not required merely
for rational pencil coefficients. An exact stronger test for all critical
points is that q splits over R and theta(r)>2|r| at every q root.
The full rational degree-twelve existence question and singular spectral
strata remain open.

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
