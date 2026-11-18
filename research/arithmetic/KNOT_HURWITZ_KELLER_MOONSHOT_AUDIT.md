# Seven knot, alternating-cover and Keller research targets

The first uncovered cases now yield a prime bridge-five augmentation gap,
a complete degree-ten Hurwitz-orbit calculation, one fixed polynomial Keller
map with infinitely many A10 fibres, and a uniform arithmetic Hasse family.
There are also uniform obstructions to several proposed construction methods.
Rational paired-critical existence in degree 26 and generic alternating
polynomial Keller realization remain open.

This document connects the current proofs. Mathematical status and the proof
dependency DAG remain in [MATH_STATUS.json](../MATH_STATUS.json); an exact
computation is not promoted to formal verification, external review, or a
historical priority claim.

| Target | Smallest new case tested | Result and next unresolved case |
| --- | --- | --- |
| Prime augmentation gaps | Bridge five | One prime gap proved; infinite prime bridge-five families and one prime gap at bridge six remain open. |
| Rational square pencils | Degree 26 | Exact symmetric coefficient criterion and necessary conditions; no rational degree-26 model yet. |
| Hurwitz classification | Degree ten | One absolute full-A10 orbit proved by exhaustive finite calculation; degree twelve and inner refinements remain open. |
| One fixed Keller map | Degree-ten pencil | Solved uniformly for a reparametrized pencil; infinitely many complete A10 field fibres of one fixed A3 map. |
| Generic alternating Keller group | Natural A4/A5 actions and the knot-derived degree-eight/ten targets | Elementary group tests pass in those specified actions; the existing compiler and Laurent architectures are obstructed. Universal plane exclusion remains open. |
| Hasse algebras | A10 subgroup covers including ramified groups | Infinite linearly disjoint rank-301 examples; minimal rank 301 for the whole abstract double-transposition local-group domain, using a pinned subgroup classification. |
| Extremal real nodes | Rational degree eight beyond the known degree ten | Infinite rational full-A8 family with four real nodes; uniform real-root restrictions. Degree twelve and complete real-component classification remain open. |

## Prime knots and augmentation rank

The [protected-tail theorem and prime example](../meridional-rank/representations/UNBOUNDED_PRIME_AUGMENTATION_GAP_AUDIT.md)
give a uniform augmentation mechanism. For a specified subgroup of tails,
the old paired root/covector contradiction survives every number of added
strands. A literal five-strand braid in that subgroup closes to a knot with

\[
 \operatorname{ar}_{\mathbb C}\le4<\operatorname{mr}
                    =\operatorname{br}=5.
\]

Its full A10 meridional quotient gives the bridge lower bound, and a rigorous
interval hyperbolicity calculation proves primeness. The natural elementary
chain in this subgroup is composite, so the all-bridge augmentation mechanism
does not by itself prove the requested prime-family theorem.

There is also an exact obstruction to reusing arbitrary powers of the old
five-strand braid: a full S6 transposition quotient exists precisely when
2 divides the power or 147 divides it. For the corresponding known prime
knots it forces augmentation rank five. The first old power left undecided
by this test is the third power. The first missing bridge number for a prime
gap is six.

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
necessary-and-sufficient equations for that symmetric locus, rather than
an existence proof. The next Hurwitz connectivity case is degree twelve.

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
universal exclusion here. A genuinely different polynomial construction or
a universal plane obstruction is still required.

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
family with more restricted actual ramification. The complete algebras become
Keller fibres of generic-S301 compilers; a single fixed map carrying this
entire higher-rank product family remains open.

## Extremal real geometry

The [real-geometry theorem](ALTERNATING_HURWITZ_UNIFORMITY_AND_REAL_NODES.md)
gives an infinite rational degree-eight family with geometric A8, eight real
critical points, four real paired branch values and four real collinear
ordinary nodes. Its arithmetic group is S8, with constant quadratic
resolvent Q(sqrt(-7)); a literal rational discriminant square in degree eight
is excluded by the known leading-coefficient obstruction.

Uniformly, every regular real fibre of a paired-simple-critical even-degree
pencil has real-root count congruent to two modulo four. If H is even and
n is two modulo four, its count is at most n-4. Thus the symmetric
degree-ten example has spectrum {2,6}, despite all ten critical points being
real. The degree-eight real family also has spectrum {2,6}. These jump
sequences give necessary real-component invariants; no completeness theorem
for those invariants or real connected components is proved.

The next uncovered real paired-node degree is twelve. For literal squares
over Q the next admissible degree remains 26. These are distinct rational
existence questions.
