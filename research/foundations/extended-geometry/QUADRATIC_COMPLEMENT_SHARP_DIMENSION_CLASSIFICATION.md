# Sharp Mathieu dimension classification for quadratic norm complements

Ordinary internal mathematical deduction. This is one written proof, with no
invented checker or computational certification. All assurance flags remain
false. The low-dimensional proof and the even absolute-nonpower argument are
included in full. The main higher-dimensional input is the uniform theorem identified precisely
in Section8. Section7 retains a separately supported even-dimensional route
as an independent comparison; the sharp classification does not depend on it.

Let K be ANY field of characteristic zero, q a NONDEGENERATE quadratic form on
a positive-dimensional K-vector space V of dimension N, and
A_q=K[V,q^{-1}]. Let Lambda_q be the normalized Reynolds functional for the
geometrically transitive connected reductive action SO(q) x G_m.
Define E(q) to mean:

  there exist regular F,H in A_q such that
  (a) F is absolutely non-proper-power EVEN UP TO A UNIT: over Kbar, F is
      not U G^e for any unit U, any G in A_q tensor_K Kbar, and integer e>=2;
  (b) Lambda_q(F^m)=0 for EVERY integer m>=1;
  (c) Lambda_q(H F^m)!=0 for EVERY integer m>=1, with ONE H fixed for all m.

**Theorem SharpQuadraticDimension.** For every such K,V,q, the following are
equivalent:

  ker(Lambda_q) is a Mathieu--Zhao subspace;
  N<=2;
  E(q) is false.

Equivalently, ker(Lambda_q) is not Mathieu iff N>=3 iff E(q). The last
existential assertion is obtained by the uniform all-dimensional construction in
this geometric category, not by a general conversion of every unsafe phase
into a nonpower phase or an all-orders witness.

Both higher-dimensional constructions satisfy the strong unit-times-power
exclusion in E(q). In particular they exclude ordinary proper powers. The same
classification therefore also holds if E(q) is defined using only the weaker
ordinary absolute non-proper-power condition: that weaker witness would still
contradict Mathieu, and dimension>=3 already supplies the stronger one.

The theorem is not about degenerate forms of small rank in larger dimension,
nor about arbitrary connected reductive homogeneous varieties. Nondegeneracy,
positive dimension and the stated normalized trace are part of every clause.

## 1. Exact statement and category

Let K be ANY field of characteristic zero. Let V be a K-vector space of positive
dimension N, and let q be a NONDEGENERATE quadratic form on V. Put

    X_q = Spec K[V,q^{-1}],
    G_q = SO(q) x G_m,

with rotations and scalar multiplication, and let Lambda_q be the unique
normalized G_q-invariant Reynolds functional on O(X_q). Say that a phase f is
pure when Lambda_q(f^m)=0 for every integer m>=1. Say that f is safe when, for
every FIXED h in O(X_q), some cutoff M(h) satisfies
Lambda_q(h f^m)=0 for all m>=M(h). The kernel is a Mathieu--Zhao subspace exactly
when every pure phase is safe.

THEOREM. If N=1 or N=2, ker(Lambda_q) is a Mathieu--Zhao subspace over K.
This holds for every such K and q, without an ordering of K, an embedding of K
in C, or a splitness or compactness hypothesis on K-points.

More precisely, after the explicit splitting extension and Laurent coordinates
below, a NONZERO f is pure iff zero is outside its Newton polytope. This implies
safety against every fixed regular marker. For N=1 it says exactly that f has
strictly one-sided Laurent support. For a nonsplit binary form it says exactly
that all nonzero homogeneous scalar weights of f have ONE strict sign.

The statement concerns nondegenerate q in dimensions one and two. It does NOT
assert anything about degenerate forms of rank<=2 in a larger ambient space,
whose symmetry group and geometry belong to a different category. N=0 gives
q=0 and an empty principal open; a normalized functional Lambda(1)=1 is not
part of this category.

## 2. Scalar extension: the exact implications that are valid

Let A be any commutative unital K-algebra and Lambda:A->K any K-linear
functional. For a field extension L/K write A_L=A tensor_K L and
Lambda_L=Lambda tensor id_L. Fix f IN A, so the phase itself descends.

(2a) f is pure over K iff f_L is pure over L.

Indeed Lambda_L(f_L^m)=Lambda(f^m) tensor1, and K->L is injective.

(2b) f is safe over K iff f_L is safe over L.

For the forward direction, every fixed upstairs marker H in A_L is a FINITE
sum H=sum_{i=1}^r c_i h_i, with c_i in L and h_i in A. If every h_i has a
cutoff, the maximum of the finitely many cutoffs annihilates
Lambda_L(H f_L^m)=sum_i c_i Lambda(h_i f^m). For the reverse direction, use
upstairs safety on h tensor1 for each given h in A, and descend its zero values
by injectivity. Both cutoffs may depend on the given marker.

Thus unsafety also holds iff for this GIVEN phase. Explicitly, if a fixed
upstairs marker survives infinitely often, its finite expansion above has at
least one downstairs h_i surviving infinitely often; otherwise all terms, and
hence H, vanish eventually. A downstairs unsafe marker also pulls up directly.

(2c) Global Mathieu property DESCENDS: if ker(Lambda_L) is Mathieu in A_L,
then ker(Lambda) is Mathieu in A. Apply (2a) to a downstairs pure f, upstairs
Mathieu to f_L, and (2b) back down.

No converse asserting global Mathieu ascent follows: a new pure phase in A_L
need not come from A. In this note, the split and nonsplit torus cases over all
fields are proved directly, not by making that invalid ascent inference.

## 3. Laurent constant terms over every characteristic-zero field

LEMMA. Let F be any characteristic-zero field and B=F[z_1^{+/-1},...,z_d^{+/-1}]
with CT the zero Laurent coefficient. If f is nonzero, the following are
equivalent:

  (i) CT(f^m)=0 for all m>=1;
  (ii) 0 is outside the Newton polytope of f;
  (iii) there is an integral covector ell:Z^d->Z strictly positive on the finite
        support of f.

Consequently ker(CT) is Mathieu over F, and every pure f is safe against each
fixed Laurent marker.

Proof of (i)=>(ii), including the all-field point. The FINITELY MANY nonzero
coefficients of f generate a finitely generated field F_0 over Q inside F.
There is an embedding F_0->C: embed a finite transcendence basis as algebraically
independent complex numbers and extend over the finite algebraic extension,
using algebraic closedness and characteristic zero. This embeds only F_0; it
need not and does not presume an embedding of the whole F. All zero moment
identities still hold after this embedding. Duistermaat--van der Kallen's
complex Laurent theorem puts zero outside the Newton polytope. The support
has not changed, since the coefficient-field embedding is injective. Therefore
the SAME finite Newton polytope over F excludes zero.

For (ii)=>(iii), strictly separate the compact convex hull of finitely many
integral points from zero. A nearby rational covector still gives strict
positive values on this finite set; clearing denominators gives integral ell.
For (iii)=>(i), every exponent in f^m has positive ell-weight, so cannot be zero.

For safety, let delta=min ell(support f)>0, and for a nonzero FIXED h let
b=min ell(support h). Every exponent in h f^m has ell-weight at least b+m*delta.
For all sufficiently large m this is positive, so CT(h f^m)=0. If h=0 or f=0,
the conclusion is immediate. This proves all directions and the stated
marker-dependent cutoff. It uses the elementary sufficiency mechanism as well
as the imported necessity theorem.

The source is the retained actual primary PDF
research/foundations/extended-geometry/sources/dvk-powers-1998.pdf,
SHA256 bb002d12b539a1ff3fdd099f343d3accd921ccf60f8a149f29b07cd0f7e87e1b:
J. J. Duistermaat and W. van der Kallen, Constant terms in powers of a Laurent
polynomial, Indag. Math.9(1998),221--231, DOI10.1016/S0019-3577(98)80020-7.
Actual primary introductory precise Newton-polytope claim on printed221--222,
Section3 minimal-face reduction and Theorems4--5, printed223--224, and the full
proofs in Sections3--5, printed224--230 were read for this note. This is the
ARBITRARY-RANK theorem, not merely the one-variable Theorem2. Field transfer
from F_0 and the marker weight argument are the separate deductions above.
The existing EXT-DVDK import is sufficient; no new primary file or formal
axiom-gap discharge is asserted.

## 4. Dimension one

Write q=a x^2 with a in K^*. Then

    O(X_q)=K[x,(a x^2)^{-1}]=K[x,x^{-1}].

The special orthogonal group in dimension one is trivial and the scalar G_m
acts with x-weight one. Invariance forces every nonzero Laurent weight to have
zero functional value, while normalization gives Lambda_q(1)=1. Thus
Lambda_q=CT_x, and Section3 proves the theorem. This also identifies its scalar
extension with CT over every extension field. In dimension one pure phases are
exactly zero or Laurent polynomials with all exponents strictly positive or
all strictly negative.

## 5. Dimension two: explicit splitting and normalized trace

Diagonalize over K:

    q=a x^2+b y^2,       a,b in K^*.

Such diagonalization is valid in characteristic zero and transports the
normalized trace by the equivariant linear change of coordinates. Choose

    delta^2=-b/a,
    L=K(delta),
    u=x+delta y,         v=x-delta y.

Here [L:K]<=2, delta!=0, and

    q=a u v,
    x=(u+v)/2,           y=(u-v)/(2delta).

Inverting uv makes u and v individually units, so

    O(X_q) tensor_K L = L[u,u^{-1},v,v^{-1}].

Over L, rotations in SO(q) act as (u,v)->(t u,t^{-1}v), while the scalar torus
acts as (u,v)->(s u,s v). The character of the monomial u^i v^j is

    t^{i-j} s^{i+j}.

It is trivial iff i-j=i+j=0, hence iff i=j=0. The unique normalized invariant
functional is therefore exactly CT_{u,v}.

This proves trace existence and scalar compatibility without a pointwise
transitivity assumption on G_q(K). If L=K we are done directly. If [L:K]=2,
its nontrivial automorphism exchanges u and v and conjugates coefficients;
CT_{u,v} is Galois-stable. For any downstairs element its CT belongs to K.
It therefore descends to a K-linear functional, invariant by faithful scalar
descent and normalized at1. Its L-extension is CT. Any other normalized
invariant functional extends to CT by the character calculation, proving
uniqueness. Thus this is precisely Lambda_q.

Section3 proves the split Laurent kernel is Mathieu over L, and Section2(c)
descends this property to K. This proves the binary assertion for every K and
q. Equivalently, purity for a nonzero downstairs f is exactly exclusion of
zero from its split Newton polygon: each direction follows from Sections2(a)
and3. Every fixed downstairs marker is safe by Sections2(b) and3.

As a geometric description, with E=K[T]/(T^2+b/a), the binary principal open is
the torus Res_{E/K} G_m: q=a*N_{E/K}(x+yT), and norm invertibility is equivalent
to invertibility in the quadratic etale algebra E. The multiplication map
SO(q) x G_m -> Res_{E/K} G_m is a degree-two central isogeny, with kernel mu_2;
over a splitting field it is (t,s)->(st,s/t). Its possible failure to be onto
on K-points causes no problem. The coordinate-ring and weight proof above uses
only the scheme action and its rational characters.

## 6. Extra exact description in the nonsplit binary case

Suppose -b/a is nonsquare, so L/K is quadratic. For a K-defined f, its Laurent
support S in the u,v coordinates is stable under (i,j)->(j,i): applying the
nontrivial field automorphism swaps u,v and sends each nonzero coefficient to
a nonzero conjugate coefficient.

If f is nonzero and pure, Section3 supplies a real/integral separating covector
ell(i,j)=p i+q j positive on S. Its swapped covector is also positive on S.
Their sum is

    (p+q)(i+j),

strictly positive on every point of S. In particular p+q!=0, and every scalar
weight i+j occurring in f has the same strict sign. Conversely, if all these
weights have one strict sign, every power has nonzero scalar weight and zero
trace; every fixed marker is safe by the same finite weight-bound argument.

Thus for the nonsplit binary torus the exact equivalence is:

    pure f  <=>  f=0 or all its nonzero scalar homogeneous components have
                  strictly positive degrees, or all strictly negative degrees.

The statement is independent of an ordering on K. The covector separation is
ordinary real convex geometry of the INTEGRAL character lattice; it does not
order coefficients or K-points. This extra equivalence is not asserted in the
split binary case: there the two-dimensional Newton criterion is the precise
one, allowing any separating character direction.


## 7. Independent even-dimensional comparison and its absolute nonpower

Suppose N=2r>=4. Diagonalization over K gives an orthogonal decomposition
q=A+B with A nondegenerate of rank2 and B nondegenerate of rank2r-2.
The existing all-field even-quadratic theorem supplies the REGULAR functions

    F=(q-1)-(q-1)^2 B^(r-1)/q^r,       H=q^{-1},

and proves, for EVERY positive integer m,

    Lambda_q(F^m)=0,
    Lambda_q(H F^m)=(-1)^(m-1)/(m+1) !=0.             (7.1)

The exact imported mathematical proposition is
FMQ-EVERY-EVEN-QUADRATIC-COMPLEMENT-ALL-FIELDS-FAILURE,
with the explicit phase in FMQ-EVEN-QUADRATIC-ACTUAL-REGULAR-DESCENDING-PHASE,
in foundations/extended-geometry/FIELD_MATHIEU_PRODUCT_AND_QUADRATIC_COUNTEREXAMPLES.md.
This theorem holds for every characteristic-zero K and every nondegenerate
q of the stated even dimension, without splitness or K-points on {q=1}.
Its proof uses the regular polynomial beta CDF, an all-order formal-uniform
identity and equivariant trace restriction. Those are the existing exact
mathematical input; no finite moment experiment is imported as proof.

Here is the entire additional absolute-nonpower argument. Over Kbar the
polynomial q-1 is irreducible. Indeed its homogenization q-T^2 is a quadratic
form of rank N+1>=5, while a reducible quadratic is a product of two linear
forms and has rank<=2 in characteristic zero. Thus q-1 is prime in the
polynomial UFD Kbar[V]. It remains a nonunit prime after inverting q, since
q and q-1 are coprime. In A_q tensor Kbar,

    F=(q-1) * [1-(q-1) B^(r-1)/q^r].

The bracket reduces to1 modulo the prime q-1. Therefore the associated
height-one prime valuation satisfies

    ord_(q-1)(F)=1.                                  (7.2)

Every unit has valuationzero. A unit times an e-th power has valuation divisible
by e>=2, which cannot equal1. Hence F is absolutely non-proper-power even up to
a unit. This is a sufficient simple-zero argument; no converse from nonpower
to the existence of a simple zero is asserted. Equations(7.1)--(7.2) establish
E(q) in every even dimension N>=4.

The extra argument is internal algebra, included here in full; it introduces
no new classical source or companion-file dependency. The fixed marker H in
(7.1) is independent of m, and its displayed rational values are nonzero in
every characteristic-zero K.

## 8. Uniform higher-dimensional construction: the main exact input

For every characteristic-zero K, every N>=3 and every nondegenerate q, use

    QUN-ALL-FIELDS-ALL-DIMENSIONS-ABSOLUTE-NONPOWER-COUNTEREXAMPLE

in foundations/extended-geometry/QUADRATIC_COMPLEMENT_UNIFORM_NORMALIZER_COUNTEREXAMPLES.md.
Its written source is SRC-QUADRATIC-UNIFORM-NORMALIZER-WRITTEN-PROOF and its
principal is QUADRATIC-COMPLEMENT-UNIFORM-ABSOLUTE-NONPOWER-COUNTEREXAMPLES.
This one theorem covers BOTH parities with no old even-family premise.

More concretely, write M=N-2>=1 and choose q=A+B with ranks2,M. Its regular
K-defined functions are

    P=(q+1)[q-(q+1)^2 B]^2/q^3,
    F=P^M(1+P),       H=q^M.

The uniform theorem proves that F is not even a unit times any proper power
over Kbar, Lambda_q(F^m)=0 for EVERY m>=1, and the ONE fixed H has

    Lambda_q(H F^m)=sum_{r=0}^m binom(m,r) J_(M,Mm+r) !=0,
    J_(M,k)=M int_0^1 z^(M-1)(1-z^2)^(2k) dz.

Each J_(M,k) is a strictly positive RATIONAL number, for example by expanding
the finite polynomial integrand and by its positive real integral. Hence the
finite sum is a nonzero rational and remains nonzero in every char0 field;
this argument does not order K. These are exact all-order identities, not an
inference from a bounded list of moments.

The cited proof derives beta moments directly on X_q from Lie derivations and
scalar weights, and proves the absolute unit-nonpower conclusion by a simple
zero of1+P at a smooth geometric point with B=0,q=-1/2,P=-1. It does not need a
finite-cover theorem, a regular square root, isotropy, K-points on a norm-one
quadric, a real ordering, or a whole-field complex embedding. We import only
its stated exact mathematical conclusion E(q). Its fuller circuit identities
and threshold refinements are not needed for the classification below.

This is a separately written mathematical input. It is not assumed from a
conjectural extension of another family, a historical classification label,
or finite probes. The optional even construction of Section7 is an independent
comparison, and withdrawing that old even input does not remove the uniform
all-dimensional conclusion or the main classification.

## 9. Every implication in the sharp classification

First, if N<=2, Sections3--5 prove ker(Lambda_q) is Mathieu. These arguments
are independent of the uniform higher-dimensional construction and of the
optional even-dimensional comparison.

Second, if N>=3, the uniform theorem in Section8 establishes E(q) directly
for every such dimension and form. This step has no parity split and no
old even-family premise. Section7 remains an independent sufficient route
for the even-dimensional subfamily.

Third, E(q) implies failure of the Mathieu property directly from the
definition: F is pure, but its ONE fixed marker H has a nonzero marked moment
at every positive order, hence no eventual-zero cutoff. This implication does
not need the nonpower condition; that condition strengthens the supplied
witness, rather than substitutes for purity or the fixed marker.

These three independent steps give each desired converse as follows.

* If the kernel is Mathieu, N>=3 is impossible by the second and third steps,
  so N<=2. Conversely N<=2 implies Mathieu by the first step.
* If N<=2, E(q) is impossible by the first and third steps. Conversely, if
  E(q) is false, N>=3 is impossible by the second step, so N<=2.
* If the kernel is Mathieu, E(q) is false by the definition and the third step.
  Conversely E(q) false implies N<=2 and hence Mathieu, using the separate
  dimension arguments just given.
* If the kernel is not Mathieu, N<=2 is impossible by the first step, so N>=3;
  conversely N>=3 gives E(q) and non-Mathieu by the second and third steps.
* Non-Mathieu implies E(q) through N>=3 and the uniform high-dimensional construction.
  This is a geometric dimension theorem, not a universal operation turning
  an arbitrary pure unsafe phase into an absolutely nonpower all-order one.

There is also an exact GLOBAL scalar-extension corollary WITHIN THIS FAMILY.
For EVERY field extension L/K, the form q_L remains nondegenerate of the same
positive dimension N. Applying the dimension theorem independently over both
fields gives

    ker(Lambda_q) Mathieu  <=> N<=2 <=> ker(Lambda_(q_L)) Mathieu,
    E(q)                  <=> N>=3 <=> E(q_L).

The canonical normalized traces commute with scalar extension: the extension
of Lambda_q is invariant and normalized, hence equals Lambda_(q_L) by uniqueness.
Non-Mathieu is therefore also invariant within this family. The two witness
existence directions use the dimension criterion and the explicit construction
OVER THE FIELD IN QUESTION. They do not descend an arbitrary new upstairs
phase or marker. This special quadratic-family consequence does not assert
Mathieu ascent, failure descent, or witness descent for general algebras or
homogeneous spaces.

Thus both directions of all stated equivalences are proved. No literal phase
purity/safety predicate is available from the universal theorem alone; an
actual witness comes only from the applicable specified construction.

## 10. Exact sources and scope of this one-proof corollary

The MAIN classification uses exactly two mathematical inputs:

* existing EXT-DVDK, the ordinary arbitrary-rank Laurent Newton-polytope
  theorem, with its primary PDF and coefficient-field transfer in Section3;
* QUN-ALL-FIELDS-ALL-DIMENSIONS-ABSOLUTE-NONPOWER-COUNTEREXAMPLE for the uniform
  all-field, all-N>=3 strong fixed-marker witness in Section8.

The independent optional Section7 comparison additionally uses existing
FMQ-EVERY-EVEN-QUADRATIC-COMPLEMENT-ALL-FIELDS-FAILURE for its original even
pair(7.1). That source is a premise of the comparison statement, not of the
main sharp dimension classification.

The low-dimensional field/trace/torus arguments and even simple-zero proof are
fully written here. EXT-DVDK retains its existing external origin and actual
primary resource; the counterexample inputs are INTERNAL mathematical
propositions with their own written proofs. No fake external bibliography
record is created for them, and no historical priority is asserted.

There is no checker for this corollary. Its registry checker is null; the
proof does not invent a computational test to fill that field. Hashes identify
these exact written sources, not a proof-verification or external-review
status. All such assurance flags remain false.

The earlier restricted degree-nine fixed-height theorem is independent of
this classification and is not a premise here. The family-specific GLOBAL
field-extension equivalence in Section9 follows from the sharp dimension
criterion; no such invariance for general algebras or varieties is asserted. The descent and
fixed-phase equivalences used in the low-dimensional proof are exactly those
proved in Section2. Existing unrelated results, publication states and
dependency objects are preserved by the proposed append-only integration.
