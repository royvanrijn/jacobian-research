# Boxed universal Mathieu validity and constructible counterexample restrictions

20 November 2026. Internal written proof of a terminating theoretical
algorithm. The effective Gaussian whole-source and pure-prefix results are
retained premises. No coefficient-space run in unrestricted degree, practical
complexity bound, rational-point algorithm, formal verification, external
review, or historical-priority assertion is included.

## 1. Inputs, outputs, and the two moment settings

The input is a fixed finite polynomial phase family F_c with coefficients
polynomial in parameters c=(c_1,...,c_r) over Q, together with a specified
constructible parameter restriction defined by finitely many rational
polynomial equalities and inequalities. A specified embedded number field
can replace Q with the same effective arithmetic. Every coefficient field
below is then understood to be a characteristic-zero EXTENSION of that
coefficient base, with the specified base embedding retained. For Q this
is simply every characteristic-zero field. The constructible set is
presented as a finite union of pieces

    V(J) intersect D(g),   D(g)={c:g(c)!=0}.              (B1)

An empty list of equalities is allowed, as is g=1. General finite Boolean
polynomial conditions can effectively be put in this form by finite
Boolean expansion; several simultaneous inequalities become one product.
A condition that some member of a finite list is nonzero is a UNION of
principal opens, not their product.

There are two settings:

* Gaussian: F_c is a polynomial of bounded degree in N fixed variables,
  Gamma_N is the product standard Gaussian Wick functional, and the mixed
  moment for a marker Q is Gamma_N(Q F_c^m).
* Image: F_c is a polynomial of bounded total degree in n fixed contraction
  pairs (xi,s), E_n(xi^alpha s^beta)=partial_s^alpha s^beta, and the mixed
  moment is the ENTIRE polynomial E_n(Q F_c^m) in the output coordinates.
  No balance or homogeneity assumption is imposed.

In either setting pure means that every positive unit moment is zero.
Safe means that for every ORIGINAL polynomial marker Q there is a cutoff
beyond which every mixed moment is zero. The cutoff may depend on Q.
For Image moments, every output coefficient must vanish. The zero marker
is harmless; its degree can be omitted from degree assertions.

**Theorem BoxedUniversalMathieuDecision.** There is a terminating theoretical
algorithm deciding whether EVERY complex pure phase in the input
constructible family is safe. It has the following alternative outputs:

1. A generic counterexample over an explicitly presented finitely generated
   characteristic-zero field K, an ORIGINAL monomial marker, and an
   existence proof of a faithful embedding K into C realizing a complex
   pure unsafe phase in the restriction. The marker has infinitely many
   nonzero mixed moments after that ONE embedding.
2. A nonnegative rational number C_family such that EVERY pure phase in
   the restriction, over EVERY characteristic-zero extension of the coefficient base, obeys

       mixed_moment(Q,F)=0 for every integer m>C_family*deg(Q)             (B2)

   for EVERY original polynomial marker Q over that field.

In particular the algorithm decides complex existence of a pure unsafe
phase in any fixed degree box with any specified constructible coefficient
restriction. The affirmative-validity output is stronger than a pointwise
safety assertion: it supplies a uniform linear marker cutoff for that box
or restriction. The counterexample output is a symbolic generic field model,
not an algorithm producing rational or algebraic complex coefficients.
No statement that the whole safety locus is constructible is needed.

## 2. Precise retained algorithmic premises

Let R=Q[c] and write the original pure moments and mixed moments coefficientwise
in R. Gaussian moments are scalars; Image moments have finitely many output
coefficients at each order. They are rational polynomials in c by literal
Wick contraction or differentiation. There are no coefficient denominators
in these ORIGINAL moments.

The following retained results supply all non-elementary algorithmic steps.

1. GAUSSIAN_COMPUTABLE_UNIFORM_PURITY_PREFIX.md proves a computable finite
   pure prefix B_G(N,D) and, through the shifted Gaussian identity, a
   computable B_Image(n,D)=B_G(2n,D) for arbitrary Image phases. Hence one
   can compute a FINITELY GENERATED coefficient ideal I whose vanishing,
   over every characteristic-zero field, is exactly purity. For Image this
   ideal contains all output coefficients of every pure moment in that
   finite prefix. The theorem asserts radical finite-prefix stabilization,
   not equality of ordinary infinite moment ideals.
2. GAUSSIAN_FORMAL_COMPLETION_AND_POLYNOMIAL_HOM.md and
   WHOLE_GAUSSIAN_SOURCE_HOM_AND_EFFECTIVE_CUTOFFS.md give a terminating
   all-original-marker safety decision over effective characteristic-zero
   fields, using the exact WHOLE source and its polynomial Hom. At a pure
   safe phase they compute a sufficient rational pure-linear constant C.
   A zero Hom alone is not a purity test; purity is imposed separately by I.
3. SHIFTED_GAUSSIAN_IMAGE_EQUIVALENCE.md transfers those field safety
   decisions, original monomial marker witnesses, and pure-linear cutoffs
   to unrestricted Image inputs. The inverse affine shift preserves marker
   degree, and the finite original-monomial test selects an original bad
   marker when safety fails.

The required generic fields are K=Frac(R/P) for computed Q-primes P.
Arithmetic, equality, Weyl Groebner and syzygy operations are effective in
these presentations. Scalar integer roots are also effective: clear field
and parameter denominators, reduce polynomial representatives modulo P,
and express a nonzero coefficient polynomial h(M) in the finitely many
standard parameter monomials actually occurring. Select one nonzero rational
coordinate h_b(M). Every integer root of h is among the finitely many integer
roots of h_b; test those candidates in K. This is the generic-field step
proved in the retained pure-prefix argument. The theoretical source/Hom
algorithms apply to these computable subfields of C. The executed older
producer interface over Q and simple number fields is not being identified
with an implemented generic-field backend.

## 3. The specialization lemma, with the quantifiers retained

**Lemma GenericSafetySpecializes.** Let P be a prime coefficient ideal
containing I, and let F_eta be the phase over K=Frac(R/P). If F_eta is safe,
every point of V(P) in EVERY characteristic-zero extension of the coefficient
base is safe. If a
pure-linear marker constant C works for F_eta, exactly the same C works
for every such point.

Fix one original monomial marker Q. Generic safety says that each original
mixed moment has every output coefficient zero in K for all m beyond a
cutoff for Q. These coefficients belong to R, so they lie in P. They
therefore vanish at every coefficient specialization whose kernel contains
P. Repeat this argument for each fixed original monomial. This proves
safety at each point. There is no exchange of infinitely many cutoffs:
for an arbitrary polynomial marker, take its finite monomial expansion and
the maximum of its finitely many cutoffs.

If the generic estimate is m>C deg(Q), the same argument gives that exact
bound for each monomial. For a polynomial marker all its monomial degrees
are at most deg(Q), and C>=0, so the same bound handles their entire sum.
This proves the uniform assertion. The generic source computation may use
coefficient denominators which vanish at some specialization. We never
specialize that computation: we specialize the ORIGINAL polynomial moment
identities already known to belong to P. Thus no nonvanishing denominator
guard is needed in this lemma.

**Lemma GenericUnsafetyGivesAComplexPoint.** If the generic phase over
K=Frac(R/P) is unsafe, it supplies a complex pure unsafe point on V(P).
If g is not in P the point also lies in D(g).

The retained field algorithm can select one original monomial marker with
nonpolynomial period, hence infinitely many nonzero original mixed moments.
Alternatively, a polynomial unsafe marker has a finite monomial expansion;
if every monomial tail were zero, the whole marker tail would be zero.
Thus some original monomial is bad. In the Image case infinitely many mixed
moments are nonzero OUTPUT POLYNOMIALS, rather than merely a selected scalar
coefficient.

Every finitely generated characteristic-zero field over Q embeds faithfully
in C: choose an algebraically independent complex tuple for a transcendence
basis, then extend the resulting embedding through the finite algebraic
extension into algebraically closed C. The embedding fixes a specified
embedded number field when that is the base. Its restriction to R gives
a complex point with kernel EXACTLY P. All pure moments vanish there.
Every nonzero coefficient of every nonzero mixed output remains nonzero
under the SAME injective field embedding, so infinitely many mixed moments
survive. If g is not in P, its image remains nonzero as well.

Generic unsafety is NOT asserted to specialize to every point of V(P).
For example lambda times a fixed pure unsafe phase is generically unsafe,
but its specialization lambda=0 is safe. This prevents an invalid reversal
of GenericSafetySpecializes.

## 4. The finite component algorithm and uniform constant

For each piece (B1), compute the minimal Q-primes P of I+J. Discard precisely
those primes containing g. If I+J is the unit ideal, there are no primes
and that piece is empty. No absolute prime decomposition is necessary.
For every retained P run the field all-marker safety algorithm on F_eta.

If any generic phase is unsafe, the preceding lemma gives output 1 of the
theorem, including the declared inequalities. Otherwise compute its retained
pure-linear constant C_P for each generic safe phase and output

    C_family=max({0} union {C_P over every retained component}).          (B3)

All operations are finite. The prefix computation and each field decision
terminate by the retained theorems; prime decomposition is a terminating
exact polynomial algorithm. This proves effective termination even though
the cost can be enormous and no numerical bound is supplied here.

For completeness let c be ANY allowed pure coefficient tuple in a
characteristic-zero extension of the coefficient base, and let p be the
kernel of evaluation R->K_c.
It contains I+J for some input piece and does not contain g. It therefore
contains a minimal prime P of I+J, and this P cannot contain g. The retained
component was tested. GenericSafetySpecializes proves (B2), with C_P<=C_family.
Every complex point is covered, and so is every point over every other
characteristic-zero extension of the coefficient base. Conversely every unsafe generic component gives
an actual complex restricted point by the embedding lemma. This proves
both directions of the decision theorem, including empty and reducible
restrictions.

A Q-prime need not remain prime over C. This causes no gap. The polynomial
identities in P specialize to ALL of its complex points and geometric
components. A faithful embedding of its fraction field realizes one complex
generic point and preserves every nonzero element. For instance a^2+b^2 is
prime over Q but splits over C; the argument does not assume it stays prime.

The output does not describe where safety holds on an unsafe component.
It decides the universal assertion and counterexample existence on the
specified constructible restriction. It does not decide rational points,
Qbar points, least cutoffs, or an unrestricted unbounded-degree problem.

## 5. Radial-primitive and proper-power-free balanced restrictions

Fix n and a positive balanced degree d, and let V_d be the coefficient
space of all polynomials bihomogeneous of bidegree (d,d) in (xi,s). Put
R_0=sum_i xi_i s_i. The forms divisible by R_0 are the image of the fixed
linear injection V_(d-1)->V_d, H->R_0 H. Their coefficient space is a
computable linear closed subspace, so radial primitivity is a constructible
open guard. This definition includes its absolute meaning after coefficient
field extension, since solving the fixed linear divisibility equations
commutes with field extension.

**Lemma PolynomialPowerImageClosed.** For e,k>=1, the map

    V_e -> V_(ke), H -> H^k

has closed image on complex coefficient points, and its image ideal is
computable by polynomial elimination over Q.

Here is an elementary finite-integral proof. Choose finitely many rational
evaluation points p_j for which evaluation is an invertible linear map
on V_e. Such points can be found effectively by selecting independent rows
from the finite grid {0,...,e}^(2n): a polynomial of individual degrees at
most e vanishing on that grid is zero, by repeated univariate interpolation.
Write B=Q[h_1,...,h_M] for the input coefficient ring and A for its subring
generated by all coefficients of the universal H^k. For each selected point,
L_j=H(p_j) satisfies the monic equation

    L_j^k-(H^k)(p_j)=0 with (H^k)(p_j) in A.             (B4)

Every h_i is a rational linear combination of the L_j. Thus B=A[L_1,...,L_M]
is finite integral over A. The kernel of Q[c_output]->B is computed by
eliminating h from c_output-coefficients(H^k). After base change to C this
remains a finite integral extension onto its image ring. Lying-over gives
a maximal ideal upstairs over every maximal ideal downstairs; its residue
field is finite over C and therefore C. Hence every complex point of the
elimination variety is an actual kth power. This proves closedness, including
zero, without treating an affine polynomial image as automatically closed.

If a nonzero bihomogeneous F of degree (d,d) is a scalar proper power
c H^k, the maximum and minimum xi-degrees in H^k are k times those in H;
thus H is xi-homogeneous, and its degree divides d by k. The identical
argument in s shows that H has bidegree (d/k,d/k). A prime divisor p of k
already exhibits F as a pth power after absorbing a scalar root over C.
Consequently the scalar proper-power locus is the FINITE union of the
closed images V_(d/p)->V_d, H->H^p, over primes p dividing d. Its complement
is an effectively specified open guard. At a prime balanced degree d,
only the dth-power image of bilinear forms is needed.

Apply Theorem BoxedUniversalMathieuDecision to this degree box, with the
radial-primitive and proper-power-free guards. It decides complex existence
of a pure unsafe radial-primitive nonpower phase in the box. In the unsafe
branch the generic field phase is absolutely radial-primitive and nonpower,
because these closed geometric bad loci have been excluded. In the safe
branch (B2) holds for every restricted pure phase.

Equivalently, if one works directly with minimal pure components, discard
those wholly contained in the union of bad closed loci. Primality and
irreducibility show that a component contained in a finite union of closed
sets is contained in one of them. Remaining components meet the common
open guard, so their generic phases are legitimate tests. The general
principal-open decomposition in Section 1 gives a mechanical alternative.

For n=2 and d=5 the phase coefficient space has dimension36, the R_0-divisible
linear subspace has dimension25, and the scalar fifth-power cone is the
image of the four-dimensional bilinear coefficient space. The algorithm
therefore applies to the primitive quintic question without asserting its
answer or supplying a feasible computation. Nothing here replaces the
all-order witness obligation by a bounded experimental scan.

## 6. Finite controls and assurance

The portable exact checker verifies specialization bookkeeping on explicitly
proved families, including a reducible parameter restriction, empty guards,
and the generically unsafe family lambda*F_4 whose zero member is safe. It
checks a Q-prime which splits geometrically, original marker cutoffs through
the isotropic Gaussian identity, finite rational interpolation matrices and
the integral equations (B4), and the dimensions and prime-divisor guard logic.
The retained two-pair F_4 supplies the all-order pure unsafe antecedent;
finite contractions are supplementary controls only. No general minimal-
prime/Hom loop is executed by this checker.

Exact support: [source ledger](boxed_universal_mathieu_decision_source_ledger.json), [finite inputs](boxed_universal_mathieu_decision_inputs.json), [checker](verify_boxed_universal_mathieu_decision.py), and [receipt](boxed_universal_mathieu_decision_verification.json). All mathematical reviews and replays recorded here are internal.
