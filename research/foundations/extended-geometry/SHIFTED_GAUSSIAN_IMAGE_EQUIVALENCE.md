# Shifted Gaussian recognition of unrestricted Image safety

20 November 2026. Internal written proof and exact finite controls. The
classical Gaussian integration and polynomial-Hom algorithms are imported
with their stated exact-source and effective-field hypotheses. No new
implementation over a transcendental coefficient field is asserted.

## 1. The shifted contraction identity

Let K be a characteristic-zero field, n>=1, and
B=K[xi_1,...,xi_n,s_1,...,s_n]. Define the K-linear map

    E_n(xi^alpha s^beta)=partial_s^alpha s^beta in K[s],
    M_n=sum_i (partial_(s_i)-xi_i) B.

Then M_n=ker E_n. Indeed E_n(xi_i H)=partial_(s_i)E_n(H)
=E_n(partial_(s_i)H). Conversely xi_i H is congruent to
partial_(s_i)H modulo M_n. Induction on xi-degree gives
H-E_n(H) in M_n. This is a vector-space reduction, not an assertion
that M_n is an ideal.

Adjoin one square root i of -1 if necessary, take algebraically independent
parameters z_1,...,z_n, and put L=K(i)(z_1,...,z_n). Let Gamma be the
normalized standard Gaussian Wick functional on L[X_1,...,X_n,Y_1,...,Y_n]:
independent variables have zero odd moments and
Gamma(X_j^(2r))=Gamma(Y_j^(2r))=(2r-1)!!.
Set

    V_j=(X_j-iY_j)/2,   W_j=X_j+iY_j,
    S(H)=H(V,z+W).

Here the s_j are the original coordinate variables; the z_j are coefficient
parameters. The map S:L[xi,s] -> L[X,Y] is an affine L-algebra isomorphism,
with inverse substitutions

    X_j=xi_j+(s_j-z_j)/2,
    Y_j=((s_j-z_j)-2xi_j)/(2i).                         (S1)

**Theorem ShiftedGaussianImageIdentity.** For every polynomial H in B,
without a homogeneity or balance assumption,

    Gamma(S(H))=E_n(H)(z) in K[z] subset L.             (S2)

To check it, the formal exponential generating function of (V,W) is
Gamma exp(a.V+b.W)=exp(sum a_j b_j). Thus
Gamma(V^alpha W^gamma) is zero unless alpha=gamma and then is alpha!.
Expanding (z+W)^beta yields zero when some alpha_j>beta_j and otherwise
beta!/(beta-alpha)! z^(beta-alpha), exactly the displayed definition of E_n.
The generating-function calculation is coefficientwise polynomial Wick
arithmetic, so it uses no analytic convergence. Equivalently it follows
from the covariances Gamma(V_jV_k)=Gamma(W_jW_k)=0 and
Gamma(V_jW_k)=delta_jk. The identity also holds for H with coefficients in L
when z is treated as a coefficient parameter.

## 2. Pure phases and arbitrary fixed multipliers

For F in B put q_F=S(F). A phase is called Image-pure if E_n(F^m)=0 for
every m>=1. Gaussian-pure means Gamma(q_F^m)=0 for every m>=1.

**Corollary ShiftedGaussianPurityEquivalence.** F is Image-pure if and only
if q_F is Gaussian-pure over L. More generally, for every fixed Q in B and
every m>=0,

    E_n(Q F^m)(z)=Gamma(S(Q) q_F^m).                    (S3)

In particular the original Q tail vanishes if and only if the Gaussian
S(Q) tail vanishes. The assertion includes multipliers of mixed or unequal
bidegrees and preserves cancellations between all their terms. Faithfulness
of K[s] -> L, s_j -> z_j, proves both zero implications. Consequently F^m
belongs to M_n if and only if its corresponding Gaussian moment is zero;
the same holds for QF^m. No comparison of merely constant contractions is
being used.

## 3. All-multiplier safety, in both directions

Call F Image-safe if for every Q in B there is an integer M_Q such that
E_n(QF^m)=0 for every m>=M_Q. Define Gaussian safety of q_F by the identical
condition for every A in L[X,Y], with Gamma in place of E_n.

**Theorem ShiftedGaussianAllMarkerSafetyEquivalence.**

    F is Image-safe over K
      iff q_F is Gaussian-safe over L.                (S4)

The reverse implication follows by choosing A=S(Q) and using (S3).
For the forward implication take any A in L[X,Y]. The inverse (S1) expresses
it as S(Q_A) with Q_A in L[xi,s]. Write the finite sum

    Q_A=sum_(alpha,beta) c_(alpha,beta)(z) xi^alpha s^beta,
    c_(alpha,beta)(z) in L.

Each monomial xi^alpha s^beta is an allowed original marker over K, so its
original tail is zero beyond its own cutoff. The maximum of the finitely
many cutoffs annihilates their L-linear combination, hence annihilates
Gamma(A q_F^m). This argument permits rational functions of z and the
adjoined i as coefficients; neither a common cutoff for infinitely many
markers nor specialization of a rational parameter is required.

The same finite-combination argument proves stability and reflection of
these all-marker assertions under arbitrary coefficient-field extensions.
In particular K(i) causes no additional field-descent obligation. The
universal quantifier on the Gaussian side really covers all polynomial
Gaussian markers, because S is the full affine polynomial-ring isomorphism.
These equivalences do not assume that F is pure. A pure, unsafe F remains
pure and unsafe after this shifted Gaussian transformation.

The independent output parameters are essential. Setting them all to zero would send F=s_1, whose contractions s_1^m never vanish and which is Image-unsafe already for marker1, to W_1. That isotropic Gaussian linear phase is pure and safe for every polynomial marker: each factor W_1 must pair with a V_1 from the marker, so its mixed moments vanish once m exceeds marker degree. Thus scalar contraction at zero cannot replace the faithful generic shift.

## 4. Theoretical pointwise decision and original cutoff transfer

The following are existing proof inputs, not new Gaussian algorithms:

* GAUSSIAN-SAFETY-HOM, proved in
  research/moments/vanishing/GAUSSIAN_FORMAL_COMPLETION_AND_POLYNOMIAL_HOM.md
  Sections 1 and 3: for a polynomial phase q over a characteristic-zero
  field, the exact full integration source H_q satisfies
  Hom_D(H_q,K[[t]])=K Phi_q and all-marker safety is equivalent to
  Hom_D(H_q,K[t])!=0. The whole source, with all amplitude maps, is essential.
* GAUSSIAN-EFFECTIVE-FIELD-POINTWISE-SAFETY-DECISION, proved in
  research/moments/vanishing/WHOLE_GAUSSIAN_SOURCE_HOM_AND_EFFECTIVE_CUTOFFS.md
  Sections 1--5: exact full integration, polynomial Hom and original-amplitude
  certificates give a terminating theoretical decision over effective
  characteristic-zero fields with exact integration and scalar integer-root
  operations. Sections 4--5 also yield computable rational constants C,E>=0
  such that every Gaussian marker A has zero tail for m>C deg_(X,Y)(A)+E
  when q is safe. In the pure safe case a possibly enlarged C gives the
  sufficient bound m>C deg_(X,Y)(A).

**Corollary UnrestrictedImageSafetyDecision.** For F over Q or a specified
number field, Image all-marker safety is theoretically decidable by applying
those retained theorems to q_F over L. There is no restriction on n, degree,
homogeneity or balance. This is a pointwise decision for a specified F,
not a claim that every pure F is safe or a classification of safe forms.

Here is why the effective-field hypothesis is met in this theoretical use.
A finite algebraic extension K(i) is effective, and rational functions in
finitely many indeterminates over it have exact field arithmetic, equality,
polynomial factorization and Groebner algorithms. The classical integration
and polynomial-solution algorithms operate over this computable field.
For the scalar integer-root step, take a nonzero h(N) in L[N], clear
parameter denominators, and express its coefficients in a fixed Q-basis of
K(i), then in parameter monomials. This gives finitely many rational
coordinate polynomials h_j(N), at least one nonzero. Every nonnegative
integer root of h is a root of that selected nonzero h_j. Enumerate the
finite candidates by an exact rational root bound and test them in L.
This also supplies integer roots for the integration and indicial steps.
The effective field embeds in C after choosing algebraically independent
complex parameters; no numeric embedding is used by the arithmetic.

A zero polynomial Hom certifies unsafety, not failure of purity. Pure phases
with unsafe multipliers exist. Purity can be decided separately by the
retained Gaussian actual-unit-period annihilator and finite origin-selector
method; (S2) transfers that answer. Merely finding some polynomial solution
of a scalar annihilator is not a safety test.

**Corollary OriginalAffineMarkerCutoff.** When F is safe, the retained
Gaussian constants C,E give, for every original polynomial marker Q,

    E_n(QF^m)=0 whenever m>C deg_(xi,s)(Q)+E.             (S5)

Indeed A=S(Q) has Gaussian-variable degree at most the total degree of Q:
the translation parameters have degree zero over L and all substitutions
are affine linear. Apply the Gaussian estimate and (S3). At a pure safe
point the sufficient bound can be written m>C deg(Q), using the retained
pure-linear enlargement of C. These constants concern the fixed input F,
not a uniform numerical bound for coefficient families.

**Corollary FiniteOriginalMonomialMarkerTest.** For each specified F over
an effective field as above, a computable integer B_F>=0 suffices to test
all-marker Image safety using only the finitely many original monomial
markers xi^alpha s^beta of total degree at most B_F. In the unsafe case at
least one such original monomial marker has infinitely many nonzero mixed
moments. The assertion is a bound for this fixed F, not an explicit uniform
function of its degree.

To prove this, use the maximal integer root of the whole Gaussian
integration b-function, or the retained verified recurrence, to obtain
D_t-generating amplitudes X^gamma Y^delta with |gamma|+|delta|<=B_F for
H_(q_F). The exact full-source construction explicitly retains all these
amplitude maps. A negative root threshold would give a zero integration
module, impossible because its formal unit period has constant term 1;
using max(0,B_F) is harmless. If all the indicated original monomial tails
vanish, (S1) expresses every Gaussian generating amplitude as a finite
L-linear combination of their shifts. Its actual period is therefore a
polynomial in t. Every other source element is a polynomial differential-
operator combination of the generating amplitudes, so its period is
polynomial too. This proves safety by (S4). The converse is immediate.
If safety fails, some generating amplitude has a nonpolynomial period.
At least one original monomial in its finite inverse expansion must then
have a nonpolynomial period, and hence infinitely many nonzero coefficients.
Applying exact fixed-marker period tests to the finite original list can
select one. The test never treats an inverse marker with parameter-dependent
coefficients as an original K-marker: it uses the finite monomial expansion.

The retained producer's executed interface supports Q and simple number
fields. It does not implement the transcendental-parameter backend L used
here. The theorem establishes the theoretical decision and cutoff transfer;
it does not assert that the existing command accepts this new input field.
No claim of practical efficiency, formal verification, external review,
classification, or a decision of an unrestricted existential coefficient
problem is made.

## 5. Exact controls and evidence limits

The accompanying standard-library checker independently computes original
contraction by falling factorials and computes the shifted expression by
literal sparse substitution followed by ordinary Gaussian Wick expectation.
It includes arbitrary unbalanced monomials, mixed-degree phases and markers,
multiple pairs, the affine inverse, and countercontrols for an omitted
shift, a wrong factor of two and an incorrect complex sign.

The retained two-pair degree-four balanced witness from
research/foundations/extended-geometry/TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md
is included as a finite original-variable control. Its all-order purity and
unsafe marker Q=Z are imported from that proof's Theorem 1.1; finite checker
moments do not establish those tails. This supplies a concrete check that
(S4) accommodates pure unsafe phases and does not confuse purity with safety.
The checker is finite exact evidence for its declared fixtures; the proofs
above establish the unrestricted equivalences.

Proof inputs and replay: [primary/retained source ledger](shifted_gaussian_primary_source_ledger.json), [exact inputs](shifted_gaussian_image_equivalence_inputs.json), [checker](verify_shifted_gaussian_image_equivalence.py), [receipt](shifted_gaussian_image_equivalence_verification.json). Primary sources actually inspected: [Oaku whole integration](sources/oaku_holonomic_integration_2013.pdf), PDFpp14--15, and [Tsai--Walther](sources/tsai_walther_hom_0007139v1.pdf), Introduction and Algorithm3.3. All audits and replays are internal; no external or formal verification is asserted.
