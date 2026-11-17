# Algebraic resolvent splitting and the zero-value logarithmic obstruction

16 November 2026. Internal written proofs; no novelty, independent replay,
external review or formal verification is claimed. See the
[checker](verify_algebraic_resolvent_splitting.py),
[receipt](algebraic_resolvent_splitting_verification.json) and
[source audit](algebraic_resolvent_splitting_source_audit.json).

For EVERY pure Gaussian phase, in every dimension and degree, all-marker
safety is equivalent to algebraicity of ALL its ORDINARY marked resolvents.
The same statement holds in the full intrinsic Haar setting, including both
exceptional quotients. Thus every pure failure requires a TRANSCENDENTAL
ordinary marked resolvent, sharpening the
[previous nonrational obstruction](POLYNOMIAL_RECOMPOSITION_REFLECTION_AND_UNIVERSAL_FAILURE.md).

For compact Haar, regularity gives a more local conclusion. Every pure failure
requires an integer-exponent logarithmic term at PHASE VALUE ZERO in some
continued original marked Cauchy period. Semisimplicity at eigenvalue1 of the
COMPLEMENT/Cauchy Gauss-Manin monodromy at zero is sufficient for full marked
safety. Monodromy of the generic fibre of P is a different object and cannot
be substituted. Neither hypothesis is proved for every exceptional pure phase;
both full MZ cases remain OPEN.

## 1. Exact Borel conversion and the time-zero source

Use the setting and notation of
[MVRECOMPREFLECT20261116](POLYNOMIAL_RECOMPOSITION_REFLECTION_AND_UNIVERSAL_FAILURE.md).
K has characteristic zero, h(1)=1, and the actual formal marked image

    V_P=sum_Q K[t] E_Q(t), E_Q=sum_n h(QP^n)t^n/n!,

has dim(V_P/tV_P)=1, also after algebraic closure. The retained exact Gaussian
and intrinsic Haar sources establish this for ALL their phases. Every nonzero
K[t]-submodule of K[[t]] has a nonzero quotient by t, so V_P is indecomposable.
If P is pure, V_P contains O=K[t]. No generic source-rank assertion is used.

Define the coefficientwise linear isomorphism

    B: K[[t]] -> W=K((s^-1))/K[s],
    B(sum_n a_n t^n/n!)=[sum_n a_n s^(-n-1)].                 (1)

Every class in W has a unique strictly proper Laurent representative at
infinity. Direct coefficient calculations give

    B(tF)=-partial_s B(F),
    B(partial_t F)=s B(F),                                 (2)

where the second equality is in W: before quotienting, its discrepancy is
the constant F(0). Thus(1) is a Fourier-compatible Weyl-module isomorphism,
with t=-partial_s and partial_t=s. It is NOT a claim about convergence or
an analytic Laplace integral. In particular

    B(O)=D_0=span_K{s^-1,s^-2,...} subset W.                (3)

D_0 is the point module at s=0. Keeping the polynomial quotient in W is
essential for multiplication by s to preserve the displayed module.

## 2. Normalized algebraic trace supplies a global retraction

Let F=K(s), and let L_alg be the field of elements of K((s^-1)) algebraic
over F. There is a canonical F-linear map

    T:L_alg -> F, T(f)=Tr_(L/F)(f)/[L:F],                   (4)

where L is ANY finite extension of F in L_alg containing f. This does not
depend on L: trace transitivity gives the same normalized value after any
larger finite extension, and common finite composita prove additivity.
It fixes F. It commutes with partial_s, since the derivation extends uniquely
to a finite separable extension and its embeddings commute with that
derivation; differentiate the sum of conjugates defining trace. No finite
generation of L_alg is needed.

For rational g(s), let PP_0(g) be its finite principal part at zero. Modulo
K[s], this projection commutes with s and partial_s: the analytic Taylor part
at zero contributes no negative powers under differentiation; multiplying
a principal part by s can add only a constant, killed in W. It kills K[s].
Consequently

    rho=PP_0 o T : L_alg/K[s] -> D_0                        (5)

is Weyl-linear and restricts to the identity on D_0.

**MVALGORDSAFE20261116.** For a pure P satisfying the exact time-zero source
hypothesis, the following are equivalent:

1. P is safe for every original marker.
2. Every ordinary C_Q(s)=sum_n h(QP^n)s^(-n-1) is algebraic over K(s).

Indeed safety makes each C_Q a finite Laurent polynomial. Conversely, under
the algebraicity premise B(V_P) lies in L_alg/K[s]: algebraic functions are
closed under differentiation, rational multiplication and finite sums.
Purity puts D_0 inside B(V_P). Restriction of(5) is a retraction onto D_0.
After(1), it splits off O from V_P. Indecomposability forces V_P=O, proving
all-original-marker safety. This uses the exact actual image, not an
unrelated algebraic solution of a larger scalar equation.

If a finite set of original marker images E_(Q_i) generates V_P as a Weyl
module, only those finitely many ordinary C_(Q_i) need be algebraic.
Such a finite set EXISTS for the retained holonomic Gaussian and intrinsic
Haar sources: choose finite source generators, expand their polynomial-time
amplitudes into finitely many original markers, and map to V_P. This is a
finite criterion with a supplied exact source, not an implemented producer
of that source or a new numerical prefix criterion.

No statement with EXPONENTIAL marked series in place of C_Q is asserted.
The retained factorial failure has a rational nonpolynomial exponential
marked series, whereas its ordinary factorial series is divergent and
transcendental. The normalization in(1) is substantive.

Purity is also essential. On the intrinsic Haar model G_m with h=CT_z and
P=z+z^-1, every Laurent marker has an algebraic ordinary resolvent:
for Q=z^k it is r(s)^abs(k)/sqrt(s^2-4), with
r(s)=(s-sqrt(s^2-4))/2 and the branch sqrt(s^2-4)~s at infinity.
Factoring the quadratic denominator in the Laurent resolvent proves this
formula. Yet h(P^(2n))=binomial(2n,n) is nonzero for every n. Thus the
algebraic criterion cannot replace the rational criterion of the previous
UNCONDITIONAL scalar-centre reflection theorem.

## 3. Transcendental witnesses survive every polynomial recomposition

**MVTRANSCENDBLOCK20261116.** Every pure unsafe P in either setting has some
fixed original marker Q with C_Q transcendental over K(s). For every
nonconstant H of degree d, at least one marker in

    Q, QP, ..., QP^(d-1)                                   (6)

has a TRANSCENDENTAL ordinary marked resolvent for H(P). In particular its
marked tail is infinite. The same Q is chosen once for P, independently of
H. Its degree block is bounded by degQ+(d-1)degP when a degree filtration
is available. The same single marker need not work alone.

Existence follows from section2 by contraposition. Use the exact retained
divided-difference reconstruction

    C_Q^P(s)=sum_(j<d) r_j(s) C_(QP^j)^(H(P))(H(s)),
    sum_j r_j(s)T^j=(H(s)-H(T))/(s-T).

If every term on the right were algebraic, the left would be algebraic,
a contradiction. Algebraicity is preserved by this polynomial substitution
and finite sums. This strengthens the previous finite bad-marker block,
without changing its bounds or claiming eventual nonzero values at every index.

The fixed transcendence also survives every coefficient field extension.
If a series over K satisfies an algebraic equation of bounded bidegree over
an extension, its finitely many unknown polynomial coefficients solve a
homogeneous linear system with entries in K. The infinite system has finite
row rank because it has finitely many columns, and a finite row basis over K
gives a nonzero solution already over K. Hence algebraicity descends.
This is not a finite-prefix test without a supplied relation-degree bound.

The proof of section2 also applies to an eventually pure P: its nonzero
polynomial unit period generates O under differentiation, so O still lies
in V_P. The same transcendental-witness conclusion then follows for an
eventually pure failure.

## 4. Compact Haar failures need a zero-value integer logarithm

Return to complex compact Haar in the intrinsic top-cohomology-one setting.
The C_Q are convergent at infinity and are the actual compact resolvents

    C_Q(s)=integral_gamma Q Omega/(s-P).

The retained [MVCOMPACTREGULAR20261116](COMPLEX_HOMOGENEOUS_GROWTH_AND_BOREL_TRICHOTOMY.md)
makes them periods in a regular-singular system, INCLUDING s=0. Use ONE
compatible continuation of the entire family from infinity to a sector at
zero. Locally its periods have finite power/logarithm expansions

    sum_(alpha mod Z) sum_(k=0)^M s^alpha (log s)^k
                                         f_(alpha,k)(s),   (7)

with meromorphic convergent f_(alpha,k) and finitely many exponent classes.
The classical regular-singular normal form was audited in the retained proof;
no new comparison theorem is imported here. Choose alpha=0 for the integer
class. Only this class matters for the following assertion.

**GHAARZEROLOG20261116.** A pure phase is safe iff ALL the continued original
C_Q have NO logarithmic terms in their INTEGER exponent class. Absence along
one compatible continuation suffices; for a failure every such choice has
an original Q with a nonzero integer-class logarithm.

For the converse, the same absence holds throughout B(V_P), by(2).
Discard all noninteger classes in(7), and take the finite negative-power
part of f_(0,0). This defines a map B(V_P)->D_0. It kills polynomials,
commutes with partial_s and s modulo polynomials, and is the identity on
D_0. Noninteger classes remain noninteger under these operations, even
when they contain logarithms. The integer class is log-free by hypothesis,
so its Taylor part creates no negative-power derivative. Thus this is a
retraction and the same indecomposability argument gives safety.
The forward direction is immediate from the finite Laurent-polynomial C_Q.

A sufficient geometric condition is therefore: the local monodromy at s=0
on the COMPLEMENT/Cauchy Gauss-Manin connection has no nontrivial unipotent
part on its generalized eigenvalue1 subspace. Its fundamental power/log
matrix then has no integer-class logarithms for any rational cohomology
section and compact-cycle period. Every actual full exceptional failure
must violate this condition and exhibit an actual original marked logarithm.
Semisimplicity on other eigenvalues is unnecessary. The condition on the
full complement source is sufficient, not necessary: the compact period
can kill source components.

This connection is the one for the smooth generic OPEN family

    U_s=X minus {P=s}, with forms Q Omega/(s-P),             (8)

not merely the middle cohomology of the generic fibre {P=s}, not the
generic connection of the ordinary top direct image, and not the exponential
source localized at t. Forgetting their point/extension data loses exactly
the obstruction being tested.

An exact countercontrol makes the distinction concrete. In the retained
entire nonsplit module take f(t)=(exp(t)-1)/t and include1. Its ordinary
resolvent is log(s/(s-1)), and

    (1-s)*partial_s B(f)=1/s.                               (9)

In the point quotient (1), 1/s generates the point module at zero.
Localizing that ordinary direct-image module away from the point kills
the right side, leaving a trivial generic connection. Nevertheless the
actual Cauchy function in(9) has a genuine log s and nontrivial unipotent
monodromy at zero. Thus generic-fibre/source semisimplicity does not give
the hypothesis for the COMPLEMENT/Cauchy system in(8).

## 5. Controls and the next full exceptional implication

The checker verifies the Borel/Fourier boundary constant, normalized traces
in two small separable extensions and a tower, principal-part commutation,
the failure of that commutation for log s, the harmlessness of a noninteger
logarithm, and(9). The universal trace/retraction proofs, global algebraicity
criterion and full local-log implication are written arguments, not
extrapolations of these finite checks. No general source generator, marked
algebraicity solver or exceptional monodromy computation is implemented.

The next complete exceptional gate is now more precise: exclude or realize
the actual eigenvalue1 logarithmic extension at phase value0 for EVERY pure
mixed-isotypic phase with image interior, strict leading nullcone and the
retained orbit/stabilizer conditions. A finite-branch algebraic replacement
cannot model a full failure. Proving the complement monodromy condition for
all such phases would close MZ, but no universal proof of it is claimed.
The earlier sharper factorial witnesses, recomposition/ring normalization,
complex homogeneous growth and all protected frontiers remain intact.
External review and outreach remain deferred under closure-first priority.
