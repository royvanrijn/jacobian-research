# Full marked-safe Haar families and exact torsion normalization

16 November 2026. Internal written all-degree/all-ring deductions and an
implemented family adapter. Independent replay, formal verification and
external review remain false. Both FULL exceptional MZ cases remain
**OPEN**. Review is deferred and is not the next task.

For each full bounded coefficient space on either exceptional quotient,
ALL-marker eventual safety, WITHOUT a pure premise, is controlled by one
closed field-safe ideal. Over ANY coefficient ring, safety is exactly
nilpotence of that ideal's image. The exact normalization ideal of the
ORIGINAL Haar exponential source is its coefficient torsion, Gamma_I(R).
The full ordinary affine cutoff topology is its formal completion.

The [family adapter](exceptional_haar_family_producer.py) now transfers
the separately implemented Gaussian family solver to polynomial Haar
families of arbitrary degree and support. It keeps the ORDINARY parameter
ideal, tests invariance modulo that ideal, and returns the original safe
ideal and exact normalization ideal. This family route inherits trusted
CAS holonomic operations; it is distinct from the
[pointwise original-operator certificates](EFFECTIVE_POINTWISE_EXCEPTIONAL_HAAR_CERTIFICATES.md).

See the [checker](verify_exceptional_haar_family.py),
[receipt](haar_family_verification.json), and
[source audit](haar_family_source_audit.json).

## 1. Full safety transfers before imposing any pure equations

Let A_Gamma=C[SL2/Gamma], Gamma=2O or2I. Fix an even matrix filtration
bound D and a universal phase P_c in its ENTIRE finite-dimensional
filtered space. The retained
[radial bridge](GAUSSIAN_RADIAL_TRANSFER_AND_EXCEPTIONAL_UNIFORMITY.md)
identifies it linearly with a homogeneous degree-D invariant polynomial
Ptilde_c in four quaternion variables. For every homogeneous invariant
marker Qtilde of even degree e and EVERY coefficient ring R,

    E_G[Qtilde*Ptilde_c^m]
       =2^((Dm+e)/2)*((Dm+e)/2+1)! *h_R(Q*P_c^m).     (1)

All displayed scalar factors are rational units. Gaussian markers can
be Reynolds-averaged over Gamma; homogeneity and every marked moment
are preserved, and odd markers have zero moments. Conversely every
original invariant marker has a lift at its ORIGINAL filtration degree.
Thus full eventual marked safety is equivalent in the two models,
WITHOUT assuming any pure moments vanish. The proof is polynomial in
phase and marker coefficients, so it applies to nonreduced and
non-Noetherian R, not only to field points.

Let J_G be the universal full Gaussian field-safe radical ideal from the
retained [Gaussian mixed-completion theorem](GENERAL_GAUSSIAN_FAMILY_LOCUS_AND_NORMALIZATION.md).
For the linear lift ell define

    J_H=sqrt(ell^*J_G).                                (2)

This radical is taken in the regular finite coefficient polynomial
ring, not asserted to commute with arbitrary base change. It defines
the full Haar field-safe locus. For an R-valued phase let I_H=J_H(P)R.
The Gaussian nilpotence criterion and finite generation give, in ALL rings,

    P safe against EVERY original marker
      iff I_H NILPOTENT
      iff P safe over R_red
      iff every prime-residue phase is safe.           (3)

For clarity, replacing a finitely generated ideal by its radical before
evaluation preserves this nilpotence test: some fixed power of each
radical generator belongs to the original ideal; finite generation then
makes the radical image nilpotent whenever the original image is.
No ordinary ideal equality or radical base-change identity is inferred.

There are common field-safe affine constants C,E at fixed D. If I_H^L=0,
one sufficient all-marker bound inherited from the Gaussian theorem is

    m>C*deg_matrix(Q)+E+(CD+1)*(qL-1),
    q=binomial(D+4,4).                                 (4)

It uses the four-variable regular universal Gaussian coefficient space.
The same Gaussian constants work after the linear lift. This is an
existence/transfer statement for C,E, not a computation of numerical
universal constants by the new family adapter. A supplied pointwise
certificate can still give numerical constants for that supplied phase.

## 2. Ordinary affine cutoff ideals and formal completion

Before imposing ANY pure equations, define the Haar coefficient ideals

    K_N^H=( h(Q*P_c^m): m>N*(deg_matrix(Q)+1) ), N>=0.  (5)

Include Q=1. All original markers are included. The same notation K_N^G
uses all Gaussian polynomial markers and their usual degrees. Then

    K_N^H=ell^*(K_N^G) as ORDINARY ideals.              (6)

For one containment, lift an original marker at its minimum filtration
degree and use (1). For the other, average a homogeneous Gaussian marker
of degree e. Its restriction is an original marker of degree at most e;
m>N(e+1) is therefore also a forbidden original index. Divide by the
nonzero scalar in (1). Linear combinations of homogeneous markers prove
both containments. The group and quaternion coefficient field can be
adjoined faithfully flat when needed, then the identities contract.
For a presented coefficient family the equality is in its ORDINARY
quotient ring, where the phase is actually invariant.

The retained Gaussian cofinal containments imply that K_N^H are cofinal
with powers of J_H. If the literal pulled-back Gaussian ideal is not
radical, its powers and J_H powers are cofinal by finite generation.
Consequently the full marked-safe functor on arbitrary discrete
coefficient algebras is the formal completion along V(J_H). There is
NO ordinary stabilization claim and NO thickness-independent bound.

The constant-phase family P=u is a useful exact boundary. Q=1 gives

    K_N^H=(u^(N+1)) for integer N>=0.                  (7)

Other generators are multiples of this one. These powers do not stabilize
in C[u]; over a thickening they retain its exact ordinary nilpotents.
The older zero-intercept pure-quotient ideals are a different object.

## 3. A direct formal-Hom proof for the ORIGINAL Haar source

This section identifies the normalization ideal on the original quotient,
not merely its Gaussian analogue. Let R be ANY C-algebra, P in R tensor
A_Gamma. Choose the three left Lie-algebra vector fields D_i, divergence-
free for the original invariant volume. The full amplitude source is

    H_P^H=R[t] tensor A_Gamma /
           sum_i (D_i+t*D_i(P)) (R[t] tensor A_Gamma),
    partial_t acts on amplitudes as partial_t+P.        (8)

The operators commute with this time connection and give the actual
full top integration source. Every original regular amplitude is kept;
this is not the cyclic module generated by the unit's scalar operator.

A differential map to R[[t]] has values F_Q satisfying
partial_t F_Q=F_(PQ). Hence for its initial functional lambda(Q)=F_Q(0),

    F_Q=sum_m lambda(Q P^m)*t^m/m!.

The spatial relations at t=0 give lambda(D_i Q)=0. The left
representation A_Gamma is a direct sum of finite-dimensional irreducible
SL2 modules; every nontrivial summand equals its Lie-algebra span, while
the trivial summand is C*1 because the left action is transitive.
Thus its coinvariants are exactly C, with quotient the normalized h.
Tensoring this explicit decomposition with R preserves it. Therefore
lambda=a*h_R for a UNIQUE a in R. Conversely h_R-invariance checks all
relations in (8), coefficient by coefficient. We obtain, for EVERY R,

    Hom_Dt(H_P^H,R[[t]])=R*Phi_P,
    Phi_P(Q)=h_R(Q exp(tP)).                           (9)

This is an explicit initial-functional argument, not relative-Hom base
change. Formal series are sufficient; no convergence of a Gaussian
exponential is assumed. The ORIGINAL compact Haar series are entire
on complex field fibres, but that fact is not required for (9).

Let N_P^H be the ideal of scalars a for which a*Phi_P takes every
amplitude to R[t]. Reduction modulo Ann_R(a) says exactly that P is
marked-safe over R/Ann_R(a). Applying (3) gives the exact formula

    N_P^H=Gamma_(I_H)(R)
         =union_L Ann_R(I_H^L)=H^0_(I_H)(R).            (10)

All original global relations are kept via (8)-(9); the actual period
image kernel need not be computed. Safety means 1 belongs to this ideal,
not merely that the ideal is nonzero. For a finitely presented family
R=B/I0, its lift is the ORDINARY saturation

    (I0 : J_family^infinity) in B.                     (11)

Nilpotent directions of I0 must remain. For example R=C[u,v]/(u²v²),
P=u has I=(u) and normalization ideal (v²), nonzero but not containing1.
Over C[u]/u³, P=u is safe with normalization R, although its first
positive moments are nonzero. These are scalar-phase controls, not
new exceptional pure failures.

## 4. General family adapter and its exact input contract

Input supplies `group`, a positive `parameter_count`, sparse rational
`phase_terms` whose exponent lists contain four spatial coordinates
followed by the parameters, and an optional ORDINARY `parameter_ideal`.
The phase represents a polynomial family on r²=1. The adapter reduces
spatially modulo r²-1 and coefficientwise modulo I0, without taking a
radical. It rejects nonzero odd parts and checks the even lift under
exact group generators MODULO I0. For example u*q0² over R=C[u]/u²
is not invariant merely because its reductions at field points vanish;
the adapter correctly rejects it.

Algebraic coefficient fields can be encoded by parameter equations. When
the SAME parameter also defines the group-field square root, specify
`group_root_parameter`. Its square must be2 for2O or5 for2I in the
ORDINARY coefficient algebra, not just after reduction. The group action
then uses that very parameter. This covers families over the natural
quadratic group field and both of its conjugate models consistently.
Otherwise the quadratic group field is used only for the exact invariance
check of the rationally presented family. No algebraic disjointness is
assumed. Fixing a complex embedding gives the usual original quotient.

The homogeneous degree D is fixed for the WHOLE family. It need not equal
the smaller degree at a special fibre. Equation(1) remains valid with
this fixed D even on such fibres. Inactive spatial Gaussian variables
can be removed using the independent Wick factorization proved in the
pointwise adapter. All coefficients and all parameter derivatives are
then sent to the UNCHANGED Gaussian joint-family backend.

`universal_family(group,d)` constructs a surjective linear Reynolds
parametrization of the ENTIRE original filtered phase space, using
all homogeneous monomials of degree2*floor(d/2). The last parameter is
the quadratic group-field root. Parameters are deliberately redundant;
no independence of this spanning set is claimed. Every invariant
homogeneous form is its own Reynolds average, so this parametrization
omits NO phase at the requested degree. The constructor prepares inputs
only; it does not launch a degree-by-degree calculation.

On completion the adapter returns the original field-safe radical ideal,
the exact ordinary polynomial-Hom normalization ideal, its annihilation
exponent and colon-stability check, and coefficient-ring safety. It does
NOT test the pure locus or report an MZ family decision. For finitely
presented inputs all values are obtained from the retained Gaussian
family operations and (1)-(11), with no new Gaussian algorithm.

## 5. Assurance, interface and executed controls

The family backend computes joint holonomic operations, ambient local
cohomology and conormal MODULE multiplicities in Macaulay2. Those are
trusted CAS outputs. Unlike pointwise IBP/jet certificates, the family
output is NOT an independent original-operator proof. The adapter keeps
that boundary explicitly. For a reused result it checks the exact phase,
parameter count, precision and ORDINARY ideal, core hash, generated script,
stdout and artifact hashes; reported ideals must match that stdout.
It does not promote trusted CAS functors to formal assurance.

    python exceptional_haar_family_producer.py input.json result.json
    python exceptional_haar_family_producer.py input.json result.json --prepare-only
    python exceptional_haar_family_producer.py input.json result.json \
      --gaussian-family-certificate retained-result.json --backend-workdir retained-cas

Default cap30seconds gives UNKNOWN; `--timeout 0` disables the cap and
inherits the backend's theoretical termination assumptions. Errors fail
closed. Numerical universal cutoffs, least cutoffs and individual ordinary
cutoff-ideal producers are not implemented by this adapter.

The checker transfers FOUR already-completed family results (zero and
scalar-parameter cusp, on both groups), checks eight ENTIRE filtered-
space constructor inputs at small bounds, and prepares a nonconstant
algebraic octahedral degree8 family WITHOUT a holonomic family run.
It rejects nilpotent-only reduced invariance, a group-root relation true
only after radicalization, altered family ideals and provenance, and
other malformed inputs. The CLI retained-result path passes. There were
ZERO new holonomic family CAS constructions and no new nonconstant
spatial family safety-locus computations. The first preparation attempt
had an int0/SymPy0 implementation error; it was corrected before any
backend invocation and produced no mathematical decision.

## 6. The FULL exceptional target remains

This closes the executable family-safe-locus/ordinary-normalization
transfer and the full all-ring Haar normalization formula, WITHOUT a
pure premise. It does not show that every pure phase belongs to that
safe locus. At each D the positive target is still inclusion of the
ENTIRE pure locus in V(J_H,D), equivalently the established intrinsic
unit-extension criterion. A failure requires an actual original pure
phase outside this locus with a fixed original bad marker.

Writing I_pure,D for the ordinary ideal of ALL positive pure moments,
the fixed-degree universal inclusion is exactly

    J_H,D subset sqrt(I_pure,D).                       (12)

Here I_pure,D uses ALL positive orders; we do not identify it with a
closed eventually-pure locus at the same D. Across ALL degrees the tests
are equivalent: for an eventually-pure P choose N with P^N all-positive
pure, then use the finitely many original markers Q*P^j, 0<=j<N.

All degrees D are required for full MZ; neither ideal nor this containment
has been computed on the entire remaining exceptional locus here. The
nilpotent criterion shows that a full complex-field proof would already
extend to every coefficient C-algebra; nonreduced rings are not a separate
full obstruction. This credits the retained radial nilpotent bridge and
makes its ordinary normalization explicit.

No finite list of safe degrees proves the unbounded statement. The
four/six all-ring theorem, five/seven failure-support lower bounds,
all-support capacity inequality and full actual marked-source relation
requirements remain in force. No new coefficient/degree/phase campaign
or Gaussian-backend duplication was introduced. Both FULL MZ cases remain
OPEN; direct exceptional ownership remains, and review stays deferred.
