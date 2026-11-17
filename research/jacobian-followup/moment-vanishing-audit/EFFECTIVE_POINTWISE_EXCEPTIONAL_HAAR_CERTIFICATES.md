# Effective pointwise certificates on both exceptional Haar quotients

16 November 2026. Implemented adapter with written all-degree transfer
proof. Independent replay, formal verification and external review remain
false. Both FULL exceptional MZ cases remain **OPEN**; review deferred.

The general Gaussian pointwise producer is now an executable input to the
ENTIRE exceptional Haar problem. The adapter accepts a supplied exact
phase on either fixed quaternion model, in ANY matrix-entry degree and
with arbitrary weight support. It decides eventual purity and all-original-
marker safety when the backend completes. A negative safety certificate
returns ONE fixed regular invariant original marker. Numerical safe
cutoffs transfer with the ORIGINAL quotient filtration.

This is a pointwise producer, not a proof that the intrinsic extension
vanishes on the full pure locus. No phase enumeration, new nonconstant
safety campaign or new Gaussian CAS construction was run. The Gaussian
backend is separately owned and is imported unchanged.

Implementation: [adapter](exceptional_haar_pointwise_producer.py),
[checker](verify_exceptional_haar_pointwise.py),
[receipt](haar_pointwise_verification.json),
[source audit](haar_pointwise_source_audit.json).

## 1. Exact input and canonical quotient lift

Use quaternion coordinates q=(a,b,c,d), identified with

    [[a+i*b,c+i*d],[-c+i*d,a-i*b]],
    Delta=a^2+b^2+c^2+d^2=r^2.

The original compact cycle is r^2=1; its complexification is SL2.
Input specifies group2O or2I and a sparse polynomial in q over QQ or
a specified simple number field. These are the explicit standard finite
quaternion groups in the retained exact group checker. Any algebraic
phase on these models is representable after a finite coefficient-field
extension. The root convention for a simple number field is exactly that
of the Gaussian backend (its specified CRootOf index0); no numerical
embedding arithmetic is used. Other algebraic conjugate models require
transport to this fixed representation, not a silent change of group.

The adapter first reduces modulo r^2-1 by replacing d^2 with
1-a^2-b^2-c^2. This is a Groebner normal form with d-degree at most1;
reduction does not increase total degree. Its total degree is therefore
the ORIGINAL quotient filtration degree: if a lower-degree representative
existed, its normal form would have lower degree too. Odd terms in this
normal form are rejected. They cannot represent a function invariant
under the binary centre. Odd multiples of r^2-1 have already reduced to
zero and are correctly harmless.

For even filtration degree D and canonical phase sum P_j, set

    Ptilde=sum_(j even) (r^2)^((D-j)/2)*P_j.            (1)

This is a homogeneous degree-D polynomial restricting to the SAME original
phase. It is independent of the original representative. A homogeneous
polynomial vanishing on r^2=1 is zero (by radial scaling on real S3 and
polynomial identity), so the prescribed-degree lift is unique.

The adapter checks Ptilde under exact generating quaternions of the
specified group. Group coefficients lie in QQ(sqrt2) or QQ(sqrt5); it
constructs the exact compositum with the supplied phase field, without
assuming linear disjointness. Failure of invariance is rejected: phases
are NEVER averaged to manufacture purity. Orthogonality preserves r^2,
so invariance of the homogeneous lift is equivalent to invariance of the
original regular quotient function.

## 2. Whole-marker transfer and harmless inactive coordinates

The retained [radial bridge](GAUSSIAN_RADIAL_TRANSFER_AND_EXCEPTIONAL_UNIFORMITY.md)
proves, for a homogeneous quaternion marker R of even degree e,

    E_G[R*Ptilde^m]
       =c_(D*m+e)*h(R|_(r²=1)*P^m),
    c_L=2^(L/2)*(L/2+1)! !=0.                          (2)

For an invariant original marker, its canonical homogeneous lift has
its original filtration degree. Conversely any ambient Gaussian marker
can be averaged over Gamma without changing its moments against the
invariant Ptilde. Averaging preserves homogeneous degree, and odd degrees
have zero moments because -I belongs to Gamma. Decomposing an arbitrary
marker into finitely many homogeneous pieces proves exact transfer of
ALL-original-marker eventual safety. This is not a single-unit or finite
marker-prefix inference.

If Ptilde omits some quaternion variables, the adapter compresses to its
active variables before invoking the backend. Every omitted Gaussian
variable is independent. A general marked monomial period factors as an
active marked period times the omitted variables' Wick moment. Thus
all-marker safety is preserved AND reflected; a bad active monomial
embeds by inserting zero exponents in the inactive slots. The unit period
is unchanged. If the phase is constant, one dummy variable suffices.
This is why the retained one-variable zero/constant certificates can be
reused without creating new CAS sources.

## 3. A homogeneous bad monomial gives ONE original quotient marker

The [Gaussian whole-source producer](WHOLE_GAUSSIAN_SOURCE_HOM_AND_EFFECTIVE_CUTOFFS.md)
returns a bad amplitude from its checked basis of ALL original cutoff
monomials, not an arbitrary sum of untracked amplitudes. Embed it into
four quaternion variables and write it R=q^alpha, of degree e. A bad
monomial must have even degree: odd R has identically zero Gaussian
moments against every power of the even Ptilde. Define

    Q(q)=1/|Gamma| sum_(gamma in Gamma) R(q*gamma),
    Q_original=Q|_(r²=1).                              (3)

The adapter returns the exact finite orbit-average presentation: the
monomial, full rational-pair quaternion group, its quadratic radicand,
normalization and coefficient field. This is ONE explicit polynomial
expression for a fixed regular invariant original marker, in every
order; expansion is unnecessary for its definition or replay.
Its original filtration degree is at most e, possibly strictly smaller.

Gaussian invariance and orthogonality give

    E_G[Q*Ptilde^m]=E_G[R*Ptilde^m],
    h(Q_original*P^m)=E_G[R*Ptilde^m]/c_(D*m+e).       (4)

The SAME nonzero radial factor applies to the entire averaged marker,
since averaging preserves homogeneity. No unspecified choice among
homogeneous components is needed. Infinite nonvanishing transfers in
both directions. In particular Q_original cannot be zero, since that
would make every Gaussian marked coefficient zero. The negative finite
backend certificate therefore supplies a genuine infinite ORIGINAL
marked tail, not just a nontrivial abstract integration module.

Averaging (3) is permitted for the marker because the phase was already
verified invariant. It is not a norm of a noninvariant phase, a weighted
Haar identity or an attempted averaging of a pure counterexample phase.

## 4. Eventual purity is tested separately from safety

The backend provides a nonzero scalar annihilator of each original
basis amplitude, original polynomial IBP witnesses, an infinity-indicial
bound for EVERY possible polynomial solution, and an origin-indicial
bound making a finite actual Wick jet injective. It also supplies a
checked degree-lowering recurrence covering EVERY marker. The adapter
replays those certificates through the retained backend verifier.
Completeness of the CAS integration relation matrix is not needed for
the decision certificate; the actual full period kernel is not computed.

In particular it reads the selected ACTUAL unit period, even when the
whole-source polynomial Hom is zero. If the selected Gaussian unit
period is a polynomial, (2) makes the original unit period polynomial,
so P is eventually pure. If it is exactly1, P is pure in ALL positive
orders. If it is nonpolynomial, P is not eventually pure. A zero Hom
by itself is NEVER used as a purity decision.

The output reports separately

    eventual_purity,
    all_positive_pure,
    all_original_marker_eventual_safety,
    actual_MZ_failure = eventual_purity AND NOT safety. (5)

The last uses the retained eventual-radical formulation of MZ. An unsafe
nonpure phase is not a counterexample; the constant1 control exercises
this distinction. When the selected unit period is polynomial, the
adapter divides its coefficients by the nonzero factors in (2), retaining
the finite original Haar moment polynomial as well.

No claim that a finite raw prefix decides a tail is made: the injective
origin bound and complete polynomial degree bound are essential. Nor is
a nonzero Hom on an unrelated larger source accepted as the actual map.

## 5. Numerical cutoffs retain the ORIGINAL filtration

A safe backend result carries rational C,E with

    E_G[R*Ptilde^m]=0 when m>C*deg(R)+E

for EVERY original Gaussian marker R, from the all-amplitude recurrence.
Lift an arbitrary original invariant Q at its matrix-entry filtration
degree e (the same as quaternion filtration degree under the linear
coordinate change). Equation(2) gives

    h(Q*P^m)=0 when m>C*deg_matrix(Q)+E.               (6)

Thus the SAME numbers transfer. Averaging an ambient witness does not
increase its degree. When the unit period is1, the backend's pure linear
slope also transfers without an intercept. These are certified sufficient
cutoffs, not least cutoffs or a universal bound over the full phase family.

## 6. Interface and fail-closed execution

A JSON input has `group`, sparse `phase_terms` in four quaternion
variables, and optional `field` in the Gaussian backend's exact format.
For example the zero phase is

    {"group":"2O","phase_terms":[]}

Run, from this directory,

    python exceptional_haar_pointwise_producer.py input.json result.json
    python exceptional_haar_pointwise_producer.py input.json result.json --prepare-only
    python exceptional_haar_pointwise_producer.py input.json result.json --gaussian-certificate retained-result.json

The default CAS stage cap is30seconds; timeout returns UNKNOWN. The
backend can be run uncapped with `--timeout 0`; theoretical termination
uses its retained holonomic integration/finite certificate theorem.
No degree or weight-support cap is built into the adapter. Exact
validation errors fail closed. A supplied Gaussian certificate must match
the PRECISE normalized homogeneous phase, active variables and field,
and must pass finite original-IBP/indicial/jet/coverage replay.

`verify(result)` reconstructs the quotient lift, exact group action,
radial and inactive-coordinate transfer, backend decision, fixed original
marker, unit test and cutoff. It checks source hashes and rejects changed
decision fields. This is finite algebraic certificate replay within the
same internal workflow, not independent specialist or formal assurance.

## 7. Executed evidence and the FULL remaining obligation

The checker reuses the already retained zero/constant Gaussian certificates
on BOTH groups, including a nonzero polynomial representative of zero on
the quadric. Six positive/negative quotient-level results are replayed.
The degree8 octahedral highest-weight phase is prepared and its exact
invariance/lift are checked; its Gaussian safety source is NOT computed.
Homogeneous marker averaging checks Reynolds(a²)=r²/4 on BOTH groups,
including the nontrivial radial factor4. The CLI certificate path passes,
and ten false phase/field/decision/marker/cutoff/timeout mutations are rejected.
UNKNOWN propagates without a mathematical decision.

There were ZERO new CAS constructions and NO new nonconstant safety
runs. These controls validate the adapter; they are not evidence of a
new pure phase, a degree-by-degree MZ campaign or full closure.

The old theoretical pointwise Haar route now has an executable route
through the separately implemented Gaussian backend. General coefficient-
family stratification, specialization-valid ordinary-ideal algorithms
and arbitrary transcendental coefficients are NOT implemented here.
They must not be inferred by specializing pointwise certificates.

The FULL target remains: eventual purity forces the corrected local unit
component, equivalently the intrinsic epsilon_f obstruction, to vanish
on the ENTIRE exceptional algebra, or an actual invariant pure phase with
one fixed original bad marker. A pointwise procedure cannot prove that
unbounded universal statement by listing safe examples. The sparse
four/six theorem and all-support capacity constraints remain in force;
their residual class is not replaced by a smaller implementation target.
No degree-reduction theorem for arbitrary functions is deduced from a
theorem about automorphism tuples. Both FULL exceptional cases stay OPEN.
Review is deferred and is not the next task.
