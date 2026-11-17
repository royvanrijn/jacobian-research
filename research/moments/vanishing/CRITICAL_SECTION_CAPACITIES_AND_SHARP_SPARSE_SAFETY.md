# Critical section capacities and sharp sparse safety

16 November 2026. Written all-degree/all-ring proof extending the retained
logarithmic capacity criterion. Both FULL exceptional Mathieu-Zhao cases
remain **OPEN**. Independent replay, formal verification and external
review are not claimed; review remains deferred. No new degree search,
CAS construction or pointwise safety run is used.

The new general input is a section estimate at capacity k-1 instead of
k. Together with strict leading instability this improves the all-support
capacity obstruction, not just a finite family: every mixed eventually
pure phase must have one entire charged side of capacity at most k-2.
It also closes every phase with at most FIVE weights on 2O and SEVEN
on 2I over arbitrary coefficient algebras. The mixed-nonpurity threshold
is sharp: actual pure SAFE mixed phases occur at six/eight weights.

## 1. At a fibre section the strict charged term saves one derivative

Use the notation and proved logarithmic estimates of
[the capacity theorem](LOGARITHMIC_JETS_AND_SPARSE_ALL_RING_SAFETY.md).
For k distinct effective charges, delta=t partial_t extracts each
component by a polynomial of degree k-1. At a valuation above the
diagonal infinity base B put

    b=v(R)>0, c=v(t)>0 at the zero section,
    A=log discrepancy, d=A-b-c.

If the ORIGINAL volume had a pole, its exact order A-1-2b-c gives
d<=b. Suppose f has finite value there, and choose a NEGATIVE component
f_i with i=-a<0, of minimum degree 2D and capacity at least k-1.
With h=v(A_i) and m the leading binary root multiplicity at the centre,
the retained extraction inequalities give the following contradictions.

If m=0, h=0 and

    (k-1)d >= D b+a c.

Capacity at least k-1 implies D>=k-1. Since d<=b and a c>0,
the displayed inequality is impossible, including D=k-1.

If m>=2, the root-jet bound h<=m d gives

    (k-1+m)d >= D b+a c.

Here cap>=k-1 implies D-M>=k-1 for the maximum multiplicity M,
so k-1+m<=D. Again d<=b and a c>0 give a contradiction.

If m=1, the stronger simple-root logarithmic coordinate gives

    (k-1)d >= D b+a c+(k-2)h.

We have k>=2 and D>=k-1, so the right side is strictly greater
than (k-1)b. This is also impossible. It applies even when another
root is multiple. All cases preserve the original volume and frame.

Thus at the ZERO SECTION a negative component of capacity >=k-1
excludes every finite-phase volume pole. At the INFINITY SECTION a
positive component of that capacity does the same in the inverse frame.
The saving uses c>0; it cannot be used directly off the two sections.
The retained finite-base theorem already treats b=0 for every mixed
phase, including all collisions, without a capacity condition.

## 2. Strict leading roots localize the remaining possible poles

Let the minimum degree of a nonconstant f be 2D. On PGL2 compactified
as P3 with determinant boundary B0, write

    f=H/Delta^D,     s=H|B0.

The boundary is P1_left times P1_right. The retained
[strict-leading theorem](STRICT_LEADING_NULLCONES_AND_FULL_THIN_CLOSURE.md)
says that EVERY nonzero eventually pure exceptional phase has a common
LEFT root of its leading tuple of multiplicity greater than D.
Such a root is unique: each left coefficient has degree2D.

We use also the precise adjunction statement already proved in
[the single-isotypic/leading audit](SEMISTABLE_LOG_PERIODS_AND_SINGLE_ISOTYPIC_SAFETY.md),
Sections2--3. At coefficient 1/D, the boundary pair is klt away from
common left fibres of multiplicity at least D. The exceptional right
root-orbit bound is essential here. If a unique common root has
multiplicity >D, every other common root has multiplicity <D, so
this klt statement applies away from that ONE left fibre.

Consequently all finite-phase original-volume poles must be centred
over that left fibre. Explicitly, near every other point inversion of
adjunction gives

    A(v)>v(Delta)+(1/D)v(H).

Finite value of f requires v(H)>=D v(Delta), so A(v)>2v(Delta)
and the original volume order A(v)-1-2v(Delta) is nonnegative.
These are the retained klt/plt inputs, not a new application with an
unverified lc equality case.

If the unique root is at left coordinate0 or infinity, its inverse
image in the global torus bundle lies on the corresponding fibre
section. In the chart [t:t xi:1:eta], at the diagonal eta=xi the
left coordinate is [t:1]. Thus a valuation above that root is either
over finite base (already regular) or has b>0 and section order c>0,
exactly the situation of Section1. A common graph may be used to
compare the two compactifications; divisorial volume order is intrinsic.

This proves the ORIENTED original-volume criterion:

    if the leading tuple has a common left root >D at0,
    and SOME negative component has cap>=k-1,
    then f has a proper model with ORIGINAL holomorphic volume.

At infinity replace negative by positive. Off the chosen section the
adjunction argument removes the poles; on it Section1 removes them.
The retained valuative model criterion then gives the smooth projective
phase model. EVERY nonconstant polynomial in f is not eventually pure.
No unit-splitting assertion for the full marked source follows merely
from this proper-volume conclusion.

## 3. The full capacity obstruction improves from k to k-1

Suppose f is mixed and has a positive and a negative component whose
capacities are BOTH at least k-1. Assume, for contradiction, that SOME
nonconstant polynomial P(f) is eventually pure.

The leading form of P(f), with e=deg P, is a nonzero scalar times s^e.
Order of a common left factor multiplies by e, because the coefficient
algebra on the right is a domain. Strict leading instability of P(f)
therefore implies that s itself has a common left root >D. This does
NOT reflect purity from P(f) to f.

Expand the leading tuple in any basis of right invariant binary forms.
Each left coefficient has at most k nonzero monomials, since its
exponents come from the original torus support. A nonzero polynomial
with l monomials has multiplicity at most l-1 at a nonzero finite
point: its first l logarithmic derivatives give an invertible
Vandermonde matrix in the distinct exponents. Hence a common root
away from0 and infinity has multiplicity at most k-1.

The capacity assumption gives D>=k-1, since each selected capacity
is bounded by that component's half-degree, which is at most D.
The strict common root >D cannot be away from the two axes. It is
therefore at0 or infinity. The corresponding selected opposite-sign
component satisfies the criterion of Section2. This gives a proper
original-volume model for f, which excludes P(f), a contradiction.

We have proved for EVERY support size and EVERY degree:

    if max_positive cap >= k-1 AND max_negative cap >= k-1,
    EVERY nonconstant P(f) is NOT eventually pure.        (1)

In particular a mixed eventually pure f must satisfy

    min(max_positive cap,max_negative cap) <= k-2.        (2)

In a torus gauge that already aligns its strict leading root to an
axis, Section2 says more: it is specifically the ENTIRE OPPOSITE
charged side whose capacities are <=k-2. Lower opposite terms are
retained; their disappearance or intrinsic splitting is not proved.

The coarse original matrix-degree bounds for one entire charged side
are consequently

    degree <=3(k-2) on2O,
    degree <=(12/5)(k-2) on2I,                            (3)

by the retained root-orbit bounds. These are bounds on ONE side,
not on the phase's largest matrix degree. They apply afresh after
every allowed automorphism and choice of maximal torus.

## 4. Five/seven weights are closed over every coefficient algebra

Every nonzero charged component has capacity at least4 on2O and
at least6 on2I. Thus (1) applies to EVERY mixed phase with

    k<=5 on2O,     k<=7 on2I.                            (4)

Every nonconstant polynomial recomposition of every such mixed phase
is NOT eventually pure. For a one-sided phase, the neutral projection
and the retained zero neutral radical give exactly the old field
classification: eventual purity iff the neutral part is0. Such phases
are pure and safe by weight selection. For constant neutral part c,
the pure recompositions are exactly P(c)=0 and are safe; a nonconstant
neutral part admits no nonzero eventually pure recomposition.

The arbitrary-ring proof from the retained sparse theorem now applies
WITHOUT CHANGE to the enlarged ranges (4). For an R-valued phase in
fixed component bases put

    J = I_0 + I_+*I_-,

where I_0 is the ideal of neutral coefficients and I_+,I_- are the two
charged coefficient ideals. Over EVERY commutative C-algebra R,

    f eventually pure iff J is nilpotent;
    then f is safe for EVERY fixed original marker.      (5)

For clarity about the transfer: work first in the finitely generated
C-subalgebra of phase coefficients. Every closed-point phase is
classified by the field proof above, so Nullstellensatz makes J
nilpotent. Conversely the same weight-counting argument proves safety
directly, including all nonreduced structure. It does not choose a
global sign over R. If J^nu=0, w is the smallest absolute nonzero phase
charge, W the largest, and B_Q the largest absolute marker charge, the
retained explicit sufficient bound is

    w m > (W+w)(nu-1)+B_Q.                               (6)

Purely neutral and zero phases have the evident nilpotent bound instead.
No new arbitrary-ring topological argument or numerical universal
cutoff computation is needed.

## 5. Pure SAFE mixed phases show the nonpurity threshold is sharp

Let n=8 on2O or n=12 on2I and fix the nonzero right-invariant binary
form F_n from the retained Haar-shear/single-isotypic construction.
It gives an injective copy V_n -> C[SL2]^Gamma of left modules.
Choose in that copy the left binary vector

    v_n=X^(n/2-1) (X+Y)^(n/2+1).                         (7)

Its nonzero monomials have indices j=0,...,n/2+1, so it has EXACTLY
n/2+2 distinct weights n-2j. Both signs and zero occur. Injectivity
of the irreducible matrix-coefficient embedding preserves that exact
weight support; this is not a bounded numerical moment observation.

For an explicit regular formula, put m=n/2+1 and

    E_j = [u^(n-j) v^j] F_n(u x+v z,u y+v w)/binomial(n,j),
    P_n = sum_(j=0)^m binomial(m,j) E_j.

The coefficient extraction is a polynomial in the ORIGINAL matrix
entries and is right-Gamma-invariant. The normalization identifies
X^(n-j)Y^j with E_j equivariantly, so P_n is exactly the vector (7).

The root at X+Y=0 has multiplicity n/2+1, strictly greater than n/2.
The entire vector is in the left nullcone. A left rotation makes its
weights strictly one-sided, proving purity and safety against EVERY
original marker. These are actual regular invariant phases, with
SIX weights on2O and EIGHT on2I, in degrees8 and12 respectively.

Thus the assertion "every mixed phase is nonpure" cannot extend one
weight past (4). These controls are SAFE, not MZ failures. Full MZ
could still hold beyond (4); its proof must allow mixed pure phases.

## 6. Updated pruning of the FULL unresolved locus

Any genuine failure has at least SIX weights on2O or EIGHT on2I
after EVERY allowed automorphism and torus choice. Degree-gap pruning
now makes one entire charged side low PW8/12 for k<=7 on2O, or
PW12 for k<=11 on2I. This is (2) with the audited higher-capacity
thresholds6 and10; it does not remove that side.

Combine (1) with
[composition amplification](CAPACITY_AMPLIFICATION_UNDER_POLYNOMIAL_COMPOSITION.md).
For f=P(r), d=deg P>=2 and r mixed, the actual extreme components
have capacities at least3d on2O and5d on2I. Pure f must now obey

    k>=3d+2 on2O,     k>=5d+2 on2I.                      (8)

In the2O sign-root branch f=G(q0^2), d=deg G>=1, the bound is
k>=4d+2. Consequently every mixed pure phase with k<=7 on2O or
k<=11 on2I is primitive on the ORIGINAL quotient. The primitive2O
phase may still be quadratic upstairs; that sign branch is retained.

The whole-quotient target remains pure-to-safe inclusion in EVERY
degree, equivalently splitting of the intrinsic unit extension on
every eventual-pure phase. The oriented low side, primitive character,
capacity and automorphism constraints are INPUTS. The new sharp
sparse closure is not a solution of that full implication. No actual
pure bad phase and fixed original bad marker has been constructed.
The ordinary conormal amplification mechanism remains unproved.
