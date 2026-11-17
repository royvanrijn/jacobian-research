# The all-degree obligation as ordinary moment-ideal amplification

Status: an algebraic reformulation of the retained family-transfer theorem,
with a written proof below. It does **not** prove the required amplification,
split the intrinsic extension, or produce an exceptional counterexample.
Both full Mathieu-Zhao cases remain **OPEN**. External review is deferred.

The input is the full, no-pure-premise, arbitrary-coefficient-ring theorem in
[the family-transfer proof](FULL_MIXED_HAAR_FAMILIES_AND_TORSION_NORMALIZATION.md).
The [capacity constraints](LOGARITHMIC_JETS_AND_SPARSE_ALL_RING_SAFETY.md)
and the retained automorphism constraints remain inputs. No further
degree campaign, CAS construction, or literature theorem is used here.

## Setup and exact equivalence across all degrees

Fix either original algebra A = C[SL2]^Gamma, Gamma = 2O or 2I, and its
normalized Haar functional h. Let V_D be its full matrix-degree filtration,
choose a basis e_1,...,e_q, and put

    B_D = C[c_1,...,c_q],    F_D = sum c_i e_i,
    mu_{D,m} = h(F_D^m),    I_D = (mu_{D,m} : m >= 1).

Here I_D is the **ordinary** all-positive pure-moment ideal, without taking
its radical or imposing eventual purity. It is finitely generated because
B_D is Noetherian. Let J_D be the radical full marked-safe ideal supplied
by the family theorem; V(J_D) is the safe locus.

The following statements are equivalent, with EVERY D quantified:

1. Every original all-positive pure phase is safe for every fixed original
   regular marker.
2. J_D is contained in sqrt(I_D) for every D.
3. For every D and every integer r >= 1, there exists N(D,r) such that
   mu_{D,m} belongs to I_D^r for every m > N(D,r).
4. For every D, there exists N(D) such that mu_{D,m} belongs to I_D^2 for
   every m > N(D).

Thus the all-degree closure target can be stated as eventual vanishing of
the distinguished moment classes [mu_{D,m}] in the ordinary conormal
module I_D/I_D^2, in EVERY degree. This is a uniform identity problem,
not a consequence of computing that module in finitely many degrees.
No numerical N(D) is asserted.

### Proof

The equivalence of 1 and 2 is the definition of the two loci over C.
Assume 2. Noetherianity gives an integer a_D with J_D^{a_D} contained in
I_D. For fixed r, work over the ordinary ring R = B_D/I_D^r. The ideal
J_D R is nilpotent, since its (a_D r)-th power is zero. The retained
full family theorem therefore makes F_D safe over R, without any purity
premise on this nonreduced ring. In particular its unit-marked moments
vanish eventually in R. This is precisely 3. Clearly 3 implies 4.

Assume 4. Fix ONE actual all-positive pure P in A and ONE actual original
marker Q in A. Choose D large enough that BOTH P and Q belong to V_D;
this enlargement is essential. Write p for P's coefficient point and
v_Q for Q's constant coefficient direction in V_D. For every m >= 1,

    (partial_{v_Q} mu_{D,m})(p) = m h(Q P^{m-1}).

Every element of I_D vanishes at p. By the product rule every first
derivative of every element of I_D^2 vanishes at p. Statement 4 therefore
gives h(Q P^{m-1}) = 0 for every m > N(D). P and the same fixed Q were
arbitrary, proving 1. This argument never changes the marker with m.

This proves the equivalence. By the retained source criterion, statement 1
is also equivalent to splitting the intrinsic unit extension on every
pure phase. The result concerns ALL degrees: a square-tail identity in
one V_D only controls markers in that V_D by this first-derivative
argument. It must not be called full safety at that single degree.

Eventual purity reduces to all-positive purity ACROSS degrees: if h(P^m)
is zero for all large m, choose k with P^k all-positive pure. Apply 1 to
P^k with the finitely many FIXED original markers QP^j, 0 <= j < k.
Their eventual vanishings cover all residue classes of powers of P.
The full field conclusion then implies the retained all-ring conclusion
by the family theorem. This does not identify the two pure loci at a
fixed filtration degree.

## What an actual failure would have to do

If P is an actual pure exceptional phase and Q is ONE fixed original bad
marker, then for arbitrarily large m,

    (partial_{v_Q} mu_{D,m})(p) != 0

in any fixed V_D containing both P and Q. Those mu_{D,m} cannot belong to
I_D^2, or even vanish to second order at p. Conversely, such an actual
point p and fixed direction Q with infinitely many nonzero derivatives
are already the required original counterexample. An isolated nonzero
conormal class, a finite jet, or failures in degrees with varying points
and markers do not supply that counterexample.

The general equivalence also says that if square-tail amplification
fails somewhere across all degrees, the full field statement cannot hold.
That logical conclusion is not an explicit phase or a fixed-marker
certificate; those must still be constructed and checked in the original
algebra.

For comparison, the retained nonexceptional Long example has

    P = (1+c)(ad+b),    Q = -c

in Long's matrix convention [[a,c],[b,d]]. Its exact identities give
h(P^n)=0 and h(QP^n)=(-1)^{n-1}/(n+1), n>=1. Hence for m>=2,

    d/du h((P+uQ)^m) at u=0 = (-1)^m.

Every one of these scalar moment equations vanishes at P but has a
nonzero derivative in the SAME regular marker direction. Thus scalar
purity alone cannot imply conormal amplification in a general Haar
algebra. This control is already retained in
[the source-model audit](CONSTANT_VOLUME_MODELS_AND_ACTUAL_MARKER_DEFECTS.md);
it is not an invariant phase on either exceptional quotient.

## Why ordinary finiteness and existing constraints do not finish this

Even homogeneity and a finitely generated pure ideal do not imply 4:
in C[x,y], the sequence a_m = x y^{m-1} is homogeneous of degree m and
generates (x), but no a_m lies in (x)^2. Finite generation of I_D/I_D^2
allows distinguished classes of arbitrarily high degree to remain
nonzero by multiplication with coefficient functions. This example is
an algebraic scope control, not an exceptional Haar construction.

The full family theorem proves cofinality for the FULL marked-safe
ideal J_D. Transferring that cofinality to powers of I_D is exactly the
missing implication; it cannot be assumed in the preceding proof.

Capacity and automorphism restrictions can be used to eliminate possible
points supporting a persistent first derivative. A remaining failure
still needs at least five effective left-torus weights for 2O, or seven
for 2I, after every allowed automorphism and torus choice, and must obey
the all-support capacity inequality. These are restrictions on P, not a
proof that the whole conormal tail vanishes. No new reduction theorem
for arbitrary invariant functions is supplied here.

In particular the published tame degree reduction for automorphism
tuples cannot be applied directly to an arbitrary invariant phase.
The independent scope audit is recorded in
[the topological-coupling note](../../jacobian/cross-programme-connections/BREAKTHROUGH_SEARCH_TOPOLOGICAL_COUPLING_2026_11_16.tex).
Images of locally nilpotent derivations would also require BOTH an
appropriate image theorem and a proof that scalar-pure powers enter
that image; neither bridge has been established here.

The next mathematical work remains a proof of this amplification for
EVERY D, a direct proof of intrinsic splitting, or an actual pure phase
and fixed original marker violating it. No new decidability milestone
or finite list of degrees closes that obligation.

## A sufficient finite-generation mechanism, and an actual obstruction

Define the scalar moment subalgebra

    S_D = C[mu_{D,1},mu_{D,2},...] subset B_D.

It is connected and graded by coefficient degree: mu_{D,m} has degree m.
Finite generation of S_D for EVERY D would prove both full statements.
This is a sufficient mechanism, NOT an established theorem about the
exceptional algebras, and NOT asserted to be equivalent to MZ.

Indeed, finite generation permits choosing a finite subset of the moments
as generators: every proposed finite generator is already a polynomial
in finitely many moments. Let M_D bound their indices. For m > M_D,
write mu_{D,m} as a polynomial in those generators and retain only its
homogeneous coefficient-degree-m part. It has neither a constant term
nor a term linear in a generator. Every term is a product of at least
two positive moments. Thus

    mu_{D,m} belongs to (S_{D,+})^2 subset I_D^2
        whenever m > M_D.

More generally mu_{D,m} belongs to I_D^r whenever
m > (r-1) M_D, since a product of fewer than r generators has degree
at most (r-1) M_D. The preceding all-degree equivalence applies. This
mechanism would supply all the required powers directly, without a new
decidability calculation or marker-by-marker search.
Its sufficient direction also has a direct proof independent of the
family theorem: at a pure P, differentiate the product expression for
mu_{D,m} in the direction of ONE fixed Q in V_D. Each product contains
at least two factors that vanish at P, so m h(QP^{m-1})=0 for every
m>M_D. Enlarge D once to contain any prescribed original P and Q.

There is also an EXACT weaker formulation, without hypothesizing finite
generation of S_D itself. Consider its image in B_D/I_D^2. Its positive
part has square zero, and is precisely the C-linear span of the classes
[mu_{D,m}]. Because these classes have distinct degrees, this image is
finite dimensional over C iff [mu_{D,m}]=0 for all sufficiently large m.
Consequently finite dimensionality of these images for EVERY D is
equivalent to full closure. Coherence of I_D/I_D^2 over B_D/I_D does
not imply this finite dimensionality over C.

### Why Hilbert finite generation does not prove this mechanism

The moments are left-SL2-invariant, so

    S_D subset R_D := B_D^(SL2).

Hilbert finite generation of R_D does not give finite generation of its
subalgebra S_D. A common stronger route, proving R_D integral or finite
over S_D, is actually impossible here in every D >= 72 for 2O and
every D >= 156 for 2I. This obstruction is proved using actual PURE SAFE
phases; it is not evidence of an MZ failure.

Use the retained
[Haar-shear construction](HAAR_SHEARS_AND_GENERAL_NULLCONE_CLASSIFICATION.md).
With n=8 or12, it gives P_t in V_{n(n+1)}, t!=0, which is pure and safe,
and a positive-degree homogeneous left invariant A_D in R_D with

    A_D(P_t) != 0.

Explicitly A_D(f)=h((pi_n f)^(n+1) pi_{n(n+1)} f), of coefficient
degree n+2, and A_D(P_t)=t^n h((pq)^(n+1)) !=0 in the retained
notation. The same formula works on every larger V_D via the canonical
isotypic projections. No new nullcone classification is assumed.

Restrict coefficient polynomials to the affine line f=lambda P_t.
EVERY positive scalar moment is zero on this entire line. Therefore
restriction sends S_D to constants. It sends A_D to the nonconstant
polynomial lambda^(n+2) A_D(P_t).

If A_D were integral over S_D, its monic integral equation would restrict
to a monic polynomial equation over C for this nonconstant polynomial
in lambda. That is impossible: the degree of the leading power cannot
cancel with any of the lower powers having constant coefficients.
Thus A_D is NOT integral over S_D. In particular R_D is neither integral
nor module finite over S_D in these degrees.

This rules out using the full left-invariant coefficient ring as a finite
module to force scalar-moment finite generation. It does NOT rule out
finite generation of S_D itself. The two finite-generation assertions
must remain separate. Likewise the original Haar Reynolds map acts in
the A-variable of A tensor B_D and has image all of B_D; coefficient
Reynolds projection has image R_D. Neither is a projection onto S_D.

### A fixed bad marker forces non-finite generation on one actual pencil

Suppose P is an actual all-positive pure phase and Q is ONE fixed original
bad marker. In C[a,b] define

    nu_m(a,b)=h((aP+bQ)^m),    S_{P,Q}=C[nu_1,nu_2,...].

Every nu_m is homogeneous of degree m and vanishes at b=0. If this
algebra were finitely generated, the homogeneous argument above would
put every sufficiently high nu_m in the square of its positive ideal,
hence in (b^2). But

    partial_b nu_m(a,0) = m a^(m-1) h(QP^(m-1))

is nonzero for arbitrarily large m. Contradiction. Therefore every
actual pure fixed-marker failure yields a non-finitely generated scalar
moment algebra on a TWO-dimensional pencil. Non-finite generation of
an arbitrary pencil alone is not asserted to construct a pure bad pair.

For the retained Long P,Q above, this gives an explicit all-order
algebraic control. Modulo b^2,

    nu_1=0,    nu_m=(-1)^m b a^(m-1), m>=2.

The image of S_{P,Q} is exactly

    C direct_sum b a C[a],    (b a C[a])^2=0.

Its square-zero positive part is infinite dimensional, so this quotient
and S_{P,Q} are not finitely generated C-algebras. This uses the SAME
P and Q at every order. It is a proof from the retained exact moment
identities, not a finite collection of degree checks.

For the exceptional quotients, finite generation of S_D, or just finite
dimensionality of its image in B_D/I_D^2, remains UNPROVED across all D.
The safe-shear obstruction explains why ordinary reductive invariant
theory cannot be substituted for that proof. Capacity and automorphism
constraints are retained as restrictions on any genuine failure.
