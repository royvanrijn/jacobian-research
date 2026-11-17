# Integral source products and complete Hamiltonian replicas

2026-11-16. HC4MYLV0/1.

FULL THEOREM TARGET: prove HC4/JC2 or construct an actual counterexample.
CURRENT LAST UNSOLVED IMPLICATION: ACTUAL polynomial plane source forces
Tr(Lambda_L*Lambda_L) subset B, or full marked-energy integrality; OPEN.
External review/outreach is DEFERRED. The previous goal turn made PROGRESS
with LU0/1 and a nonadopted parity-source audit. No parity theorem is used here.

The trace-dual question requires a global constraint. This note determines what
regularity of a SOURCE COORDINATE PRODUCT would give, and separates complete
holomorphic dynamics from locally finite algebraic dynamics. The latter has an
exhaustive obstruction; the former survives on proper exact whole-plane controls.
Neither the needed product integrality nor ordinary trace pairing is automatic.

## 1. LV0: product integrality upgrades the primitive and pins the boundary

Let X=Spec R=A2_(f,g) be the WHOLE plane, openly embedded in a normal AFFINE
irreducible Z=Spec S with the same function field. Assume a GLOBAL REGULAR
KAHLER one-form theta on Z restricts to a primitive of Omega=df wedge dg.
Write div_Z(Omega)=sum a_i E_i on omitted prime divisors. Necessarily a_i>=0,
because dtheta is regular on the smooth locus. Set h=f*g.

If h belongs to S, then BOTH conclusions hold:

* f*dg itself extends to a GLOBAL REGULAR KAHLER one-form on Z;
* v_(E_i)(h)>0 for EVERY omitted prime: the ENTIRE boundary lies over h=0.

Proof. On X polynomial Poincare gives

    f*dg-theta|X=dF, F in R.                         (1)

At a generic boundary point write f=u^m A, g=u^n C with A,C units. Since h
is regular, m+n>=0. Hence f*dg=u^(m+n-1) A(n C du+u dC) has at worst a
SIMPLE LOGARITHMIC pole. If F had a pole of order q>0, the du coefficient of
dF would have order -q-1<=-2 in characteristic zero, contradicting (1).
Thus F is regular at every boundary divisor. It is regular on X; normal Krull
intersection gives F in S. Therefore theta+dF is the asserted GLOBAL KAHLER
extension of f*dg, including singular points.

At a boundary valuation f,g cannot both be regular: that would give a centre
in X, contradicting separatedness of Z. Consequently one of m,n is negative
and the other positive. If m+n=0, the du/u coefficient of f*dg is n*A*C,
whose residue is nonzero. This contradicts its regularity. Thus m+n>0.
The boundary is supported on its divisorial components: otherwise normality
at a point outside those components would make both f,g regular there and
give the same forbidden source centre. Hence the entire boundary is over zero.
The argument covers horizontal boundary as well; no action or fibration
classification is used. It does not conclude that f or g separately belongs
to S. In particular product integrality is not interchangeable with balanced
derivation preservation.

## 2. The EXACT hyperbolic-replica threshold

Continue with h in S and put r_i=v_(E_i)(h)>0. If Z is proper, define

    M=max_i ceil((a_i+1)/r_i).                       (2)

For ANY nonconstant p(T) in C[T], the source Hamiltonian derivation is

    delta_(p(h))=p'(h)*(f partial_f-g partial_g).

The exhaustive LT0 Hamiltonian criterion gives

    delta_(p(h))(S) subset S
       iff p(T)-p(0) belongs to T^M*C[T].            (3)

Indeed the first nonzero exponent m of p(T)-p(0) has exact valuation m*r_i
at E_i, and LT0 requires m*r_i>=a_i+1. At a_i=0 this inequality is automatic
for m,r_i>=1, consistently with LT0's regularity-only condition. Normality
handles every singular finite point. Moreover M>=2: otherwise delta_h
preserves S, and LU0's balanced-action theorem forces Z=X.

Conversely, if delta_(p(h)) preserves S for any nonconstant p, then LT0 puts
p(h) in S. The monic equation obtained from p(T)-p(h) and normality put h
in S. Thus product regularity is EXACTLY equivalent to existence of such a
nonzero preserving hyperbolic replica. For a proper Z all preserving replicas
have p'(0)=0, and sufficiently high powers ALWAYS preserve S.

For an ACTUAL plane Keller pair, the full finite normalization supplies
theta=J^-1 P*dQ. Here h in S means precisely that the source product is
INTEGRAL over B. Replace a_i+1 by the actual ramification index e_i in (2).
LT0 also makes (3) equivalent to preservation of the FINITE trace-dual Lambda.
This is a conditional criterion, not product integrality or inversion proved.

## 3. Every preserving replica has an ENTIRE holomorphic flow on Z

The conclusion is stronger than source completeness alone. Consider generally
a locally finite derivation D of R, an invariant alpha in S with D(alpha)=0,
and delta=alpha*D preserving S. Then delta extends to a COMPLETE HOLOMORPHIC
vector field on Z, even if D does not preserve S. No smoothness of Z is needed.

For s in S choose a finite D-invariant complex vector space containing it.
Some monic polynomial q(T)=T^m+sum_(k<m)c_k T^k annihilates s. Invariance
of alpha gives the MONIC relation

    delta^m(s)+sum_(k<m)c_k*alpha^(m-k)*delta^k(s)=0. (4)

Every entry of F_s=(s,delta(s),...,delta^(m-1)(s)) is in S. Equation (4)
writes its derivative as C_s(alpha)*F_s with a companion matrix polynomial
in the REGULAR invariant alpha. Therefore

    F_s(t)=exp(t*C_s(alpha))*F_s(0)                  (5)

is jointly holomorphic on C_t times Z. Apply (5) to finitely many algebra
generators of S. On dense X it agrees with the known entire source flow:
for a fixed starting point alpha is constant along the D-orbit, so it is
the D-flow at time t*alpha. All equations defining Z and all group-law
identities hold on X and hence on Z. These coordinate functions thus define
holomorphic automorphisms Phi_t of Z with inverse Phi_(-t), for EVERY complex
time. This is a finite matrix argument, with NO division by alpha and NO
claim that the holomorphic functions in (5) are algebraic.

Apply it to D=f partial_f-g partial_g and alpha=p'(h). The flow on X is

    (f,g) -> (f*exp(t*p'(h)), g*exp(-t*p'(h))).      (6)

Thus EVERY replica in (3) is holomorphically complete on the WHOLE Z.
For p(T)=T^N, N>=2, it is NOT locally finite on the polynomial source:

    delta^k(f)=N^k*h^(k*(N-1))*f, k>=0,             (7)

are linearly independent. Completeness is not algebraic local finiteness.
For fixed s, (4) does give a finite C[alpha]-span of its entire delta-orbit;
if S is finite over B that span is contained in a finite B-module. This is
not a finite COMPLEX orbit, and does not divide away alpha at the boundary.

## 4. LV1: all divergence-free LOCALLY FINITE stabilizers under Keller failure

Let S be the FULL finite normalization of an actual plane Keller pair and
suppose its inverse does not exist. Then EVERY locally finite divergence-free
polynomial source derivation preserving S is LOCALLY NILPOTENT. The same is
true with the FINITE trace-dual Lambda in place of S. Consequently the entire
set of such fields is either {0}, or EXACTLY

    p_min(f)*C[f]*partial_g                         (8)

in LO0/1's unique integral source-coordinate pencil. In particular all these
fields commute. This includes arbitrary complex eigenvalues, arbitrary
polynomial frames and every finite-dimensional algebraic symplectic source
action preserving S. No target-field descent is assumed.

Proof. KX0's written Jordan argument applies to any invariant subspace S,
not just the target B: both Jordan parts of a locally finite delta preserve
S by finite orbit projectors. The nilpotent part delta_n is a source LND
and hence divergence-free. If delta is divergence-free, so is delta_s.

If delta_s is nonzero, the eigenvalues of finitely many source generators
give the finitely generated free weight group Gamma, exactly as in KX0.
Spectral interpolation makes S Gamma-graded. Its source volume has Gamma
weight zero, since delta_s(Omega)=0 and Gamma is an actual subgroup of C.
Choose a nontrivial one-parameter subtorus. Use the ALREADY RETAINED Gutwirth
plane linearization import from PF2GRADE1/KX0. In its polynomial coordinates
the integer weights satisfy a+b=0 and are nonzero, so its effective generator
is the balanced field. It preserves S, contradicting LU1. Thus delta_s=0.
LO0/1 now gives (8), including the no-integral-coordinate case. LT0 gives
the S/Lambda equivalence. This is a consequence of retained classification
inputs and the new LU criterion, not a new independent linearization theorem.

Non-divergence-free locally finite fields and arbitrary unbalanced tori are
not excluded. Holomorphically complete fields outside the locally finite
class need not satisfy (8), as the next all-parameter controls show.

## 5. New complete dynamics on RETAINED proper whole-plane controls

Use PF2RULED40 Section 6's EXISTING smooth affine proper completions Z_r,
for EVERY r>=2. Their geometry is credited and NOT replayed here. On its
cyclic root cover write

    B_r=C[z,v,v*(v^r-1)^j/z^j : 1<=j<=r-1],
    y=z^r, x=v/z^(r-1),
    (z,v)->(zeta*z,zeta^(r-1)*v).

The quotient contains the WHOLE source A2_(x,y), with sole omitted line E,
canonical volume order a=1 and regular primitive x*dy. Boundary charts have
v=c+z*t, c^r=1. Source translation partial_x extends algebraically to Z_r.
Those facts and affineness are retained PF2RULED40 evidence, not new tests.

The NEW observation is the globally REGULAR invariant

    h=x*y=z*v, v_E(h)=1.                            (9)

Thus M=2 for EVERY r. For EVERY N>=2 the field delta_(h^N) preserves Z_r
and is complete holomorphically on it by Sections 2–3. On a root chart its
NEW explicit expressions are

    delta(z)=-(N/r)*z^N*v^(N-1),
    delta(v)=(N/r)*z^(N-1)*v^N,
    delta(t)=(N/r)*z^(N-2)*v^(N-1)*(c+2*z*t).        (10)

For N=2 it translates the boundary line, with constant speed 2*c^2/r on
each etale root chart. For N>=3 it fixes that ENTIRE line pointwise. All
N>=2 fields are NON-locally-finite by (7), and preserve the geometric
canonical coefficient module Gamma(omega_Zr)/Omega. This latter finite
S-module is not claimed to be a finite B trace-dual of an unsupplied target.

There are even TWO GENERICALLY INDEPENDENT NONCOMMUTING complete Hamiltonian
fields with DISTINCT invariant rational fields. For any c0!=0 use the source
translation, which already extends to Z_r, and set

    h_0=x*y, h_c=(x+c0)*y,
    H=h_0^N, K=h_c^N, N>=2.

Conjugation gives completeness and preservation for both. Their Poisson
bracket is

    {H,K}=N^2*c0*y*h_0^(N-1)*h_c^(N-1).             (11)

It is nonconstant and nonzero, so their fields do not commute and are
generically independent. Since {h_0,h_c}=c0*y!=0, their invariant rational
fields are distinct. Both preserve the source line y=0; no tangent spanning
of the ENTIRE source is asserted. The bracket in (11) is not a unit, so
LT0's Poisson-trace obstruction is respected.

These are proper normal AFFINE whole-plane completions, with exact volume,
positive nontorsion canonical class, and complete dynamics on the entire
completion. They refute extending LO's TWO-LND inverse criterion to TWO
HOLOMORPHICALLY COMPLETE fields from these geometric hypotheses alone.
They do NOT give an actual finite Keller target, an actual finite-B Lambda,
or a counterexample to HC4/JC2. Their canonical order is not promoted to a
ramification index of an unconstructed target. PF2RULED40's historical
quartic queue is retired by the retained PUBLISHED Domrina theorem; this
use of its all-r geometric controls does not reopen that queue.

## 6. Evidence and next global constraint

The new exact controls verify (4), (7), the new root-chart Hamiltonian
transcriptions (10), both-field brackets and invariant distinctions. They
do not replay the registered PF2RULED40 geometry or any other registered
checker. Retained LU, LT, LO, KX, KU, Tang and quartic sources are byte-pinned.
No external theorem is newly adopted and no primary PDF is newly cached.
The classification and analytic extension are written universal proofs;
finite controls are not independent whole-proof review or formalization.

This closes the divergence-free LOCALLY FINITE stabilizer branch under
failure. It also determines exactly why passing to complete replicas can
evade that obstruction: integral source products and a coefficient vanishing
on the boundary give entire flows without dividing that coefficient in the
finite lattice. Product integrality is UNSUPPLIED for actual Keller pairs.
Even if supplied, its implication to inversion remains OPEN here.

The next effort must use the ACTUAL finite B-algebra/trace-dual compatibility,
not generic completeness or an affine exact-volume recognition claim. The
ordinary trace pairing/full marked-energy implication and unrestricted
HC4/JC2 remain OPEN. External review stays DEFERRED; all assurance flags
remain false and the user's FULL goal stays active.
