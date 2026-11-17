# The full normalization's unique additive source pencil

2026-11-16. HC4MYLO0/1.

FULL THEOREM TARGET: prove HC4/JC2 or construct an actual counterexample.
CURRENT LAST UNSOLVED IMPLICATION: ACTUAL polynomial plane source forces
Tr(Lambda_L*Lambda_L) subset B, or full marked-energy integrality; OPEN.
External review/outreach is DEFERRED. The previous goal turn made PROGRESS
by proving the contracting finite-lattice criterion LN0/1.

This note keeps the FULL finite normalization, including RAMIFIED omitted
primes. It classifies EVERY locally nilpotent SOURCE derivation extending
to that normalization, rather than only testing a chosen coordinate's
replicas. Two distinct invariant source fields extending to it force the
Keller inverse. This is stronger than LM0/1's two-action conclusion for
the maximal etale completion, which only removes unramified added primes.
Existence of even one extending action is NOT proved automatic.

## 1. LO0: a proper affine enlargement has at most one regular coordinate pencil

Let X=Spec R=A2_C be an open subvariety of a NORMAL AFFINE surface
Z=Spec S with the same function field, so S subset R. Assume Z!=X.
The complement has a prime divisor: otherwise normal Hartogs gives
S=R and the open immersion is an isomorphism.

Suppose an ACTUAL polynomial source coordinate f belongs to S, and
write R=C[f,g]. At any omitted prime E, g must have a pole. Indeed if
both f and g belonged to the valuation ring O_(Z,E), its centre would
give a map Spec O_(Z,E)->X. Composing with X->Z and comparing with
the canonical valuation map gives the same generic map; separatedness
of Z makes them equal. Its closed point would then lie in X, a
contradiction. Thus v_E(g)<0 and v_E(f)>=0.

Apply the ALREADY PROVED PF2RULED38 valuation theorem, with its retained
Jung--van der Kulk normal-form import. It says that the ONLY polynomial
source coordinates nonnegative for this valuation are

    a*f+b, a in C*, b in C.

Consequently the coordinates in S are EXACTLY that single affine pencil.
This is a global consequence of the retained pole theorem, not a new
automorphism-classification import. It requires no degree bound, smooth
normalization, finite Keller target, or hypothesis that two individual
coordinates jointly complete a coordinate system.

In particular, for a plane Keller inclusion B subset R with full finite
normalization S, two source coordinates integral over B and not in one
affine pencil FORCE the inverse. If Z=X, finiteness of S/B gives a
finite Keller map and the retained affine-space etale criterion applies.
If Z!=X, the preceding uniqueness excludes the two coordinates.

### Extending any source additive action produces such a coordinate

Let 0!=delta be a locally nilpotent derivation of R with delta(S) subset S.
The restriction to S is nonzero since Frac(S)=Frac(R). The retained
plane Rentschler input, as used and credited in KX0, writes its source
kernel as C[f] for an actual source coordinate f.

We justify f in S rather than assume it. A local slice for the nonzero
LND on S gives

    Frac(ker(delta|S)) = Frac(S)^delta = C(f).

For completeness choose s in S with delta(s)!=0, let a=delta^m(s)!=0
be its last nonzero derivative, and put t=delta^(m-1)(s)/a. Then
delta(t)=1 and S_a=(ker delta)_a[t], by subtracting the highest powers
of t successively. Taking fractions proves the displayed equality.
Hence ker(delta|S) contains a nonconstant h(f). Since S subset R,
h is a polynomial. The equation h(T)-h(f)=0, divided by its nonzero
leading scalar, is monic over S. Normality and f in Frac(S) give f in S.

Thus if Z is proper, ALL extending nonzero source LNDs have the SAME
kernel C[f]. If no source coordinate belongs to S, there are NONE.
Two extending source LNDs with distinct invariant fields force Z=X.
These statements concern actions on the CHOSEN whole source plane.
An arbitrary LND intrinsic to S is not asserted to preserve that plane
open or extend to R.

## 2. LO1: the entire source additive stabilizer and the trace-dual

Return to an ACTUAL complex plane Keller map, R=C[f,g] when a regular
source coordinate exists, and let S be the FULL finite normalization.
For all omitted primes E_i let

    f|E_i=c_i, r_i=v_(E_i)(f-c_i)>0,
    e_i=transverse ramification index over the target.

Verticality and the COMPLETE replica threshold are already proved in
KX1, crediting PF2RULED99. Define the retained polynomial

    p_min(f)=product_c(f-c)^a_c,
    a_c=max_(c_i=c) max(0,ceil(e_i/r_i)-1).             (1)

LO0 now makes that one-coordinate calculation EXHAUSTIVE. Under
nonautomorphy the set of ALL source LNDs preserving S, including zero, is

    LND_R(S) = p_min(f)*C[f]*partial_g.                (2)

If there is no integral source coordinate, it is {0}. Every nonzero
field has kernel C[f], all fields in (2) COMMUTE, and their finite
exponentials form the additive group of translations

    (f,g) -> (f,g+p_min(f)*h(f)), h in C[f].

No field in (2) preserves the target FIELD under nonautomorphy, by KX0.
Thus extension to S is still strictly weaker than descent to B or F.
The formula closes the search over every coordinate system and arbitrary
polynomial replicas; it does not supply an integral coordinate.

### The same stabilizer acts on the actual finite trace-dual lattice

Let Lambda=Hom_B(S,B) be embedded in L by trace as in KU1. Finite duality
and KN0 identify regular reflexive canonical forms on Z with

    omega_S = Lambda*Omega,
    Omega=dP wedge dQ=J*df wedge dg, J in C*.

This is the retained canonical-module identification, not a new duality
import. Its multiplier ring is EXACTLY S:

    (Lambda:Lambda)={a in L:a*Lambda subset Lambda}=S.   (3)

Indeed at every height-one prime Lambda is a nonzero principal fractional
ideal of its DVR, so its multiplier ring is that DVR. Intersect the DVRs
and use normality of S. The reverse containment is the S-module property.

If a source LND delta preserves Lambda, then for s in S and lambda in
Lambda,

    delta(s)*lambda=delta(s*lambda)-s*delta(lambda)

belongs to Lambda. Equation(3) makes delta(s) belong to S. Conversely,
if delta preserves S, its polynomial Ga flow acts on Z and on X.
Its source Jacobian is1, so it preserves Omega. Automorphisms preserve
regular canonical forms, including their reflexive extension on normal
Z. Thus they preserve Lambda in omega_S=Lambda*Omega. Differentiating
gives preservation by delta. We obtain the EXACT equality

    LND_R(Lambda)=LND_R(S).                            (4)

This DOES NOT assume that the target algebra B is preserved by the flow.
It uses its constant-Jacobian volume, not an averaged target velocity.

Likewise any source polynomial automorphism sigma preserves S iff it
preserves Lambda: one direction uses (3), and the other uses
sigma(Omega)=constant*Omega. If (2) is nonzero, this full source stabilizer
preserves the unique integral coordinate pencil and therefore has

    sigma(f)=a*f+b, sigma(g)=c*g+h(f), a*c!=0.          (5)

Only those triangular transformations actually extending to S are meant;
no sufficiency of arbitrary a,b,c,h is asserted.

In particular TWO source LNDs with DISTINCT invariant fields both
preserving Lambda, or both preserving S, FORCE THE KELLER INVERSE in
EVERY degree. Their existence is not automatic. Two independent actions
on A=B[Lambda] alone still only give LM0/1's b0=0; (4) concerns the
FINITE MODULE Lambda and the FULL S, and removes the RAMIFIED boundary
as well. Neither the intrinsic ML invariant of S nor the extension of
all intrinsic actions to R has been classified here.

## 3. Sharp whole-plane controls in every positive canonical boundary order

The one-action branch cannot be discarded from source geometry, exact
volume, or a nontorsion canonical boundary class alone. Here is a uniform
normal affine WHOLE-PLANE enlargement with those properties for EVERY
a>=2 and q>=2. It generalizes PF2RULED38's retained a=1 chart construction;
the new controls use a>=2, not a replay of that checker.

Start in R=C[x,f]. Put

    G=x*f^q-1, H=x*G,
    N=ceil((q+a)/q), M=q*N-q-a (so 0<=M<q),
    T=x^N*f^M*G, S_(a,q)=C[f,G,H,T] subset R.          (6)

On a second chart use u,v and

    f=u, x=u^a*v+u^(-q),
    G=u^(q+a)*v,
    H=u^a*v*(1+u^(q+a)*v),
    T=v*(1+u^(q+a)*v)^N.                             (7)

These identities prove the ENTIRE affine gluing. The two principal
opens D(G) and D(G+1) cover Spec S because (G+1)-G=1.
On D(G), x=H/G recovers the source chart C[x,f,1/G].
On D(G+1), v=T/(G+1)^N and u=f recover
C[u,v,1/(1+u^(q+a)*v)]. All four generators have the polynomial
expressions (7), so these are equalities of localized rings.
Consequently Z_(a,q) is smooth and normal, and

    U=D(G) union D(f) = A2_(x,f)

is the WHOLE source plane, not a punctured plane. Its complement is
EXACTLY E={u=0}=A1_v on the second chart. The retained source line
f=0 lies in D(G), where G=-1, and is disjoint from E.

The source volume extends with divisor

    Omega=dx wedge df = u^a*dv wedge du,
    div_Z(Omega)=a*E.                                (8)

It is EXACT for q>=2. The regular primitives on the two charts are

    theta_U=x*df-d(x*f/(1-q)),
    theta_W=u^a*v*du-d(u^(a+1)*v/(1-q)).               (9)

Their difference on the overlap is
u^(-q)du-d(u^(1-q)/(1-q))=0. Their differentials equal (8).
The retained free-boundary class-group proof gives Cl(Z)=Z[E], so
its canonical class a[E] is nonzero and nontorsion, just as required
for a hypothetical full Keller normalization. This is COMPATIBILITY,
not a finite target morphism.

The exact source additive stabilizer is

    f^a*C[f]*partial_x.                              (10)

Indeed p(f)*partial_x=(p(u)/u^a)*partial_v on W; regularity is exactly
ord_0(p)>=a. Necessity is witnessed by its action on T, whose v-derivative
is a unit generically on E. Sufficiency gives a global action on both
charts. LO0 excludes any other invariant source-coordinate pencil.
The minimal field f^a*partial_x fixes the retained source line f=0
pointwise but translates the ADDED line E. One action, exact volume and
the full positive canonical boundary are therefore compatible in every
positive order a>=2. The a=1 source geometry is retained prior work.

NO finite map Z_(a,q)->A2 with the Keller restriction on U is supplied.
The parameter a+1 is a canonical-volume analogue of a transverse
ramification index, not an asserted ramification index of an unconstructed
target map. These controls do NOT contradict any actual Keller theorem,
and do not reopen the published full generic-degree-four plane case.

## 4. Scope, evidence and the remaining implication

PF2RULED38 supplies the all-word valuation theorem; KX0 supplies the
retained Rentschler source normal form and field-descent exclusion;
KX1/PF2RULED99 supply the fixed-coordinate multiplier threshold.
KN0 and KU1 supply the actual normalization, volume and trace-dual
canonical-module descriptions. LM0/1 and LN0/1 supply the weaker algebra
and stronger contracting finite-module comparisons. No external theorem
is newly adopted, and no registered mathematical checker is rerun.

The new checker verifies the NEW all-(a,q) ring/chart transcriptions,
primitive compatibility, volume order and exact additive defect witness,
plus retained byte pins and scope flags. Universal statements are proved
above, not inferred from sample charts. No formal verification,
independent whole-proof replay or external review is claimed.

The most direct remaining source-symmetry landing would be TWO distinct
source additive invariant fields preserving the ACTUAL finite Lambda or S.
Only one may exist under failure, with the exhaustive formula (2), or
none. There is no proof that the source supplies either pair. The main
ordinary trace-dual multiplication/full marked-energy gap and full
HC4/JC2 remain OPEN; the user's goal remains active.
