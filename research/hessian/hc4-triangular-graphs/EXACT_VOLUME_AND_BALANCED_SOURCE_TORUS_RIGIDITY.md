# Exact volume excludes an extending balanced source torus

2026-11-16. HC4MYLU0/1.

FULL THEOREM TARGET: prove HC4/JC2 or construct an actual counterexample.
CURRENT LAST UNSOLVED IMPLICATION: ACTUAL polynomial plane source forces
Tr(Lambda_L*Lambda_L) subset B, or full marked-energy integrality; OPEN.
External review/outreach is DEFERRED. LT0/1 made PROGRESS by classifying
all preserving divergence-free source fields and showing that finite
source jets cannot supply a global action. This note proves a GLOBAL
rigidity statement for one balanced torus action, in EVERY degree.

## 1. LU0: an exact affine enlargement cannot preserve the balanced action

Let X=A2_C=Spec C[f,g] be a dense open subvariety of an irreducible
AFFINE variety Z=Spec S. Thus S is a finitely generated subalgebra of
R=C[f,g] with the same function field. Suppose

* the source action (f,g)->(t*f,t^-1*g) extends ALGEBRAICALLY to Z;
* there is a GLOBAL REGULAR KAHLER ONE-FORM theta on Z whose
  differential restricts to df wedge dg on X.

Then Z=X. Normality, smoothness of the added boundary, a finite target
map and nonvanishing of d theta on the boundary are NOT needed for
this general rigidity theorem. Affineness and the global primitive
are needed, as the explicit controls below show.

Write delta=f partial_f-g partial_g. Take the weight-zero part theta_0
of theta under the torus action. This is a regular KAHler form: every
element of the algebraic module Omega^1_S has a finite weight
decomposition. Exterior differentiation commutes with the action,
so on X

    d theta_0=df wedge dg,   Lie_delta theta_0=0.

Cartan's identity on X gives

    d( i_delta theta_0+f*g )=0.

The function in parentheses is constant. Evaluate at the SOURCE
origin f=g=0, which is fixed, to see that the constant is zero. Hence

    h=f*g=-i_delta theta_0 belongs to S,            (1)
    h(z)=0 for EVERY fixed point z of Z.           (2)

The second assertion includes singular fixed points: contraction of
sum a_j d b_j is sum a_j delta(b_j), and each delta(b_j) vanishes at
a fixed point. No nonsingular tangent-space assertion is required.
Restriction to X is injective on functions; equality(1) is therefore
an equality in S, not just an equality of rational functions on a
chosen chart.

The invariant algebra is EXACTLY

    S^Gm=C[h].                                    (3)

Indeed S^Gm is contained in R^Gm=C[fg], while(1) gives the reverse
containment. We now give the required affine orbit arguments explicitly
rather than import a quotient theorem.

First, two disjoint CLOSED invariant orbits in an affine variety are
separated by an invariant regular function. Their ideals are comaximal,
so choose a function with values0 and1 on the two orbits. Its
weight-zero component still has these values, because both orbit
ideals are homogeneous. It is the desired invariant function.

Second, every nonempty closed invariant subset of an affine variety
with a Gm action contains a closed orbit. One elementary verification
uses homogeneous algebra generators. On an orbit they restrict to
constant multiples of t^w. If nonzero restrictions occur in both
positive and negative weights, their generated semigroup gives a
Laurent coordinate algebra and the orbit is closed. If only one sign
occurs, its closure has a fixed endpoint at t=0 or t=infinity. If all
weights are zero, the starting point is already fixed. This argument
also handles finite stabilizers and does not assume normality.

For c!=0 the source fibre h=c is the SINGLE Gm orbit

    O_c={(f,g):fg=c}.

Its closure in Z can add only finitely many points, which must be
fixed because Gm is connected. All those points would still have
h=c, contradicting(2). Thus O_c is CLOSED IN Z. For c=0 the source
origin O_0={(0,0)} is already a closed fixed orbit. Consequently every
h-fibre contains a closed orbit lying in X.

If Z minus X were nonempty, it would be a closed invariant subset,
so would contain a closed orbit O. Say h(O)=c. It is disjoint from
O_c in X. The first orbit argument separates O and O_c by an
invariant. But by(3) every invariant is a polynomial in h and takes
the same value on both orbits. Contradiction. This proves Z=X.

The source-plane orbit structure and the exact primitive both enter
the proof. It is not an argument about completeness of an arbitrary
polynomial vector field, nor an automatic linearization theorem.

## 2. LU1: one balanced field preserving the finite trace dual forces inversion

Let K=(P,Q):A2->A2 be an ACTUAL complex plane Keller map. Put
B=C[P,Q], R=C[f,g] in ANY actual polynomial source coordinates,
S the FULL finite normalization of B in Frac R, and Lambda its finite
trace dual. The retained source open embedding is X subset Z=Spec S.
Write J=Jac_(f,g)(P,Q) in C*. The GLOBAL regular form

    theta=J^-1*P*dQ in Omega^1_S

satisfies d theta|X=df wedge dg. It is crucial that P and Q are in
S: this is a primitive ON THE ENTIRE affine normalization.

Suppose the SINGLE source derivation

    delta=f partial_f-g partial_g                 (4)

preserves S. Every polynomial source element has finitely many delta
weights. Spectral projection in delta recovers each weight component
of every s in S inside S. Thus S is graded and the balanced Gm action
extends algebraically. LU0 gives S=R. Since S is B-finite, the original
Keller map is finite, hence invertible by the retained finite-etale
affine-plane criterion.

By LT0/LO1, preserving the actual finite Lambda is EQUIVALENT to
preserving S for this divergence-free Hamiltonian field. Alternatively
the implication Lambda-stable=>S-stable follows directly from
(Lambda:Lambda)=S and the Leibniz identity. We have proved

    K inverse
      iff SOME polynomial source coordinate frame (f,g) has
          (f partial_f-g partial_g)(Lambda) subset Lambda
      iff SOME such frame has delta(S) subset S.    (5)

The inverse case works in EVERY source coordinate frame. Under
failure EVERY balanced source torus in EVERY polynomial frame fails
to preserve Lambda and S. Source translations and arbitrary source
polynomial changes are included in this quantification.

No preservation of B or of the target FIELD is assumed. This is a
new exact-volume exclusion of the BALANCED mixed-sign branch; LN0's
contracting positive-weight argument does not apply to(4). It does
not exclude unbalanced mixed-sign actions or all locally finite fields.
No general source-torus linearization result is imported or needed.

For the FINITE-TYPE algebra A=B[Lambda], LU0 similarly shows that
balanced stability would make its whole affine etale completion equal
to X. That only removes the unramified added boundary; it does NOT
remove ramified primes outside this completion. Criterion(5) uses the
FULL finite S or FINITE MODULE Lambda and therefore forces the full
inverse. These domains must not be interchanged.

## 3. New sharp whole-plane controls: balance alone does not suffice

For EVERY a>=2 define inside R=C[x,f]

    G=x*f-1, H=x*G, T=x^(a+1)*G,
    S_a=C[f,G,H,T].                               (6)

This is the NEW q=1 member of the chart pattern whose q>=2 members
were proved in LO1. We prove its geometry here; no registered checker
or old q>=2 chart calculation is replayed.

On the second chart use

    f=u, x=u^a*v+u^-1,
    G=u^(a+1)*v,
    H=u^a*v*(1+u^(a+1)*v),
    T=v*(1+u^(a+1)*v)^(a+1).                      (7)

The principal opens D(G) and D(G+1) cover Spec S_a. On D(G), x=H/G
recovers exactly C[x,f,1/G]. On D(G+1), u=f and
v=T/(G+1)^(a+1) recover exactly
C[u,v,1/(1+u^(a+1)*v)]. Formula(7) gives both inclusions in these
localized ring equalities. Thus Z_a=Spec S_a is smooth and normal.

The open D(G) union D(f) is the WHOLE source A2_(x,f): on D(f),
x=(G+1)/f; in R the ideal(G,f) is the unit ideal. The complement is
EXACTLY E={u=0}=A1_v. In particular Z_a is a proper AFFINE enlargement.
The source volume extends regularly with positive order

    Omega=dx wedge df=u^a*dv wedge du, div Omega=aE. (8)

The balanced source action x->t*x,f->t^-1*f preserves S_a:
the four generators f,G,H,T have weights -1,0,1,a+1. In the second
chart u has weight -1 and v weight a+1. There are TWO distinct fixed
points: the source origin has moment x*f=G+1=0; the boundary point
u=v=0 has moment G+1=1. Their different moment values prove that
Omega has NO global regular Kahler primitive, by the averaging and
Cartan argument in Section1.

The failed primitive gluing is also explicit. On X and on the second
chart take

    theta_X=x*df,   theta_W=u^a*v*du.

Both differentiate to Omega, but on the overlap

    theta_X-theta_W=du/u.                          (9)

Their invariant moment constants are respectively0 and1. No claim
of exactness is made for these NEW q=1 controls. This pinpoints why
the q>=2 exactness proof in LO1 could not be used at q=1.

Conversely the RETAINED q>=2 LO1 controls have exact volume and an
unbalanced mixed-sign action with source weights(-q,1). Their source
volume has weight1-q!=0. They show that the BALANCE condition cannot
be omitted. No finite Keller target is supplied for either control
family, and a+1 is not called the ramification index of an
unconstructed target.

## 4. Affineness is also essential

For completeness, there is an elementary NONAFFINE exact control.
Glue two affine planes (q,p) and (v,w) by

    v=1/q,   w=-q^2*p.

This is the total cotangent line bundle T*P1. Its two one-forms
-p*dq and -w*dv agree and give a global primitive of the volume
dq wedge dp=dv wedge dw. The balanced action on the first source
plane extends with weights(-1,1) on (v,w). It is a proper enlargement
of the WHOLE first A2 chart, with an added A1 and exact volume.
It is NOT AFFINE: its zero section p=w=0 is a closed P1. Both zero-
section fixed points have moment value0; affine invariant separation
is precisely the unavailable step. Passing to its affine coordinate
ring contracts this zero section and loses the given source open
embedding, so does not produce an affine counterexample to LU0.

## 5. What this closes and what remains

All degrees are allowed in the Keller map and in the supplied
polynomial source frame. The proof uses algebraic weight projection,
Cartan's identity, an elementary orbit-closure argument and invariant
separation by Chinese remainders. Standard affine invariant-theory
facts are proved in the special form needed here; no external
linearization, classification or new proved theorem is imported.

This closes the balanced mixed-sign source-action branch WITHOUT
target-field descent. It gives a concrete sufficient global landing
for the actual finite trace dual beyond positive contracting fields.
It does NOT show that ANY balanced field preserves that lattice.
LT1's arbitrary finite source jets can still imitate such a field at
each finite order without producing the exact global field(4).

The remaining whole-plane implication must force this landing, the
original ordinary trace pairing, or another genuinely universal
global condition. Full HC4/JC2 and the marked-energy gate remain
OPEN. All assurance flags and the published full quartic retirement
are unchanged. External review/outreach remains DEFERRED.
