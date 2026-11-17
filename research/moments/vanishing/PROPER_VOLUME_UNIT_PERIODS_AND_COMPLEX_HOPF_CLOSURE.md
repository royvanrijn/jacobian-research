# Proper-volume unit periods and complex Hopf closure

16 November 2026. Internal written proofs with retained named Hodge inputs.
Independent replay, formal verification and external review are not claimed.
Both COMPLETE exceptional Mathieu--Zhao cases remain **OPEN**; review deferred.

The complex-height question has an all-degree answer for the ENTIRE
torus-equivariant Hopf architecture: no such original exceptional phase
is even eventually pure. No reality assumption or degree bound is used.

Two general results give the proof. A quadric phase with a projective
model carrying its ORIGINAL holomorphic volume has a locally split unit
in its actual UNIT Cauchy period module. Separately, eventual Hopf purity
is equivalent to a polynomial-density law for the even-height moments.
Every nonzero such density has a nonsplit unit logarithm at the endpoint,
which the proper model forbids. The former uniform-even law is the
constant-density special case.

This closes the torus-equivariant complex-height gate, not the full
exceptional pure locus. Other homogeneous section weights, mixed-weight
sections and general mixed-isotypic phases remain in the full problem.
No statement that purity of H(P) reflects to P is made.

See the [checker](verify_haar_proper_unit.py),
[receipt](haar_proper_unit_verification.json) and
[source audit](haar_proper_unit_source_audit.json).

## 1. A proper model gives an auxiliary source with no unit Tate logarithm

Let X=Q_N be a smooth complex affine quadric of dimension N>=2. Let M
be its compact middle cycle, and normalize an algebraic N-form Omega by
integral_M Omega=1. Suppose f:X -> A1 extends through an OPEN inclusion
j:X -> Y to a PROJECTIVE morphism bar_f:Y -> A1 with Y smooth, and
Omega extends to a HOLOMORPHIC algebraic form Omega_Y on Y.

The rank-one ambient MHS is

    H^N(X)=L=Q(-ceil(N/2)),
    weight(L)=N for N even, N+1 for N odd.              (1)

For completeness this follows from the Gysin sequence of the projective
quadric and its smooth infinity quadric. In even dimension the middle
cohomology of the completion has two Tate generators and the boundary
Gysin image has rank1; in odd dimension the kernel of the next Gysin
map has rank1 and an extra twist. The usual quadric cell decomposition
gives exactly (1).

Choose one common connected algebraic regular base S for the projective
fibres and for the complement variations. Put

    B_s=H^N(Y minus Y_s),
    K=ker(H^N(Y) -> H^N(X)).                            (2)

Restriction in (2) is onto: [Omega_Y] restricts to a nonzero generator
by the compact period, so the rational rank-one image is all of L.
K is a CONSTANT Hodge subobject. Localization gives

    H^(N-2)(Y_s)(-1) -> H^N(Y) -> B_s
       -> H^(N-1)(Y_s)(-1) -> H^(N+1)(Y).              (3)

The first Gysin image I lies in K. Indeed for sufficiently large
physical s, M avoids X_s. The map H^N(X) -> H^N(X minus X_s) is then
nonzero, hence injective by (1). Naturality of localization shows that
I restricts to zero in H^N(X). The Gysin image is a flat subobject of
the constant ambient cohomology on S, so this holds on the full base.
No injectivity of H^N(Y) -> B_s has been assumed.

Quotient B by the IMAGE of K. Its ambient term is
(H^N(Y)/I)/(K/I)=L. Thus there is an actual geometric admissible quotient

    0 -> L -> E_Y -> P -> 0,
    P subset H^(N-1)(Y_s)(-1).                         (4)

The projective smooth fibres make P pure polarizable of weight N+1,
hence GLOBAL semisimple. Disconnected fibres and all components are
allowed. This source is auxiliary to the UNIT period; it is not an
assertion that every original marker form on X extends to Y.

For N even, take a finite unipotent local cover at ANY P1 puncture.
Strictness of the actual geometric limit MHS and the nilpotent diagram
give the connecting map

    delta:ker N_P -> L(-1).                            (5)

For pure P centred at weight N+1, ker N_P has limit weights at most
N+1. The target in (5) has weight N+2. Hence delta=0: the unit splits
LOCALLY in E_Y. Genus monodromy need not be semisimple and M=W is
not assumed. Averaging over the finite cover restores the original
monodromy. Infinity is included by geometric admissibility.

For N odd, both terms in (4) have weight N+1. E_Y is therefore a
PURE polarizable geometric variation and is GLOBAL semisimple. Its
unit splits GLOBALLY, a stronger conclusion.

## 2. Descend to the actual unit Cauchy function, through its real kernel

We spell out the descent because a local functional on an arbitrary
source need not survive a quotient. For a global sequence

    0 -> L -> E -> P -> 0, P GLOBAL semisimple,

let q:E -> Q be ANY GLOBAL local-system quotient retaining L. Its
kernel D projects injectively into P. A global semisimple projection
onto its image, followed by the inverse D -> image(D), retracts E
onto D and kills L. Therefore

    E=D direct_sum E', L subset E', E' ~= Q.           (6)

Any local unit retraction on E restricts to E', hence one on Q.
No Hodge compatibility of D or local semisimplicity is needed. This
is the elementary quotient lemma already retained in the
[graph theorem](TATE_GRAPH_NECESSITY_AND_SEMISIMPLE_BOUNDARY_QUOTIENTS.md).

Use the ACTUAL compact-cycle orbit in H_N(Y minus Y_s) to define its
period quotient of B. It kills the image of K on EVERY continuation:
its image in the fixed ambient H_N(Y) does not change under monodromy,
and its original pairing with K factors through H^N(X), where K is zero.
It retains L because its normalized Omega_Y period is1. Thus this is
a genuine GLOBAL unit-retaining quotient of (4), not a chosen local
relation module. By (6) it has the local unit retraction at every
puncture; for odd N it has a global one.

The actual unit Cauchy function

    A_f(s)=integral_M Omega/(s-f)                       (7)

is a period in this quotient: Omega_Y/(s-bar_f) is a REGULAR closed
N-form on Y minus Y_s. The literal function1 is the period of Omega_Y.
Consequently the differential submodule generated by A_f and1 inherits
the local unit retraction. Regular singularity turns the invariant
covector into a differential morphism to convergent meromorphic germs,
fixing1. This is the FUNCTION module before any polynomial quotient.

We have proved, BEFORE ANY PURITY, in EVERY dimension N>=2:

    projective model + ORIGINAL holomorphic volume
      => the actual unit Cauchy module has a LOCAL
         unit retraction at EVERY puncture;
         for odd N, a GLOBAL one.                     (8)

Every polynomial reparametrization of f has the same property: compose
the projective map with the finite polynomial map, retaining the same
Omega_Y. This does not infer purity reflection or extension of all
original marker forms to Y. The model hypothesis is substantive.

## 3. Eventual Hopf purity is exactly a polynomial-density even-height law

Let A be ANY commutative C-algebra with a rational Cstar action and
normalized invariant h:A -> C. Suppose

    uv+tau^2=1, wt(u)=a, wt(v)=-a, wt(tau)=0, a!=0,
    P=(1+u)(v-(2+u)tau^2).                            (9)

The retained weight-zero calculation gives

    h(P^m)=ell(C_m), ell(R(T))=h(R(tau)),
    C_m(T)=sum_(j=0)^m (-1)^j binom(m,j)
                              binom(m+2j,m) T^(2j).   (10)

There is a universal Rodrigues identity

    C_m(T)=1/m! * d^m/dT^m [T^m(1-T^2)^m].            (11)

It follows immediately by binomial expansion. Integration by parts m
times, with all boundary terms zero at0 and1, gives

    integral_0^1 T^k C_m(T)dT=0 if k<m,
    integral_0^1 T^m C_m(T)dT
       =(-1)^m integral_0^1 T^m(1-T^2)^m dT !=0.       (12)

Thus the pairing of 1,T,...,T^(M-1) with C_0,...,C_(M-1) is triangular
and invertible. If h(P^m)=0 for m>=M, there is a UNIQUE polynomial
q(T) of degree<M matching ell(C_m) for m<M. Equation(12) makes its
integral against C_m zero also for all m>=M. Since the C_m are a
basis of the even polynomials, we obtain the full equivalence

    P EVENTUALLY PURE
      iff there is a polynomial q with integral_0^1 q=1
          and h(tau^(2j))=integral_0^1 q(T)T^(2j)dT
          for EVERY j>=0.                            (13)

Conversely a degree-K density makes h(P^m)=0 for all m>K. The density
is unique among ALL polynomials: if a nonzero difference has degree K,
its pairing with C_K has the nonzero diagonal term in (12).
The density may be COMPLEX. No positivity of q, reality of tau or
finite-prefix inference is used. The all-positive-pure special case
is q=1, recovering the earlier uniform-even criterion.

There is also an explicit eventual marked tail. Let K=deg q, let k
be its first nonzero coefficient index, and write that coefficient q_k.
For the SAME original Hopf phase in (9), the marker u^(k+1) satisfies

    h(u^(k+1) P^m)=q_k integral_0^1 T^k(1-T^2)^m dT
                         !=0 for EVERY m>K.           (14)

To verify it, the weight-zero polynomial for u^s P^m, m>=s, is

    B_(m,s)(T)=T^-s/(m-s)! *
                 d^(m-s)/dT^(m-s)[T^m(1-T^2)^m].     (15)

For m>K, integrate q(T) B_(m,s)(T) by parts m-s times.
The polynomial part of q/T^s differentiates to zero. A term q_j with
j<s contributes

    q_j binom(m-j-1,s-j-1)
                       integral_0^1 T^j(1-T^2)^m dT.

The boundary terms still vanish at0: the factors have a remaining
positive power of T, even for the Laurent terms in q/T^s. They vanish
at1 by the order-m zero. Set s=k+1; lower coefficients are zero and
the displayed binomial is1. This proves (14), including complex q_k.
Thus every abstract eventually pure Hopf phase is actually nonsafe,
with a fixed original power of u. For q=1 it gives the retained
all-order tail 4^m(m!)^2/(2m+1)! for the marker u.
If the first moments of P are nonzero, P^(K+1) is nevertheless pure in
EVERY positive order and has the same fixed bad marker by (14). Thus
this statement also gives the usual all-positive-pure MZ failure on
any abstract algebra realizing (13); eventual purity is not silently
substituted for the definition.

## 4. Every nonzero polynomial density has a unit logarithmic obstruction

Let q be ANY nonzero complex polynomial. Its even-moment Cauchy function is

    A_q(s)=integral_0^1 q(T)/(s-T^2)dT.                (16)

Split q(T)=q_e(T^2)+T q_o(T^2). Polynomial division gives

    A_q(s)=q_e(s) A_0(s)+q_o(s) A_1(s)+R(s),
    A_0(s)=1/(2sqrt(s))*log((sqrt(s)+1)/(sqrt(s)-1)),
    A_1(s)=1/2*log(s/(s-1)), R(s) polynomial.          (17)

Work over s=w^2, on the branch w=1 above s=1. This cover is unramified
there. Put L_1=log((w+1)/(w-1)), L_2=log(w^2/(w^2-1)). Then

    A_q(w^2)-R(w^2)=a(w)L_1+b(w)L_2,
    a=q_e(w^2)/(2w), b=q_o(w^2)/2.                    (18)

The differential module generated by A_q and literal1 contains a
logarithm with derivative of NONZERO residue at w=1. If just one of
a,b is nonzero, divide by that rational coefficient to obtain L_1
or L_2 modulo a rational function. If both are nonzero, b/a is a
nonzero ODD rational function, so it is not constant. Differentiate

    (A_q-R)/a=L_1+(b/a)L_2.

Subtract the rational terms L_1'+(b/a)L_2' and divide by (b/a)' to
obtain L_2 in the same differential module. These are operations over
C(w), not an assertion that a torus character extends to an unrelated
group. At w=1,

    residue L_1'=-1, residue L_2'=-1.                  (19)

A differential unit retraction to meromorphic germs would send this
logarithm to a meromorphic primitive of its rational derivative, which
is impossible: the derivative of a meromorphic germ has zero residue.
Therefore the module in (16) has a NONZERO local unit obstruction at1
for EVERY nonzero q. Zeros of q at T=1 do not remove it; division and
differentiation in (18) remain legitimate meromorphic operations.

In particular (16) cannot be the unit Cauchy function of a phase with
the proper-volume model in Section1. Even eventual equality of the
moments is impossible: finitely many discrepancies only add a Laurent
polynomial in1/s and do not change the obstruction at1. Constants tau
are excluded as well, since their Cauchy function 1/(s-tau^2) is rational
and cannot have (19).

For comparison, the ordinary noninvariant sphere height has q=1 and
exactly A_0 as the Cauchy function of its square. Its affine fibres have
punctures; the original volume does NOT extend to a projective model
over that height. Thus the model hypothesis in (8) cannot be omitted.

## 5. The entire exceptional complex-height Hopf gate is CLOSED

For Gamma=2O or2I, the full original left-torus-fixed algebra is

    A_Gamma^T=C[Q_2]^(barGamma), barGamma=O orI.        (20)

The retained
[volume theorem](VALUATIVE_BOUNDARY_VOLUMES_AND_TORUS_INVARIANT_CLOSURE.md),
Sections2--3, supplies a projective holomorphic-volume model for EVERY
nonconstant element tau of (20), including all complex coefficients and
lower terms. This is stronger input than merely saying tau is not pure.
The orbit-size estimates6/12 make every finite-boundary volume order
nonnegative. The model is an isomorphism on the original affine quadric,
so the actual normalized sphere period is the original Haar functional.
No replacement by a weighted period occurs.

Suppose u,v,tau are ORIGINAL invariant regular functions with the weights
and relation (9). Eventual purity of P forces (13). If tau is nonconstant,
compose its projective model with tau -> tau^2. Apply Sections1--2 to
the ORIGINAL unit Cauchy function of tau^2. It has a local unit retraction,
whereas (13),(16)--(19) say it cannot. If tau is constant, the rational
Cauchy contradiction after (19) applies directly. Hence

    on BOTH exceptional quotients, EVERY Hopf triple
    of opposite nonzero left weights and invariant
    weight-zero height gives P in (9) that is NOT
    EVEN EVENTUALLY PURE, in ANY degree.              (21)

No assumption that tau is real on the compact core remains. Algebraic
automorphism transport gives the same statement for the corresponding
automorphically conjugated left torus, using the retained intrinsic
Haar invariance. This does not assert that every full phase has (9).

In particular (21) applies to EVERY regular Bezout parameter b of weight
L=n+m, for which the earlier Hopf phase P_b has

    u_b=-2(U+Fb)(V+Gb), v=2FG,
    tau_b=(U+Fb)G+(V+Gb)F.

It closes the entire formerly open complex-height torus-equivariant
branch, including eventual purity and all polynomial degrees of b.
It also excludes the stronger all-positive and eventual uniform-even
height laws themselves. The earlier real-density proof remains valid
as a simpler special case, with its old reality restriction superseded.

The other homogeneous section weights and genuinely mixed-weight b
do not automatically have (9). They are NOT excluded here. Nor does
eventual purity of a polynomial H(P_b) imply eventual purity of P_b:
the subsequence or composition conditions do not force the full
triangular tail in (10). No such purity reflection or all-recomposition
exclusion is claimed.

## 6. Keep the FULL exceptional target

The complete problem remains: EVERY actual exceptional eventually pure
phase must have ZERO reduced local component, equivalently ALL corrected
actual torus character orders at0 must vanish, or one must produce an
actual invariant pure failure with an original bad marker. That whole
locus is larger than (9). The new result proves an actual all-degree
purity exclusion using exceptional invariance and unweighted moments;
it is not an equivalent criterion presented as full closure.

The all-b pullback nonzero-component controls and their uniform marker
remain. Their primitive Long phase is different from the Hopf P_b here.
The proper-volume hypothesis is not inferred from invariance and purity
of a general phase. Indeed the known nonzero bottom-row invariant weight
phases are pure and safe; the retained
[projective-volume period theorem](SEMISTABLE_LOG_PERIODS_AND_SINGLE_ISOTYPIC_SAFETY.md)
already precludes a regular-original-volume projective model for those
phases. Here the model is supplied for the AUXILIARY torus-fixed height
by (20), not assumed for the mixed-weight Hopf phase itself. A regular
dualizing generator on a singular graph is likewise not substituted
for a regular top form on the smooth model required in Section1.
No blanket exclusion of all b or all H(g_b), new degree/support campaign,
coefficient-ring promotion, universal machine proof or review gate is
introduced. Both COMPLETE exceptional MZ cases remain OPEN.
