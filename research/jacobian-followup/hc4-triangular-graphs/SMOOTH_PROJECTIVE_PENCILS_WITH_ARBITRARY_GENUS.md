# Smooth whole-plane projective pencils do not force rationality

16 November 2026. HC4MYMD0. All integers n>=1, over C. This is a
written parameter proof with elementary symbolic transcription checks.
HC4/JC2 and the full codifferent trace-pairing implication remain OPEN.
External review/outreach remains DEFERRED.

## 1. The new scope, and the prior examples

The retained [PP2 control](../hc4-pointed-pencil/PROOF.md), Section5,
already realizes every genus with a fixed singleton-zero ideal and a
fixed nonconstant Jacobian. The retained [LP1 controls](SUBMERSION_GAUSS_MANIN_AND_MARKED_VOLUME_OBSTRUCTION.md)
already show why unmarked whole-plane submersion arguments do not give a
Keller mate. Neither is superseded here.

Here the simultaneous conditions are stronger and different: the source
is the WHOLE A2, the polynomial pair has FINITE FIBERS, its image is
EXACTLY A2 minus ONE point, and its associated morphism to the COMPLETE
base P1 is smooth, surjective, with EVERY fiber irreducible. Nevertheless
the generic projective genus is arbitrarily large. Thus even all those
geometric hypotheses together do not imply rationality. The pair has a
NONCONSTANT Jacobian; this is NOT an actual plane Keller counterexample.

The degree-two seed is known. Jelonek, *A number of points in the set
C2 minus F(C2)*, Example9, gives

    H(a,b)=(a(1+ab)-b,ab),

as a quasi-finite map missing (0,-1). The author's
[uploaded primary text](https://www.researchgate.net/publication/266363557_A_number_of_points_in_the_set_C_2_FC_2)
was consulted for that example, not adopted as a new theorem. Our seed
F1 below satisfies H(iy,ix)=(-i Q1,-P1-1); it is affine equivalent to
that published example. The all-n construction is its finite source
power pullback. No novelty claim is made for the seed or general method.
All properties used here are proved directly below.

## 2. The polynomial pair and its exact image

Put u=x^n and define

    P_n=u*y-1,        Q_n=y*P_n+u=u*y^2-y+u.

There is an exact unit-ideal identity

    1=y*Q_n-(y^2+1)*P_n.                                  (1)

Thus (0,0) is omitted. Conversely, a target fiber (p,q) is described by

    u=q-p*y,       p*y^2-q*y+p+1=0,       x^n=u.           (2)

These equations are equivalent to P_n=p,Q_n=q: the quadratic says
(q-p*y)*y=p+1, after which both output identities follow. If p!=0,
the quadratic has at least one complex root and only finitely many;
each gives finitely many x. If p=0,q!=0, then u=q and y=1/q,
giving exactly n distinct x. If p=q=0, the middle equation is 1=0.
Consequently ALL fibers are finite and

    F_n(A2)=A2\{(0,0)}.                                    (3)

For p!=0,-1 and q^2-4p(p+1)!=0, the quadratic has two distinct roots,
both with u!=0. Each gives n distinct x; thus the generic degree is 2n.
The determinant is

    J_n=n*x^(n-1)*(y*P_n-x^n).                             (4)

It is nonconstant and vanishes at the source origin for EVERY n>=1.
In particular n=2 is only a NON-KELLER degree-four control. It does not
reopen the published quartic plane Keller exclusion or make any claim
about proper quartic intermediate models.

## 3. A smooth surjective pencil, including every special fiber

The unit identity makes pi_n=[P_n:Q_n]:A2->P1 a morphism. Write

    alpha_n=Q_n*dP_n-P_n*dQ_n=A_n dx+B_n dy,
    A_n=n*x^(n-1),
    B_n=x^(2n)-P_n^2=-1+2*x^n*y+x^(2n)*(1-y^2).

Its coefficient ideal is the whole ring, since

    ((2*x*y+x^(n+1)*(1-y^2))/n)*A_n-B_n=1.                (5)

On the Q_n chart, d(P_n/Q_n)=alpha_n/Q_n^2; the other chart is analogous.
Thus pi_n is smooth EVERYWHERE, including where J_n=0. It is surjective:
every nonzero linear combination aP_n+bQ_n is nonconstant, hence has a
zero over C, and there are no simultaneous zeros of P_n,Q_n.

All its fibers are irreducible. For t=P_n/Q_n finite and t!=0, the equation is

    x^n A_t(y)-B_t(y)=0,
    A_t(y)=y-t*y^2-t,       B_t(y)=1-t*y.                 (6)

At the only zero y=1/t of B_t, one has A_t(1/t)=-t!=0.
Thus A_t,B_t are coprime. In the DVR of C(y) at y=1/t,
X^n-B_t/A_t is Eisenstein: its constant term has valuation one.
It is irreducible over C(y), and primitivity gives irreducibility of
(6) in C[x,y]. This also handles t=1/2 and t=-1/2, where A_t has a
double root. For t=0, P_n=0 has coordinate ring C[x,x^-1] with y=x^-n.
For t=infinity, Q_n=0 gives x^n=y/(y^2+1), and the same Eisenstein
argument at y=0 applies. All these fibers are smooth by (5).

## 4. A uniform genus and puncture computation

For t!=0,1/2,-1/2, the normalization of the projective completion of
(6) is the cyclic degree-n cover of P1_y

    x^n=(1-t*y)/(y-t*y^2-t).                              (7)

The right side has four distinct divisor points: a simple zero at
y=1/t, a simple zero at infinity, and two simple poles at the roots
of A_t. Its denominator discriminant is 1-4t^2, and numerator and
denominator have no common root. At each of these four points the
cover has one point with ramification index n; elsewhere it is unramified.
For n=1 the same statements mean index one. Riemann-Hurwitz gives

    2g-2=-2n+4(n-1)=2n-4,       g=n-1.                    (8)

Exactly three points of this smooth compact completion are omitted
from the affine fiber: the two denominator poles (where x has a pole),
and the point above y=infinity. At y=1/t the unique point with x=0
is retained. On every other finite y with A_t!=0 the finite cyclic
model is already the affine fiber; at roots of A_t there is no affine
solution because B_t!=0. There are no further punctures. Hence

    r=3,       b1=2g+r-1=2n,       chi=1-2n.              (9)

This is a parameter proof, not a bounded genus census. Smoothness over
the nonproper base does not force the topology to be constant at all
values: behavior at infinity remains essential.

## 5. The exact differential mark that the controls lack

For an ACTUAL Keller pair K=(P,Q) with constant J!=0, suppose a target
point b=(p0,q0) is missing and write P'=P-p0,Q'=Q-q0. Then (P',Q')
is the unit ideal and pi_b=[P':Q']:A2->P1 is smooth and surjective.
Indeed dP',dQ' form a polynomial coframe and no linear combination
can be constant. Define the polynomial lifted target derivatives

    delta_p=(Q_y partial_x-Q_x partial_y)/J,
    delta_q=(-P_y partial_x+P_x partial_y)/J,
    E_b=P'*delta_p+Q'*delta_q.

They satisfy E_b(P')=P', E_b(Q')=Q'. The invertible coframe implies
E_b is zero-free exactly when b is missing. Its divergence with respect
to dx^dy is 2: the two lifted derivatives have divergence zero, and
the product rule contributes delta_p(P')+delta_q(Q')=2. Moreover

    E_b(P'/Q')=0,
    alpha_b=Q'*dP'-P'*dQ',
    i_(E_b)(dx^dy)=-alpha_b/J,       d alpha_b=-2J dx^dy. (10)

These are necessary consequences, not a theorem excluding such a field.
For the controls, (5) gives a zero-free polynomial tangent field
V_n=-B_n partial_x+A_n partial_y, but

    div(V_n)=2J_n,       d alpha_n=-2J_n dx^dy,            (11)

with J_n nonconstant. No automatic constant-divergence normalization,
polynomial lifted coframe, complete flow or local finiteness follows.
Nor is irreducibility of every pencil fiber proved for every actual
Keller pair by the preceding argument. The stronger control shows
that even adding it would not make the unmarked hypotheses sufficient.

The retained Chau2011 pencil theorem requires rationality as well as
irreducibility of every member. Its proof on printed p52 defines a
regular value to have a NONEMPTY fiber without critical points; the
empty fiber here does not satisfy that convention. The existing cached
[primary PDF](../hc4-pointed-pencil/sources/chau-pencil-2011.pdf) and PP1
are preserved. The all-rational n=1 control contradicts no published
theorem. For n>=2 the rationality hypothesis itself fails. No new
external theorem is imported, and no new primary PDF is cached.

## 6. Consequence for the research route

The smooth-surjective-projective-pencil-to-rationality shortcut is
RETIRED even with all-fiber irreducibility, finite fibers of the pair,
unit ideal, punctured-plane image and whole-plane source. Previous PP2
already retired the weaker ideal/row-operation shortcut; the new scope
is this simultaneous global control. The actual nonzero constant
Jacobian and its differential mark must be used in any stronger argument.

Excluding missing values, if established by another theorem, would
still require its precise connection to the full goal; it is not
silently substituted for inversion or full trace-pairing integrality.
The original unsolved implication remains

    ACTUAL complex polynomial plane source ==> Tr(Lambda*Lambda) subset B,

or the full marked-energy assertion. Natural alternating/symmetric and
binary affine monodromy survivors and their sheet-block quotients are
not excluded here. No intermediate model is assumed to be A2. No
independent replay, formal verification or external review is claimed.
Registered mathematical checkers are NOT replayed.
