# A submersion lift and the exact literal regular sector of an angular example

Let K have characteristic zero, and let E be the original normalized
Gaussian Wick functional with nonsingular covariance Gamma. For a
polynomial phase P write

    Phi_P[Q]=E[Q exp(tP)], V_P=Phi_P[K[X]], M_P=K[t]V_P.

The distinction between fixed spatial markers and polynomial-time
amplitudes matters even for a nonconstant angular core. This note gives
a general sufficient condition for literal equality V=M, then computes
the embedded global module and exact polynomial-marker filtration for

    S=ZW, P0=S^2, P=S^2+Z,                          (L1)

where E[ZW]=1. An independent ordinary Gaussian Y is allowed throughout.
The phase does not depend on Y; its precise filtered effect is recorded
below. The [general angular-core theorem](GAUSSIAN_ANGULAR_CORE_FULL_IMAGE_INVARIANCE.md)
already identifies M_P=M_P0. The proofs here also verify that equality
directly for this example and determine more of its original image.

## 1. A polynomial submersion makes the literal image time-closed

Suppose there are polynomial A_i such that

    sum_i A_i partial_i P=1.

For every ORIGINAL polynomial marker Q define the fixed spatial operator

    TQ=sum_i[((Gamma^-1 X)_i A_i Q)-partial_i(A_iQ)]. (L2)

Gaussian integration by parts, coefficientwise in t, gives

    Phi_P[TQ]=t Phi_P[Q].                            (L3)

Indeed apply E[partial_i F]=E[(Gamma^-1X)_i F] to
F=A_iQ exp(tP), and sum. Thus V_P is closed under multiplication by t,
and

    V_P=M_P.                                         (L4)

This holds in every Gaussian dimension, without purity, safety, finite
generation or a rank assumption. It is sufficient, not necessary: for
P=Z^2 the gradient ideal is proper, but
Phi[W^(2n)]=(2n)!t^n/n!, so V=M=K[t].

For P in (L1) one has the explicit polynomial identity

    (1-2ZW^2) P_Z+2W^3 P_W=1.

Its lift is

    TQ=(W-4W^2)Q-(1-2ZW^2)Q_Z-2W^3Q_W.             (L5)

The highest fourth-degree multiplication terms cancel. In particular
deg TQ<=deg Q+2. No time-dependent marker substitution has been used.

## 2. The radial core and its explicit literal image

Define

    F_j(t)=E[S^j exp(tS^2)]
          =sum_(n>=0)(2n+j)!t^n/n!,
    f=F_0, c_n=(2n)!/n!.

The original core literal image is exactly

    V_P0={sum_(n>=0)c_n q(n)t^n:q in K[n]}.           (L6)

To prove this, an unbalanced Z,W marker has zero period, and

    [t^n]F_j=c_n product_(i=1)^j(2n+i).

These multiplier polynomials form a triangular basis of K[n]. Ordinary
Y markers only multiply by their scalar Gaussian moments. A polynomial
period would make q(n) vanish for all sufficiently large integers, so

    V_P0 intersect K[t]={0}.                         (L7)

There is also a useful Laurent-truncation description. Set

    q_0=f,
    q_j=t^-j(f-sum_(n=0)^(j-1)c_n t^n), j>=1.

Their coefficient multipliers c_(n+j)/c_n are polynomials of degree j,
with nonzero leading coefficient. Hence

    V_P0=span_K{q_j:j>=0}.                           (L8)

This is a literal K-span, not a K[t]-span. In particular (L8) does not
put the constant period one in the core literal image.

## 3. The exact embedded global module

The original Gamma moments E[S^j]=j! give the coefficientwise identities

    2tF_1-F_0=-1,
    jF_(j-1)+2tF_(j+1)-F_j=0, j>=1.                 (L9)

For example these follow by applying the Gamma integration-by-parts
identity to S^j exp(tS^2), with the endpoint term retained for j=0.
Every coefficient is an ordinary polynomial moment, so no convergent
integral of exp(tS^2) is assumed.

The exact global actual module is

    M_P0=(K[t,t^-1]f+K[t,t^-1]) intersect K[[t]].     (L10)

Here the intersection is taken inside K((t)). For one inclusion, (L9)
expresses every F_j as a Laurent-polynomial combination of f and one,
and all original F_j are regular power series. Conversely (L8) puts
every regularized negative shift q_j in M_P0, and (L9) puts one there.
For any Laurent polynomials A,B with Af+B regular, write

    A=A_+(t)+sum_(j>=1) a_-j t^-j.

Subtract A_+f and the finite sum of a_-j q_j. The remainder is a
regular Laurent polynomial, hence belongs to K[t]. This proves (L10)
globally; no unspecified rational denominators or localization at an
ordinary parameter have been substituted for this lattice.

The series f is not rational. Over C its coefficient ratio is
c_(n+1)/c_n=2(2n+1), so it has radius zero, whereas a rational series
regular at zero has positive radius. A putative rational expression
over arbitrary characteristic-zero K uses finitely many coefficients;
embedding their generated field into C gives the same contradiction.
Thus one and f are linearly independent over K(t). Since both lie in
(L10), its ACTUAL generic rank is exactly two.

The module is not finite over K[t]. A finite list of Laurent pairs
(A,B) has a uniform lower bound on its t exponents, preserved by
polynomial-time combinations. But q_j has unique f coefficient t^-j,
with no uniform lower bound. Uniqueness follows from nonrationality of f.

## 4. Adding Z makes the entire global module literal

For nonnegative a,b,k, angular balance gives exactly

    Phi_P[Z^aW^bY^k]=0 if a>b or k is odd;
    Phi_P[Z^aW^bY^k]
      =(k-1)!! t^(b-a)F_b/(b-a)! otherwise.          (L11)

When k=0, the factor (-1)!! is one. In particular every F_b is a
literal period, from the balanced marker S^b, and all periods in (L11)
belong to M_P0. The submersion lift (L4) therefore proves

    V_P=M_P=M_P0.                                    (L12)

This is also the particular embedded equality supplied by the general
angular-core theorem. The core literal image is strictly smaller:

    Phi_P[1-2W]=f-2tF_1=1,

while (L7) excludes every nonzero polynomial from V_P0.

## 5. Exact polynomial sector at every original marker degree

Let V_P,D be the image of original markers of total spatial degree at
most D, including Y degree. Then

    V_P,0 intersect K[t]={0};
    V_P,D intersect K[t]=K[t]_(degree<=floor((D-1)/2)), D>=1. (L13)

In particular the least marker degree of a period of exact polynomial
degree n is 2n+1, even for the constant target one.

For the upper bound, define polynomials H_j,K_j by

    H_0=1,H_1=1, K_0=0,K_1=1,
    R_(j+1)=R_j-2jt R_(j-1) for both sequences.

The recurrence (L9) gives

    F_j=(2t)^-j(H_j f-K_j),
    deg H_j<=floor(j/2),
    deg K_j<=floor((j-1)/2) for j>=1.               (L14)

The initial values and recurrence prove these bounds by induction.
For b>=a, the corresponding Y-free monomial in (L11) is

    2^-b t^-a(H_b f-K_b)/(b-a)!.                    (L15)

Its Laurent-polynomial K_b part has upper degree at most
floor((b-1)/2)-a<=floor((D-1)/2) whenever a+b<=D.
For b=0 that part is zero. A linear combination giving a polynomial
period must have zero f coefficient, because f is not rational.
Its remaining Laurent polynomial is regular and has the same upper
degree bound. Even Y powers only rescale periods already realized by
lower-degree Y-free markers, and odd powers give zero. Thus Y cannot
improve the bound.

For the lower bound set R_0(W)=1-2W and use (L5):

    R_(n+1)=(W-4W^2)R_n-2W^3 R_n',
    Phi_P[R_n]=t^n.

Its degree is exactly 2n+1: the leading coefficient is multiplied at
each step by -2(2n+3), which is nonzero in characteristic zero.
These markers produce every target in (L13) by a finite linear
combination and complete the proof at all degrees.

## Evidence and scope

The companion checker
`verify_gaussian_angular_core_literal_regular_sector.py` compares bounded
original Wick contractions, the submersion identity/lift, the radial
recurrences, regularized Laurent shifts, and the filtered polynomial
sector. Its receipt pins this proof and the angular-core import. Finite
checks support the displayed algebra; the unrestricted image and
filtration statements are the written proofs above. No independent
replay, external review, formal verification or priority claim is made.
