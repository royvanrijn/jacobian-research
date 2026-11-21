# Polynomial Hopf shears give primitive Image witnesses in every prime balanced degree at least eleven

This internally reviewed research proof deduces a degree-propagation theorem
from two retained constructions. No historical-priority claim is asserted.

Work over any characteristic-zero field K. Balanced degree d means bidegree
(d,d) in (xi1,xi2;z1,z2), hence ordinary total degree2d.

## 1. Exact retained inputs

Use the quadric algebra

    A=K[x,y,t]/(t^2+2xy-1)

and its normalized polynomial moment functional

    Lambda(f)=CT_x (1/2) integral_(-1)^1
                    f(x,(1-t^2)/(2x),t) dt.                  (S1)

Integration is formal coefficientwise integration over Q, extended to K.

The retained HYPERBOLIC_ODD_HEIGHT_IMAGE_WITNESSES.md, Sections1--4,
gives, for every integer n>=1, with lambda=1,

    U=1+x,  S_n=((1+x)^(2n)-1)/x,
    C_n=U^(n+1),
    D_n=-U^(n+3) S_n^2/(4n^2),
    p_n=t(C_n+t^2 D_n)/x.                                 (S2)

This is a polynomial class in A, since C_n(0)+D_n(0)=0.
It has least degree d_n=5n+3. Put P_n=4n^2 p_n. For all m>=1,

    Lambda(P_n^m)=0,
    Lambda(x P_n^m)=0                                   (m odd),
    Lambda(x P_n^(2r))=(4n^2)^(2r) kappa_(2r) mu_(n,r)
                         !=0                            (r>=1), (S3)

where

    kappa_m=2^m m!/product_(j=0)^m(m+2j+1)               (m even),
    mu_(n,r)=(-n^2)^r (1/(2n))_r(1-1/(2n))_r
                         /((1/2)_r r!).

All factors in the last nonvanishing assertion are nonzero rational numbers.
The retained theorem also already gives radial-primitive nonproper-power
original lifts for n=1,2,3, of degrees8,13,18.

The other retained input is the polynomial moment-preserving shear in
QUADRATIC_HEIGHT_HOPF_CLASSIFICATION.md, Section4. The complete argument
is reproduced below, so no unpublished automorphism-invariance assertion
is needed.

## 2. A polynomial automorphism preserves every polynomial moment

Let delta=x partial_t-t partial_y. Then delta(x)=0 and
delta(t^2+2xy)=0. It is locally nilpotent: delta^2(t)=0 and
delta^3(y)=0, and the derivation rule proves local nilpotence on every
polynomial. For any h(x) in K[x], D_h=h(x)delta is also locally nilpotent,
because delta(h)=0 and D_h^j=h^j delta^j.

Its finite exponential is the quadric automorphism

    phi_h(x)=x,
    phi_h(t)=t+xh(x),
    phi_h(y)=y-th(x)-x h(x)^2/2.                          (S4)

Its inverse is phi_(-h). For any polynomial representative f, the localized
action of D_h on the chart in (S1) is xh(x) partial_t. Thus

    Lambda(D_h f)
      =CT_x [xh(x)/2 * (f(x,0,1)-f(x,0,-1))]=0.          (S5)

Both boundary values are POLYNOMIALS in x; the additional factor x kills
their constant terms. Applying (S5) to D_h^(j-1)f and summing the finite
exponential proves

    Lambda(phi_h f)=Lambda(f)                            (S6)

for EVERY polynomial class in A. This also applies to every marked
polynomial; it is not merely preservation of a selected finite moment
list. Since phi_h(x)=x, the marker x stays fixed.

## 3. An explicit integral polynomial representative

For n>=1 define integer polynomials

    C_n^*=4n^2(1+x)^(n+1),
    D_n^*=-(1+x)^(n+3) S_n(x)^2,
    B_n=(C_n^*+D_n^*)/x.                                 (S7)

The quotient B_n is a literal integer polynomial: its numerator has
constant term4n^2-(2n)^2=0. In A,

    P_n=t B_n-2ty D_n^*.                                 (S8)

For each integer k>=1 set h=x^k and t_k=t+x^(k+1). Define the ordinary
integer polynomial

    G_(n,k)
      =t_k [B_n-2yD_n^*+2t x^k D_n^*
                         +x^(2k+1)D_n^*].                (S9)

This is exactly phi_(x^k)(P_n). Its localized expression, useful only
for the proof, is

    G_(n,k)= (4n^2/x) t_k [C_n+t_k^2 D_n].                (S10)

No localization is used to define G: (S9) is the literal representative.
By (S3) and (S6), every positive pure moment of G vanishes and the ONE
fixed marker x has exactly the same every-even nonzero marked moments.

## 4. Exact least degree and a literal original lift

One has deg(D_n^*)=5n+1 with leading coefficient-1, and deg(B_n)=5n.
For k>=1, the unique highest total-degree term of (S9) is

    -x^D,    D=5n+3+3k.                                 (S11)

Indeed it comes from x^(k+1)*x^(2k+1)*D_n^*. The other three bracket
contributions have degrees at most5n+k+1,5n+k+3,5n+2k+3,
strictly less than D. There are no other terms of degree D.

The leading quadric t^2+2xy cannot divide x^D: this quadric is an
irreducible polynomial involving y, whereas x is not one of its factors.
Consequently no use of t^2+2xy=1 lowers the least degree in (S11).

Introduce the actual original bilinears

    R=xi1*z1+xi2*z2, Z=xi1*z2, W=2xi2*z1,
    T=xi1*z1-xi2*z2,  T^2+2ZW=R^2.                      (S12)

Write G=sum c_(i,j,l) x^i y^j t^l. Define the LITERAL polynomial lift

    F_(n,k)=sum c_(i,j,l) R^(D-i-j-l) Z^i W^j T^l.        (S13)

All exponents of R are nonnegative by (S11), and every coefficient is an
integer. Thus F has balanced degree(D,D), dehomogenizes to G and is
defined over every characteristic-zero field.

At the original point

    (xi1,xi2,z1,z2)=(1,0,0,1)

one has (R,Z,W,T)=(0,1,0,0). By the UNIQUE top term (S11),

    F_(n,k)(1,0,0,1)=-1.                                (S14)

Hence R does not divide F, with no generic-parameter or field-ordering
assumption. This explicitly proves radial primitivity.

## 5. Proper powers are excluded in the original polynomial ring

Take the original one-variable slice

    (xi1,xi2,z1,z2)=(1,0,1,z).

Then (R,Z,W,T)=(1,z,0,1), so by (S10)

    F_slice(z)=(4n^2/z) delta_k(z)
                 [C_n(z)+delta_k(z)^2 D_n(z)],
    delta_k(z)=1+z^(k+1).                               (S15)

This is a polynomial because (S9) was polynomial. Choose in the algebraic
closure a root zeta of delta_k with zeta!=-1. Such a root exists:
k+1>=2, and all roots of delta_k are nonzero and simple in characteristic
zero, with at most one equal to-1. Then

    C_n(zeta)=(1+zeta)^(n+1)!=0,
    F_slice'(zeta)=4n^2 delta_k'(zeta) C_n(zeta)/zeta!=0.  (S16)

The slice therefore has a SIMPLE zero. If F were a nonzero scalar times
H^r for any r>=2 over ANY coefficient-field extension, its nonzero
one-variable slice would be a scalar rth power and could not have that
simple algebraic zero. This contradiction excludes every scalar proper
power in the original four-variable polynomial ring. It does not merely
exclude powers in a selected invariant subring.

## 6. All-order original contractions and the unchanged marker

The retained original contraction/quadric identity says, for a balanced
homogeneous polynomial H of degree e,

    E2(H)=(e+1)! Lambda(H|_(R=1)).                       (S17)

Apply this to F^m and ZF^m. Equations(S3),(S6),(S13) give

    E2(F_(n,k)^m)=0                                     (all m>=1),
    E2(Z F_(n,k)^m)=0                                   (m odd),
    E2(Z F_(n,k)^(2r))
       =(2rD+2)! (4n^2)^(2r) kappa_(2r) mu_(n,r)!=0      (r>=1). (S18)

The factorial is (2rD+2)!, with the NEW degree D. The clearing factor
(4n^2)^(2r) is retained. Z is one fixed original polynomial independent
of r; the result is an actual two-pair Image/Mathieu failure over every
characteristic-zero field.

## 7. Degree coverage

**Theorem ShearPrimeDegreeImage.** Every prime balanced degree p>=11
admits an explicit rational/integer two-pair Image witness that is
radial-primitive, is not a scalar proper power over any extension, and
has one fixed marker Z surviving every positive even exponent.

For p=2 mod3 take n=1 and k=(p-8)/3>=1 in (S9)--(S13).
For p=1 mod3 and p>=19 take n=2 and k=(p-13)/3>=2.
The remaining p=13 is the already retained unsheared n=2 witness.
The only prime divisible by3 is3, outside this range. This proves the
EVERY-prime statement, rather than density one.

**Corollary AllLargeBalancedDegreesImage.** Every integer D>=16 is
covered by the sheared n=1,2,3 seeds and their already retained
unsheared degree13/18 members. Use n=2 for D=1mod3, n=1 for D=2mod3,
and n=3 for D=0mod3; k=(D-(5n+3))/3 is nonnegative, and the only
needed k=0 case in this range is the retained degree18.

Combining the prime theorem with the RETAINED all-even and all-odd-
composite primitive nonpower families gives a witness at EVERY balanced
degree D>=8 (and also degrees4,6). This does NOT settle degree5 or7,
or prove any arbitrary-phase minimum prime. In this combined all-degree
statement the one fixed marker may depend on the witness; the retained odd-
composite family need not use Z. The new sheared family and every prime
p>=11 DO use Z. The minimum-prime13 theorem
for the retained regular TWO-HEIGHT class is unchanged: the sheared
phase expands into height degrees0,1,2,3 and is outside that normal form.
Balanced degree11 here means original total degree22.

**Corollary DegreeCoverageEquivalence.** Under the retained all-even and
all-odd-composite primitive nonpower constructions, the following predicates
are equivalent: every balanced integer degree D>=8 has such an Image witness
with SOME fixed polynomial marker; every prime balanced degree p>=11 has
such a witness with SOME fixed polynomial marker. The first implies the
second by restriction. For the reverse, an integer D>=8 is even, an odd
composite (then D>=9), or an odd prime (then D>=11); apply the corresponding
input. This equivalence does not assert the same marker for every degree
and does not recover any technical shear identity from a degree-existence
statement. The stronger fixed-Z prime assertion was proved above.

## 8. Antecedents and verification boundary

The unipotent moment-preserving automorphism itself is already proved
in QUADRATIC_HEIGHT_HOPF_CLASSIFICATION.md Section4 and
TWO_PAIR_IMAGE_CUBIC_STRATIFICATION.md, “Safe mixed cubics outside the
orthogonal nullcone.” It is NOT a newly discovered automorphism or
invariance theorem. Equations(S2)--(S3) are the retained hyperbolic seed
and recurrence, not a new proof of the original all-order family.

PRIME_IMAGE_DEGREE_DENSITY.md already gives density ONE of prime degrees,
with its a1 quarter-density corollary. The present transport deduction
strengthens coverage to each prime>=11; it does not rediscover or claim
priority for the retained density result. The fixed-coordinate quotient
classification remains valid within its fixed-coordinate class; the
shears change that coordinate/height architecture.

The argument above supplies the all-n,k,m assertions. Finite original
contractions are useful controls only; no finite vanishing prefix
establishes those universal identities. No worldwide priority,
external review or formal verification is claimed. The all-quantifier proof
is the written argument; exact controls and source review remain internal.

