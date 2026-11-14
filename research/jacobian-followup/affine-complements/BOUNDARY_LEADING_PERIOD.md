# A leading boundary period restricts every gcd-one and gcd-two pole pair

14 November 2026. PF2RULED11 proves a necessary rational-primitive
condition for a prescribed pole-three volume. It gives all-height root
restrictions, including the complete exclusion of the square-free
cubic pole-four case left by [PF2RULED10](POLE_FOUR_NONSQUARE_EXCLUSION.md).
The condition is on a proposed coordinate's leading symbol; it is
stronger than merely asking that the global volume be exact.

Under the compatible F1 hypothesis, let m>=2 be the pole of P. Suppose
its leading symbol is h^(m/d), where d is1 or2, d divides m, and
h is a section of L^d, L=O_B(B)=O_P1(3). If a function Q with only
boundary poles satisfies dP wedge dQ=Omega, div(Omega)=E-3B, then

 eta_h = e(w)/h(w)^(2/d) dw is a rational exact one-form on B.       (1)

Here e is the leading coefficient of Omega at B, with its unique
zero divisor4q. In a coordinate on B with q=infinity, h is a polynomial
of degree at most3d. Up to nonzero scalar, its complete possible forms
are

 d=1: h=constant, or (w-a)^D, D=1,2,3;
 d=2: h=constant, or (w-a)^D, D=2,3,4,5,6.                         (2)

This is necessary, not sufficient for global functions or a finite map.
The core condition and root classification require no finiteness or
ordinary-degree theorem. Their height-eight application below uses
the already imported2022 degree alternative and retains its72,108 case.

## The intrinsic one-form

Choose local coordinates w,z with B=(z=0). Write

 Omega=(e(w)+O(z))*z^-3 dz wedge dw,
 P=h(w)^(m/d)*z^-m+(lower normal poles).

The leading volume coefficient is a section of

 (K_S+3B)|B = K_B tensor L^2.

Its divisor is E|B=4q. Under a transverse change z_new=u(w)z+O(z^2),
h_new=u^d*h and e_new=u^2*e, so the ratio in(1) is unchanged.
The usual differential factor handles changes of w. Thus eta_h is
an intrinsic meromorphic one-form, not a chart-dependent coefficient.
Its divisor is4q-(2/d)*div(h), of degree minus two.

## Formal commutation down to the first allowed remainder

Work in C(B)((z)). A leading term rho*z^-j in a remainder R has
leading wedge coefficient

 j*p'*rho-m*p*rho', at pole m+j+1,
 p=h^(m/d).

If j>2-m, that pole exceeds three and the coefficient must vanish.
It follows that rho=constant*h^(j/d).

For d=1 there is a formal m-th root V=P^(1/m), of pole one, and every
such rational homogeneous term is removed by a constant multiple of
V^j. This is a finite descending process through all j>2-m; it may
use negative powers. For m=4, in particular, it must retain the V^-1
term. Omitting that term would impose an unjustified stronger condition.

For d=2, if h is nonsquare in C(B), a nonzero rational rho exists
only for even j. The formal root W=P^(2/m), of pole two, removes
those terms by W^(j/2); odd j are forced absent. If h is square,
a formal m-th root of pole one instead removes all j. Either way,
subtracting a finite Laurent polynomial F in a formal root gives

 Q-F=z^(m-2)*rho+O(z^(m-1)),                              (3)

with rho in C(B). Since F commutes with P, the leading pole-three
coefficient is still e. At the first unforced order j=2-m it is

 e=-m*h^(2/d)*(h^((m-2)/d)*rho)'.                        (4)

The exponent (m-2)/d is an integer for d=1,2. Dividing(4) by
h^(2/d) proves(1). The primitive factor h^((m-2)/d)*rho is invariant:
rho transforms by u^(-(m-2)), canceling the change of h. It is an
actual rational function on B. The formal roots, their negative powers and F
need not be global regular functions or target coordinates. No
intermediate plane or field factorization has been inferred.

For any two pole orders m,n with gcd(m,n)=d in{1,2}, their leading
wedge vanishes in the present range m+n>2. The established leading
symbol theorem gives p=h^(m/d), q_symbol=h^(n/d), up to target scalings,
with h in H0(L^d). Therefore(1) applies to every such pair. It also
applies whenever a leading symbol has the indicated power structure,
without requiring that d be the gcd of the actual two poles.

## Complete rational-primitive classification

Choose an abstract coordinate w on B with q=infinity and trivialize
L in the standard degree-three frame. The section e then becomes a
nonzero constant. This coordinate choice is for the one-form calculation;
it is not an asserted automorphism of S, the boundary pair or the source.
Put a0=2/d, so eta_h is a nonzero constant times dw/h(w)^a0.

If h is constant, a linear primitive exists. Otherwise let D=deg(h)
and let r be its number of distinct finite roots, with multiplicities
k_1,...,k_r. A root with a0*k_i=1 would be a simple pole with nonzero
residue, incompatible with a rational derivative. Exclude that case.
If a rational primitive G exists, its pole at the i-th root has order
a0*k_i-1. Because a0*D>=2, G has a finite value at infinity and
G-G(infinity) vanishes there to order a0*D-1. Its entire pole degree is

 deg(G)=sum(a0*k_i-1)=a0*D-r.

A zero of a rational function cannot exceed its pole degree. Thus

 a0*D-1 <= a0*D-r, hence r<=1.                            (5)

This is an elementary degree argument; no classification theorem for
rational maps is imported. Conversely, a single finite root gives
an explicit primitive

 h=c(w-a)^D:  G=c^(-a0)*(w-a)^(1-a0*D)/(1-a0*D),

unless a0*D=1. This proves the full lists(2), including infinity and
all contact multiplicities. Equivalently, writing p0!=q:

* Cubic h: div(h)=D*p0+(3-D)*q, D=0,1,2,3.
* Sextic h: div(h)=D*p0+(6-D)*q, D=0,2,3,4,5,6.

For D=0 all roots are at q and p0 is unused. In particular there is
at most one root away from q. A mixed cubic divisor must contain q,
and a square-free cubic can never pass this period test.

The allowable sextic5*p0+q has a simple root at q; exactness alone
does not exclude it. Confusing a simple root at q with one away from
q would incorrectly strengthen the result.

## All-height consequence for pole four

PF2RULED9 excludes odd companions to a pole-four coordinate and lowers
divisible companions. PF2RULED10 gives its fourth-power symbol tau^4,
with tau cubic. Apply(1) with m=4,d=1. The square-free tau case, or
sextic partition(2,2,2), is now completely excluded at every companion
height. The remaining possibilities for tau are exactly

 3*p0, 2*p0+q, p0+2*q, 3*q, p0!=q.                      (6)

Thus its sextic square has only partition(6) or(4,2); the mixed case
necessarily includes q. This does not close either of those cases.
The minimum target-coordinate pole remains four.

The formal fourth-root congruence used here is explicitly

 Q=sum_{j=-1}^n c_j*V^j+O(z^2), V=P^(1/4),

in the meromorphic completion. Positive exponents divisible by four
and the constant can be removed by polynomial target shears. General
odd terms and the negative term c_-1*V^-1 must still be kept. The
period criterion follows from the next remainder coefficient, without
assuming that any of those terms vanishes.

## Sharpening the source-degree transfer at the first three contact orders

The [source-coordinate proof](SOURCE_COORDINATES.md) gives

 deg(F|U)<=max(16(m-1),16m-4r), r=ord_q(symbol(F)).

Its exact-degree statement at r=0 extends to r=1,2,3:

 0<=r<=3 implies deg(F|U)=16m-4r.                        (7)

To see this at every height, normalize the four V1 symbols to
1,xi,xi^2,xi^3 in a local cubic frame at q. Their top ordinary terms
have the form H^(4-i)/b^(3-i), with a nonzero sign. The product for
any binary monomial of total order r has degree16m-4r and a nonzero
coefficient independent of its chosen factorization. In the first
source chart its sign is(-1)^r, and in the second(-1)^m. The first
nonzero binary coefficient of a symbol of order r therefore supplies
that exact homogeneous term. Higher binary orders have smaller degree;
the lower-pole remainder has degree at most16m-16. For r<4 it cannot
cancel the term. At r>=4 this argument does not assert exactness.

## The first height-eight pairs after the period and degree gates

The restricted map is a plane Keller map of field degree four, as
proved in PF2RULED7. Use only the pinned
[2022 Theorem2.1](https://arxiv.org/html/2204.14178v1): its maximum
ordinary degree is at least125 or its degrees are72,108 up to exchange.
The [source pins](sources/degree-transfer-manifest.json) already retain
that precise theorem. No conditional2026 removal of its exception is used.

| Pole pair | Leading common symbol away from contact | Remaining contact alternative |
| --- | --- | --- |
| (5,8) | tau has divisor3*p0, p0!=q; ordinary degrees(80,128) | none |
| (6,8) | sigma has divisor6*p0, p0!=q; ordinary degrees(96,128) | div(sigma)=(6-k)*p0+k*q, k=2,3,4,6, and ordinary degrees(72,108) |
| (7,8) | tau has divisor3*p0, p0!=q; ordinary degrees(112,128) | div(tau)=(3-k)*p0+k*q, k=1,2,3, and ordinary degrees(72,108) |

In the noncontact column the relevant symbol is nonzero at q and the
source degree is exactly16 times its pole. The period classification
makes all its zeros the single point p0.

For(5,8), if tau has a zero at q, the degree bounds are at most64,112.
These satisfy neither imported alternative, so that case is excluded.
For(7,8) with a contact root they are at most96,112, leaving only
(72,108) in the displayed order.

For(6,8), let k=ord_q(sigma). The period test allows k=0,1,2,3,4,6.
If k=1, the P symbol sigma^3 has order three, so(7) gives ordinary
P degree84 exactly, while Q has degree at most112. This cannot be
the72,108 exception and is below125. For k>=2 the upper bounds are
80,112, leaving precisely the stated72,108 alternative. The absent
k=5 case would have a simple root away from q and already fails(1).

No height-eight pair is closed in full. The numerical72,108 alternative
is not assigned to a particular standard Newton-corner case from an
unproved identification of its ordinary top term. That additional
comparison, and the remaining global equations, are unresolved gates.

## Evidence and limits

The [checker](verify_boundary_leading_period.py) controls the universal
leading differential equation and primitive identity, normal-frame
invariance, every bounded cubic/sextic partition and explicit primitive,
nonexact controls, the source top-term calculation and height-eight
numeric handoff. The [receipt](boundary_leading_period_receipt.json)
separates these exact controls from the written formal descending
argument, global line-bundle interpretation and rational-degree proof.
No formal verification, independent full reproof, external review or
novelty priority is claimed.

The unrestricted finite-map problem on the marked Wright model, the
remaining mixed/triple pole-four layers, all three height-eight pairs,
actual ruled compatibility and full quartic packets remain open.

The subsequent [PF2RULED12 full(5,8) exclusion](POLE_FIVE_EIGHT_EXCLUSION.md)
uses the restricted triple-root case above and two finite-fibre pole
contradictions; it also excludes the noncontact(7,8) branch. The table
retains PF2RULED11's necessary period/degree
rows before that extra exclusion. The current first open pairs are
(6,8),(7,8), with the retained contact72,108 alternatives.

The subsequent [PF2RULED14 global cusp-jet theorem](GLOBAL_CUSP_LEADING_OBSTRUCTION.md)
excludes every one-root noncontact symbol, including all noncontact(6,8)
profiles. The current height-eight pairs(6,8),(7,8) are both contact
cases with ordinary source degrees72,108; the older tables retain their
source-specific necessary inputs.

The subsequent [PF2RULED16 full height-eight exclusion](POLE_HEIGHT_EIGHT_EXCLUSION.md)
closes both remaining pairs(6,8),(7,8), including every contact stratum.
The earlier height-eight handoffs in this note are historical. The current
minimum height is9; first reduced pairs are(5,9),(6,9),(7,9),(8,9).
