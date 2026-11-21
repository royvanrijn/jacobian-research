# Mixed cubic Hopf phases: zero center and a conditional off-center reduction

This note proves two restricted reductions inside the full mixed cubic
height box. It does not classify all cubic-height phases and does not settle
the global primitive balanced degrees 5 or 7.

Work over a field K of characteristic zero in
A=K[x,y,t]/(t^2+2xy-1), embedded in K[x,x^-1,t] by
y=(1-t^2)/(2x). The functional is

    Lambda(f)=CT_x (1/2) integral[-1,1] f(x,(1-t^2)/(2x),t) dt.

Here integration means the rational polynomial functional. A phase is pure
if Lambda(f^m)=0 for every positive integer m. It is safe if, for every fixed
g in A, Lambda(g f^m)=0 for all sufficiently large m.

Consider a regular cubic-height expression

    f=x^-1 P(x,t),   P=C+tB+t^2 D+t^3 E,

where C,B,D,E are polynomials in x. Regularity is exactly

    C(0)+D(0)=0,   B(0)+E(0)=0.                       (M1)

The associated graded quadric shows that its least polynomial degree is

    max(deg C-1, deg B, deg D+1, deg E+2),             (M2)

with the convention deg 0=-infinity. Here are direct justifications.
Modulo x the relation is t^2=1, so (M1) follows by the linear
independence of 1,t over K[y]. Conversely the regular representative is

    2C0 y+2B0 yt+(C-C0)/x+t(B-B0)/x
      +t^2(D+C0)/x+t^3(E+B0)/x.

It has no xy monomial. Under a degree-compatible monomial order with
leading term xy in t^2+2xy-1, this is a reduced normal form. A nonzero
top homogeneous part of a reduced polynomial cannot be divisible by
t^2+2xy, whose leading monomial is xy. Thus the sphere quotient cannot
lower its degree. Its displayed monomials give exactly (M2), including
the constant terms represented by y and yt.

## 1. Classical interval rigidity and the complete zero-center branch

A classical result of Francois--Pakovich--Yomdin--Zhao, Corollary 4.1
(2011, printed p.27), says that a complex polynomial whose interval power
moments vanish eventually must be zero. Muger--Tuset, Theorem 2.1
(2019, printed p.2), proves the stronger positive-limsup statement for
nonzero polynomials on [0,1]. The affine substitution t=2z-1 supplies
the same vanishing conclusion on [-1,1].

This implication holds over every characteristic-zero field K: the
finitely many coefficients of a putative nonzero polynomial generate
a finitely generated extension of Q, and that field embeds in C.
Apply the complex theorem to the embedded polynomial. One need not,
and does not, assume that the whole field K embeds in C.

Consequently, for EVERY f in K[x,t], with no height bound,

    f is pure  iff  f is divisible by x in K[x,t].     (M3)

Indeed let h(t)=f(0,t). Nonnegative angular weights give
Lambda(f^m)=(1/2) integral[-1,1] h(t)^m dt. Interval rigidity forces h=0.
Conversely, a phase divisible by x has strictly positive angular weight.
Every fixed g in A has bounded negative Laurent x support, so g f^m
has no constant x term once m exceeds that bound. Thus every pure
phase in K[x,t] is safe.

In the cubic presentation, assume the entire center vanishes:

    C(0)=B(0)=D(0)=E(0)=0.

Then f lies in K[x,t], and (M3) becomes

    f pure iff x^2 divides C,B,D,E; every such f is safe.       (M4)

In particular this entire branch supplies no unsafe degree-5 or degree-7
phase. This conclusion is an application of classical interval rigidity,
not a new claim about that classical theorem.

## 2. A direct first-four-moment certificate for height at most three

The implication in this section is purely rational and works directly
over every characteristic-zero field. If deg h<=3 and

    (1/2) integral[-1,1] h(t)^m dt=0,  m=1,2,3,4,

then h=0. This is an exact calculation, independent of the analytic
theorem in Section 1.

If the cubic coefficient is nonzero, divide h by it. The first moment
then puts h in the form

    h=t^3+b t+c(t^2-1/3).

Write

    q=105b^2+126b+45,
    P=56b^2+105b+45,
    R=28028b^4+84084b^3+93275b^2+46410b+8955.

The next moments are

    mu2=(q+28c^2)/315,
    mu3=4c(63b^2+108b+4c^2+45)/945,

and

    135135 mu4 =
      27027b^4+77220b^3+37752b^2 c^2+90090b^2
      +61776b c^2+49140b+2288c^4+26520c^2+10395.

Modulo mu2=0 these reduce to

    mu3=8c P/2205,    mu4=-16 R/735735.               (M5)

If c=0, then q=R=0. Polynomial division gives

    R mod q = -84(117b-110)/125.

Hence b=110/117, whereas q(110/117)=1169375/4563 is nonzero.

If c is nonzero, (M5) gives P=R=0. This time

    R mod P = -315(299b+159)/128.

Thus b=-159/299, but P(-159/299)=446976/89401 is nonzero.
These contradictions exclude a nonzero cubic coefficient.

For a nonzero quadratic coefficient, normalize to
h=t^2+b t-1/3. Then

    mu2=b^2/3+4/45,   mu3=4b^2/15+16/945,
    mu3-(4/5)mu2=-256/4725.

A linear polynomial is killed by its first two moments, and a constant
by its first. This proves the assertion.

The number four cannot be replaced by three for all complex cubics:
take c=0 and any complex root b of q. The first and third moments vanish
by oddness and the second vanishes by q=0, although h is nonzero.
The fourth moment is nonzero by the preceding coprimality calculation.
This note makes no historical-priority claim for this finite certificate.

## 3. A separate depressed-cubic invariant restriction

The remainder imposes an ADDITIONAL polynomial identity which has not
been derived from purity in the unrestricted mixed cubic box:

    I=27 C E^2-9 B D E+2D^3=0.                        (M6)

Assume E(0) is nonzero. Put, in K(x),

    s=-D/(3E),   h^2=D^2/(3E^2)-B/E.

Then (M6) is exactly the identity

    P=E[(t-s)^3-h^2(t-s)].                            (M7)

Regularity gives, at x=0,

    I(0)=2 C(0)(9 B(0)^2-C(0)^2).

The following argument concerns only the two OFF-CENTER cases

    C(0)=3B(0) or C(0)=-3B(0),    B(0) nonzero.        (M8)

It does not treat the remaining C(0)=0, E(0) nonzero case.
Reflection t -> -t reduces (M8) to C(0)=3B(0), which has

    s(0)=-1,   h(0)^2=4.

Choose the formal square root h with h(0)=2 and set

    A0=-E h^3,    L=h/2.

All of these are units at zero. Substituting v=(t-s)/h turns the
numerator into A0 v(1-v^2). Its transformed endpoints are 0 and 1
modulo x. An antiderivative of [v(1-v^2)]^m changes by a multiple of
x^(m+1) when either endpoint is perturbed by O(x): its derivative
has a zero of order m at each endpoint. Therefore

    (1/2) integral[-1,1] P(x,t)^m dt
       = kappa_m A0^m L  modulo x^(m+1),              (M9)
    kappa_m=integral[0,1] [v(1-v^2)]^m dv.

Every kappa_m is a nonzero rational number; for example it is
(1/2) Beta((m+1)/2,m+1), or directly a positive real integral with
rational polynomial integrand. Thus f is pure precisely when

    [x^m] A0^m L=0  for every m>=1.                   (M10)

For a prescribed unit A0 and L(0)=1, these equations determine each
successive coefficient of L uniquely, since its coefficient in the
mth equation is A0(0)^m. Their unique solution is

    L=1-x A0'/A0,    hence h/2=1-x(Eh^3)'/(Eh^3).      (M11)

Indeed the coefficient in question is
[x^m](A0^m-x(A0^m)'/m)=0. This also proves the converse.

## 4. Rationality, poles, and low-degree consequences of (M6)--(M8)

Equation (M11) makes h rational: h^2 and h'/h are rational, as is E'/E,
so its right side is rational. At a nonzero finite zero of h of order j,
the rational function A0=-Eh^3 has positive order ord(E)+3j.
The right side of (M11) would have a nonzero simple pole there, contrary
to the left side vanishing. There is no zero at zero. Thus h has no
finite zero, and can be written

    h=2/Q,  Q in K[x], Q(0)=1.                        (M12)

Every pole of h is simple, because a rational logarithmic derivative
has at most a simple pole and the pole is not at zero. Hence Q is
squarefree. If q=deg Q>=1, taking the value of (M11) at infinity gives

    deg E=3q+1.                                      (M13)

If Q=1, (M11) instead gives E constant.

Suppose the least degree of f is at most 7. By (M2), deg E<=5, so
(M13) leaves only the following possibilities.

(a) Q=1. Here E is constant, h=2 and D polynomial makes s polynomial.
Since s(0)=-1, write s=-1+xH. Then f is the Hopf shear t -> t-xH,
x fixed, of

    E[(t+1)^3-4(t+1)]/x = -2E y(t+3).                (M14)

The unsheared phase has strictly negative angular weight, and is safe.
For completeness, the needed arbitrary-polynomial shear invariance has
a direct proof. Put J=x partial_t-t partial_y. For H in K[x] and g in A,
the Laurent chart gives Jg=x partial_t g, and integration gives

    Lambda(H Jg)=CT_x {xH(x)[g(x,0,1)-g(x,0,-1)]/2}=0.

Both boundary restrictions are polynomials in x and xH has strictly
positive x order (unless H=0, when the assertion is immediate).
The derivation HJ is locally nilpotent; its finite exponential fixes x
and sends t to t+xH. Every positive derivative term has Lambda zero
by the same calculation. Thus exp(HJ) preserves Lambda, and its
inverse transports each fixed marker. Applying it with -H proves that
f remains safe. No classification of Aut(A) is used.

(b) q=1 and deg E=4. An E-zero away from Q would make the right side
of (M11) have a pole while h is a nonzero finite value, impossible.
Therefore E is a nonzero constant times Q^4.

The only possible finite pole of s=-D/(3E) is the root of Q. If its
order were j>=2, polynomiality of
B=E(3s^2-h^2) would first force j<=2, while polynomiality of
C=E(-s^3+h^2s) would then fail at j=2. Hence

    s=p+a/Q,    p in K[x], a in K.                    (M15)

If p has degree at most zero, (M2) gives exact least degree 6 because
deg E+2=6 and every other entry is at most 6. If j=deg p>=1, the
degree is

    max(3j+3,2j+4,j+5,6)=3j+3.                       (M16)

For j=1 this is 6; for j>=2 it is at least 9. Thus case (b) has no
least degree 5 or 7.

We have proved: inside the off-center restricted identity branch
(M6)--(M8), every pure phase of least degree 5 or 7 is safe. Thus there is no
unsafe phase of either least degree 5 or 7 in this branch. Purity has NOT been
proved to imply (M6). The generic I!=0 branch, the E(0)=0 nonzero-center
cases, and other height directions remain outside this conclusion.

## 5. Purity does not imply the additional identity, even in degree five

The restriction (M6) is not merely unproved as a necessary condition:
it is FALSE as a necessary purity condition inside the actual mixed
quintic box. For any alpha in K, start with the safe pure phase
y(t+alpha), whose angular weight is strictly negative. Apply the
polynomial shear phi_1=exp(xJ). Its literal regular image is

    f_alpha=(y-xt-x^3/2)(t+x^2+alpha).                    (M17)

Section4(a)'s direct shear argument proves that f_alpha is pure AND
safe for every alpha. Its cubic-height numerator coefficients are

    C=(alpha+x^2-alpha*x^4-x^6)/2,
    B=(1-2alpha*x^2-3x^4)/2,
    D=-(alpha+3x^2)/2,   E=-1/2.                       (M18)

They satisfy (M1). The unique degree-five top term of (M17) is
-x^5/2; it cannot be removed by the leading quadric. Equivalently
(M2) gives exact least degree5. Its third derivation is
J^3 f_alpha=-3x^2, so its derivation height is exactly3.

Direct substitution in the cubic invariant gives

    27CE^2-9BDE+2D^3=alpha(9-alpha^2)/4.                (M19)

In particular alpha=1 gives I=2, not zero. This is a safe pure
quintic, not an unsafe Image witness. It rules out using purity
alone to impose (M6), without answering the unsafe degree5/7 problem.

## 6. Exact context and reproduction

The two classical interval sources were read in their actual primary
versions, not inferred from a search-result summary:

- Francois, Pakovich, Yomdin, Zhao, *Moment vanishing problem and
  positivity: Some examples*, Bull. Sci. math. 135 (2011), Corollary 4.1,
  printed p.27. Author-hosted publisher PDF:
  https://www.math.bgu.ac.il/~pakovich/Publications/FPYZ.pdf
- Muger, Tuset, *On the moments of a polynomial in one variable*,
  arXiv:1910.12148v1, Theorem 2.1 and Remark 2.2, printed p.2:
  https://arxiv.org/pdf/1910.12148

Their roles here are exactly the interval vanishing implication in
Section 1. Neither is cited as a finite first-four-moment bound.
Muger--Tuset Conjecture 2.3 concerns an exact limsup formula; it is
not a conjecture about that finite bound.

The elementary Hopf shear, full quadratic classification, and least-degree
normal forms are already retained in the research repository. The
conditional calculations in Sections 3--4 are given literally above;
no all-automorphism classification or global Image conclusion is used.

The standard-library checker
verify_mixed_cubic_hopf.py reconstructs the four
interval moments, performs the two exact remainder contradictions,
checks a nonzero first-three-moment control over Q[b]/(q), verifies
the universal depressed-cubic and endpoint-profile identities, and checks the
safe-quintic counterexample (M17)--(M19). Its
finite controls authenticate the displayed identities; the all-order
endpoint and valuation arguments are the written proof above.

No historical priority, external review, formal verification or independent replay
is asserted.

