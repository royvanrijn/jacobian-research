# Pole reduction excludes the first two quartic layers

14 November 2026. PF2RULED2 is a written consequence of the compatible
surface in [PF2RULED1](PROOF.md). It excludes every finite degree-four map
in that setting which has any target-coordinate presentation of maximum
pole order at most three along B. The exclusion holds for the entire
allowed contact-four boundary family.

The new tools are triangular reduction of divisible pole orders and an
odd-degree obstruction to the first nondivisible pair (2,3). The original
ruled-compatibility hypothesis remains essential. This does not exclude
the quartic packets outside that hypothesis or maps of higher pole order.

## The common hypotheses and the filtration

Let S=F1, B=C0+2f, E=C0+3f be smooth irreducible sections, with
B.E=4 supported at one point. Put N=S\B and U=S\(B union E)=A2.
Suppose pi=(P,Q):N -> A2 is finite of degree four and has ramification
divisor R=E intersect N, with coefficient one. Then

\[
 \Omega=dP\wedge dQ,\qquad \operatorname{div}_S(\Omega)=E-3B.
 \tag{1}
\]

For a nonconstant regular function F on N, write
delta(F)=-ord_B(F)>0. A nonconstant function has positive pole order,
since a function with no pole would extend to projective S. Write

\[
 V_m=H^0(S,\mathcal O_S(mB)).
\]

Its elements are regular functions on N with delta at most m. The leading
symbol of a function of exact pole order m is a nonzero section of
L^m, where L=O_B(B) has degree three. In particular

\[
 V_1/\mathbb C=H^0(B,L).
 \tag{2}
\]

Surjectivity in (2) follows from H1(S,O_S)=0 in the restriction exact
sequence. For the P1-bundle F1 -> P1, the pushforward of O is O and its
first higher direct image vanishes; hence this H1 equals H1(P1,O)=0. Also dim V_m=(m+1)(3m+2)/2, giving dimensions 5,12,22,35 at
m=1,2,3,4.

Use z as a local equation for B and w as a parameter along it. If P,Q
have pole orders m,n and symbols p,q, the coefficient of
z^(-m-n-1) dz wedge dw in Omega is

\[
 n p'q-m p q'.
 \tag{3}
\]

When m+n>2, (1) makes this coefficient zero. These are statements on the
function field of B; ratios of sections of the same line bundle are
globally defined rational functions. Their derivative is zero precisely
when they are constant, in characteristic zero.

## Divisible pole orders admit a target shear

Suppose n=k m with k>=1. Formula (3) says

\[
 (q/p^k)'=0,\qquad q=c p^k,\quad c\in\mathbb C^*.
\]

Consequently Q-cP^k has strictly smaller pole order than Q. The target
change (P,Q) -> (P,Q-cP^k) is a polynomial automorphism with Jacobian one.
It preserves finiteness, generic degree, the subalgebra C[P,Q] and the
differential divisor. Interchanging P,Q is harmless as well.

Iterating these steps terminates: the sum of the two positive pole orders
strictly decreases. A transformed coordinate cannot become constant,
since Omega is nonzero. A terminal pair therefore either has orders (1,1)
or has two nondivisible orders, each at least two.

If the two orders are at most two, the only terminal possibility is
(1,1). Two general coordinate fibres then belong to |B|, so their affine
intersection number is at most B^2=3. This contradicts generic degree four.

This proves the complete m=2 exclusion, including all target-automorphic
presentations with that pole bound. It does not rely on a fixed example
of B,E or a search over base-point configurations.

## The first nondivisible pair has cubic leading sections

If the maximum pole order is at most three, the only other terminal pair
is (2,3), after interchange. Let p in H0(L^2) and q in H0(L^3) be its
symbols. Equation (3) gives q^2=c p^3.

After multiplying P,Q by nonzero constants, there is a section
tau in H0(L) such that p=tau^2 and q=tau^3. Indeed the equality of zero
divisors 2 div(q)=3 div(p) implies div(p)=2D, div(q)=3D for an effective
degree-three divisor D. On B=P1, D is the zero divisor of a section of
L. The remaining ratios are constants.

Lift tau to T in V1 using (2), and put A=P-T^2 in V1. Set

\[
 D_2=Q-T^3-\tfrac32 T A\in V_2.
\]

A direct differentiation gives

\[
 dP\wedge dQ=-\tfrac32 A\,dT\wedge dA+dP\wedge dD_2.
 \tag{4}
\]

The first term has pole order at most four. If D2 has exact pole two
with symbol d, vanishing of the pole-five term in the second term forces
d/tau^2 to be constant. Subtracting that constant times P from Q lowers
D2 into V1. If D2 already has pole at most one, no subtraction is needed.

Thus an affine target change puts the pair in the global normal form

\[
 \boxed{P=T^2+A,\qquad
 Q=T^3+\tfrac32 T A+C,\qquad T,A,C\in V_1.}
 \tag{5}
\]

Here A or C may be constant; their leading symbols are then zero. Write
a,c in H0(L) for their symbols, allowing that value.

## A conic in a line bundle of odd degree forces a common zero

For (5),

\[
 dP\wedge dQ=
 -\tfrac32 A\,dT\wedge dA+
 2T\,dT\wedge dC+dA\wedge dC.
 \tag{6}
\]

The pole-four coefficient is

\[
 -\tfrac32 a^2\tau'
 +\tfrac32\tau a a'
 +2\tau\tau'c-2\tau^2c'.
\]

It is -2 tau^3 times the derivative of
c/tau-(3/8)(a/tau)^2. Its vanishing is equivalent to

\[
 \boxed{8\tau c-3a^2=\lambda\tau^2,\qquad\lambda\in\mathbb C.}
 \tag{7}
\]

The three cubic sections tau,a,c must have a common zero on B.
To prove this, choose a zero of tau of odd multiplicity l; one exists
because deg L=3 is odd. Equation (7) implies tau divides a^2 locally,
so ord(a)>=(l+1)/2. Solving (7) for c shows ord(c)>=1 there. The
argument also covers a=0 identically. Thus tau,a,c all vanish at that
point.

Equivalently, a base-point-free triple satisfying the smooth conic
relation (7) would define a morphism P1 -> a conic whose pullback of
O(1) had degree three. That degree would have to be even. The local
odd-multiplicity proof avoids any nonconstant-image qualification.

Oddness is substantive: binary quadratic sections
tau=s^2, a=s*t, c=(3/8)t^2 satisfy (7) with no common zero.

## The common zero creates a forbidden maximal base point

At a common zero x of tau,a,c, let z=0 define B. The local regular
functions

\[
 t=zT,\qquad \alpha=zA,\qquad \gamma=zC
\]

all vanish at x. The three sections defining [P:Q:1] in |3B| have local
expressions

\[
 z^3P=z t^2+z^2\alpha,\quad
 z^3Q=t^3+\tfrac32 zt\alpha+z^2\gamma,\quad z^3.
\]

They all vanish to order at least three. The last has order exactly
three, so x is a proper base point of multiplicity exactly three.

That is impossible under (1). More generally, for a presentation in
|mB|, a proper base point of multiplicity m gives the first exceptional
curve coefficient m-m=0 in the pullback of target infinity. The rational
map is defined at its generic point and maps it into affine target A2.
But the pullback of Omega to that blowup has order

\[
 \operatorname{mult}_x(E)-3+1
 =\operatorname{mult}_x(E)-2<0,
\]

since E is smooth and its multiplicity is either zero or one. A morphism
into A2 cannot pull back du wedge dv with a pole. Later point blowups
cannot change this divisorial contradiction.

This excludes the pair (2,3), and together with divisible-order reduction
proves

\[
 \boxed{\text{No compatible quartic finite map has any
 target-coordinate presentation with }\max(\delta(P),\delta(Q))\leq3.}
 \tag{8}
\]

The proof uses finiteness and degree four in the terminal (1,1) intersection
bound. The (2,3) contradiction itself applies to any dominant regular map
satisfying (1). Smoothness of E is essential in the maximal-base-point
step; replacing it by a singular or nonreduced zero divisor would change
that step.

## Exact equations at the next pole pair

After the reductions and exclusions above, the first possible pair is
(3,4). The subsequent [complete marked atlas](POLE_FOUR_EXCLUSION.md),
PF2RULED4, now excludes that pair as well. The following matrices remain
finite calibrations of the necessary differential equations.

The subsequent [four-section normal form](POLE_FOUR_FRONTIER.md),
PF2RULED3, reduces this pair to21 scalar coefficients and excludes one
explicit cubic-symbol lift family. PF2RULED4 supplies the complete
leading-data coverage and exclusion.

Every B in our completion class can be normalized to the strict transform
of XZ-Y^2=0 on the blowup of P2 at [0:0:1]: projective transformations act transitively on
smooth conics with a marked point. In the affine chart Z=1, set b=x-y^2.
A basis for V_m consists of the quotients

\[
 x^i y^j/b^m,\qquad m\leq i+j\leq2m.
\]

This basis has the dimension above. For calculations modulo constants,
replace x^m by b^m and then remove that constant.

For P=A/b^m and Q=C/b^n the numerator of their differential is

\[
 W=b(A_xC_y-A_yC_x)
 -nC(A_xb_y-A_yb_x)-mA(b_xC_y-b_yC_x),
\quad dP\wedge dQ=(W/b^{m+n+1})\,dx\wedge dy.
\]

Pole order at most three is exactly the divisibility b^(m+n-2) | W
at the generic point of B. It is a finite set of bilinear coefficient
equations. With two different pole bounds it becomes linear on the
coefficient tensor; decomposability means that tensor has rank one.

The retained rational matrices and exact kernels give:

| Pole orders | Section dimensions | Coefficient columns | Linear constraint rank | Kernel dimension |
| --- | ---: | ---: | ---: | ---: |
| (2,2) | 12,12 | 55 exterior coordinates | 19 | 36 |
| (2,3) | 12,22 | 231 tensor coordinates | 33 | 198 |
| (3,4) | 22,35 | 714 tensor coordinates | 70 | 644 |

The first two are calibrations of the layers now excluded by the written
argument. For (3,4), the kernel supplies necessary linear conditions;
rank-one tensors, the exact zero divisor E, contact four, degree four and
finiteness remain to be solved. A linear kernel is not evidence that such
a map exists.

The [symbolic checker](verify_pole_reduction.py) verifies the differential
identities, target shears, representative root partitions, the sharp
even-degree control, a nonzero degenerate-volume example and the saved
linear certificates. Written arguments above supply all-degree and
all-family implications; finite controls do not certify them formally.

## What this changes for actual compatibility

A quartic packet that admits the compatible F1 completion now has pole
height at least five in every polynomially changed target-coordinate
presentation, by PF2RULED4. The first unexcluded reduced possibilities
are now (6,8),(7,8), after the
[source-coordinate transfer](POLE_SEVEN_EXCLUSION.md), the
[complete (3,5) proof](POLE_THREE_FIVE_EXCLUSION.md) and the subsequent
[all-height coordinate bound](LOW_COORDINATE_POLE.md) removes every (2,n).
The old eight-simple-
center m=2 pattern is now excluded by the differential argument, even
before testing proximity or constructing its sections.

The larger goal remains to control higher pole pairs and prove or obstruct
a compatible embedding of the actual finite normalization. The
[Wright connection audit](WRIGHT_CONNECTION.md) identifies a relevant
older affine-line-bundle problem and a concrete multiple-fibre pitfall.
No affine-line-bundle structure preserving the actual omitted prime has
been inferred from an abstract surface isomorphism.


The subsequent [all-height pole-three theorem](POLE_THREE_ALL_HEIGHTS.md)
excludes every(3,n) presentation, independently of the2022 degree bound.
The current first open reduced pairs are (6,8),(7,8).

The adopted frontier(6,8),(7,8) incorporates the subsequent
[PF2RULED12 full(5,8) exclusion](POLE_FIVE_EIGHT_EXCLUSION.md); it is
not an additional exclusion proved by this earlier note.
