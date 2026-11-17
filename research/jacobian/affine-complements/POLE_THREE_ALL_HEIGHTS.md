# No compatible quartic coordinate has pole order three

14 November 2026. PF2RULED8 excludes every target coordinate of B-pole
three at any companion height on the compatible F1 model. Combined
with PF2RULED5, every coordinate in every polynomial target presentation
has pole at least four. The first open reduced height-eight pairs are
(6,8),(7,8).

The argument is intrinsic to the full contact-four boundary family. It
uses the earlier intrinsic low-pole exclusions, not the2022 ordinary-
degree theorem or PF2RULED7. The key is a weighted local differential
criterion, followed by two complete generic-fibre degree checks.

## Initial reduction

Suppose delta_B(P)=3. Repeated divisible-pole target shears lower a
companion pole divisible by three. A constant remainder is impossible
by dominance. The intrinsic [low-pole proof](POLE_REDUCTION.md),
[full(3,4) proof](POLE_FOUR_EXCLUSION.md) and
[full(3,5) proof](POLE_THREE_FIVE_EXCLUSION.md) handle every remaining
companion pole at most five. The next nondivisible pole is at least
seven. Thus assume

n=delta_B(Q)>=7, gcd(3,n)=1.

The leading differential equation gives symbols tau^3,tau^n, after
nonzero target scalings, for a nonzero cubic tau in H0(L),
L=O_B(B)=O_P1(3). Lift tau to T in V1 and put P=T^3+A2, A2 in V2.
All zeros of tau, including infinity and the marked point, are retained.

## A weighted local obstruction

Choose smooth coordinates w,z near a point of B, with B=(z=0), and
suppose ord_w(tau)=k. Give w,z positive integer weights a,b and put
c=b-ak>0. These are local valuation weights, not an assumed descent
of the global source grading. Assume P has weighted initial form

P_in=t^(-3c)*h(Z), under w=t^a,z=t^b*Z,

where h has a pole of order three at Z=0 and has the form

 h(Z)=(A+B*Z+C*Z^2)/Z^3,

with A!=0 and nonconstant numerator. In particular h has a zero whose
multiplicity is one or two, so it is not a cube of a rational function.

The leading B-pole term of Q contains z^-n*w^(kn) with nonzero
coefficient. Its weighted pole K is therefore at least cn. Its initial
form is t^-K*g(Z), where g is a nonzero polynomial in Z^-1 with pole
order at most n. No positive Z powers occur at a negative weighted
value, because Q lies in C[[w,z]][z^-n].

Write the prescribed volume as U*ell*z^-3 dz wedge dw, with U a unit
and ell regular. Its weighted value is at least a-2b. If

\[
 (n+1)c>a(2k-1),\tag{1}
\]

then a-2b>-(n+3)c>=-3c-K. The leading wedge must consequently vanish.
Differentiating the two initial forms gives

\[
 3c\,h g'-K h'g=0,
 \qquad g^{3c}=\lambda h^K.\tag{2}
\]

At Z=0, equation(2) makes the pole order of g equal to K/c. Since it
is at most n whereas K>=cn, we have K=cn. Then g^3/h^n is constant.
At a zero of h of multiplicity one or two, its valuation would give
3*ord(g)=n*ord(h), contradicting gcd(3,n)=1. Thus the stated initial
configuration is impossible.

This is a formal local argument with positive weights; t is an auxiliary
variable in the field extension w=t^a, not a change of source coordinates
or a claim of finiteness. Condition(1) is checked below on the complete
finite multiplicity list, with n>=7 handled by monotonicity.

## First consequence: tau divides the quadratic symbol

Let a2 be the symbol of A2. If h=ord(a2)<k at a zero of tau, use

 a=1, b=3k-h, c=2k-h.

The two leading terms of P are tau^3/z^3 and a2/z^2. All z^-1 and
regular terms have strictly larger weighted value. The initial h(Z)
has a nonconstant linear numerator. Since k<=3 and h<k,
condition(1) holds already at n=7. The local criterion excludes every
such defect. Hence tau divides a2 as a global section.

Lift a2/(3tau) to X in V1 and write the remaining V1 function as3Y:

\[
 P=T^3+3TX+3Y,\qquad T,X,Y\in V1.\tag{3}
\]

Let tau,u,v be their cubic symbols.

## Second consequence: u is proportional to tau

If u were nonzero at any zero of tau, use a=1,b=2k,c=k. The leading
terms in(3) are tau^3/z^3 and3tau*u/z^2, with a nonconstant linear
numerator h(Z). Condition(1) holds. Thus u vanishes at every zero of tau.

If u is not proportional to tau, cubic degrees leave only:

| tau root type | k at its multiple root | h=ord(u) |
| --- | --- | --- |
| double plus simple | 2 | 1 |
| triple | 3 | 1 or2 |

At that multiple root, write the fibre numerator as

F=z^3(P-c0)=tau^3+z*A(w)+z^2*B(w)+z^3*C(w,z,c0),

where ord(A)=k+h and let j=ord(B), including infinity if B=0.
Regular terms in the T lift do not change these leading orders.
Choose a=1,b=2k-h,c=k-h. If j+k-2h>=0, the initial numerator of h(Z)
is linear or quadratic with nonzero linear coefficient. Condition(1)
holds at n=7 for all three rows, so all these cases are excluded.
The only leftover is k=3,h=2,j=0.

## First exceptional fibre: triple tau, u of order two, unit B

For that leftover, a generic P=c0 fibre has local equation with leading
terms w^9+A0*w^5*z+B0*z^2+C0(c0)*z^3, with A0,B0 nonzero and generic
C0(c0) nonzero. Its local Weierstrass degree in z is two, and its unique
normalized branch has

w=t^2, z~t^9, ord_t(F_z)=9.

All points at infinity lie over the sole triple root of tau, so this is
the entire infinity list of the completed generic fibre. The volume gives

\[
 \frac{dQ}{dt}=\frac{U\ell\,w'(t)}{F_z}.\tag{4}
\]

If this boundary point is not the contact point q, the derivative has
order -8 and Q has pole degree seven, contradicting degree four.
If it is q, contact four makes ord_t(ell)=8, since w=t^2 and z~t^9.
The derivative then has order zero and Q is regular at a missing point,
contrary to finiteness. This excludes the last nonproportional-u case.

## Depress the proportional case

We now have X=mu*T+chi. Replacing T by T+mu and adjusting Y, all in V1,
puts P in the form

P=T^3+3Y.

This is a change in the auxiliary lift, not in the actual map. The cubic
symbol tau is unchanged. If v is proportional to tau, Y=nu*T+constant
and P is a polynomial in T, excluded by the
[polar-centralizer lemma](LOW_COORDINATE_POLE.md).

Otherwise some zero of tau has h=ord(v)<k. Use

 a=2, b=3k-h, c=k-h.

The initial form is h(Z)=(A+C*Z^2)/Z^3 with A,C nonzero. The intermediate
terms from T's regular lift have larger weighted value. Condition(1)
holds at n=7 for every k<=3,h<k except k=3,h=2. Thus only that configuration
requires a further fibre calculation.

## Second exceptional fibre: depressed cubic, triple tau, v of order two

Now tau has its only zero at this point. The generic fibre numerator has
leading terms w^9+B0*w^2*z^2+C0(c0)*z^3, with B0 nonzero. It has exactly
two normalized branches:

| Branch | w,z orders in its parameter | ord(F_z) | Q pole away from q |
| --- | --- | --- | --- |
| high | 2,7 | 11 | 9 |
| low | 1,2 | 4 | 3 |

The possible z term from the T lift has coefficient of order at least
six and lies above the Newton faces; it changes neither branch nor the
orders. The low branch is simple for generic c0. Equation(4) gives the
pole orders in the table. Away from q their total is twelve, impossible.
At q, contact four has orders seven and two on the respective branches,
so the Q poles become two and one. Their total is three, again not four.
There are no other infinity points, and the dominant contact terms are
nonzero, so neither numerical conclusion omits a cancellation branch.

## Finiteness and coverage

A generic P fibre is finite over the Q-line of total degree four.
Normalize and complete it, summing degrees over components if necessary.
Every missing point is a pole of Q by finiteness; the sum of pole orders
is four. These facts justify both exceptional-fibre contradictions and
make no assumption of geometric irreducibility of the generic fibre.
The local root types exhaust all cubic divisors, and all arguments are
intrinsic at their points, including the marked point and infinity.

Therefore no compatible quartic target coordinate has pole three at any
companion height. Combining the pole-one/two exclusions gives

\[
 \boxed{\delta_B(P)\ge4,\qquad\delta_B(Q)\ge4
 \text{ in every polynomial target presentation}.}
\]

This is a coordinate bound, not a prohibition on arbitrary low-pole
elements of the target algebra. Alongside the separate complete
height-seven exclusion, the first open reduced height-eight pairs are
(6,8),(7,8); divisible(4,8) lowers to the closed region.
The whole higher-pole and actual ruled-compatibility problem stays open.

The [checker](verify_pole_three_all_heights.py) verifies the universal
weighted ODE, the complete finite cubic inequality list, non-cube
multiplicities, both exceptional Newton models and contact pole totals.
The local-initial coverage, all-height argument, generic-fibre geometry
and finiteness interpretation are written proofs, without formal
verification, independent whole-proof implementation, external review
or a novelty claim.

The adopted frontier(6,8),(7,8) incorporates the subsequent
[PF2RULED12 full(5,8) exclusion](POLE_FIVE_EIGHT_EXCLUSION.md); it is
not an additional exclusion proved by this earlier note.
