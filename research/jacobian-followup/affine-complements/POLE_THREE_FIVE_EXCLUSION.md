# Complete exclusion of the quartic pole pair three and five

14 November 2026. PF2RULED6 excludes (3,5) under the compatible F1
completion of PF2RULED1. This uses the actual finite degree four and
volume divisor E-3B, not merely an abstract A2 complement.
The subsequent [source-coordinate transfer](POLE_SEVEN_EXCLUSION.md)
excludes every height through seven using an imported degree theorem.

Write L=O_B(B)=O_P1(3), V_m=H0(S,O(mB)), and let q be the unique
contact-four point E cap B. The marked point p0 of the plane blowup
model is a separate point; no equality q=p0 is assumed.

## A local formula and the finiteness requirement

Take a generic fibre P=c. Every point added in its smooth projective
completion must be a pole of Q: the affine fibre is finite over A1_Q,
so its completion adds points only over infinity. This holds component
by component and does not assume irreducibility. The sum of Q's pole
orders on the completed fibre is four.

Near B use z=0 and u along B, and write

\[
 dP\wedge dQ=U\ell z^{-3}dz\wedge du,
 \qquad F=z^3(P-c),
\]

where U is a unit and ell is an equation of E. On a normalized branch
u=u(t), z=z(t), implicit differentiation gives the particularly useful
identity

\[
 \boxed{\frac{dQ}{dt}=\frac{U\ell\,u'(t)}{F_z}.}\tag{1}
\]

Put a=ord_t(u), s=ord_t(F_z). Away from q, a pole of Q has order
s-a. At q subtract ord_t(ell) from that number. A resulting number
zero would give a nonzero simple-pole differential, impossible for
a meromorphic derivative. A negative number gives regular Q, contrary
to finiteness at a missing point. If z has order e, contact four gives
ord(ell)=min(e,4a) unless e=4a, when the order can only increase.

For a Weierstrass polynomial of local degree D in z, the sum of the
s-a values is ord_u(discriminant)-D. We use this only where each
individual value is positive. Discriminant orders also follow by
summing twice the valuations of the pairwise root differences.
All uses below are in characteristic zero and at generic c.

## The first global normal form

Scale leading symbols to tau^3,tau^5, with tau a nonzero cubic section.
Choose T in V1 lifting tau and put P=T^3+A2, A2 in V2, with symbol a2.
The first truncated power has the exact differential

\[
 dP\wedge d(T^5+\tfrac53 T^2A2)
 =-\tfrac{10}3TA2\,dT\wedge dA2.
\]

The remainder in Q has pole at most four. Its pole-four symbol is
kappa*tau^4, since its prospective pole-eight differential must vanish.
Thus Q=T^5+(5/3)T^2A2+kappa*T^4+D3 with D3 in V3.
The pole-seven equation integrates to

\[
 \overline{D3}=\frac{5a2^2}{9\tau}
       +\frac{4\kappa}3\tau a2+\lambda\tau^3.\tag{2}
\]

For example its differential coefficient is

\[
 -3\tau^6\left(
 \frac{\overline{D3}}{\tau^3}
 -\frac{5a2^2}{9\tau^4}-\frac{4\kappa a2}{3\tau^2}
 \right)'.
\]

Subtract lambda*P from Q. Regularity in (2) forces tau to divide a2^2.

In fact tau divides a2. The only failures of that stronger divisibility
are a double zero k=2 of tau with h=ord(a2)=1, or a triple zero k=3
with h=2. These follow from k<=2h and the fact that tau is cubic.
Locally the fibre numerator has the form

\[
 F=\tau(u)^3+z\,a2(u)+z^2 b(u)+z^3 C(u,z,c),\tag{3}
\]

up to the regular terms from the lift of T. Those terms have higher
orders than the displayed balances. C(0,0,c) is nonzero for generic c.
Both failures have a smooth branch z~u^(3k-h), with
ord(F_z)=h. By (1), away from q its differential order is -h.
For h=1 this is a forbidden residue. At q, contact adds order four,
so Q is regular, contradicting finiteness. This excludes the double case.

For the triple case h=2, that branch gives a simple pole of Q away
from q and regular Q at q. If q is elsewhere, put j=ord(b).
For j=0 the other branch has z~u^2 and a simple pole of Q; for j>=1
there are two other smooth branches z~u, each with a simple pole.
Their leading quadratic has distinct roots for generic c, because its
discriminant depends nontrivially on c. Thus the total degree is two or
three, not four. This proves the stronger divisibility.

Lift a2/(3*tau) to X in V1, and write the remaining V1 function as 3Y.
Every candidate therefore has

\[
 P=T^3+3TX+3Y.
\]

Let tau,u,v be the cubic symbols of T,X,Y. Define

\[
 \begin{aligned}
 G5&=T^5+5T^3X+5T^2Y+5TX^2+10XY,\\
 G4&=T^4+4T^2X+4TY+2X^2.
 \end{aligned}
\]

Equation (2) now puts Q=G5+kappa*G4+C2, C2 in V2.
The pole-six equation integrates to

\[
 \overline{C2}=-\frac{5u^3}{3\tau}+\lambda\tau^2.\tag{4}
\]

Indeed the highest exterior coefficient of dP wedge dG5 is
15*T*X^2*dT wedge dX; together with dP wedge dC2 it gives
-3*tau^5 times the derivative of C2bar/tau^2+5u^3/(3tau^3).
Thus tau divides u^3. Choose any W in V2 with symbol u^3/tau.
Surjectivity of V2 -> H0(L^2), with kernel V1, gives the complete form

\[
 \boxed{Q=G5+\kappa G4-\tfrac53 W+\lambda T^2+C1,
 \qquad C1\in V1.}\tag{5}
\]

## Proportional u, including squarefree tau

If u=mu*tau, then X=mu*T+x0. Replacing T by T+mu and adjusting Y
puts P=T^3+3Y with T,Y in V1. Reconstructing (5) in this presentation
has u=0 and gives

\[
 Q=T^5+5T^2Y+\kappa(T^4+4TY)+\lambda T^2+C1.
\]

Its pole-five equation is C1bar=5v^2/tau+rho*tau, so tau divides v^2.
If v is proportional to tau, P is a polynomial in T, excluded by the
[polar-centralizer lemma](LOW_COORDINATE_POLE.md). Otherwise the cubic
multiplicities give just two possibilities:

| tau roots | ord(v) at the multiple root | Q poles away from q |
| --- | --- | --- |
| double plus simple | 1 | 5,1 at the double root; 1,1,1 at the simple root |
| triple | 2 | 9,3 |

These are Newton balances of tau^3+3v*z^2+C*z^3. In the first row,
the high branch has u=t^2,z~t^5 and the low branch z~u; in the second,
the high branch has u=t^2,z~t^7 and the low branch z~u^2.
The generic fibre degrees are nine and twelve when q is elsewhere.
Putting q at the double or simple root in the first row gives a forbidden
residue on a low branch. Putting q at the triple root in the second row
subtracts seven and two from the pole orders nine and three, leaving
total degree three. Every case is excluded.

In particular squarefree tau is impossible: (4) already makes u
proportional to tau. The zero section u is included above.

## Triple-root nonproportional u

Suppose tau has a triple root. Divisibility tau|u^3 makes u vanish
there, and nonproportionality leaves h=ord(u)=1 or2. In (3) the
coefficient of z has order 3+h. Put j=ord(b), with infinity allowed.
The complete generic Newton census away from q is:

| h | j | Q pole orders at this root |
| --- | --- | --- |
| 1 | 0 | 3,3 |
| 1 | 1 | 3,3,1 |
| 1 | >=2 | 3,3,3 |
| 2 | 0 | 7 |
| 2 | 1 | total at least9 |
| 2 | 2 | 4,4,3 |
| 2 | >=3 | 4,8 |

For h=1 there is always a distinguished smooth branch z~u^5 with
pole three away from q. For h=2,j=0, the branch is u=t^2,z~t^9.
For h=2,j>=2 there is a distinguished smooth branch z~u^4 with pole four.
All other entries follow directly from the successive Newton slopes.

The potentially repeated face h=2,j=1 needs care. There is a low smooth
branch z~u, with pole one, and two roots of order at least four in the
remaining quadratic cluster. Their mutual difference has order at least
four; each differs from the low root by order one. The discriminant
order is therefore at least12, giving total pole degree at least9.
Coalescence can only raise this bound; no generic-separability assumption
on that leading face is needed.

If q is elsewhere, every row exceeds four. If q is this triple root:
for h=1 the distinguished pole three loses four and Q becomes regular;
for h=2,j=0 the pole seven loses eight and Q becomes regular;
for h=2,j=1 the low branch gives a forbidden residue;
for h=2,j>=2 the distinguished pole four loses at least four, giving
regular Q or a forbidden residue. Finiteness excludes every contact case.

## Mixed roots: the next coefficient prevents the only dangerous collision

The remaining data have

\[
 \tau=L^2M,\qquad u=LMN,
\]

with distinct L,M and N not proportional to L. The freedom T->T+c
changes u to u-c*tau when P is reconstructed, so write N=n1*M with
n1 nonzero. No normalization of n1 is needed. Every cubic v has the form

\[
 v=\delta n1^2M^3+LD,
\]

where D is an arbitrary binary quadratic. In (5) choose
W=lift(n1^2*L*M^2)*lift(n1*M^3). All choices of that lift differ by V1
and are absorbed in C1.

At the double root the fibre's leading quadratic, after z=u_local^3 Z,
has discriminant proportional to n1^2*(9-12delta). If delta is neither zero nor3/4,
there are exactly two smooth branches z~u_local^3, each with Q pole two
away from q. These statements hold with all constant lifts.

The next pole equation forces delta!=3/4. To check it over the entire
marked geometry, use the conic model and the V1 lift formula from
[PF2RULED4](POLE_FOUR_EXCLUSION.md):

\[
 \widetilde f=wf(w)/z+f_2w+2f_3w^2+f_3z+c.
\]

Automorphisms fixing p0 give precisely the three charts
(L,M)=(s,t),(t,s),(s-t,t). Keep every coefficient of D, all constants
in T,X,Y and all nonconstant coefficients of C1, and kappa,lambda.
The constant in C1 only translates Q and has zero differential. Let
F5(w) be the coefficient of
z^-5 in dP wedge dQ. Direct quotient differentiation gives

| Marked chart | Functional of F5 | Exact value |
| --- | --- | --- |
| double at p0 | coefficient w^7 | -20*n1^4*(3delta^2-3delta+1) |
| simple at p0 | coefficient w^11 | 20*n1^4*(3delta^2-3delta+1) |
| neither at p0 | F5'''(1) | -120*n1^4*(3delta^2-3delta+1) |

Thus f(delta)=3delta^2-3delta+1=0. The literal identity

\[
 16f(\delta)-(12\delta-3)(4\delta-3)=7
\]

shows delta!=3/4 over both complex roots. Also f(0)=1, so delta is
nonzero and the leading quadratic really has degree two. The
[full coefficient tables](pole35_repeated_root_equations.json) retain
all free parameters; the checker reconstructs them by a separate exterior
calculation. These are identities, not parameter samples.

If q is not a root of tau, the double root already contributes degree
four and the simple root contributes at least one, exceeding four.
If q is the double root, each of its two branches loses contact order
three from its pole two, making Q regular, contrary to finiteness.

Finally put q at the simple root. Because u vanishes there, if b is a
unit there is one branch u_local=t^2,z~t^3, with Q pole one away from q;
contact makes Q regular. If b vanishes, generic c gives three distinct
smooth branches z~u_local, each with pole one; contact gives forbidden
residues. The dependence on c of the leading cubic excludes repeated
slopes for generic c. This exhausts all positions of q.

## Status and assurance

No compatible finite quartic map has reduced pole pair (3,5). Combining
PF2RULED2, PF2RULED4 and PF2RULED5 initially left (4,5). The subsequent
[PF2RULED7 proof](POLE_SEVEN_EXCLUSION.md) excludes every height through
seven. The first open reduced pairs are (6,8),(7,8).
Higher poles and actual normalization-to-completion compatibility
remain open. Neither full quartic packet nor JC2 is excluded.

The [checker](verify_pole35.py) checks the global differential identities,
Newton-model discriminants, pole/contact tables and all three complete
mixed-chart coefficient functionals. The global normal-form coverage,
generic-fibre geometry, repeated-face lower bound and finiteness argument
are written proofs, without formal verification, independent whole-proof
implementation, external review or a novelty claim.

The subsequent [all-height pole-three theorem](POLE_THREE_ALL_HEIGHTS.md)
excludes every(3,n) presentation, independently of the2022 degree bound.
The current first open reduced pairs are (6,8),(7,8).

The adopted frontier(6,8),(7,8) incorporates the subsequent
[PF2RULED12 full(5,8) exclusion](POLE_FIVE_EIGHT_EXCLUSION.md); it is
not an additional exclusion proved by this earlier note.
