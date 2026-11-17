# The nonsquare sextic case of pole four is excluded at every height

14 November 2026. PF2RULED10 excludes the entire (3,3) leading-sextic
case left by [PF2RULED9](POLE_FOUR_COMPANIONS.md), at every companion
height. Consequently a remaining pole-four coordinate has, after the
previous reductions, poles (4,2r), r>=5 odd, and symbols

 tau^4,tau^(2r), tau in H0(L), L=O_B(B)=O_P1(3).

The common sextic is necessarily a cubic square. Its remaining root
partitions are (6),(4,2),(2,2,2). This does **not** exclude those cases,
raise the minimum coordinate pole above four, close height eight or
prove an actual compatible embedding. The first open height-eight
pairs remain (6,8),(7,8).

The proof combines a new finite formal-jet constraint with the actual
quartic generic-fibre degree. It is intrinsic for r>=5; the reduction
to those r uses the earlier height-seven base and its pinned2022
ordinary-degree theorem. No new external theorem is imported.

## Setup and a small-multiplicity improvement

Use the full compatible F1/contact-four geometry, including all marked
positions. Write the prescribed volume locally as

 Omega=U*ell*z^-3 dz wedge dw,

where B=(z=0), U is a unit, and ell is a unit away from the unique
contact point q. The previous theorem leaves symbols sigma^2,sigma^r,
with sigma a sextic having no simple zero, and

 P=S^2+H2, S,H2 in V2.                                   (1)

Here Vm denotes the pole-m space. Let eta be the symbol of H2.
We first prove that at every sigma root of multiplicity k<=3,

 ord(eta)>=k.                                           (2)

Suppose h=ord(eta)<k. In local Laurent expansions write

 S=z^-2*(sigma(w)+z*beta(w)+z^2*gamma(w,z)),
 H2=z^-2*(eta(w)+z*b(w)+z^2*c(w,z)),

with regular coefficients. Put j=ord(beta), allowing j=infinity.
If 2j>=h, use local weights a=2,b_weight=2k-h and scale
c_weight=2*b_weight-a*k=2(k-h)>0. The initial numerator of P is
(A+B Z)^2+C Z^2, C!=0, possibly with a cubic tie. It cannot be a
square: its quadratic discriminant is -4*A^2*C!=0; a nonzero cubic
term also prohibits a rational square. Lower terms cannot dominate,
since 2*ord(b)+2k-3h>=0 for k<=3,h<k. The even-power criterion of
PF2RULED9 applies already at r=5:

 (r+1)*c_weight > a*(k-1).

If 2j<h, the complete list is

 (k,h,j)=(2,1,0),(3,1,0),(3,2,0).

Here beta is a unit. The two close roots of the generic fibre
F=z^4(P-c0) are centered at z~ -sigma/beta, of order w^k.
The nonzero eta*z^2 term is the first splitting term, of w-order
2k+h; all b*z^3 terms have order at least3k>2k+h.
For h=1 there is one cusp with

 ord(w)=2, ord(z)=2k, ord(F_z)=2k+1.

The fibre formula

 dQ/dt = unit*ell*z*w'(t)/F_z                            (3)

gives a regular derivative and hence Q regular at a missing point,
contrary to finiteness. For h=2 there are two smooth branches with
orders (1,k,k+1). Their Q derivatives have a nonzero simple pole if
ell is a unit, impossible for a derivative, and are regular if ell
vanishes, again violating finiteness. Contact can only improve these
regularity obstructions. This proves(2), with every lower Laurent
term retained.

For sigma of type (3,3), or of type (2,2,2), (2) gives eta=mu*sigma
by the equal sextic degrees. Complete the square once more to obtain

 P=S^2+H, S in V2, H in V1.                              (4)

This normal form is also a proved necessary improvement for the still
open (2,2,2) case. For the rest of this proof assume sigma has type
(3,3), so it is not a square in the rational function field C(B).

## Formal centralizer constraint through normal order one

Work in the completion C(B)((z)), using a transverse rational parameter
regular at the boundary point under examination. Since the leading
coefficient of P is sigma^2, there is a formal square root

 W=sqrt(P)=sigma*z^-2+(lower normal poles).

If a Laurent series R has leading term rho*z^-j, the leading
coefficient of dP wedge dR is

 j*(sigma^2)'*rho - 4*sigma^2*rho', at pole j+5.

Whenever j>=-1 this pole exceeds three, so it must vanish if the
wedge has pole at most three. The equation gives

 rho=constant*sigma^(j/2).

Because sigma has an odd zero of order three, a nonzero rational rho
can occur only when j is even. For even j>=0 subtract the corresponding
constant multiple of W^(j/2), which commutes with P exactly, and lower
the leading normal pole. At j=0 subtract a constant. At j=-1 no nonzero
rational leading coefficient is possible. The finite descending
process proves that for the actual Q there is a polynomial f with
constant coefficients such that

 Q-f(W) is in z^2*C(B)[[z]].                             (5)

Its leading term is W^r after the usual nonzero coordinate scaling;
thus f has odd degree r. This is a formal congruence, not a claim that
W is a global regular function, an intermediate plane or a target
coordinate. Its square class and the parity argument are independent
of rational trivializations of the normal bundle.

The symbol v of H is nonzero. Otherwise H is constant and P is a
polynomial in S. The same leading equation forces every positive pole
of Q to be even; subtracting actual powers of S then lowers it while
preserving dP wedge dQ, eventually reaching a constant. This contradicts
the nonzero volume. This is the nonsquare-symbol version of the polar
centralizer argument; it uses no unproved factorization assumption.

## The first forced coefficient puts zeros of v at both sextic roots

Set r=2l+1, l>=2. The formal binomial identity is

 W^r=(S^2+H)^(l+1/2)
     =sum_{i>=0} binomial(l+1/2,i)*S^(2l+1-2i)*H^i.       (6)

The sum through i=l is an actual polynomial in S,H, whose Laurent
coefficients are regular at every point of B. The first negative-S
term, i=l+1, begins at normal exponent z^(1-l), with coefficient

 binomial(l+1/2,l+1)*v^(l+1)/sigma.                      (7)

The subsequent negative-S terms start three normal orders later.
Lower odd powers in f(W) first acquire a negative-S term at strictly
larger normal exponent; even powers are polynomials in P. Since
1-l<=-1<2, (5) and regularity of Q force (7) to be regular. Thus
sigma divides v^(l+1). In particular v vanishes at both triple roots.
As v is a nonzero cubic, its orders at those two roots are each one
or two. All constants in(7) are nonzero in characteristic zero.

## At an away-contact root, quartic finiteness forces beta to be a unit

At least one of the two distinct sigma roots is away from q. Choose
one and normalize its leading coefficient so sigma=w^3 times a unit.
There v has order d=1 or2. Suppose beta has order j>=1. In the actual
local fibre

 F=(sigma+z*beta+z^2*gamma)^2+v*z^3+(h-c0)*z^4,

all gamma,h are regular; put gamma(0,0)=g. The complete Newton census
is determined by j=1 versus j>=2, and d=1 versus d=2. Models retaining
the constant regular lift g are

 F0=(w^3+a*w^j*z+g*z^2)^2+V*w^d*z^3+xi*z^4,

where a!=0 when j=1, V!=0, and xi is a generic fibre parameter.
For j>=2 the j=2 model suffices; larger j only removes higher terms.
The local Weierstrass degree is four and the constant term has w-order
six. The discriminant in z has the following orders and leading terms:

| j,d | disc order | leading discriminant coefficient | forced signed Q pole total |
| --- | --- | --- | --- |
| 1,1 | 15 | -4 V a^3*(-V^2-4 V a g+4 a^2 xi) | 5 |
| 1,2 | 16 | -16 a^4 xi*(V a-xi) | 6 |
| >=2,1 | 16 | -27 V^4 | 6 |
| >=2,2 | 18 | 256 xi^2*(g^2+xi) | 8 |

These coefficients are nonzero for generic c0, including arbitrary g.
For the first two rows the face centered at z~-w^2/a is a double
cluster. Its first splitting comes respectively from V*w*z^3, or
from V*w^2*z^3 together with the free xi*z^4. Higher terms of sigma,
beta,gamma,h change the center at higher orders and cannot lower
these splitting orders. The low branches are simple for generic xi.
For the last two rows the ordinary simple Newton faces give the table
directly. Thus this is coverage of the full local family, not an
inference from selected numerical models. Further cancellations would
raise the discriminant order and do not weaken the contradiction.

To see the pole formula without assuming irreducibility, sum(3) over
all normalized local branches. If (a_i,b_i,s_i) are the orders of
w,z,F_z, their signed pole sum is

 sum(s_i-a_i-b_i)=ord_w(discriminant)-4-6.

Here sum(a_i)=4 by the Weierstrass degree, sum(b_i)=6 by the constant
term, and sum(s_i) is the discriminant order by its norm formula.
The actual sum of positive Q pole orders is at least this signed sum;
a zero signed pole would instead be a forbidden derivative residue.
Every row therefore already contributes more than the entire quartic
degree four. This is impossible at a point where ell is a unit.
Consequently beta is a unit at the chosen away-contact root.

No contact-case classification was used here. At the contact root,
vanishing of ell can reduce these pole totals; applying the away table
there would be invalid. The existence of the other, away root is enough.

## The normal-order-one coefficient contradicts that unit

At the chosen root, sigma has order three, beta is a unit, and v has
order d=1 or2. We prove the following uniform coefficient statement:
for every odd power 2l+1, l>=0, the z^1 coefficient of W^(2l+1) has
w-valuation

 (l+1)*(d-3)<0,                                        (8)

with nonzero leading contribution

 binomial(l+1/2,l+1)*(-1)^l * beta^l*v^(l+1)/sigma^(l+1).
                                                               (9)

Subtract the actual polynomial part i<=l of(6); it contributes no
negative w-valuations. Index a remaining term by i=l+1+k, k>=0.
It has the form

 binomial(l+1/2,l+1+k) *
 z^(1-l+3k) * (v+z*h)^(l+1+k) /
                  (sigma+beta*z+z^2*gamma)^(2k+1).

To contribute to z^1 it needs l-3k further normal orders, so k<=l/3.
The smallest w-valuation comes from using only beta*z in the denominator
and the leading v in the numerator. Its value is exactly

 (l+1)*(d-3)+k*(d+3).

For k>0 this exceeds(8). Replacing a numerator v by a higher normal
coefficient costs at least 3-d>0 in w-valuation compared with the
all-beta choice. Using a denominator term of normal order b>=2
instead of b first-order beta terms costs at least 3*(b-1)>0.
All regular higher coefficients are thus covered. For k=0 the unique
lowest contribution is(9), and its binomial coefficient is nonzero.
This proves(8), uniformly for unbounded l.

Lower odd powers in f(W) have strictly larger w-valuation in their
z^1 coefficients; even powers are actual polynomials in P with regular
coefficients. Therefore the highest odd power W^r has an uncancellable
pole in the z^1 coefficient of f(W). But (5) identifies that coefficient
with the z^1 coefficient of Q, which is regular because Q has only
boundary poles. This contradiction excludes sigma of type(3,3).

## Scope, evidence and remaining problem

The [checker](verify_pole_four_nonsquare.py) retains exact controls for
the formal leading ODE, the complete small-multiplicity atlas, the
universal valuation comparisons, symbolic quartic discriminants and
the sextic partition handoff. Its
[receipt](pole_four_nonsquare_receipt.json) separates exact controls
from the written all-height coverage, local expansions, discriminant
interpretation and finiteness proof. No formal verification, independent
whole-proof implementation, external review or novelty priority is claimed.

The surviving common sextic is a cubic square, with root partitions
(6),(4,2),(2,2,2). In the last case the new small-multiplicity lemma
also gives P=S^2+H, H in V1, with S's leading symbol the square of a
square-free cubic. The remaining cases allow a formal fourth root of
P and hence additional odd terms in the formal centralizer; the
nonsquare parity argument does not transfer to them. The first open
height-eight pairs, unrestricted mixed-weight finite map and actual
ruled-compatibility gate remain unresolved.

The subsequent [PF2RULED11 leading-period theorem](BOUNDARY_LEADING_PERIOD.md)
excludes the entire(2,2,2) case and requires q to be a root in the mixed
(4,2) case. This note and its receipt retain the PF2RULED10 handoff
before that refinement. Current pole-four sextic types are only(6),(4,2),
with the stated contact restriction; the full pole-four layer stays open.

The adopted frontier(6,8),(7,8) incorporates the subsequent
[PF2RULED12 full(5,8) exclusion](POLE_FIVE_EIGHT_EXCLUSION.md); it is
not an additional exclusion proved by this earlier note.
