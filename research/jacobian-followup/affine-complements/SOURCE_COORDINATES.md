# Explicit affine-plane coordinates and the pole-to-degree transfer

14 November 2026. PF2COORD1 gives polynomial coordinates on U=S\(B union E)
for every forced F1 boundary pair: B=C0+2f, E=C0+3f, with contact four
at one point q. It proves a uniform bound deg(F|U)<=16m for F in
V_m=H0(S,O(mB)), and a sharper bound detecting vanishing at q.
The construction concerns this entire boundary family, not just its
previous positive control. It supplies no finite morphism N=S\B -> A2.

## Two normal forms cover the family

Use Cox coordinates s,t,alpha,beta, of respective classes
f,f,h,h-f on F1, where h=C0+f. The sections are represented by columns
(b1,b2), (e1,e2) of degrees (1,2), (2,3), with no common zero in each
column. Their equations are

Bcal=alpha*b1+beta*b2, Ecal=alpha*e1+beta*e2.

Their determinant vanishes to order four at q. Change base coordinates
so q is s=0 and scale the E column to make the determinant s^4.
The linear b1 is nonzero: otherwise the degree-two b2 would have a zero.
A bundle automorphism adds b1 times a linear form to b2; it reduces
b2 to a nonzero square modulo b1. Scaling gives exactly two cases:

| Position of q | B column | E column |
| --- | --- | --- |
| q is the zero of b1 | (s,t^2) | (s(as+bt), s^3+t^2(as+bt)) |
| q is not the zero of b1 | (t,s^2) | (-s^2+t(as+bt), s^2(as+bt)) |

Here a is arbitrary and b is nonzero. Solving the determinant equation
shows these are all E columns. In either row b!=0 is exactly the
no-common-zero condition for the E column: at s=0 it is necessary,
and at any other zero of as+bt the remaining entry is nonzero.
The first row has q=B cap C0; the second has q different from B cap C0.
Thus these are precisely the two positions relative to the marked
blowup point, with all remaining parameters retained.

## The Cox slice is the actual open complement

The irrelevant locus is (s,t)=(0,0) or (alpha,beta)=(0,0). It does not
meet Bcal=Ecal=1. The Cox torus acts freely here, and the classes
[Bcal]=h+f and [Ecal]=h+2f form a lattice basis. Consequently every
orbit in the inverse image of U has a unique representative with
Bcal=Ecal=1, regularly in both directions. Hence

U is isomorphic to {Bcal=1,Ecal=1} in affine four-space.

This is an isomorphism to the specified open inside S. When S contains
the actual finite normalization compatibly, it preserves that actual U;
it is not a replacement of N or of its boundary valuations.

## Polynomial coordinates in the first chart

Use independent affine coordinates x,y and set

\[
 \begin{aligned}
 s&=x,&\beta&=b^2+xy,& A&=a+\beta x^2,\\
 t&=(1-xA)/b,&
 \alpha&=[\beta(2A-xA^2)-y]/b^2.
 \end{aligned}\tag{1}
\]

All expressions are polynomials over C[a,b,b^-1]. Direct substitution
gives Bcal=Ecal=1. Conversely, on that slice,

\[
 x=s,\qquad
 y=\beta[2(a+\beta s^2)-s(a+\beta s^2)^2]-b^2\alpha.
 \tag{2}
\]

Indeed Ecal-(as+bt)Bcal=beta*s^3, so bt=1-as-beta*s^3.
The Bcal equation then gives beta-b^2=s*y. These identities recover
(1), including s=0; no inversion of s is used. The retained polynomial
identities certify both directions with only b inverted.
The degrees of s,t,beta,alpha are at most1,5,2,11 respectively.

## Polynomial coordinates in the second chart

Set

\[
 \begin{aligned}
 s&=x,&\alpha&=b+abx+x^2y,\\
 t&=(1-ax+\alpha x^2)/b,&
 \beta&=[a^2b-(1-ax)y-\alpha^2]/b.
 \end{aligned}\tag{3}
\]

Again Bcal=Ecal=1. The inverse is

\[
 x=s,\qquad y=a^2\alpha-(1+as)(\alpha^2+b\beta).\tag{4}
\]

Here Ecal-(as+bt)Bcal=-alpha*s^2. After eliminating t, the single
relation is

(1-as)*alpha+s^2*(alpha^2+b*beta)-b=0.

It implies alpha-b=as*alpha-s^2*(alpha^2+b*beta), and then
alpha-b-ab*s=s^2*y. This proves the inverse across s=0 as well.
The degrees of s,t,alpha,beta are at most1,5,3,6.

## Four section generators with controlled degree and contact order

In the first chart put

(f0,f1,f2,f3)=(alpha*t,beta*t^2,beta*s*t,beta*s^2).

In the second put

(f0,f1,f2,f3)=(beta*t^2,beta*s*t,alpha*t,alpha*s).

Together with1 these are a basis of V1 after division by Bcal; on the
slice Bcal=1 this division is already built in. Write H=x^3*y.
For every complex a and b!=0 their exact highest homogeneous terms are:

| Chart | f0 | f1 | f2 | f3 |
| --- | --- | --- | --- | --- |
| first | H^4/b^3 | H^3/b^2 | -H^2/b | H |
| second | -H^4/b^3 | -H^3/b^2 | H^2/b | H |

Their degrees are therefore16,12,8,4 uniformly, without exceptional
parameter specializations. On B their restrictions have vanishing
orders0,1,2,3 at q. More explicitly, with r=s/t and a local frame f0,
the ratios are (1,-r,-r^2,-r^3) and (1,r,-r^2,-r^3).
Their global cubic symbols form a full binary-cubic basis.

## The all-height filtration bound

Every section of mB is a sum of Cox monomials

alpha^(m-j)*beta^j*s^i*t^(m+j-i), 0<=j<=m, 0<=i<=m+j.

Each factors into m sections of B: give one base variable to each alpha
factor and two to each beta factor. Thus Sym^m(V1) surjects onto V_m.
After restriction to U, every product has ordinary degree at most16m:

\[
 \boxed{F\in V_m\quad\Longrightarrow\quad\deg_{x,y}(F|U)\le16m.}\tag{5}
\]

Let F have exact pole m>=1 and leading symbol p in H0(B,L^m),
L=O_B(B)=O_P1(3). If r=ord_q(p), use the basis of binary monomials of
degree3m to lift p by products of the four f_i, each with total
vanishing order at least r and degree at most16m-4r. The difference
between F and that lift lies in V_(m-1). Consequently

\[
 \boxed{\deg(F|U)\le
 \max\{16(m-1),16m-4r\}.}\tag{6}
\]

If r=0, the coefficient of f0^m is nonzero. It is the unique binary
basis term with degree16m; every other lift has degree at most16m-4,
and the lower-pole remainder has degree at most16m-16. Therefore

\[
 \boxed{p(q)\ne0\quad\Longrightarrow\quad\deg(F|U)=16m.}\tag{7}
\]

This is an all-height argument, not an extrapolation from finite tests.
In particular, the ordinary degree detects the actual contact point q,
including when q is the marked point.

## Evidence and limits

The [producer](compute_source_coordinates.py) retains universal
[coordinate and inverse certificates](source_coordinate_certificates.json)
over Q[a,b,b^-1]. The [checker](verify_source_coordinates.py) checks both
Cox equations, both inverse directions, the full section polynomials,
the exact top terms, the contact orders and bounded factorization controls.
Normal-form coverage, the Cox-slice interpretation and all-height
factorization are written arguments, without formal verification,
independent whole-proof implementation or external review.

The [degree-bound application](POLE_SEVEN_EXCLUSION.md) imports the
published plane degree theorem to exclude all compatible quartic pole
heights through seven. Neither this coordinate construction nor that
bounded-height exclusion proves actual ruled compatibility or excludes
the unrestricted quartic packets. The subsequent
[volume-period test](VOLUME_PERIOD_OBSTRUCTION.md) rules out most of this
boundary family at all heights and leaves the
[explicitly marked Wright model](WRIGHT_MARKED_MODEL.md).
