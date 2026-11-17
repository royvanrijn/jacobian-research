# Exact volume forces a regular source coordinate

14 November 2026. PF2RULED45 closes the no-regular-coordinate gap for a
sole generic-index-two omitted prime, in every finite field degree.

**Geometric bridge.** Let N be a smooth affine complex surface containing
the marked U=A2_(x,y), with sole complement R=A1. If div(dx wedge dy)=R
and this two-form is globally algebraically exact on N, then some actual
polynomial source coordinate extends regularly to N.

A finite Keller normalization with one index-two omitted prime meets
all these hypotheses. Thus the regular-coordinate assumption in
PF2RULED38-44 is now automatic there. In degree four, PF2RULED44 then
excludes the whole single-index-two-prime normalization packet. In all
finite degrees the full regular-coordinate affine classification,
marked end, radial recognition and boundary-infinity constraints apply.
Higher-index or multiple-prime packets and higher finite field degrees
are not thereby excluded. The full objective remains active.

The proof is a global canonical/key-form argument. It contains no
search for deeper Laurent coefficients. A nonexact smooth affine
family in Section5 shows why exactness cannot be omitted.

## 1. A canonical bound for every ramification index

More generally let N be normal affine with U=A2 and sole prime R,
and suppose the rational source form has divisor (e-1)R, e>=1. Put
v=ord_R and delta=-v. The normal intersection description PF2RULED36 is

    A=O(N)=C[x,y] intersect {v>=0}.                      (1)

A nonzero element of the height-one ideal of R has delta<0. Mondal's
published sign criterion, Michigan Math.J.66(2017), Theorem1.4(1),
therefore says that the last key value of delta is negative.

Assume every source coordinate has a pole. Mondal's normal-form theorem
(arXiv1307.5577v4, Theorem4.6(1), Definition4.5) supplies coordinates
with a nonweighted key sequence (w_0,...,w_(n+1)), n>=1, such that

    w_0>w_1>0, w_1/w_0 notin {1/k:k>=1}.

The weighted case is impossible: its positive last value would give
nonnegative delta on every polynomial. Write

    w_0=g*a, w_1=g*b, a>b>=2, gcd(a,b)=1.

Set d_i=gcd(|w_0|,...,|w_i|) and alpha_i=d_(i-1)/d_i. The signed key
inequalities and canonical formula, Definition2.8 and Theorem7.2 of
the same primary manuscript, give

    w_(i+1)<alpha_i*w_i,
    e=-w_0-w_(n+1)+sum_(i=1)^n (alpha_i-1)*w_i.         (2)

The canonical formula applies to ANY divisorial semidegree; it does
not require a primitive compactification or positive key values.
Here e is one plus the actual coefficient of R in the source volume
divisor, exactly the convention of that theorem.

Define C_0=-w_0, C_i=C_(i-1)+(alpha_i-1)*w_i and
z_i=w_i-C_(i-1). Then

    z_(i+1)<z_i, z_(n+1)=-e,
    C_1=g*(a*b-a-b)>=g>=1.                             (3)

Let j be the first nonpositive key value. If j=n+1, then C_n=w_(n+1)+e<e.
If j<n+1, then w_j-C_(j-1)>-e, so C_(j-1)<e+w_j<=e.
All earlier values are positive, hence their contributions to C are
nonnegative. In either case

    g*(a*b-a-b)<=e-1.                                  (4)

This is an unbounded index-relative reduction. Since

    a*b-a-b >= a-2,  a*b-a-b >= b-1,

it also gives

    w_0<=3*(e-1), w_1<=2*(e-1).                        (5)

Consequently index one always has a regular source coordinate, without
an exactness assumption. At index two the only possible all-polar
normal-form coordinate pair is (w_0,w_1)=(3,2). At index three the only
possible initial pairs are (3,2) and (6,4). These are necessary initial
pairs, not higher-index exclusions. No estimate on the degrees of
unknown Keller polynomials is being asserted.

## 2. The index-two exception is a complete three-parameter normal form

Return to smooth N and e=2, still assuming all source coordinates polar.
Equality in (4) forces g=1, a=3, b=2 and C_1=1. Thereafter all alpha_i
are one, so C_i=1 and the last key value is -1. Equivalently the essential
key sequence is (3,2,-1); every remaining intermediate value is nonnegative.

Use x for the pole-three coordinate and y for the pole-two coordinate.
BGV Lemma11.1 extends x to P1 with infinity fiber 3R. Normalize its
cyclic base change x=z^3 as in PF2RULED37. There are three lines above R.
On each entire single-boundary piece, the reduced infinity ruling is an
affine-line torsor of degree -5. Keeping the actual source coordinates
z,y, its whole infinity plane has coordinates u=1/z,t and transition

    x=u^-3, y=u^5*t+h(u).                              (6)

The lowest exponent of h is -2, since v(y)=-2. A polynomial shift of t
removes all terms of order >=5; thus h is a Laurent polynomial supported
in [-2,4]. Rescale y to make its leading coefficient one. Affine source
coordinate changes eliminate the coefficients of u^-1,u^0,u^1:

* x -> x+A*y changes the u^-1 coefficient by -2A/3;
* a translation of y removes the constant;
* a translation x -> x+B changes the u^1 coefficient by -2B/3.

Later operations do not reintroduce the previously removed coefficients.
These are genuine source automorphisms, preserving exactness up to a
nonzero constant. The whole normalized chart is therefore

    x=u^-3,
    y=u^-2+a*u^2+b*u^3+c*u^4+u^5*t, a,b,c in C.        (7)

The letters a,b,c here are new coefficients, not the coprime integers
of Section1. Whole-chart existence is inherited from the affine torsor
construction, rather than deduced from a formal expansion alone.

## 3. Exactness and affineness impose opposite conditions on b

In (7) the source volume is -3u du wedge dt. Regular local primitives on
the source plane A2_(z,y) and the infinity plane A2_(u,t) are respectively

    theta_U=-y dx,
    theta_W=3u*t du.

Their difference on Gm x A1 is

    theta_U-theta_W=3*(u^-6+a*u^-2+b*u^-1+c) du.        (8)

If the global form on N is exact, pull back its primitive to the cyclic
cover and restrict to these two planes. The differences from theta_U
and theta_W are closed polynomial one-forms on affine planes, hence
polynomial differentials. Thus (8) is a Laurent differential and has
zero residue. Its residue is 3b, so

    global exactness implies b=0.                     (9)

Now suppose b=0 and define the actual source polynomial

    G=y^3-x^2-3a*y-3c.

Its expansion is 3u*t+terms of order at least two, so v(G)=1. Every
source polynomial has a unique finite expansion

    sum d_(i,j,k)*x^i*y^j*G^k, i,k>=0, 0<=j<3.         (10)

This follows by division by the monic degree-three polynomial G in y.
The leading valuation of a monomial is -3i-2j+k, with leading coefficient
(3t)^k. For a fixed value and a fixed k, the range 0<=j<3 makes i,j
unique. Consequently different terms at the least value cannot cancel:
they have distinct powers of the independent boundary coordinate t.
Hence (1) and (10) give EXACTLY

    A=C[G, X=x*G^3, Y=y*G^2],
    Y^3-X^2-3a*G^4*Y-3c*G^6-G^7=0.                   (11)

In particular every member of A is constant on the whole source curve
G=0. That nonempty affine plane curve contains infinitely many points.
But regular functions of an affine N must separate the distinct source
points of its open embedded plane. Thus N cannot be affine with that
source embedding. This is an actual source-curve contraction, not a
local valuation, smoothness or boundary-jet test.

We conclude

    affineness of the marked plane extension implies b!=0.   (12)

Combining (9) and (12) contradicts the all-coordinate-polar assumption.
The geometric bridge follows, in every field degree. It does not say
that smoothness and K=R alone force a regular coordinate.

## 4. Consequences for actual finite Keller normalizations

For a finite normalization pi:N->A2 of a plane Keller map with sole
index-two omitted prime, PF2OMIT7 gives smooth N and R=A1. The rational
source form equals c0^-1*dP wedge dQ, where P,Q are the two pulled-back
target coordinates, hence it is the global differential

    d(c0^-1*P*dQ).

Its divisor is R. The bridge supplies a regular source coordinate.
After shifting its boundary value, the source has y regular,
v_R(y)=r>=1, and a complementary coordinate with pole q.
Thus PF2RULED38-43 now give the complete geometric normal forms for
EVERY such actual normalization, with no remaining coordinate-negative
case. The exact-volume condition and actual finite target algebra
remain separate: the geometric models are not constructed Keller maps.

The source admits an actual nontrivial additive action preserving the
marked open and R. On U use delta=y*partial_x. For r=1 its boundary
expression is partial_t in the simple-zero chart. For r>=2 it is
z^(2r-2)*partial_t in the finite-value root chart. It is regular,
descends, and is locally nilpotent on A as the restriction of this
locally nilpotent derivation on C[x,y]. Its kernel is C[y]. No
invariance of the actual target algebra C[P,Q] follows.

All these sources have the primitive odd marked end homology of
PF2RULED39/43 and the actual radial identification of PF2RULED44.
Therefore target torus/unknot infinity, coprime target boundary degrees
and every prime target-coordinate boundary degree are excluded in all
finite field degrees in this whole one-index-two-prime setting.
These are degrees on the target boundary, not the source q,r.

In field degree four, PF2RULED44 excludes every regular-coordinate
model. Therefore

    no actual quartic normalization has exactly one omitted prime
    of generic index two.                             (13)

For several omitted primes, BGV Corollary11.4(2) says that U union one
chosen prime R_i is an affine open. At an index-two prime this open is
smooth by PF2OMIT7, and restriction of the source form is exact with
divisor R_i. The bridge gives a regular source coordinate on that
individual affine open. At an index-one prime Section1 gives the same
coordinate-existence conclusion without exactness. The coordinates
may depend on i. Removing the other primes generally loses finiteness,
so (13) cannot be applied to each component to exclude all packets.

The next quartic gates are higher-index or multiple-prime packets.
Higher finite field degrees with one index-two prime remain open,
now entirely inside the regular-coordinate classification and its
all-degree fiber/end restrictions. The original objective is not closed.

## 5. Actual nonexact affine controls for the exceptional normal form

Here is a converse geometric construction for every a,c in C and b!=0.
It confirms that the discarded exactness hypothesis is essential and
that affineness in (12) is an actual global distinction.

On U=A2_(x,y) put

    f=y^3-x^2-3a*y-3c, T=f^2-9b^2*y, X=27b^3*x.

Then U is isomorphic to the finite chart of the hyperelliptic family

    X^2=(f^2-T)^3-243a*b^4*(f^2-T)-2187c*b^6-729b^6*f.   (14)

The inverse is x=X/(27b^3), y=(f^2-T)/(9b^2). Compactify it as a finite
flat double cover of P1_f x A1_T. Its chart at infinity, s=1/f,V=X/f^3,
is

    V^2=(1-T*s^2)^3-243a*b^4*s^4*(1-T*s^2)
                         -2187c*b^6*s^6-729b^6*s^5.     (15)

At s=0 there are two disjoint sections R_+,R_- with V=+1,-1. The total
surface is smooth: on the finite chart it is U, while 2V is a unit at
the two sections. Every projective fiber is integral. Indeed the
right side of (14) is a monic sextic in f with zero f^5 and f^3 terms
and nonzero f term. A square cubic would have zero f^2 and constant
coefficients in its square root, forcing its f coefficient to be zero.
So that sextic is never a square, over any extension of C(T) either.
The generic fiber is smooth of genus two: a finite-chart critical point
of T has f!=0 and x=0, and then satisfies the degree-five equation
6*(y^3-3a*y-3c)*(y^2-a)-9b^2=0. There are only finitely many such points;
the two infinity sections are smooth. The generic squarefree double
cover of P1 has six branch points.

Let N_(a,b,c) be this proper family with R_- removed. The divisor R_-
has degree one on every integral proper fiber. Stacks tags 0B5X and
0D2S make its line bundle relatively ample. A high relatively very
ample power, including a power of the defining section among its
generators, embeds its complement as a closed subscheme of affine
space over A1_T. Thus N_(a,b,c) is affine. It is smooth and has

    N_(a,b,c)=U disjoint-union R_+, R_+=A1.

The exact rational identity

    dx wedge dy=(3b/2)*df wedge dT/X
               =-(3b/2)*(s/V)*ds wedge dT             (16)

shows that its divisor is R_+. On its normalized x-root chart the
expansion is (7), with t|R=T/(18b). Hence its volume residue is 3b!=0,
and it has no global regular primitive. Its source pole orders are
(3,2). Jung--van der Kulk normal form shows that every source coordinate
still has a pole: the initial triangular factor has unequal positive
poles max(3,2*degree p),2, every later nonlinear alternating factor
strictly dominates, and the final affine factor retains one of the
two positive poles. There is no cancellation since 3 is not an even
multiple of 2. Thus these are actual coordinate-negative K=R controls.

Conversely Section2 puts every smooth affine coordinate-negative K=R
plane extension into (7), and Section3 forces b!=0. The construction
has the same marked source valuation. Equation(1) then identifies the
affine surface itself with N_(a,b,c). This is a complete parametrized
normal form with chosen normalizing source coordinates; no uniqueness
of the triple modulo all source automorphisms is claimed.

## 6. Primary inputs and assurance

The new primary key-form imports are Mondal's published sign criterion
Theorem1.4(1), and his 1307.5577v4 signed key inequalities, normal-form
Theorem4.6(1) and general canonical formula Theorem7.2. Complete relevant
PDF pages were rendered and visually inspected. These statements apply
to general divisorial semidegrees; no invalid passage from affine
partial extension to primitive compactification is made. Retained BGV
and PF2RULED37 provide the whole cyclic torsor charts. The small
relative-ampleness imports are used only to prove the actual nonexact
control family is affine. Precise bytes and limits are in
[sources/coordinate-negative-canonical/manifest.json](sources/coordinate-negative-canonical/manifest.json).

[verify_exact_coordinate_bridge.py](verify_exact_coordinate_bridge.py)
writes [exact_coordinate_bridge_receipt.json](exact_coordinate_bridge_receipt.json).
Its symbolic checks support the general defect inequality, full
three-parameter residue and contraction identities, whole affine
control family, smoothness certificate, inverse substitutions and
additive-action boundary formulas. No bounded experiment proves the
unbounded theorem. Primary inputs are imported; formal verification,
independent full reproof, external review and novelty priority remain
false. There is no finite Keller construction or full-goal completion.
