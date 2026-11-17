# Disjoint images reduce to two nodal fibres of one degree-ten polynomial

14 November 2026. PF2OMIT16 strengthens the intersection obligation
of PF2OMIT15 beyond odd degree parameters. It determines a coordinate
ring basis for every known gap21 image and reduces any disjoint
additional image to a specific, still unresolved critical-value gate.
It does not construct a Keller normalization or prove properness.

## 1. Source scope and the elementary disjoint-curve lemma

The [versioned Jelonek text](https://arxiv.org/html/2011.03472v3),
Theorem1.2, assumes that the entire nonproperness set is smooth.
Its connected-set Euler conclusion retains that hypothesis. Our image
union is singular, so that conclusion cannot be applied here. The
[version audit](https://arxiv.org/abs/2011.03472) rejects v4/v5; the
selected v3 and history are [pinned](sources/disjoint_image_manifest.json).
No assertion from the rejected versions is used.

The proof of Theorem1.1 in v3 contains an elementary observation which
does apply. Let A,B be distinct reduced irreducible plane curves with
polynomial normalizations A1, and let H_A,H_B be defining equations.
If A and B are disjoint, H_B restricted to A's normalization is a
nowhere-zero polynomial, hence a nonzero constant c. Thus H_A divides
H_B-c. Reversing A,B proves both degree inequalities. Consequently

    deg A=deg B, H_B=alpha*H_A+beta,
    alpha!=0, beta!=0.                              (D1)

Conversely two distinct such fibres are disjoint whenever both are
irreducible curves. This lemma requires no smoothness or immersion.

For the known image C of degree10 and any additional polynomial image
B of degree5k, disjointness therefore forces k=2. In monic V^6
normalization it forces H_B=H_C-lambda, lambda!=0. Every other k,
including even k>=4, forces an affine intersection. The same rule
applies to a distinct unramified image H.

## 2. The complete affine coordinate ring and its pole degrees

Write A_C=C[p(t),q(t)] inside C[t]. The infinity calculation of
[PF2R7TOP1](../jc2-r7-complement/PROOF.md) gives delta_infinity24.
The rational degree10 projective curve has arithmetic genus36, so

    dim_C(C[t]/A_C)=delta_affine=12.                (D2)

The same source, PF2R7B5CONTACT1, gives an explicit polynomial G_C(U,V)
whose restriction g(t)=G_C(p(t),q(t)) has degree9 at every coefficient
point in both complete charts, with nonzero leading coefficient.
Consider the six elements

    1, q, q^2, g, q*g, q^2*g.                      (D3)

Their degrees are0,10,20,9,19,29, all distinct modulo6. Their C[p]
span M is therefore free with this basis, and all of its displayed
polynomial multiples have distinct leading degrees. Its codimension
in C[t] is12: the minimum degrees in residues0,...,5 modulo6 are
0,19,20,9,10,29, whose sum is87, and (87-15)/6=12. Since
M is contained in A_C and both have codimension12, M=A_C.

Thus (D3) is a C[p]-basis of the full affine coordinate ring, at every
coefficient point, and its exact pole semigroup at infinity is

    S_affine=<6,9,10>, conductor24,
    gaps={1,2,3,4,5,7,8,11,13,14,17,23}.           (D4)

This is a semigroup of polynomial-function pole orders on the affine
curve. It is different from the local infinity-branch semigroup
<4,10,41>. In particular every nonconstant restriction of a target
polynomial to C has degree at least6. Uniqueness of the C[p] normal
form also gives an exact interface for restricting a proposed image
equation without interpreting repeated roots as extra boundary points.

For any actual B of degrees(3k,5k), the norm argument of PF2OMIT15
gives deg H_B(p,q)<=30k for all k; its upper-bound part needs no parity.
Combining with (D1)--(D4), either

    k=2 and H_B(p,q) is a nonzero constant, or
    deg H_B(p,q) in S_affine intersect[6,30k].      (D5)

If k is odd, the stronger9k lower bound of PF2OMIT15 remains valid.
For k=1 the later [PF2R7B5RAT1](../jc2-r7-quintic-normalizations/PROOF.md)
improves it to10 and excludes contact41; together with the exact Newton
space its necessary degree list is{10,12,16,18,20,22,24,26,28,30}.
These imports are not replayed as new coefficient classifications here.

## 3. Every fibre H_C=lambda has the same infinity branch

Fix a coefficient point and put h=H_C, monic of V-degree6. Over U=s^6
its six roots are q(t(zeta*s)), zeta^6=1. For a fixed root, four root
differences have highest s-exponent10; the remaining one, paired by
zeta=-1, has highest exponent-11 with nonzero coefficient. Therefore

    h_V(s^6,q(t(s)))=unit*s^29*(1+O(s^-1)).        (D6)

To solve h=lambda near that root, substitute V=q(t(s))+s^-29*epsilon.
The term linear in epsilon has a nonzero constant leading coefficient.
Every term of degree at least two in epsilon has negative s-exponent:
for the quadratic term the largest is40-58=-18, and the higher terms
are smaller. The formal implicit function theorem therefore gives

    V_lambda(s)=q(t(s))+O(s^-29).                 (D7)

The same holds at each of the six roots. These are algebraic convergent
Puiseux series and depend analytically on lambda near any bounded
lambda value. Their first odd exponent remains-11. The cyclic action
on all six branches is still transitive, since gcd(6,10,11)=1. Hence
every affine fibre h=lambda is reduced irreducible with one place at
infinity, of characteristic(4;10,31) and delta24. Any additional global
component would supply another infinity branch, which is impossible.

The same expansion is uniform over a small coefficient-parameter
neighbourhood. The marked Whitney calculation of PF2R7TOP1 applies
at this unchanged infinity branch: its coordinates start r^4,r^10
and the first odd characteristic coefficient remains a unit. Thus
there is no bifurcation at infinity at any finite lambda. The gradient
has only finitely many affine zeros, since a common curve of its two
partials would force a repeated fibre in characteristic zero.

A generic fibre is smooth of genus12 with one puncture, so chi=-23.
Euler integration over lambda, with the local isolated plane-curve
Milnor formula at each singular fibre, now gives

    sum_(critical points y) mu_y(h)=24.           (D8)

This is the usual polynomial-fibration Euler calculation. Smooth
infinity comparison is needed here; no positivity statement about a
singular connected nonproperness set is being substituted for it.

## 4. A disjoint polynomial image forces two twelve-node fibres

If h=lambda is a polynomial image B, its projective normalization is
P1, so (D7) gives delta_affine(B)=12. For any reduced plane curve,
mu_y=2delta_y-r_y+1. Since delta_y>=r_y-1, with equality at a singular
point only for a node, its total affine Milnor number is at least12,
with equality precisely for twelve ordinary nodes.

Both C=h^-1(0) and a disjoint B=h^-1(lambda), lambda!=0, have this lower
bound. They exhaust (D8). Consequently

    C and B each have twelve ordinary nodes;
    h has exactly two critical values, 0 and lambda. (D9)

In particular **no additional polynomial image can be disjoint from
any nonnodal first-image coefficient point**. This covers all268
nonnodal points of the complete L/K catalogue, not just the immersed
ones. A disjoint B must itself be nodal, so the excess-one cusp on B
in oriented row B also rules out disjointness, regardless of the
first-image parameter. It does not assert that the two nodal fibres
in (D9) exist.

## 5. Exact critical-value traces exclude disjointness generically

Let J_h=C[U,V]/(h_U,h_V), of dimension24 by (D8). Let T_h be
multiplication by h on this finite algebra. Its eigenvalues are the
critical values, with multiplicities equal to their local Milnor
numbers; nilpotent local parts do not change traces of powers. If a
disjoint B exists, (D9) therefore forces

    tr(T_h)=12lambda, tr(T_h^2)=12lambda^2,
    M(C):=12tr(T_h^2)-tr(T_h)^2=0.                (D10)

Also tr(T_h)!=0. These are necessary conditions only. Over C, a
vanishing quadratic moment does not imply equal critical values.

The [checker](verify_disjoint_images.py) reconstructs the exact h at
d=1 on L and x=1 on K. For each it verifies the gradient reductions,
all Buchberger S-pairs, a24-element standard-monomial basis, and the
full multiplication matrix. The ideal generated by the gradient is
contained in the displayed Groebner ideal; (D8) and its checked quotient
dimension24 imply equality. Thus no unproved solver ideal membership
is needed in the reverse direction. Both exact values of M(C) are
nonzero and are recorded with the matrices in
[disjoint_image_ledger.json](disjoint_image_ledger.json).

These witnesses imply a nonempty base open in each complete family
where no disjoint additional image exists. Indeed the critical algebras
have constant dimension24; isolated-zero conservation around the
central critical points accounts for all24 nearby points with their
lengths, preventing any escape. They form a finite flat family locally.
The traces are consequently algebraic functions of the coefficient
parameter, regular near these witnesses. Viewed over C(d) or C(x),
M is a nonzero rational function. Its zeros and possible poles comprise
a finite set on A1 or Gm. Thus only finitely many nodal coefficient
parameters can still admit a disjoint image. The explicit finite
zero set has not been calculated. The two point calculations do not
exclude every nodal parameter or construct any cover.

A capped symbolic trace attempt is retained as discovery only. Its
printed dimension24 and incomplete trace output are not the basis
of the witness computations or the generic argument.

## 6. The full image-support connectedness gate

A finite union of distinct irreducible polynomially parametrized plane
curves is either connected or consists of pairwise disjoint fibres of
one polynomial, all of the same degree. To prove this, choose A,B in
different connected components. By(D1), H_B is affine linear in H_A.
Any further component D is disjoint from at least one of A,B, so(D1)
also makes H_D affine linear in H_A. Distinct such curves are distinct
fibres and hence pairwise disjoint. This proves the assertion without
assuming that a common infinity point forces an affine intersection.

In our ledger, a disconnected nonproperness support would therefore
consist of polynomial fibres of h=H_C. Each has affine Milnor number
at least12. Budget(D8) permits at most two such fibres. Hence:

* Row A with three distinct images always has connected support.
* Rows B and N always have connected support, because one image has
  an ordinary cusp and(D9) requires both disjoint images to be nodal.
* Shared-image row A and row C can have disconnected support only
  in the remaining two-fibre case(D9)--(D10): both images degree10,
  each with twelve nodes, and the trace moment zero.

This is an exhaustive connectedness reduction for the four rows.
It neither constructs the last two-fibre case nor excludes a connected
singular nonproperness set by the smooth-set theorem of Section1.

## 7. Consequences for the omitted-prime ledger

In shared-image row A and in the index-three row C, every C/B
intersection is a singleton value by PF2OMIT15's length calculation.
Such values now necessarily exist in each of these situations:

* deg B!=10, with no oddness restriction;
* C is any nonnodal coefficient point, including when deg B=10;
* C lies in the nonempty trace-certified open of either nodal family,
  again allowing deg B=10.

When an intersection occurs, its weighted total is at least6 by(D4),
at least9k for odd k, and at least10 for k1. These are not counts of
distinct values. The shared A omission packet is2+2+1, all f1, plus
one retained point; its length-one boundary factor is unramified.
Row C has omission packet2+3 plus one retained point. Ramification,
omission, target-image support and tangency multiplicity remain separate.

The entire support is connected in the cases of Section6, in particular
at any nonnodal C and on the certified nodal opens. Rows B,N still need
their cusp incidence to force a singleton; mere connectedness supplies
no such implication. The new high-contact quintic incidence results retain
their separate degree/contact hypotheses and global-cover gap.

Actual B,H equations and fibre values, infinity attachments, a connected
S6 covering and identification of the original source A2 remain missing.
The finite exceptional nodal trace gate and higher-degree matching
have not been solved. No singleton-zero properness or whole-r7 result
is claimed.
