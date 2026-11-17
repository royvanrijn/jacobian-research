# Source-coordinate transfer excludes every compatible pole height through seven

14 November 2026. PF2RULED7 excludes every target presentation of pole
height at most seven for a compatible finite quartic normalization.
This includes the remaining pair (4,5) and every pair at height six
or seven, over the full contact-four boundary family.

## The imported theorem and its exact scope

Use Guccione--Guccione--Horruitiner--Valqui,
[Increasing the degree of a possible counterexample to the Jacobian
Conjecture from 100 to 108, arXiv:2204.14178v1, Theorem2.1](https://arxiv.org/html/2204.14178v1).
For a noninvertible complex plane Keller pair, it gives this alternative:

maximum ordinary degree is at least125, or the two degrees are72,108
in either order.

The theorem is an external input. Its full Newton reduction and prior
case exclusions are not independently reimplemented here. The
[source manifest](sources/degree-transfer-manifest.json) pins the exact
TeX and PDF already retained in this repository. We do not use the
conditional2026 exclusion of(72,108) or promote the125 frontier to an
unconditional imported theorem.

## Why the degree theorem applies to the actual restricted map

Under compatible completion, N=S\B and its actual source open is
U=S\(B union E). The [explicit source coordinates](SOURCE_COORDINATES.md)
give an isomorphism phi:A2_(x,y) -> U on this specified open, in both
possible marked positions and for every boundary parameter.

Restrict the given finite map pi=(P,Q):N -> A2 and compose with phi.
The two resulting functions are polynomials. Since dP wedge dQ has
divisor E-3B, its restriction to U has no zeros or poles. Its Jacobian
in x,y is a nonzero polynomial unit, hence a nonzero constant.
Moreover the generic field degree stays four: restricting to a dense
open and changing source coordinates preserve the function field.
Thus the resulting Keller pair is not a polynomial automorphism.

This comparison uses the compatible embedding of the entire actual N
and its actual U. An abstract ruled surface birational to N would not
suffice. No general compatibility theorem is inferred from the formulas.

## Heights at most six

Write m=delta_B(P), n=delta_B(Q). Both are positive, since a nonconstant
function regular on projective S is impossible. The uniform bound gives

deg(P|U)<=16m, deg(Q|U)<=16n.

If max(m,n)<=6, their maximum ordinary degree is at most96. This
satisfies neither alternative of the imported theorem. Thus the entire
height-six layer and all lower heights are excluded.

## Height seven: retain the72,108 exception rather than assuming it away

If m=n=7, the leading differential equation makes their leading ratio
constant. Subtract a scalar multiple of one target coordinate from the
other to lower its pole strictly. This is a target automorphism preserving
finiteness, degree and volume. The lowered function remains nonconstant
by dominance. It is therefore enough to consider 1<=m<7, n=7.

Since7 is prime, gcd(m,7)=1. The leading differential equation
n*p'*q-m*p*q'=0 implies q^m=c*p^7. Divisors on B=P1 and the absence
of Picard torsion give, after nonzero target scalings,

p=tau^m, q=tau^7, with a nonzero cubic tau in H0(L).

If tau(q_contact) is nonzero, the sharp degree formula makes
ordinary deg(Q|U)=16*7=112. Meanwhile deg(P|U)<=96, so the maximum
ordinary degree is112. It is below125 and unequal to108, contradicting
the published alternative.

If tau(q_contact)=0, the leading symbol of Q vanishes there to order
at least7. The sharper filtration bound gives

\[
 \deg(Q|U)\le\max\{16\cdot6,16\cdot7-4\cdot7\}=96.
\]

Again both ordinary degrees are at most96, a contradiction. This
exhausts the whole height-seven layer, including equal initial poles.
The surviving published72,108 pair has not been excluded in general;
it cannot occur in these bounded pole presentations.

## The next actual frontier

Every compatible quartic map now has pole height at least eight in
every polynomial target presentation. The all-height pole-two exclusion
still applies. Divisible poles can be reduced by target shears, so the
first open reduced height-eight pairs are

(6,8),(7,8).

The pair(4,8) lowers to the already closed height-seven region. Equal
pole-eight presentations lower to one of the displayed possibilities
or to a previously excluded case. No height-eight layer is excluded
by this note, and no claim is made about arbitrary high poles.

The [checker](verify_source_coordinates.py) verifies the coordinate
identities, top degrees, contact filtration and the exact numeric
transfer, with a [receipt](source_coordinates_receipt.json). Normal-form
coverage, all-height section generation, the restriction-to-Keller-map
argument and the external degree theorem remain written/imported inputs,
without formal verification, independent full reproof or external review.

Actual ruled compatibility, the unrestricted finite-map realizability
problem, both quartic packets and JC2 remain open. The ambitious goal
remains active.

The subsequent [all-height pole-three theorem](POLE_THREE_ALL_HEIGHTS.md)
excludes every(3,n) presentation, independently of the2022 degree bound.
The current first open reduced pairs are (6,8),(7,8).

The adopted frontier(6,8),(7,8) incorporates the subsequent
[PF2RULED12 full(5,8) exclusion](POLE_FIVE_EIGHT_EXCLUSION.md); it is
not an additional exclusion proved by this earlier note.
