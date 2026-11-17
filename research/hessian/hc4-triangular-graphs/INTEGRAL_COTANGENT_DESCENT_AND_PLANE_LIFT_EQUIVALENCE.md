# Integral cotangent descent and exact plane lifting equivalence

Claims HC4MYKL0/1. [MATH_STATUS.json](../../MATH_STATUS.json) is the
authority. This strengthens the preceding [affine-plane specialization](SPECIALIZATION_OF_AFFINE_PLANES_AND_FROBENIUS_SHEARS.md)
from a geometric necessary condition to an exact arithmetic reduction.
It uses TWO retained results: HCTANG2 supplies a plane under hypothetical
nonautomorphy, and HC4MYHF0 makes that plane unique. Their assurance limits
remain in force. Neither HC4 nor the resulting plane lifting problem is
resolved.

## 1. The cotangent plane descends to the coefficient field: HC4MYKL0

Let K be ANY characteristic-zero field, h in K[z1,z2,z3,z4], and
det Hess(h) in K*. Suppose grad(h) is not a polynomial automorphism.
This remains true after extending the coefficient field: translate its
constant output to zero, and use uniqueness of the formal inverse to
descend any polynomial inverse from the extension field back to K.

By the finite coefficient-field/complex embedding argument in KK0,the
retained HCTANG2 alternative supplies a constant jointly affine two-plane
over an algebraic closure of K. HC4MYHF0 proves that two distinct such
planes would force the gradient to be an automorphism. Thus the plane
W is unique GEOMETRICALLY, not just unique among K-rational planes.

For every sigma in Gal(Kbar/K), sigma(W) is another such plane because
h has coefficients in K. Uniqueness makes W Galois invariant. Put a
basis matrix into its unique reduced row-echelon form; its pivot set
and all nonpivot entries are Galois invariant. Since K has characteristic
zero, the fixed field in its algebraic closure is K. The row-echelon
entries therefore lie in K. This proves W is defined over K.

Choose a K-linear source frame whose last two vectors span W. Then

    h(T(x,y,u,v))=R(x,y)+uF(x,y)+vG(x,y),
    F,G,R in K[x,y], J(F,G) in K*.                              (1)

Indeed joint affineness gives the displayed polynomial form and the
block Hessian identity gives det Hess(h composed T)=J(F,G)^2.
The latter is a nonzero constant, so J(F,G) is a nonzero constant.
There is no algebraic coefficient extension or rational function in
the source variables in this frame.

This is a descent consequence of retained existence PLUS retained
uniqueness. Existence alone would give only a geometric plane; a chosen
plane with several Galois conjugates would not justify (1).

## 2. Integral frames and good reduction: HC4MYKL1

Let O be a DVR with characteristic-zero fraction field K and arbitrary
residue field k. Let h in O[z1,z2,z3,z4] have constant Hessian determinant
in O*, and assume its residue gradient grad(h0) is NOT a polynomial
automorphism over k. The good-reduction inverse argument in KK0 excludes
the automorphism branch over K. KL0 gives its unique plane W over K.

Take L=W intersect O^4. As in KK0,L is a saturated rank-two direct
summand of O^4. Extend an O-basis of L to an O-basis of O^4. This gives

    T in GL4(O),
    h composed T=R(x,y)+uF(x,y)+vG(x,y),
    F,G,R in O[x,y], J(F,G) in O*.                              (2)

The polynomial is jointly affine in the last two coordinates over K;
the coefficients are already in O, so all higher fibre-degree terms
vanish in O as well. The Hessian block identity gives a constant square
J(F,G)^2 in O*. Over K,a polynomial whose square is constant is itself
constant;its valuation is zero because twice that valuation is zero.
Thus its constant Jacobian is an O-unit, as asserted. The frame T is
integral with unit determinant;denominator-bearing generic frames are
replaced by this saturated frame. No ramified extension is needed.

In particular h0 has a k-RATIONAL constant jointly affine two-plane,
not merely a geometric one. Its reduction in the same frame is

    h0 composed T0=R0(x,y)+uF0(x,y)+vG0(x,y),
    T0 in GL4(k), J(F0,G0) in k*.                               (3)

Moreover (F,G) in (2) is an actual constant-Jacobian lift of the exact
pair (F0,G0) in (3).

Conversely suppose a residue polynomial h0 with noninvertible gradient
admits a frame (3),and that pair admits a constant-unit-Jacobian lift
F,G over O. Lift the finitely many coefficients of R0 arbitrarily to R.
Lift the entries of T0 to a matrix T in GL4(O);its determinant is a
unit,so its inverse is integral. Then

    h=(R+uF+vG) composed T^(-1)                                 (4)

is a polynomial over O reducing to h0 and has constant unit Hessian
determinant det(T)^(-2)*J(F,G)^2. This proves the exact equivalence:

    For a fixed noninvertible-gradient residue h0, existence of ANY
    characteristic-zero constant-unit-Hessian polynomial lift over O
    is equivalent to existence of a k-linear cotangent frame of h0
    whose two plane coefficient polynomials have a Keller lift over O. (5)

There is no correction-degree bound in either direction. A prospective
lift is allowed to have arbitrary four-variable corrections;the theorem
recovers an integral constant source frame from the actual characteristic-
zero solution. Thus such corrections provide no separate lifting escape
beyond a liftable plane pair in some residue linear frame. This is a
reduction of the existence question,not a plane lifting sufficiency
conjecture or a construction of the missing pair.

### Fixed residue frame when the plane is unique

If h0 is already in (3) with T0=I and its geometric jointly affine
two-plane is uniquely the standard u,v plane,then the reduction of L
is that standard plane. Its O-basis and its complement can be chosen
congruent to the standard coordinate vectors. Thus (2) uses T=I modulo
the uniformizer,and F,G reduce to the ORIGINAL F0,G0,R reduces to R0.
For such a residue potential,arbitrary-potential lifting is equivalent
to Keller lifting of that fixed pair over the SAME O. No additional
residue frame,field extension or nonlinear coefficient map is required.

For the retained H0=uP+vQ,the affine-line-direction argument in KK1
with b=0 first forces w_x=w_y=0 from the coefficient of v and the top
monomial x^8*y^3 of Q. Conversely every direction in u,v is affine.
Hence its unique jointly affine plane is exactly span(partial_u,partial_v).
KJ0 therefore excludes a characteristic-zero lift of H0 by (5) using
the plane divided equation alone. This is a conceptual all-degree
recovery of that NONLIFTING conclusion. It does not replace KJ1's
stronger sharp finite-precision theorem:the potential has a mod4 lift
while its plane pair has no mod4 lift. Tang and uniqueness are being
used in characteristic zero,not asserted over truncated coefficient
rings,so there is no finite-precision equivalence here.

The nonlinear square-shear family from KK1 has NO jointly affine
two-plane and is excluded already before the plane-pair stage. Other
nonlinear source changes remain outside that classification.

## 3. Consequences and limits

The next arithmetic gate is now precise:find a residue linear cotangent
frame with a genuinely liftable plane Keller pair,or prove that all
such frames are obstructed. A noninvertible residue gradient without
a k-rational plane cannot lift. A geometric plane without a rational
one is also insufficient. A supplied rational plane is necessary but
its coefficient pair still has to lift. Passing finitely many divided
equations is not proved sufficient.

The theorem applies to all residue characteristics and all ramification.
It does NOT equate arbitrary finite-precision potential lifts to plane
lifts,does not give a uniform finite-precision obstruction,does not
classify all residue frames,and does not prove that every positive-
characteristic Keller pair lifts or fails to lift. There is no new
characteristic-zero example and no full HC4/JC2 decision. The separate
native critical non-CI and quotient-action questions remain open.

The [checker](verify_integral_cotangent_descent.py),
[receipt](receipt-integral-cotangent-descent.json) and
[context](sources/integral-cotangent-descent-context.json) bind the exact
block identities,source inverse/congruence factors,rank-preserving frame
controls and retained proof bytes. The Galois and saturated-lattice
arguments are written proofs;finite symbolic checks do not certify the
retained Tang classification. No new external theorem is imported.
Independent replay,formal verification and external review are not claimed.
