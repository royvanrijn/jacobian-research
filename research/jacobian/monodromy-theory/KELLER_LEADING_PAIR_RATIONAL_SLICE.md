# Rational slices of leading pairs in three variables

**Theorem KellerLeadingPairRationalSlice.**
Let \((P,Q,R)\) be a three-variable polynomial Keller map over a
characteristic-zero field \(k\), and let \(A,B\) be the highest ordinary
homogeneous parts of \(P,Q\). Suppose \(A,B\) are algebraically
independent. Then
\[
                  \det J(A,B,\phi)=\det J(P,Q,R)
\]
for some \(\phi\in k(x,y,z)\). If \(dA,dB\) are independent at the
generic points of a specified finite set of prime divisors of
\(\mathbb A^3\), one can take \(\phi\) regular at all those generic
points.

The same statement holds for highest parts with respect to any
fixed positive rational source weights. In the formal proof,
clear the weight denominators, rescale each source variable by
its corresponding power of \(\epsilon^{-1}\), and replace
the exponent \(n\) by the sum of the three component weights
minus the sum of the three source weights.

**Corollary KellerLeadingPairHigherMemberBasepointObstruction.**
In the same setting, assume \(A,B\) are coprime, with
\(\deg P=d\ge\deg Q=e\) and \(d\ge3\). Over the algebraic closure
of \(k\), there is a common point of the projective plane curves
\(A=0\) and \(B=0\) at which \(A=0\) is singular. In particular
if \(d=e\), at least one common point is singular on both curves.

**Corollary KellerLeadingPairDisjointTangentObstruction.**
Suppose \(A,B\) are coprime in the same setting, with degrees
\(d,e\). At a projective common point let \(m,n\) be their respective
multiplicities, and put
\[
    \lambda=\min(d/m,e/n),\qquad
    L=d+e-3-\lambda(m+n-2).
\]
Call the point admissible for the following branch calculation if
the two tangent cones have no common line and, when the two ratios
are unequal, the tangent cone for the larger ratio is reduced.
An admissible point with \(L=0\) is impossible. If every common
point is admissible, at least one has \(L<0\).
In particular a sextic leading form and a coprime leading form of
degree four, five or six cannot have only common points of
multiplicity at most three on each curve, with reduced tangent
cones sharing no line.

The formal method is related to the interpolation of leading-component
coordinates in the
[quintic conic-pencils argument](QUINTIC_KELLER_CONIC_PENCIL_FORMAL_SLICE.md).
Here the leading pair is arbitrary, and the coefficient of the
first nonzero determinant gives a rational primitive on its generic fibre.

## Formal rectification

Write \(d_1,d_2,d_3\) for the component degrees and let
\(P_\epsilon(X)=\epsilon^{d_1}P(X/\epsilon)\), with similar definitions
for \(Q_\epsilon,R_\epsilon\). Put \(n=d_1+d_2+d_3-3\) and
\(\kappa=\det J(P,Q,R)\). The determinant of the rescaled map is
\(\kappa\epsilon^n\).

Let \(S\) be the localization of \(k[x,y,z]\) obtained by inverting
the polynomials that are nonzero at the generic points to be protected.
For the empty set take the fraction field. Choose a linear form
\(\ell\) such that \(J(A,B,\ell)\) is invertible over \(S\). Indeed,
the vectors of the \(2\)-minors of \(J(A,B)\) are nonzero at each
of the specified generic points. Each generic point therefore
forbids at most one proper \(k\)-linear subspace of choices for
the coefficients of \(\ell\). A finite union of such subspaces
cannot fill the space over the infinite field \(k\).

Solve
\[
 (P_\epsilon,Q_\epsilon,\ell)(X_\epsilon)=(A,B,\ell)(X),
 \qquad X_\epsilon=X+O(\epsilon)
\]
recursively in \(S[[\epsilon]]\). At each order the new coefficient
vector is multiplied by the same invertible matrix \(J(A,B,\ell)\).
Differentiation preserves \(S\), so the formal chain rule gives
\[
 \det J\big(A,B,R_\epsilon(X_\epsilon)\big)
           =\kappa\epsilon^n\det J(X_\epsilon).
\]
The second determinant has constant coefficient one. The coefficient
of \(\epsilon^n\) in \(R_\epsilon(X_\epsilon)\) is the required
\(\phi\).

## The homogeneous-pair fibre at infinity

After extending constants, use algebraically independent target
values \(a,b\) for \(A,B\). The smooth affine generic fibre of
\((A,B)\) acquires a nonzero relative differential \(\omega\) by
division of \(dx\wedge dy\wedge dz\) by \(dA\wedge dB\).
The theorem gives \(d\phi=\kappa\omega\) on the generic fibre.
In particular the pullback of \(\omega\) to the smooth projective
model of each geometric component has a rational primitive.

Assume \(A=0\) is smooth at every common point at infinity.
The projective closure is given by \(A-aw^d=B-bw^e=0\) in
\(\mathbb P^3\); coprimality ensures there is no curve contained
in \(w=0\). In a projective affine chart at a common point,
choose local coordinates \(u,v\) in the plane with \(v=A\).
The equations become \(v=aw^d\) and
\[
                    H(u,w)=B(u,aw^d)-bw^e=0.
\]
If \(j\) is the local intersection multiplicity of the two
plane curves, \(H(u,0)=c u^j+O(u^{j+1})\) with \(c\ne0\).
The coefficient of \(w^e\) in \(H(0,w)\) is nonzero over
the generic coefficient field. It is \(-b\) for \(d>e\)
and \(aB_v(0,0)-b\) for \(d=e\). There are no nonzero
terms below the segment joining \(u^j\) to \(w^e\).
After \(w=\tau^j,\ u=\tau^e V\), the initial equation is
\(cV^j+\beta=0\), \(\beta\ne0\), with distinct nonzero roots.
Lifting these roots gives all branches after that ramified
base change, by Weierstrass preparation in \(u\).

Up to a local unit the relative form is
\[
                       w^{d+e-4}\frac{dw}{H_u}.
\]
In fact \(dx\wedge dy\wedge dz\) contributes \(w^{-4}\)
in the projective chart, and, modulo the defining equations,
the two fixed affine-fibre differentials contribute \(w^{-d}\)
and \(w^{-e}\). On each branch after the base change the form
therefore has order
\[
               (d+e-4)j+j-1-e(j-1)=(d-3)j+e-1\ge0 .
\]
Regularity after ramified pullback implies regularity before
pullback, since a differential of order \(h<0\) pulls back
to one of order \(r(h+1)-1<0\) under ramification index \(r\).
Thus \(\omega\) is everywhere regular on the projective model.
A nonzero regular differential on a smooth projective curve
cannot be the differential of a rational function: any pole
of a function gives a pole of strictly higher order of its
differential, and a regular function is constant. This contradicts
the slice theorem.

When \(d=e\), replace \(A\) by a generic constant combination
of \(A,B\), keeping \(B\) independent. Coprimality and the
common points are unchanged. If some member has nonzero
differential at each common point, a combination can be
chosen smooth at all of them, since the set of points is
finite. The first assertion rules this out.

## Common points with disjoint tangent cones

At an admissible common point use plane coordinates \(u,v\), and
write \(A_m,B_n\) for the initial binary forms. Along a branch of
the generic fibre above the point, let \(\rho\) be the smaller
of the orders of \(u,v\) relative to the order of \(w\).
If \(\rho<\lambda\), the leading direction must be a common
zero of \(A_m,B_n\): both initial orders are strictly below
the orders of \(aw^d,bw^e\). This contradicts disjointness.
If \(\rho>\lambda\), the equation for the smaller ratio
cannot attain its prescribed order. Thus \(\rho=\lambda\).

Make a ramified parameter change \(w=\tau^r\) with
\(r\lambda\) integral, and set \(u=\tau^{r\lambda}U,\ v=
\tau^{r\lambda}V\). After dividing out the initial powers
of \(\tau\), the residual equations at \(\tau=0\) are
\[
 \begin{cases}
 A_m(U,V)=a,\ B_n(U,V)=b, & d/m=e/n,\\
 A_m(U,V)=a,\ B_n(U,V)=0, & d/m<e/n,\\
 A_m(U,V)=0,\ B_n(U,V)=b, & d/m>e/n .
 \end{cases}
\]
They have nonempty finite solution sets with nonsingular
two-by-two Jacobian. For equal ratios the map
\((A_m,B_n)\) is finite by disjointness of the two projective
zero sets, and its generic fibre is smooth in characteristic
zero. For unequal ratios the equation with zero right-hand
side gives distinct tangent lines. On each line the other
initial form is nonzero off the origin by disjointness, and
its derivative in the radial direction is nonzero at the
prescribed nonzero value. This derivative and the transverse
derivative of the zero equation give the required rank.
Formal implicit lifting gives a branch through each
solution and accounts for all branches after ramified
extension, by the valuation argument.

In this chart the relative differential is, up to a unit,
\[
             \frac{w^{d+e-4}\,dw}
                  {\det \partial(A,B)/\partial(u,v)}.
\]
The denominator has order \(r\lambda(m+n-2)\) on each
lifted branch. The differential has order \(rL-1\).
If \(L>0\), its pullback and hence the original differential
are regular at those points. If \(L=0\), the pullback has
a simple pole with nonzero residue, which precludes an exact
differential. If every point is admissible and all \(L\) are
positive, the differential is everywhere regular and nonzero,
also precluding exactness.

In the sextic consequence, \(d=6\), \(4\le e\le6\),
and \(m,n\le3\). Since \(d+e\ge\lambda(m+n)\),
one has \(L\ge 2\lambda-3>0\) whenever \(\lambda>3/2\).
The only remaining possibility is \(\lambda=4/3\) from
\(e=4,\ n=3\), and then
\(L=7-(4/3)(m+1)\ge5/3>0\).

## Comparing two divisors

A useful additional consequence is the following. If some rational
primitive \(\phi_0\) has a pole at one protected divisor and is
regular at another, then a completion would require a correction in
the constant field of \(\det J(A,B,\cdot)\) with a pole at the first
divisor and no pole at the second. When this constant field is
\(k(p,q)\) and both divisors map dominantly to the same curve in
the \(p,q\)-plane, this is impossible. Indeed its valuations at the
two generic points restrict to positive multiples of the valuation
at that curve, because every unit at the curve pulls back to a
unit at both generic points. The same argument works with any
specified rational coordinates of the constant field.
