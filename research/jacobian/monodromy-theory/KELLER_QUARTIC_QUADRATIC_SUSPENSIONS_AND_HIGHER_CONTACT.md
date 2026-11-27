# Quadratic suspensions and higher-contact quartic boundaries

**Theorem KellerQuadraticSuspensionPlaneMate.** Let \(F=(P,Q,R)\)
be a three-variable polynomial Keller map over a characteristic-zero
field with \(R=Z^2+B(x,y)\), where \((x,y,Z)\) are polynomial source
coordinates of constant nonzero Jacobian. Then \(B\) has a polynomial
plane Jacobian mate. With weights \(1,1,k\) on \(x,y,Z\), where
\(k\) is a positive integer and \(\deg B\le2k\), suppose the
weighted degrees of \(P,Q\) in these source coordinates are at most
\(M_P,M_Q\). There is a mate of ordinary degree at most
\(M_P+M_Q-k\).

Use the exact residue identity of
[Wang–Yan, equation (3.28)](https://arxiv.org/pdf/2610.36495v1).
On the fibre \(Z^2+B=t\) the residue form is
\(dx\wedge dy/(2Z)\), with a regular polynomial primitive obtained
from a constant multiple of \(P\,dQ\). This formula means the
restriction on \(Z\ne0\), and the residue form extends on the
smooth fibre. The displayed form is odd under \(Z\mapsto -Z\).
Taking the odd part of the primitive preserves it.

Reduce even powers of \(Z\) using \(Z^2=t-B\). The odd primitive
has the form \(Z\beta_t+D_t\,dZ\), where \(\beta_t\) is a plane
one-form and \(D_t\) a plane polynomial, both polynomial in \(t\).
Subtract \(d(ZD_t)\), differentiating with \(t\) fixed on the fibre,
and replace \(\beta_t\) by \(\beta_t-dD_t\). Thus \(Z\beta_t\) is
a primitive. Since \(2Z\,dZ=-dB\) on the fibre, its exactness is
the following polynomial identity in \(x,y,t\):
\[
 2(t-B)\,d\beta_t - dB\wedge\beta_t=dx\wedge dy.
\tag{1}
\]
The equality known on fibres gives the identity on \(t\ne B(x,y)\)
and hence everywhere.

Suppose \(m=\deg_t\beta_t\ge1\). The \(t^{m+1}\)-coefficient
of (1) says that the coefficient \(\beta_m\) is closed. On the
plane it is \(dH\) for a polynomial \(H\), by polynomial integration
in characteristic zero. Subtract the exact form
\(d(Z^3t^{m-1}H)\) from the fibre primitive. Its expression as \(Z\)
times a plane one-form is
\[
 Z\,t^{m-1}\bigl((t-B)dH-\tfrac32 H\,dB\bigr).
\]
This cancels the highest coefficient of \(\beta_t\) and reduces
its degree in \(t\). After finitely many such steps, \(\beta_t\)
is independent of \(t\). Equation (1) then gives \(d\beta=0\).
Writing \(\beta=dA\) gives \(dA\wedge dB=dx\wedge dy\).

For the degree bound assign weight \(2k\) to \(t\) and the
coordinate's weight to its differential. Reduction by \(Z^2=t-B\)
does not raise weights. Both \(\beta_t\) and \(D_t\) initially have
weighted degree at most \(M_P+M_Q-k\), as does \(\beta_t-dD_t\).
If the current bound is \(M_P+M_Q-k\), the primitive \(H\) of its
top coefficient can be chosen of degree at most \(M_P+M_Q-k-2km\).
The displayed exact subtraction preserves the bound, so the
final planar primitive \(A\) satisfies it.

## The rank-one normal jet of a quartic component

**Corollary QuinticKellerQuadraticSourceCompletion.** A three-variable
Keller map of ordinary degree at most five is invertible if a
constant target component has the form
\[
 R=c z^2+a(x,y)z+b(x,y),\qquad c\ne0,\quad
 \deg a\le3,\quad \deg b\le4
\]
in some affine source frame. In particular, this applies to any
point at infinity of a quartic component where the projective
fibre equation has zero linear jet and nonzero quadratic jet
proportional to the square of the homogenizing coordinate.
This includes every rank-one jet at a plane point of multiplicity
at least three, ordinary or otherwise, where the linear jet vanishes.

Over \(\mathbb C\), complete the square, take
\(Z=\sqrt c\,(z+a/(2c))\), and put \(B=b-a^2/(4c)\).
The change of source coordinates has constant nonzero Jacobian.
With weights one on \(x,y\) and three on \(Z\), the transformed
\(P,Q\) have weighted degree at most fifteen. The theorem gives
\(B\) a planar mate of degree at most twenty-seven. The retained
[planar degree bound](../../BREAKTHROUGH_TARGET_ASSESSMENT.md#four-variable-hessian-conjecture-and-planar-degree-bounds)
makes this pair a planar automorphism.

Let \(d=\deg B\le6\). A planar automorphism and its inverse have
the same ordinary degree: the degree of the image of a general
affine line is the degree of a general linear combination of the
inverse coordinates, and is also the degree of a general linear
coordinate on the image computed in the birational parametrization
by the line. The latter is the degree of the automorphism.

If \(d=1\), choose an affine mate \(A_0\) of \(B\).
The completion \((A_0,Z,R)\) has inverse degrees at most two
for \(x,y\) and at most six for \(z\). If \(d>1\), by planar
degree reduction choose a mate \(A_0\) of degree \(e<d\).
Here \(e\) divides \(d\). Reducing the other coordinate in
\((A_0,B)\) gives \(B_1=B-H(A_0)\) of degree at most \(e\),
with \(\deg H=d/e\). These degree reductions follow from
[García Barroso–Gwoździewicz–Płoski, Proposition 3.4 and Lemma 3.5](sources/garcia-gwozdziewicz-ploski-coordinate-lines-2015.pdf)
([author copy](https://ergarcia.webs.ull.es/Semigroup-Forum.pdf)).
For the second reduction when \(e=1\), use \(A_0\) as an affine
source coordinate: \(B\) has degree one in the other variable
with constant nonzero slope, so subtracting the remaining
polynomial in \(A_0\) leaves an affine coordinate.

The inverse of \((A_0,B_1)\) has degree at most \(e\).
Inverting the completion \((A_0,Z,R)\) means applying it to
\((A_0,R-Z^2-H(A_0))\). This input has degree at most \(d/e\),
since \(d/e\ge2\). Thus \(x,y\) have inverse degree at most
\(d\le6\), and \(z=Z/\sqrt c-a(x,y)/(2c)\) has inverse degree at
most \(3d\le18\). The
[bounded coordinate-completion theorem](LOW_DEGREE_COORDINATE_KELLER_MAPS_IN_THREE_VARIABLES.md)
applies with degree product at most ninety. Complex embedding
and formal inversion descend the conclusion to characteristic zero.

For the projective-jet application, put the point at
\([0:0:1:0]\). The local degree of a homogenized term from
\(x^iy^jz^\ell\) is \(4-\ell\). The vanishing linear jet excludes
\(z^3\) terms. In the local coordinates \(u,v,w\), a quadratic
jet proportional to \(w^2\) means the coefficient of \(z^2\)
is a nonzero constant. The remaining coefficient bounds follow
from ordinary total degree at most four. At a plane point of
multiplicity at least three a rank-one quadratic form must be
a multiple of \(w^2\), because it vanishes on \(w=0\).

## A local extension at higher-contact double points

Suppose \(R\) is quartic with reduced leading curve. At a plane
double point, assume the quadratic form of the full projective
fibre equation has rank one, and the plane curve also has
nonzero quadratic part. By the parameter-preserving analytic
splitting used in the
[quartic boundary proof](KELLER_QUARTIC_COMPONENT_BOUNDARY_OBSTRUCTION.md#ordinary-double-points-of-the-plane-leading-curve),
write the germ as
\[
 u^2+f(v,w)-tw^4=0,\qquad f=\sum_{j\ge3}f_j,
\]
with \(w=0\) the plane at infinity and \(f_j\) homogeneous.
Suppose that \(f(v,0)\) has order at least five, so the ordinary
cusp and tacnode bounds do not already apply.

**Proposition KellerHigherContactQuarticLocalTests.** In this
setting the following are sufficient local tests for general \(t\).

- If \(f_3\ne0\) and its zero at \([v:w]=[1:0]\) is simple, the
  residue extends holomorphically on a resolution.
- If \(f_3(1,W)=aW^2+O(W^3)\) with \(a\ne0\), put
  \(b=[W]f_4(1,W)\) and \(c=f_5(1,0)\). The residue extends
  holomorphically when \(b^2-4ac\ne0\).
- If \(f_3(1,W)=aW^3\) with \(a\ne0\), the residue extends
  holomorphically when either \(b\ne0\) or \(c\ne0\), with the
  same definitions of \(b,c\).
- If \(f_3=0\) and \(f_4(1,W)\) has order at most two at \(W=0\),
  the fibre form has a nonzero local tube period.

Thus if every boundary point satisfies either a prior
holomorphic-extension test or one of these tests, the
hypothetical Keller component is excluded. A local tube-period
contradiction already suffices by itself.

For the three holomorphic assertions blow up the double point
first. Any multiple exceptional root of \(f_3\) away from \(w=0\)
is handled for general \(t\) by the nonzero varying quartic
coefficient, as in the prior proof. The direction \(w=0\) in
the \(v\)-chart has coordinates \(v=s,\ w=sW,\ u=sU\), and its
equation is
\[
 U^2+s f_3(1,W)+s^2\bigl(f_4(1,W)-tW^4\bigr)
       +s^3 f_5(1,W)+\cdots=0.
\]
A simple root gives a rank-three quadratic term. For a
double root the indicated binary cubic is
\(s(aW^2+bsW+cs^2)\). Under the discriminant condition its
three roots are distinct, so the next double-point blowup
leaves only ordinary double points.

For a triple root with \(b\ne0\), blow up the point whose
binary cubic is \(s^2(bW+cs)\). The simple root is immediate.
In the chart \(W=z,\ s=zS,\ U=zU_2\) the double root \(S=0\)
has binary cubic \(zS(bS+az)\), with three distinct roots.
Indeed the terms of degrees three and four of the preceding
equation give \(z(bS^2+cS^3)\) and \(z^2(aS+O(S^2))\).
All higher terms have a factor \(S\), because the whole
nonquadratic series before this blowup is divisible by \(s\).
This gives the claimed ordinary-double-point resolution.

For \(b=0,\ c\ne0\) the chart \(W=z,\ s=zS,\ U=zU_2\)
instead gives binary cubic \(az^2S\), with simple root \(S=0\)
and double root \(z=0\). Blow up that double root using
\(S=h,\ z=hT,\ U_2=hU_3\). Its binary cubic is
\(hT(aT+ch)\), again with three distinct roots. The term \(ch^2T\)
comes from \(czS^3\) in the preceding chart, and all other
terms have degree larger than three in \(h,T\). This completes
the holomorphic assertions. All point blowups used in these
resolutions are blowups of surface multiplicity two, so the
geometric residue has no new pole on their exceptional curves.

For the last assertion the normalized double cover of the
blown-up base has exceptional curve
\(U^2+f_4(V,1)-t=0\), with residue a nonzero constant times
\(dV/(2U)\). The polynomial \(f_4(V,1)\) has degree at least two
by the order hypothesis at infinity, and its finite roots
after subtracting a general \(t\) are simple. If the quartic
has no double root at infinity it gives a smooth genus-one
curve with nonzero holomorphic residue. If it has one double
root at infinity, the residue has a nonzero simple pole on
each normalized branch over it. Explicitly in the complementary
chart \(V=1/W,\ U=U_v/W^2\) the form is \(-dW/(2U_v)\) and
\(U_v\) has order one at each branch. In the latter case a
small loop about either normalization point gives a nonzero
period. The tubes over the loops are taken over smooth points
of the exceptional curve away from strict infinity, so
Stokes gives the contradiction even if the surface model is
singular at the double root itself.

## Necessary singular point in a reduced leading curve

**Corollary KellerReducedQuarticBoundaryRequiredJet.** For a quartic
component of a three-variable Keller map with reduced leading curve
\(R_4=0\), there is a singular point outside the ordinary node,
cusp, tacnode, and ordinary triple types at which \(R_3\) vanishes
and the projective fibre quadratic jet has rank at most one.
At this point none of the applicable higher-contact tests above
can succeed. For a map of degree at most five that is noninvertible,
a rank-one quadratic jet there cannot be proportional to the square
of the homogenizing coordinate. In particular, if this required point
has plane multiplicity at least three, its quadratic jet must vanish.

If \(R_3\) does not vanish at a boundary point, the normal
direction to the plane at infinity gives a nonzero first derivative
of the fibre equation. Such a point is smooth. At a point where
the full quadratic jet has rank at least two, splitting off a
nondegenerate pair gives a surface double point with a holomorphic
resolution of the form, by the
[initial quartic argument](KELLER_QUARTIC_COMPONENT_BOUNDARY_OBSTRUCTION.md#ordinary-double-points-of-the-plane-leading-curve).
This applies irrespective of the plane singularity. Its
isolatedness follows because the reduced plane leading curve has
finitely many singular points, and the fibre is smooth both on
the affine part and at all smooth leading-curve points.

Now combine the local tests with the
[ordinary triple point argument](KELLER_QUARTIC_TRIPLE_BOUNDARY_AND_LINEAR_COMPONENTS.md#local-analysis-at-an-ordinary-triple-point)
and the tests for the other three mild singularities. Any local
tube-period test gives a contradiction. If all points admit
holomorphic resolution, they do so simultaneously on a general
fibre and its global holomorphic residue gives a contradiction
as well. The stated point must therefore occur. The last claim
uses QuinticKellerQuadraticSourceCompletion and the fact that a
rank-one quadratic form restricting to zero on the plane at
infinity is a multiple of the homogenizing coordinate squared.
