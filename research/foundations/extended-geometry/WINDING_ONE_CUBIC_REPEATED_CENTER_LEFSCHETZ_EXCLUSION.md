# A Lefschetz obstruction at a repeated center

Use the [quartic-leading septic setup](WINDING_ONE_CUBIC_QUADRUPLE_ROOT_LEFSCHETZ_EXCLUSION.md)
over a characteristic-zero field with \(U=1+x\), \(v=Ut\),
\(Y=v+d(x)/3\), and
\[
 f=\frac Ux\bigl(v^3+d(x)v^2+b(x)v+c(x)\bigr)
   =\frac Ux\bigl(Y^3+p(x)Y+q(x)\bigr).
\]
Put \(u=x/(1+x)\) and \(P(u,Y)=Y^3+p(u/(1-u))Y+q(u/(1-u))\).
If \(d(0)=\pm1\), the cubic \(P(0,Y)\) has one double root
\(y_0=-2d(0)/3\) and one simple root. Let \(h(u)\) be the formal
critical branch \(P_Y(u,h(u))=0\) with \(h(0)=y_0\), and put
\[
 \phi(u)=P(u,h(u))/u,\qquad \zeta_c=\phi(0),\qquad a=\phi'(0).
\]
These coefficients exist because the initial numerator is zero and
\(P_{YY}(0,y_0)=6y_0\ne0\). In particular
\(\zeta_c=P_u(0,y_0)\) and
\(2a=(P_{uu}-P_{uY}^2/P_{YY})(0,y_0)\).

**Theorem DegreeFourSepticRepeatedCenterExclusion.** Suppose
\(d(0)^2=1\), and use any of the smooth pole models described below.
Let hypotheses 2 and 3
of [DegreeFourSepticLefschetzExclusion](WINDING_ONE_CUBIC_QUADRUPLE_ROOT_LEFSCHETZ_EXCLUSION.md)
hold, and suppose \(a\ne0\) and \(\zeta_c\ne0\), with \(\zeta_c\)
distinct from all interior critical values. Then \(\Lambda(f^m)\)
is nonzero for infinitely many positive orders. The conditions hold
for an explicit repeated-center profile of varying depressed
invariant, given below.

The allowed models at \(u=1\) are the dominant constant-term model
\(\deg p\le4,\deg q=7\), the degree-five \(p\) model with
\(\deg d\le2\), and the high cubic-shear model with \(d_3\ne0\)
and low-branch discriminant order \(k<7\). These models are smooth
over the finite value line and have the required large-circle
connectivity by the local arguments in the
[high-pole exclusion](WINDING_ONE_CUBIC_QUADRUPLE_ROOT_GENERIC_EXCLUSION.md)
and the
[low-branch collision calculation](WINDING_ONE_CUBIC_QUADRUPLE_ROOT_LOW_BRANCH_COLLISION.md#the-quadratic-at-infinity).

## 1. The resolved center fiber

Work over \(\mathbb C\) after embedding the finitely generated field
of coefficients. The notation \(h(u)\) gives an analytic branch at
zero by the implicit function theorem. Near the double root the
identity
\[
 P(u,Y)-P(u,h(u))=(Y-h(u))^2(Y+2h(u))
\]
gives an analytic coordinate \(w=(Y-h(u))\sqrt{Y+2h(u)}\), with a
local choice of square root of the unit. The pencil \(P=\zeta u\)
then has the local equation \(w^2=u(\zeta-\phi(u))\). Write
\(s=\zeta-\zeta_c\), \(r=\phi(u)-\zeta_c\), and
\(w_1=w/\sqrt{u/r}\). Because \(a\ne0\), \(r\) is a coordinate and
\(u/r\) is a regular unit. The equation is
\[
                     w_1^2=r(s-r).
\]
This is an ordinary double point of the surface at the origin.
Blow up that point. In the \(r\)-chart, write \(w_1=rW\) and
\(s=rS\). The strict equation is \(W^2=S-1\), so the resolved map
is \(s=r(1+W^2)\) on smooth coordinates \((r,W)\). Its two critical
points are ordinary nodes at \(r=0,W=\pm i\). The exceptional
conic is smooth and the only point of it missed by the \(r\)-chart
has \(r=w_1=0\) in homogeneous blowup coordinates. It is in the
\(s\)-chart, where \(r=sR,w_1=sW\) and \(W^2=R(1-R)\). Near
\(R=W=0\) this is smooth with coordinates \((s,W)\) and the
projection is a submersion. The marked center point \(u=0,Y=y_0\)
extends as the section \(R=W=0\), away from both nodes.

On a fixed small positively oriented \(u\)-circle and for \(s\) near
zero, the two inverse branches near \(y_0\) are single-valued in
\(u\) on a neighborhood of the circle. At \(s=0\) they correspond
to \(W=\pm i\) in the \(r\)-chart. Their lifts for \(s\ne0\) are
annular core circles near these two nodes: indeed \(W=\pm\sqrt{s/r-1}\)
on the circle, and \(r(u)\) winds once about zero. Thus the circles
\(\delta_+,\delta_-\) are the two vanishing cycles, with the common
positive orientation in \(r\) or \(u\). Their sum is zero in compact
homology. They bound the smoothed part near the exceptional conic,
using the same sign on both boundary circles since the local
equation for that component is \(r=0\) in both Morse charts.

## 2. The two period derivatives

Use \(\alpha=Ydu/u\) on smooth fibers. With the preceding
orientations, the lifted-circle periods \(I_\pm(s)\) extend
analytically to \(s=0\), both with value \(2\pi i y_0\). On the
central local sheets,
\[
 Y_\pm(u,0)=h(u)\pm\kappa u+O(u^2),\qquad
 \kappa^2=-a/(3y_0)\ne0.
\]
Differentiation of the inverse equation gives
\(\partial_\zeta Y=u/P_Y(u,Y)\), so differentiation on the fixed
circles yields
\[
 \frac{I_\pm'(0)}{2\pi i}
   =\operatorname{Res}_{u=0}\frac{du}{P_Y(u,Y_\pm(u,0))}
   =\pm\frac{1}{6y_0\kappa}.
\]
For the last equality use \(P_Y(u,Y)=3(Y^2-h(u)^2)\) and the
expansions above. Hence \(I_+-I_-\) is a nonzero germ.

The monodromy about \(s=0\) is the product of the two positive
Dehn twists on disjoint annuli, also in the punctured family since
the marked section misses both annuli. Their compact intersections
with a transported cycle are opposite. Thus for a cycle with zero
period germ the period variation is, up to sign,
\(\langle\gamma,\delta_+\rangle(I_+-I_-)\), which must be zero.
The nonzero derivative forces zero compact intersection with both
vanishing cycles. The same local Dehn twist and handle inputs as
in the preceding theorem apply. [Apostolakis, Piergallini and Zuddas,
§6 and Proposition 9](https://arxiv.org/pdf/1104.4536).

## 3. The compact-fiber argument

At the other center root the finite model is smooth and the map is
a submersion by the nonzero derivative with respect to \(Y\). At
the double root and all \(s\ne0\) the derivative in \(u\) is nonzero,
so the pencil and marked section are smooth there. The explicit
blowup therefore resolves the only extra finite singularity of the
total surface compared with the preceding theorem. The models over
\(u=1,\infty\) and the ordinary interior critical values are
unchanged by the blowup at \(u=0\). In particular the large-value
inverse-branch circle is still nonseparating in the compact fiber.

The proper finite model is still rational after the blowup. The
[compact-fiber spanning argument](WINDING_ONE_CUBIC_QUADRUPLE_ROOT_LEFSCHETZ_EXCLUSION.md#2-the-vanishing-cycles-span-compact-homology)
includes the two added handles, separating their critical values
locally in disjoint Morse balls for the handle count. Thus all
finite vanishing cycles span compact first homology. Eventual
purity gives a zero inverse period by the
[same prime-height criterion](WINDING_ONE_CUBIC_QUADRUPLE_ROOT_INVERSE_PERIOD.md).
The period-variation arguments at the interior nodes and the three
nodes over zero give zero intersections there. The argument in §2
gives zero intersections at the two remaining nodes, since the
value \(\zeta_c\) is disjoint from those other critical values.
This would annihilate compact first homology, contradicting the
nonseparating inverse-branch cycle.

## 4. An exact repeated-center example

Take \(d=1\), \(b=-1-4x/3\), \(c=x^7-1-4x/9\). Its depressed
profiles are
\[
 p=-\frac43(1+x),\qquad q=x^7-\frac{16}{27}.
\]
It has least sphere degree seven by the
[sphere degree formula](MIXED_CUBIC_ZERO_CENTER_AND_CONDITIONAL_OFFCENTER.md),
and \(q^2/p^3\) varies. At \(x=-1\) the initial data in \(z=1/u\)
are \(p_-(0)=0\), \(q_-(0)=-43/27\), and \(p_-'(0)=4/3\), giving
both nonvanishing conditions at that boundary.

For the critical branch put \(V=3Y/2\) and \(D=V^2-1\). The equation
\(P_Y=0\) gives \(x=D\) and \(u=D/V^2\). The restriction of \(F\)
to the critical curve is
\[
 R(V)=V^2D^6-\frac{16}{27}\frac{V^2(V^2-V+1)}{V-1}.
\]
The canceled factor \(V+1\) gives the analytic center branch near
\(V=-1\). There \(R(-1)=8/9\), \(R'(-1)=-20/9\), and
\(du/dV=-2\). Hence \(\zeta_c=8/9\) and \(a=10/9\).

In the interior the critical points are the roots of the degree-14
polynomial \(H\) given by
\[
 R'(V)=\frac{2V}{27(V-1)^2}H(V),\qquad
 H(V)=27(V-1)^2D^5(7V^2-1)-8(3V^3-6V^2+4V-2).
\]
It is nonzero at \(V=0,\pm1\), by substitution. At an interior
critical point \(F_{YY}=6Y/u\ne0\) and \(du/dV=2/V^3\ne0\), so a
simple root of \(H\) is a Morse point by elimination of the
\(Y\)-critical equation in the Hessian. The following finite
arithmetic checks the required critical values.

Reduce modulo 29, where \(H\) retains degree 14. The nonzero
Euclidean divisors for \(H,H'\), starting at the derivative, have
degrees \((13,12,11,10,6,5,4,3,2,1,0)\). The class of \(R\) in
\(\mathbb F_{29}[V]/(H)\) has coefficients in increasing degree order
\[
 (15,13,25,12,7,20,6,16,13,14,11,13,9,16).
\]
This row can be checked by multiplication by the unit \(27(V-1)\)
and reduction modulo \(H\). For multiplication by each of its first
14 powers, sum the coefficient of \(V^j\) in the remainder of
\(V^j R^k\) over \(0\le j<14\). The resulting traces are
\[
 (26,13,1,9,10,1,6,25,16,22,20,27,11,21).
\]
Newton's identities give the monic multiplication characteristic
polynomial \(L\) with coefficients in increasing degree order
\[
 (16,21,18,11,20,6,28,27,15,1,21,4,27,3,1).
\]
For \(L,L'\) the Euclidean divisor degrees are
\((13,12,11,10,9,8,7,6,5,4,3,2,1,0)\). Also \(L(0)=16\) and
\(L(8/9)=26\) modulo 29. These computations lift to the required
nonvanishing statements over characteristic zero because the leading
coefficient of \(H\) and the denominator \(27(V-1)\) are units in the
localized quotient at 29. The eigenvalues of multiplication by
\(R\) after splitting the simple roots of \(H\) are precisely the
interior critical values. They are therefore nonzero, distinct from
one another, and distinct from \(\zeta_c\).
