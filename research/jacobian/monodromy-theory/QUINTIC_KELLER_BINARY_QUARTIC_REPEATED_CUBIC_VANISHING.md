# Cubic vanishing on a repeated line of a binary quartic

**Theorem QuinticKellerBinaryQuarticRepeatedCubicVanishing.** Let \(F\)
be a noninvertible three-variable Keller map of ordinary degree at
most five in characteristic zero. Suppose a quartic component \(R\)
has geometrically nonreduced leading form \(R_4\) in two linear
source coordinates. Then \(R_3\) vanishes identically on at least
one repeated line of \(R_4\) over the algebraic closure.

The repeated support is concurrent at the binary pencil base point.
The argument tests the finitely many projective fibre singularities
under the contrary assumption that every repeated line has nonzero
cubic restriction. In particular the result applies to the primitive
degree-one case of the
[rank-two valuation identity](QUINTIC_KELLER_RANK_TWO_QUARTIC_PENCIL.md#valuations-along-a-projective-pencil).
There all pencil members are lines in two binary coordinates, and
every factor of the cancelling quartic lies in one of those members.

## Rank-zero jets and the planar case

Work over \(\mathbb C\) by embedding the finitely generated
coefficient field. Suppose \(R_3\) is nonzero on every repeated
line. A singular fibre point at infinity with zero quadratic
jet, written as the \(z\)-direction in affine source variables,
gives \(R=a(x,y)z+b(x,y)\). In fact a monomial involving \(z^k\)
has order \(4-k\) in the affine chart of the degree-four projective
fibre at that point. Vanishing of the linear and quadratic jets
removes the terms of \(z\)-degree at least two.

We will use this observation at points where the cubic leading
coefficient \(a_3\) is zero or a nonzero cube of a linear form.
If \(a\ne0\) the
[linear-source classification](KELLER_LINEAR_SOURCE_COEFFICIENT_CUBIC_CLASSIFICATION.md)
then gives either invertibility or its repeated-parallel
residual. Indeed the other two residuals have a leading cubic
coefficient with two distinct linear factors. In the
repeated-parallel residual, after an affine change,
\[
                          R-\mu=\lambda x+x^2 H .
\]
The repeated direction \(x\) divides \(R_3\) and its square divides
\(R_4\), even after undoing the translation: translating the degree-four
term divisible by \(x^2\) contributes to the cubic only terms
divisible by the homogeneous linear part of \(x\).
This contradicts the standing assumption.

If \(a=0\), \(R\) is a polynomial in two variables. Its two partial
derivatives must have no common affine zero by the Keller condition.
Choose planar linear coordinates making both partials degree three.
Their cubic homogenizations have no common curve component, since
neither is divisible by the homogenizing coordinate \(w\) and an
affine common zero is excluded. By Bézout they have total intersection
multiplicity nine. Every point of intersection at infinity is a
repeated root of the leading binary quartic by Euler's identity.

At a repeated root use an affine chart with planar coordinate \(s\)
vanishing at the root, and write \(G(s,w)\) for the homogenized
quartic in this chart. The gradient intersection is given by
\(G_s=0,\ 4G-w G_w=0\), using Euler to add a multiple of the first
equation to the second. If the root has multiplicity \(m\ge2\), the
second equation has order \(m\) along \(w=0\) and its \(w\)-derivative
at the root is \(3R_3\) at that projective direction. This is nonzero
by the standing assumption, since \(R_3\) is binary in the planar
coordinates and does not vanish on that repeated line. Thus the
second equation can be solved for \(w\) as a series in \(s\) of
order \(m\). Along this series \(G_s\) has order \(m-1\) from the
derivative of \(R_4\). The sum of all such local multiplicities
is at most three, contrary to Bézout. This excludes the planar case.

## The projective fibre

Under the standing assumption every projective fibre is smooth in
the affine part and can be singular only at the base point \(O\)
or at the finitely many cubic zeros on the repeated lines. A point
on a simple support line away from \(O\) is smooth. The fibre is
normal because it is a hypersurface with isolated singularities.
Use the
[ordinary double-point resolution](KELLER_QUARTIC_COMPONENT_BOUNDARY_OBSTRUCTION.md#ordinary-double-points-of-the-plane-leading-curve)
for a quadratic jet of rank at least two at any singular point.

At a point on a line of multiplicity two away from \(O\), the
local plane leading equation is a square times a unit. The local
calculation in the
[line-pencil residue proof](QUINTIC_KELLER_RANK_ONE_LINE_PENCIL_CUBIC_RESTRICTION.md#points-away-from-the-pencil-base-point)
applies to the simple, double and triple zeros of the nonzero
restricted cubic. The simple zero gives rank two. At a double
or triple zero the higher-contact tests either extend the fibre
residue on a resolution or contradict its exactness by a local
period. That calculation requires the local repeated line and the
cubic orders, with no assumption on the other two target components.

The other singular points are \(O\), or lie on a line of
multiplicity at least three away from \(O\). Their plane leading
equation has multiplicity at least three. A rank-one full quadratic
jet there is a multiple of \(w^2\), so the
[quadratic-source completion](KELLER_QUARTIC_QUADRATIC_SUSPENSIONS_AND_HIGHER_CONTACT.md#the-rank-one-normal-jet-of-a-quartic-component)
gives inversion. At a rank-zero point use the linear-source
argument above. At \(O\) the leading quartic is independent of its
source direction, giving \(a_3=0\). At a point on a triple line away
from \(O\), the other factor of the leading quartic does not pass
through the point. Taking the source direction at the point
gives \(a_3\) a nonzero cube from the triple line. If the line has
multiplicity four, the leading quartic is independent of its
source direction and again \(a_3=0\).

Thus for a generic fibre either one of these local arguments gives
a contradiction or the fibre residue extends holomorphically on
a resolution at every boundary point. The
[global exactness obstruction](KELLER_QUARTIC_COMPONENT_BOUNDARY_OBSTRUCTION.md#residues-on-a-resolution)
rules out the latter by the nonzero holomorphic two-form on the
complete resolution. This proves cubic vanishing on at least one
of the repeated lines.

## The discriminant of a quadratic scalar fibre

**Proposition KellerBinaryQuadraticScalarSquareDiscriminant.** Suppose
a characteristic-zero three-variable Keller triple has a component
\[
                 R=a(x)y^2+b(x)yz+c(x)z^2+d(x)y+e(x)z+f(x).
\]
If the discriminant \(\Delta=b^2-4ac\) is nonzero, it is geometrically
a square in the rational function field of \(x\).

Work over \(\mathbb C\). Suppose that \(\Delta\) is nonsquare, so
\(\mathbb C(x,\eta)\), \(\eta^2=\Delta\), is a quadratic function
field. After a constant change in \(y,z\), take \(a\ne0\).
On a nonempty open set in this double cover, translate \(y,z\)
by the rational stationary point of their quadratic to write
\[
             R=UV+K(x),\qquad K(x)\in\mathbb C(x).
\]
In detail, after the translation, put
\(U=a y+(b+\eta)z/2\) and \(V=y+(b-\eta)z/(2a)\).
Their binary Jacobian is \(-\eta\). The affine residue form on
\(R=t\) therefore pulls back up to sign to
\[
                            \frac{dx\wedge dU}{\eta U}.
\]
It has a regular primitive \(\gamma_t\) from the other two
Keller components by the
[exact fibre residue identity](KELLER_QUARTIC_COMPONENT_BOUNDARY_OBSTRUCTION.md#residues-on-a-resolution).
This uses the polynomial primitive formed by one component times
the differential of the other, scaled by the constant Jacobian.

On \(U\ne0\) set \(V=(t-K)/U\). The coefficients of \(\gamma_t\)
are Laurent polynomials in \(U\) over \(\mathbb C(x,\eta)[t]\).
Writing its \(dU\)-coefficient at \(U^{-1}\) as \(A(x,\eta,t)\),
the equality \(d\gamma_t=dx\wedge dU/(\eta U)\), with the sign
chosen accordingly, gives
\[
                         \partial_x A(x,\eta,t)=1/\eta .
\]
The \(x\)-derivative is at fixed \(t\), and the \(U\)-derivative
of a Laurent polynomial has zero coefficient at \(U^{-1}\).
The derivation on \(\mathbb C(x,\eta)\) has constants
\(\mathbb C\): any nonconstant rational function of this
curve is transcendental over \(\mathbb C\) and gives a finite
separable extension of its own rational function field.
Consequently \(A(x,\eta,t)=L(x,\eta)+H(t)\) with \(H\in\mathbb C[t]\)
and \(L_x=1/\eta\).

At a generic point of the double cover specialize \(t=K(x)\).
The line \(V=0\) in this fibre at fixed \(x,\eta\) contains the
point \(U=0\) where the primitive is regular, since the linear
change from \(y,z\) has nonzero determinant there. Its residue at
\(U=0\) is thus zero:
\[
                           L(x,\eta)+H(K(x))=0 .
\]
This is an identity of rational functions on the double cover.
It would make \(L\in\mathbb C(x)\), whereas \(L_x=1/\eta\) changes
sign under its nontrivial involution. This contradiction proves
the proposition.

**Corollary QuinticKellerBinaryQuarticDoubleLineAffineSource.** Under
the hypotheses of the cubic vanishing theorem, suppose a repeated
line on which \(R_3\) vanishes has multiplicity exactly two.
Then \(R\) is affine linear in the source direction of the binary
pencil base point. Either \(R\) is planar or, in an affine frame with
that repeated line at infinity given by \(x=0\),
\(R-\mu=\lambda x+x^2 H\) with \(\lambda\ne0\).

Choose \(x,y\) so the line at infinity is \(x=0\) and
\(R_4=x^2(a_0x^2+b_0xy+c_0y^2)\) with \(c_0\ne0\).
The free variable \(z\) is the base-point direction.
Since \(x\mid R_3\), the component \(R\) is quadratic in
\(y,z\) together. Its coefficients of \(y^2,yz,z^2\) are
\(a(x),b(x),c(x)\) with \(\deg a=2\), \(\deg b,\deg c\le1\).
If \(c(x)\) is nonconstant, \(b^2-4ac\) has degree three,
contradicting the square discriminant proposition.
A nonzero constant \(c\) gives inversion by
[quadratic-source completion](KELLER_QUARTIC_QUADRATIC_SUSPENSIONS_AND_HIGHER_CONTACT.md#the-rank-one-normal-jet-of-a-quartic-component).
Hence \(c=0\) and \(R\) is affine linear in \(z\), with
coefficient of degree at most two. If this coefficient is nonzero,
the [quadratic linear-source classification](KELLER_QUARTIC_TRIPLE_BOUNDARY_AND_LINEAR_COMPONENTS.md#invertibility-with-a-reduced-quadratic-coefficient)
forces it to be a repeated affine square. Its homogeneous
quadratic part is divisible by \(x\), since \(x\mid R_3\),
so its linear direction is \(x\). The
[normal residue classification](KELLER_LINEAR_SOURCE_COEFFICIENT_CUBIC_CLASSIFICATION.md#normal-residues-along-a-repeated-line)
gives the stated repeated-line form after translating that square.

## Excluding the planar multiplicity-two case

**Corollary QuinticKellerBinaryDoubleLinePlanarExclusion.** Under the
hypotheses of the preceding corollary, assume the two degree-five
components have binary leading forms in the same directions as
\(R_4\). Then the planar case for \(R\) is impossible. Thus
\(R=\mu+\lambda x+x^2(z+h(x,y))\) in an affine frame, after
scaling the free variable, with \(\lambda\ne0\) and with nonzero
coefficient of \(y^2\) in \(h\).

Suppose \(R\) is planar. Its degree in \(y\) is two by
\(x\mid R_3\) and the form of \(R_4\). Write
\(R=A(x)y^2+B(x)y+C(x)\), where \(A\) is quadratic,
\(\deg B\le3\), and \(\deg C\le4\). Off the roots of \(A\),
a critical point is exactly a zero of \(K'(x)\), where
\(K=C-B^2/(4A)\), obtained by taking \(y=-B/(2A)\).

If \(A\) has two distinct roots, \(B\) is nonzero at both.
Otherwise \(R_y\) vanishes along that line and \(R_x\) is
quadratic in \(y\) with leading coefficient \(A'\ne0\),
giving a critical point. Thus \(K\) has a simple pole at each
root of \(A\). The polynomial \(A^2 K'\) has a nonzero value
at each root, and has no root elsewhere by submersion.
It is therefore a nonzero constant. But the reciprocal
of \(A^2\) has a nonzero simple-pole coefficient at each
of the two roots, whereas \(K'\) is a rational derivative.
This excludes distinct roots.

Translate to \(A=a x^2\). If \(B(0)=b\ne0\), then
\(x^3 K'\) is a polynomial with nonzero value at zero
and no other zero. Hence \(K=\beta+\theta/x^2\).
Its missing simple-pole term gives \(B'(0)=0\), since
\(C\) is polynomial. An affine shift of \(y\) now absorbs
\((B-b)/(2a x^2)\), and the resulting scalar is
\(R=a x^2 y^2+b y+\beta\).
If \(B(0)=0\), then \(R_y=0\) along \(x=0\), so
\(B'(0)=0\), \(C'(0)\ne0\) by submersion there.
An affine shift of \(y\) absorbs \(B/(2a x^2)\).
Now \(K\) is polynomial and \(K'(0)\ne0\).
Since its derivative has no zero anywhere, this gives
\(R=a x^2 y^2+\lambda x+\beta\), with \(\lambda\ne0\).

For both scalars let \(v=xy\), and use the fibre parameter
\(t=R-\beta\). In the first case \(y=(t-av^2)/b\),
\(x=bv/(t-av^2)\); in the second \(x=(t-av^2)/\lambda\)
and \(y=v/x\). Thus \(\mathbb C(x,y)=\mathbb C(t,v)\).
The derivation \(D=R_y\partial_x-R_x\partial_y\) has
constants \(\mathbb C(R)\), and \(dx/R_y\) at fixed \(t\)
equals \(dv/(t-av^2)\) in the first case and
\(-dv/(t-av^2)\) in the second. Both have nonzero residues
at the two simple roots over the algebraic closure of
\(\mathbb C(t)\). Consequently \(D\) has no rational
primitive of a nonzero constant.

The full determinant is, up to sign,
\(D(P)Q_z-P_zD(Q)=\kappa\ne0\).
The absence of a primitive and the rational constants
give the three degree-independent restrictions in the
[exceptional scalar argument](QUINTIC_KELLER_CUBIC_COMPONENT.md#the-exceptional-scalar-through-degree-five):
neither component is independent of \(z\), neither has
constant leading \(z\)-coefficient, and the pair cannot
be affine in \(z\).
If their \(z\)-degrees are \(m,n\) with coefficients
\(p_m,q_n\), the highest bracket yields
\(p_m^n/q_n^m\in\mathbb C(R)\).
The leading quintics are binary, hence \(\deg p_m\le4-m\)
and \(\deg q_n\le4-n\). The ratio is constant, since a
finite zero or pole of a nonconstant rational function of
\(R\) would require a full fibre factor of degree four to
divide a coefficient of degree at most three. The fibres
are squarefree by submersion.

If the degrees agree, subtract a constant multiple of
one component from the other to lower the degree.
This preserves the binary leading forms and reduces
to unequal degrees. The larger is at most three
by the exclusion of a constant coefficient. If it is
three, then the degree comparison
\(n\deg p_3=3\deg q_n\) is impossible for \(n=1,2\).
The remaining case has degrees two and one.
Proportional leading coefficients \(p_2\propto q_1^2\)
allow subtraction of a multiple of \(Q^2\), leaving an
affine pair or a component independent of \(z\),
both impossible. This eliminates planar \(R\).
