# The double-line conic pencil in top rank one

**Theorem QuinticRankOneDoubleLineConicPencilExclusion.** Let
\(F=(P,Q,R)\) be a complex three-variable Keller map of ordinary
degree five, with leading parts
\[
 P_5=z\,G(z^2,q),\qquad Q_4=z^4,\qquad R_4=z^2q ,
\]
where \(G\) is a binary quadratic and \(q\) is a ternary quadratic
such that the pencil \(\langle z^2,q\rangle\) has general
irreducible conic. Then \(F\) is invertible.

The top-rank-one pencil case with a repeated fixed line admits
this leading form. The conclusion requires the genuine original
degree-five bound.

## Cubic vanishing along the fixed line

Suppose the map is noninvertible. Replace \(R\) by a general
\(R+\lambda Q\), preserving the displayed form with new \(q\)
and \(G\), so its conic is irreducible. If
\((R_3+\lambda Q_3)|_{z=0}\) were not identically zero, the
projective fibre would have isolated singularities. Indeed away
from the repeated line the leading curve is smooth, and on the line
singularities occur only at finitely many zeros of the cubic. Away
from the line–conic intersections these points are resolved,
or furnish a nonzero period, by the same simple, double, and
triple-zero tests as in the
[line-pencil proof](QUINTIC_KELLER_RANK_ONE_LINE_PENCIL_CUBIC_RESTRICTION.md#points-away-from-the-pencil-base-point).
That proof applies because the leading plane equation at the
point is the square of the line times a unit and the coefficient
of the homogenizing variable restricts to the given cubic.

At a line–conic intersection, a smooth fibre point or a rank-two
or rank-three quadratic jet is covered by the
[quartic boundary proof](KELLER_QUARTIC_COMPONENT_BOUNDARY_OBSTRUCTION.md#ordinary-double-points-of-the-plane-leading-curve).
The plane multiplicity is three, so a rank-one quadratic jet is
a multiple of the square of the homogenizing variable and gives
invertibility by the
[quadratic-source completion](KELLER_QUARTIC_QUADRATIC_SUSPENSIONS_AND_HIGHER_CONTACT.md#the-rank-one-normal-jet-of-a-quartic-component).
A rank-zero jet in a singular fibre gives an affine-linear
source component \(R=a(x,z)y+b(x,z)\) by taking the direction
\(y\) at that intersection. Here \(a_3=z^2 L(x,z)\), where
\(L\) is the nonzero linear part of the conic at the chosen
direction. If the conic is transverse to \(z=0\), \(L\) is
independent of \(z\). The
[cubic source-coefficient classification](KELLER_LINEAR_SOURCE_COEFFICIENT_CUBIC_CLASSIFICATION.md)
excludes the repeated-parallel residual, which has a cube for
\(a_3\). In its other two residuals \(b\bmod a\) has a representative
of degree at most three. For the line and hyperbola residual
one obtains such a representative by interpolating the two
affine restrictions with \(xy,1-xy\) in its normalized binary
frame. Thus \(b_4=a_3 c_1\) for a linear \(c_1\), possibly zero,
which would make the conic \(L(y+c_1)\) reducible.

If the conic is tangent to \(z=0\), then \(a_3\) is a nonzero
multiple of \(z^3\). Only the repeated-parallel residual has
this type of leading coefficient, and its repeated-line
coordinate has linear part proportional to \(z\).
It makes \(R-\mu=\eta s+s^2H\), where \(s\) is the affine
coordinate normal to that line and \(\eta\ne0\). Therefore
\(R_3\) vanishes on \(z=0\), contrary to the standing assumption.
This excludes all remaining bad points.

The global exactness obstruction then rules out the contrary
generic cubic restriction. By linearity we have
\[
                         Q_3=zU,\qquad R_3=zV .
\tag{1}
\]
Here \(U,V\) are homogeneous quadratics.

## The exceptional discriminant

Use the subscript zero for restriction to \(z=0\), and put
\(r=q_0,\ D=(Q_2)_0,\ W=(R_2)_0\).
For a component \(R+\lambda Q\), blow up the line \(z=w=0\)
in the projective target-fibre compactification, where \(w\)
homogenizes the affine source. On the chart \(z=w Z\), the
exceptional curve has equation
\[
 rZ^2+(V_0+\lambda U_0)Z+W+\lambda D=0 .
\tag{2}
\]
In this formula each coefficient is dehomogenized on the base
line of the blowup. The exceptional curve is a double cover
with binary quartic discriminant
\[
            \Delta_\lambda=(V_0+\lambda U_0)^2-4r(W+\lambda D).
\tag{3}
\]
If \(\Delta_\lambda\ne0\), at a general point of this curve the
surface is smooth and the residue is a nonzero constant times
\(dx/\sqrt{\Delta_\lambda}\). Indeed the ambient volume gains
a factor \(w\) from the blowup and the equation gains \(w^2\).
The affine projective residue of a quartic before this blowup
is a dualizing generator. Thus there is a simple pole along
the exceptional curve with residue given by the derivative of
(2) in \(Z\). The formula transforms regularly across the
complementary chart of the base line away from the roots.

If the discriminant has four simple roots, the residue is a
nonzero holomorphic form on the complete double-cover curve,
and has a nonzero period. If it has a double root and no higher
root, the residue has a nonzero simple pole at a normalized
branch over the double root. The corresponding small loop
has nonzero period. These arguments also work when the
exceptional curve has two components. The loops can be chosen
away from the strict infinity boundary, the other exceptional
components and the finite bad set of (2). A shrinking normal
tube contradicts exactness on the original affine fibre, as in the
[quartic residue proof](KELLER_QUARTIC_COMPONENT_BOUNDARY_OBSTRUCTION.md#residues-on-a-resolution).
This local tube only needs the blowup surface to be smooth
along the chosen loop, even when the original projective
closure is singular along the blowup centre.

Thus for general \(\lambda\), the discriminant is either
identically zero or has a triple root. In the latter case
the triple factor is fixed. Indeed a nonzero quartic with a
triple root has a unique such root. Over \(\mathbb C(\lambda)\)
that root and the remaining linear factor are rational, since
each is determined by the quartic. They define morphisms
from the parameter projective line to the root projective line.
If their degrees are \(s_1,s_2\), the coefficient map
\(\lambda\mapsto\Delta_\lambda\) to the space of binary quartics
has pullback hyperplane degree \(3s_1+s_2\). This follows
from the base-point-free multiplication map
\((L,M)\mapsto L^3M\), of bidegree \((3,1)\).
But (3) gives a coefficient map of degree at most two after
removing any common scalar factor. Thus \(s_1=0\).

If (3) is identically zero in \(\lambda\), comparing
coefficients gives \(U_0=D=0\). Otherwise let \(L^3\) be the
fixed factor. The quadratic coefficient in \(\lambda\) shows
that \(U_0=kL^2\) for a scalar \(k\), possibly zero.
If \(L\nmid r\), the linear coefficient then gives \(L^2\mid D\).
If \(L\mid r\), the constant coefficient first gives \(L\mid V_0\).
The linear coefficient now gives \(L^3\mid rD\).
For squarefree \(r\) this gives \(L^2\mid D\) again.
If \(r\) is proportional to \(L^2\), it gives \(L\mid D\).
These possibilities are distinguished in the next step.

## The remaining affine quadratic form in the cancelling component

Choose \(x=L\) when one of \(U_0,D\) is nonzero, and any affine
linear binary frame otherwise. Equations (1)–(3) put \(Q\) in
the form
\[
 Q=f(z)+A(z)x+(kz+d)x^2+(a(z)+\eta_0 x)y,\qquad
 \deg A,\deg a\le2 ,
\tag{4}
\]
where \(f\) is monic quartic. In the case \(\eta_0\ne0\), the
coefficient of \(y\) is a squarefree polynomial in \(x,z\)
of degree at most two, so the
[reduced quadratic completion](KELLER_QUARTIC_TRIPLE_BOUNDARY_AND_LINEAR_COMPONENTS.md#invertibility-with-a-reduced-quadratic-coefficient)
gives invertibility. Thus \(\eta_0=0\).

If \(a\ne0\), the
[parallel-coefficient theorem](KELLER_LINEAR_SOURCE_COEFFICIENT_CUBIC_CLASSIFICATION.md)
forces the noninvertible case to have a repeated root \(z=\xi\)
and
\[
 Q=\mu+\eta(z-\xi)+(z-\xi)^2 H,\quad \eta\ne0,\quad \deg H\le2 .
\]
Since \(Q_4=z^4\), the quadratic part of \(H\) is \(z^2\).
In particular \(k=d=0\) in (4).
If \(a=0\) and \(k\ne0\), first translate \(z\) to absorb \(d\).
If \(A(0)=0\), then \(Q_x\) vanishes identically on \(z=0\)
and \(Q_z\) is quadratic of nonzero leading coefficient in
\(x\) there, giving a critical point. If \(A(0)\ne0\), solve
\(Q_x=0\) for \(z\ne0\) by \(x=-A(z)/(2kz)\).
The remaining equation is the derivative of
\(f(z)-A(z)^2/(4kz)\). Its numerator after multiplying by
\(z^2\) has degree five and nonzero constant coefficient,
so has a nonzero root and yields a critical point.
If \(a=0,\ k=0,\ d\ne0\), the
[quadratic-source completion](KELLER_QUARTIC_QUADRATIC_SUSPENSIONS_AND_HIGHER_CONTACT.md)
gives invertibility. The only remaining case in (4) has
\(k=d=0\). When \(a=0\), either \(A\) is constant nonzero,
giving a graph coordinate, or the absence of critical points
forces a double root \(z=\xi\) of \(A\) and \(f'(\xi)\ne0\).
This is the same displayed residual form.

In all remaining cases, after affine normalization and a
quadratic shear in the source variable \(x\) alone, \(Q\) is the
exceptional cubic \(r_*=z+z^2 x\). In detail, after translating
\(z\) and scaling, \(H=z^2+\alpha x+\beta y+\gamma z+\delta\)
with \((\alpha,\beta)\ne(0,0)\), as otherwise \(Q\) is
univariate. Choose \(\alpha x+\beta y\) as the binary
coordinate called \(x\), and add the univariate polynomial
\(z^2+\gamma z+\delta\) to it. The other binary coordinate
called \(y\) is unchanged by this shear.

## Final determinant obstruction

After the shear, write the bracket along \(Q\) as
\[
 D(P)R_y-P_yD(R)=\pm\kappa,\qquad
 D=-z^2\partial_z+(1+2zx)\partial_x.
\]
The rational kernel of \(D\) is \(\mathbb C(r_*)\), since
\(1/z\) is a rational parameter with \(D(1/z)=1\).
We use the two degree-independent restrictions from the
[exceptional-cubic calculation](QUINTIC_KELLER_CUBIC_COMPONENT.md#the-exceptional-scalar-through-degree-five):
neither component can be independent of \(y\) or have constant
leading \(y\)-coefficient, and both cannot be affine in \(y\).
Those restrictions follow there from the absence of a polynomial
\(D\)-primitive of a nonzero constant.

Write \(m=\deg_y P\), \(n=\deg_y R\), with leading coefficients
\(p_m,r_n\). Since the shear changes only \(x\) by a polynomial
in \(z\), it preserves the individual \(y\)-degrees and pulls back
the leading coefficients. All ordinary degrees of the leading
coefficients below are measured in the affine frame before that
shear. The highest bracket coefficient gives
\[
                         p_m^n/r_n^m\in\mathbb C(r_*).
\tag{5}
\]
Since \(R_3\) is divisible by \(z\), \(n\le2\). If \(n=1\)
then \(q_{yy}=0\) and \(r_n\) has degree three: the leading
coefficient in \(z^2q\) is \(z^2\) times a nonzero linear form.
That form is nonzero because the irreducible conic \(q\)
cannot be independent of \(y\). If \(n=2\)
then \(r_n\) is a polynomial in \(z\) of degree at most two.

For \(m\ge2\), both leading coefficients have degree at most
three. Equation (5) must then be constant. Indeed a zero or
pole at \(c\) of a nonconstant rational function of \(r_*\)
would make the full fibre polynomial \(Q-c\) divide the
corresponding leading coefficient in the original affine
frame. This follows by a coprime numerator-denominator
identity and the fact that \(Q-c\) is squarefree, by
submersion. Its degree is four. For \(m=1,n=2\), the only
other possibility is a full fibre factor in \(p_1\).
Reducing the determinant bracket modulo \(r_*-c\) then gives
\(D(p_0)(2r_2y+r_1)=\pm\kappa\) on the reduced fibre, and
forces that same factor in \(r_2\), impossible. Thus (5) is
constant in all cases except \(m=n=1\), already excluded.

If \(n=1\) and \(m\ge2\), proportionality then gives
\(\deg p_m=m\deg r_n=3m>5-m\), impossible.
If \(n=2\), the coefficient \(r_n\) is nonconstant by the
constant-leading-coefficient restriction. The leading form
of \(P\) has \(y\)-degree at most four. If \(m=4\), its leading
coefficient has degree at most one, which is incompatible
with \(2\deg p_4=4\deg r_n\). If \(m=3\), the degree-five
part of \(P\) has \(y\)-degree at most two. In fact the only
term that could supply \(y^3\) is the \(zq^2\) part with
nonzero \(q_{yy}\), which also supplies \(y^4\). Thus
\(\deg p_3\le1\), incompatible with \(2\deg p_3=3\deg r_n\).
For \(m=2\), subtract the proportional multiple of \(R\)
from \(P\) to reduce \(m\) to at most one. At \(m=1,n=2\),
(5) makes \(r_2\) a constant multiple of \(p_1^2\).
Subtract the corresponding multiple of \(P^2\) from \(R\).
This leaves two affine components in the free variable,
or a component independent of it, both impossible for the
exceptional cubic. The contradiction proves the theorem.
