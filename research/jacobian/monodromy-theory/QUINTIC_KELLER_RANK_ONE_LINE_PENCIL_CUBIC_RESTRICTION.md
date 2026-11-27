# Cubic vanishing in a top-rank-one line pencil

**Theorem QuinticRankOneLinePencilCubicRestriction.** Suppose a
degree-five three-variable Keller map is noninvertible, has top
rank one, and has primitive cancelling-pencil degree one in
[QuinticRankOneQuarticPencilResidual](QUINTIC_KELLER_REDUCED_BOUNDARY_PENCIL_RESIDUAL.md#rank-one-pencil).
Write \(Q,R\) for two independent cancelling components. The
common quartic factor \(g\) is nonreduced and has a unique repeated
line \(\ell\). Both cubic terms \(Q_3,R_3\) vanish identically on
\(\ell\).

Work over \(\mathbb C\). In degree one, the pencil consists of
lines through a base point \(O\). The valuation statement of the
[pencil restriction](QUINTIC_KELLER_REDUCED_BOUNDARY_PENCIL_RESIDUAL.md#rank-one-pencil)
also places every factor of \(g\) in that pencil. In particular
\(\ell\) passes through \(O\). Its multiplicity in a general
cancelling leading quartic is either two or three, since
\(\deg g\le3\) and general moving lines avoid the fixed
components.

Suppose the cubic restriction to \(\ell\) is not identically
zero for a general cancelling component, which we call \(R\).
We first show that at every boundary point away from \(O\), either
there is already a contradiction or the fibre residue extends
holomorphically on its resolution for general \(t\). A point
on a simple support line is smooth for the fibre away from
\(O\). On \(\ell\) away from \(O\) the fibre can be singular
only at a zero of \(R_3|_\ell\), finitely many points.
A rank-two or rank-three quadratic jet at any such point gives
the ordinary double-point resolution, as in the
[quartic boundary proof](KELLER_QUARTIC_COMPONENT_BOUNDARY_OBSTRUCTION.md#ordinary-double-points-of-the-plane-leading-curve).

## Points away from the pencil base point

If \(\ell\) has multiplicity three, the full quadratic jet
restricts to zero on the plane at infinity. Rank one would
force invertibility by the
[quadratic-source completion](KELLER_QUARTIC_QUADRATIC_SUSPENSIONS_AND_HIGHER_CONTACT.md#the-rank-one-normal-jet-of-a-quartic-component).
At a rank-zero point \(P\ne O\), put the source \(z\)-direction
at \(P\). As in the
[quartic triple-point proof](KELLER_QUARTIC_TRIPLE_BOUNDARY_AND_LINEAR_COMPONENTS.md#local-analysis-at-an-ordinary-triple-point),
vanishing of the linear and quadratic jets makes \(R=a(x,y)z+b(x,y)\).
Here \(a_3\) is a nonzero multiple of \(\ell^3\), because
\(R_4=\ell^3 L\) and the moving line \(L\) does not contain \(P\).
The
[cubic linear-source classification](KELLER_LINEAR_SOURCE_COEFFICIENT_CUBIC_CLASSIFICATION.md)
then leaves only its repeated-parallel residual: the two other
noninvertible forms have leading coefficient \(x^2y\) with two
distinct directions. In the repeated-parallel case
\(R-\mu=\lambda x+x^2H\), where the linear part of \(x\) defines
\(\ell\), because \(a_3\) is a cube in that direction. This forces
\(R_3|_\ell=0\), contrary to the choice of \(R\).
Thus all multiplicity-three points away from \(O\) are covered.

Suppose instead that \(\ell\) has multiplicity two. Around a
singular point \(P\ne O\) take local coordinates \(u,v,w\), with
\(w=0\) the plane at infinity and \(u=0\) the repeated line in
that plane. The leading equation restricts to \(u^2\) times
a unit on \(w=0\). If its full quadratic jet has rank one,
its square has a nonzero \(u\) coefficient. After a linear
change \(u\mapsto u+\delta w\) and a parameter-preserving
analytic square completion, the surface equation becomes
\(u^2+f(v,w)-tw^4=0\). Here
\[
              f(v,0)=0,\qquad [w]f(v,w)=R_3|_\ell(v).
\tag{1}
\]
Indeed the critical shift in \(u\) vanishes on \(w=0\), so
changes the contribution of the plane equation only by a
multiple of \(w^2\). The linear change in \(u\) also changes
this plane equation only by terms divisible by \(w^2\) after
restricting the new \(u\) to zero. The coefficient of \(w\)
is then precisely the restriction to \(\ell\) of the original
cubic term. Rank zero is impossible because of the nonzero
\(u^2\) coefficient.

A simple zero of \(R_3|_\ell\) gives a rank-two quadratic jet
even before the square completion. For a double zero, the
cubic \(f_3\) has a simple root at the direction \(w=0\), by
the nonzero \(v^2w\) term in (1). The first
[higher-contact test](KELLER_QUARTIC_QUADRATIC_SUSPENSIONS_AND_HIGHER_CONTACT.md#a-local-extension-at-higher-contact-double-points)
extends the residue. For a triple zero, \(f_4\) has nonzero
\(v^3w\) coefficient and \(f_3\) is divisible by \(w^2\).
If \(f_3=0\), the last higher-contact test gives a nonzero
period. If \(f_3\) is a nonzero multiple of \(w^3\), its
triple-direction test applies using the \(v^3w\) coefficient.
Otherwise \(f_3=w^2(av+bw)\) with \(a\ne0\). In the
double-direction test the coefficient \([w]f_4(1,w)\) is
nonzero and \(f_5(1,0)=0\) by (1). Its discriminant is
therefore nonzero. These local tests use only the vanishing
orders given in (1) and the isolated singularities of the
resulting surface. They apply here although the plane
leading curve has a repeated component.

## The pencil base point and the remaining affine-linear pair

Under the standing assumption on \(R_3|_\ell\) the projective
surface has only isolated singularities: the affine part
is smooth, every point at infinity outside \(O\) and the
repeated line is smooth, and on that line the singular
points are among the finitely many zeros of \(R_3\).
It is normal by the hypersurface depth condition and
regularity in codimension one. If the fibre is smooth at \(O\),
or its quadratic jet there has rank at least two, all
singularities extend holomorphically or give a local
contradiction for general \(t\). The global
[exact fibre form obstruction](KELLER_QUARTIC_COMPONENT_BOUNDARY_OBSTRUCTION.md#residues-on-a-resolution)
then rules out the map. Since the plane leading part has
multiplicity four at \(O\), a rank-one quadratic jet there
would be a multiple of \(w^2\), giving invertibility by
quadratic-source completion. Thus both full jets at \(O\)
vanish for a general cancelling component. Being linear in
the component, they vanish for both \(Q\) and \(R\).

Take the common source \(z\)-direction at \(O\). Both \(Q,R\)
are affine linear in \(z\), and the coefficients of \(z\) have
degree at most two, because both leading quartic forms are
binary. A nonzero reduced coefficient in any combination would
give the
[reduced quadratic completion](KELLER_QUARTIC_TRIPLE_BOUNDARY_AND_LINEAR_COMPONENTS.md#invertibility-with-a-reduced-quadratic-coefficient),
so all nonzero \(z\)-coefficients of the pencil are repeated
line squares. They must be proportional: two independent squares
of affine linear polynomials have generic sum of rank two
after homogenization. If both coefficients vanish, \(Q,R\)
are a planar Keller pair by the product formula for the full
determinant, and \(P_z\) is a nonzero constant. The planar
degree bound gives invertibility.

Otherwise a target change and affine source normalization give
\[
       Q=x^2z+B(x,y),\qquad R=r(x,y),\qquad B_x(0,y)=\lambda\ne0 ,
\]
with \(B(0,y)\) constant by the
[repeated-line classification](KELLER_LINEAR_SOURCE_COEFFICIENT_CUBIC_CLASSIFICATION.md#normal-residues-along-a-repeated-line).
At \(x=0\), the determinant equation becomes
\(\lambda P_z(0,y,z)\,r_y(0,y)=\pm\kappa\).
Hence \(P_z(0,y,z)\) is a nonzero constant and \(r(0,y)\) is
affine linear with nonzero slope.

Write \(m=\deg_z P\ge1\) and \(p_m=[z^m]P\).
The degree-\(m\) coefficient of the determinant gives
\[
              2m x r_y p_m-x^2\{p_m,r\}_{x,y}=0,
              \qquad \{p_m/x^{2m},r\}_{x,y}=0.
\]
The planar polynomial \(r\) has no critical point because
it is a Keller component independent of \(z\). It cannot
be a nontrivial polynomial composition: a root of the
outer derivative has nonempty inverse image under a
nonconstant inner polynomial. By the closed-polynomial
criterion of
[Wang–Yan, Proposition 2.2 and Lemma 2.3](https://arxiv.org/pdf/2610.36495v1),
these rational constants belong to \(\mathbb C(r)\).
Thus \(p_m=x^{2m}f(r)\). A denominator of \(f\) would divide
\(x^{2m}\) after substitution, by a coprime univariate Bézout
identity, but \(r(0,y)\) is nonconstant. Hence \(f\) is a
polynomial. The degree bound \(2m\le5-m\) forces \(m=1\),
and \(\deg r=4\) then forces \(f\) constant. This gives
\(P_z(0,y,z)=0\), a contradiction.

Therefore the assumption that a general cancelling cubic
restricts nontrivially to \(\ell\) is impossible.
Vanishing of both cubic restrictions follows by linearity.
