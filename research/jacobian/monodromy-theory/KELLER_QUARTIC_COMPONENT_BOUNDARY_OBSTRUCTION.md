# Exact fibre forms obstruct mild quartic boundaries of Keller components

**Theorem KellerComponentLogarithmicFibreObstruction.** Let \(F=(P,Q,R)\)
be a complex three-variable polynomial Keller map, with \(\deg R\ge4\).
For any \(t\in\mathbb C\), the projective closure of \(R=t\) cannot be
normal with only log canonical singularities. More generally the
projective residue form on this closure cannot become a logarithmic
two-form on a smooth projective resolution that is an isomorphism
over the affine fibre.

**Corollary KellerQuarticMildBoundaryExclusion.** Let \(R\) be a quartic
component of such a map. Its projective plane leading curve
\(\{R_4=0\}\) cannot be reduced with all singularities ordinary nodes,
ordinary cusps or ordinary tacnodes. Thus the same condition excludes
a degree-four constant target-linear combination in a three-variable
quintic Keller map. These exclusions extend to characteristic zero
whenever the hypotheses hold over an algebraic closure.

Write \(\kappa=\det JF\in\mathbb C^\times\).
On the smooth affine fibre \(U_t=\{R=t\}\), the residue \(\alpha_t\) of
\(dx\wedge dy\wedge dz/(R-t)\) satisfies
\(\kappa\alpha_t=dP\wedge dQ\). In particular its restriction to
every component of \(U_t\) is nonzero and exact. This is the
[exact fibre calculation of Wang and Yan, Section 3.5.1](https://arxiv.org/pdf/2610.36495v1),
which holds in all degrees of the three components.

## Residues on a resolution

With \(d=\deg R\), let \(X_t\subset\mathbb P^3\) be the projective
hypersurface given by the homogenization \(G_t=0\) of \(R-t\).
The global adjunction section \(w^{d-4}\) of
\(\omega_{X_t}\simeq\mathcal O_{X_t}(d-4)\) gives \(\alpha_t\)
on \(w=1\), by Poincaré residue. It is regular on the smooth
locus of \(X_t\). If \(X_t\) is normal and log canonical, its
pullback to a log resolution has at most simple poles: in a
trivialization of the Cartier canonical bundle the orders of
the pullback of a local generator along exceptional divisors
are the discrepancies, which are at least \(-1\).
The additional factor from the regular section can only increase
the orders. Choose the resolution to make the exceptional curves
and the transform of \(w=0\) a simple normal crossings divisor.
A meromorphic two-form with at most simple poles along that
divisor is logarithmic.

Here is a direct obstruction to exactness on a complement of
a simple normal crossings divisor \(D\) in a smooth projective
surface \(Y\). For a component \(C\) of \(D\), the residue of a
global logarithmic two-form is a meromorphic one-form on \(C\)
with at most simple poles at the crossings. Integration over a
small normal circle bundle above a loop on the smooth part of
\(C\) gives \(2\pi i\) times the period of its residue, by the
one-variable residue formula. If the two-form were exact
on \(Y\setminus D\), each such period would vanish by Stokes.
In particular small loops about crossings kill all simple
poles. The resulting holomorphic one-form on the compact curve
has all periods zero and hence vanishes, since its integral
is then a global holomorphic function. This holds on every
component, so the two-form has no poles.

A nonzero holomorphic two-form \(\omega\) on \(Y\) also cannot
be exact on \(Y\setminus D\). Indeed the kernel of the
restriction of \(H^2(Y,\mathbb C)\) is spanned by the divisor
classes of the components of \(D\). For example, by cohomology
with supports and Poincaré duality the kernel is the image
of \(H_2(D,\mathbb C)\); the components generate that homology
since their intersections are finite. These classes pair to
zero with \([\overline\omega]\) since a \((0,2)\)-form restricts
to zero on a complex curve. But
\(\int_Y\omega\wedge\overline\omega\ne0\).
Applying this to the extension of \(\alpha_t\) proves the theorem.

## Ordinary double points of the plane leading curve

For \(d=4\), the form is a generator of the dualizing bundle
near every point of \(X_t\). A reduced plane leading curve has
finitely many singular points. Every affine point of \(X_t\)
is smooth, and at infinity every smooth point of the leading
curve is smooth on \(X_t\) by the tangential derivatives.
Thus all singularities of \(X_t\) are isolated when the
leading curve is reduced. The hypersurface is normal by the
codimension-one regularity criterion. In a projective chart
near a singular point at infinity the equation is
\(G_0(u,v,w)-tw^4=0\), where \(w=0\) is the plane at infinity.
The Hessian is independent of \(t\).

If the Hessian has rank at least two, split off a nondegenerate
pair of variables by the analytic implicit function theorem and
quadratic completion. The germ is \(ab+\psi(c)=0\); isolatedness
forces \(\psi\ne0\). Replacing \(c\) by a local parameter times a
unit gives \(ab+c^j=0\) for some finite \(j\). The residue form
extends holomorphically to a resolution. To see this, ordinary
point blowups reduce \(j\) by two in the possibly singular chart,
with smooth other charts. Each blowup of the double point has
zero canonical discrepancy: in the residue expression the
factor from the ambient three-form and the factor from the
equation are both the square of the exceptional coordinate.
This covers the nodes of the leading curve for every \(t\).

At a cusp or tacnode, if the germ is singular and the Hessian
has rank one, split off the variable \(u\) with nonzero square
term while fixing \(v,w\). This splitting is achieved by solving
the \(u\) critical equation, translating, then taking the square
root of the remaining quadratic unit. Its change of coordinates
is independent of \(t\), since the term \(-tw^4\) is independent
of \(u\). The result is
\[
                 u^2+f_t(v,w)=0,\qquad f_t=f_0-tw^4,
                 \qquad \operatorname{ord}f_0\ge3 .
\]
At \(w=0\), \(f_0(v,0)\) has order three at a cusp and order
four at a tacnode. Normality and the residue argument are
unaffected by the local analytic coordinates, which change a
dualizing generator by a unit.

Suppose the cubic term \(f_3\) of \(f_t\) is nonzero. Blow up
the surface double point once. In the charts over the
exceptional line \([v:w]\), the strict transform has equation
\[
          U^2+z f_3(V,1)+z^2 f_{4,t}(V,1)+O(z^3)=0
\]
when \(w=z,\ v=zV,\ u=zU\), and the analogous equation in the
chart with \(v=z\). Its possible singularities are at roots of
\(f_3\) on the exceptional line. A simple root produces a
nondegenerate quadratic term \(U^2+z(V-V_0)\). At a multiple
root with \(w\ne0\), the coefficient \(f_{4,t}(V_0,1)\) is
nonzero for all but one \(t\), since its coefficient of \(t\)
is \(-1\). At a root with \(w=0\) there is a nonzero quartic
coefficient in the \(v\)-chart by the tacnode hypothesis; at
a cusp this direction cannot be a root of \(f_3\). Thus at
every multiple root the Hessian has rank at least two for
general \(t\). These germs are isolated: after splitting off
\(U\) and \(z\), a root of multiplicity \(j\) gives residual
leading order \(2j\), with nonzero coefficient from completing
the square in \(z\). Hence all singularities after the first
blowup have the preceding holomorphic-resolution argument.
The first blowup was itself a double-point blowup with zero
discrepancy.

If \(f_3=0\), the plane germ is a tacnode. For general \(t\)
the binary quartic \(f_{4,t}\) has distinct roots. Indeed
\(f_{4,t}(V,1)\) has degree four with fixed nonzero leading
coefficient and varies only in its constant, so a repeated
root can occur at only finitely many values. There is no root
with \(w=0\). Blow up the base \((v,w)\) and use the normalized
double cover, so \(v=wV,\ u=w^2 U\). Along the exceptional
curve the equation is \(U^2+f_{4,t}(V,1)=0\).
The double cover of the exceptional projective line is
smooth, with four simple branch points, and the surface
model is smooth along it by the same derivative test. The
pullback residue form is a unit times
\[
                    \frac{dV\wedge dw}{2 w U}.
\]
It has a simple pole along the exceptional curve with
nonzero regular residue proportional to \(dV/(2U)\), which
extends across the simple branch points and the two points
over \(w=0\) in the complementary chart. Its periods cannot
all vanish. A tube above a loop avoiding the strict transform
of \(w=0\) and the finite branch points lies in \(U_t\).
This already contradicts exactness locally around the tacnode
resolution curve.

There are finitely many points at infinity to check and
finitely many excluded values \(t\) in the local arguments.
For any general value, either the local pole with nonzero
period gives the contradiction directly, or every singularity
has a resolution on which the form is holomorphic. In the
latter case the global holomorphic obstruction applies.
