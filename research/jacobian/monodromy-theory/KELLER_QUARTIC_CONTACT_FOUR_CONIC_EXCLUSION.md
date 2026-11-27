# The contact-four union of two conics is excluded

**Theorem KellerQuarticContactFourConicExclusion.** Let \(F=(P,Q,R)\)
be a three-variable polynomial Keller map in characteristic zero.
A quartic component \(R\) cannot have geometrically reduced leading
curve consisting of two smooth conics meeting at one point with
intersection multiplicity four. In particular, for a degree-five
Keller map with two independent degree-five leading components,
a quartic component cannot have reduced leading form equal to a
product of two irreducible conics.

We work over \(\mathbb C\) by embedding a finitely generated
coefficient field. Orient the fibre residue \(\omega_t\) by a
fixed choice of coordinate order. Its regular primitive on \(R=t\)
is a constant multiple of \(P\,dQ\), by
[Wang and Yan, equation (3.28)](https://arxiv.org/pdf/2610.36495v1).
We use both this exactness and the fact that \(R\) has no critical
point.

## The required jet at the common point

Take coordinates \(p=[0:1:0]\) at the common point, with \(z=0\)
the tangent line. After scaling the conic equations, they have
the same \(yz\) coefficient. They can have no \(y^2\) or \(xy\)
term because of the point and the tangent condition. Their
\(x^2\) coefficients coincide and are nonzero: in \(y=1\) the
quadratic coefficient of the branch \(z=z(x)\) is their negative
\(x^2\) coefficient, and if it is zero the conic factors.
The cubic coefficient of that branch similarly fixes the \(xz\)
coefficient, so these two coefficients also coincide. An affine
linear change of \(y\) removes the common \(xz\) term and absorbs
one \(z^2\) term. Scaling coordinates, and scaling \(R\) if
necessary, gives
\[
                  R_4=q(q+z^2),\qquad q=x^2+yz .
\tag{1}
\]
The two conics in (1) have no other common point or singularity.

By the [required-jet theorem](KELLER_QUARTIC_QUADRATIC_SUSPENSIONS_AND_HIGHER_CONTACT.md#necessary-singular-point-in-a-reduced-leading-curve),
\(R_3(p)=0\) and the fibre quadratic jet at \(p\) has rank at most
one. In the chart \(y=1\), with coordinates \(v=x/y,\ u=z/y,\ w\)
the normal coordinate to infinity, the restriction of this jet
to \(w=0\) is \(u^2\). Thus the jet is \((u+\beta w)^2\).
The first condition excludes a \(y^3\) term of \(R\).
The second says the coefficient of \(y^2\) is \((z+\beta)^2\).
After translating the affine \(z\) coordinate, this coefficient
is \(z^2\), with \(R_4\) unchanged.

Write \(R=z^2y^2+A(x,z)y+B(x,z)\).
Here \(A\) has degree three with top \(2x^2z+z^3\),
and \(B\) has degree four with top \(x^4+x^2z^2\).
Set
\[
 D=\frac{A(x,z)-A(x,0)}{2z},\qquad E=B-D^2 ,
 \qquad R=(zy+D)^2+A(x,0)y+E .
\tag{2}
\]
In particular \(D=x^2+z^2/2+\text{an affine polynomial}\)
and \(E_4=-z^4/4\).

For the local splitting of the fibre equation in \(u\),
write its cubic as
\[
                  w^3A(v/w,0)+u U_2(v,u,w).
\]
The notation denotes its polynomial homogenization, including at
\(w=0\). Its critical \(u\) starts at \(-U_2(v,0,w)/2\).
It follows that the residual cubic in \((v,w)\) is \(w^3 A(v/w,0)\)
and its quartic at parameter \(t\) is \(w^4 E(v/w,0)-tw^4\),
since \(U_2(v,0,w)/2=w^2D(v/w,0)\).
The residual plane function at \(w=0\) has order eight:
the two branches intersect to order four, and after splitting
the difference of their \(u\) coordinates is twice a square
root of the negative residual plane function, up to a unit.
In particular \(f_5(1,0)=0\) in the notation of the higher-contact
local tests.

If \([x^2]A(x,0)\ne0\), the residual cubic has a simple zero
at the plane direction \(w=0\). The first higher-contact test gives
holomorphic resolution, which is impossible with exact residue on
the sole singular point. We have \(A(x,0)=\lambda x+\mu\).
If this polynomial is nonzero, the residual cubic has a double
or triple zero at \(w=0\). In either case, using \(f_5(1,0)=0\),
the higher-contact tests require \([x^3]E(x,0)=0\).
If \(A(x,0)=0\), the residual cubic is zero. The last
higher-contact test then requires the residual quartic to have
order at least three in the plane direction, so
\(\deg_x E(x,0)\le1\).

## The affine critical-point tests

Suppose \(\lambda x+\mu\ne0\). The jet constraint gives
\[
 D=x^2+d_1x+d_0(z),\qquad
 E=e_2(z)x^2+e_1(z)x+e_0(z),
\tag{3}
\]
where \(d_1\) is constant, \(\deg d_0=2,\ \deg e_2\le1,\ \deg e_1\le2\),
and \(e_0(z)=-z^4/4+O(z^3)\).
The \(x^3\) term of \(E\) is absent globally because its coefficient
would be constant, and its coefficient at \(z=0\) vanishes.
On \(z\ne0\), solving \(R_y=0\) gives
\[
 y=-D/z-(\lambda x+\mu)/(2z^2).
\]
There the critical equations reduce to the two derivatives of
\[
               K(x,z)=E-\frac{(\lambda x+\mu)D}{z}
                         -\frac{(\lambda x+\mu)^2}{4z^2}.
\tag{4}
\]

### Nonconstant linear term

If \(\lambda\ne0\), translate \(x\) so that \(\mu=0\), preserving (3).
The polynomial \(K_x\) in \(x\) is quadratic with leading
coefficient \(-3\lambda/z\). The product \(N(z)\) of \(K_z\)
at its two roots, counted with multiplicity, is rational in \(z\).
It has no possible finite pole outside \(z=0\), since the quadratic
leading coefficient never vanishes there. A finite nonzero root of
\(N\) gives a critical point of \(R\).

At zero, multiplying \(K_x=0\) by \(z^2\) gives
\[
 -3\lambda z x^2+(-\lambda^2/2+O(z))x
                     -\lambda d_0(0)z+O(z^2)=0.
\]
The finite root has \(x=-2d_0(0)z/\lambda+O(z^2)\).
Thus its corresponding \(y\) in (4) is holomorphic at \(z=0\).
The value of \(K_z\) there is \(R_z\) along that branch,
and its limit cannot be zero because \(R_x=R_y=0\) also extend.
The second root is \(x=k/z+O(1)\), with \(k=-\lambda/6\).
At that root \(K_z\) has a pole of order five with coefficient
\(\lambda k^3+\lambda^2 k^2/2\ne0\).
These are exactly the leading terms from differentiating
\(-\lambda x^3/z\) and \(-\lambda^2x^2/(4z^2)\) at fixed \(x\).
Thus \(N\) has a pole at zero.

At infinity let \(\alpha=[z]e_2(z)\).
If \(\alpha\ne0\), the two quadratic roots have respective forms
\(x=k_\infty z^2+O(z)\), with \(k_\infty=2\alpha/(3\lambda)\),
and \(x=O(z)\). The corresponding values of \(K_z\)
have leading terms
\((\lambda k_\infty^3+\alpha k_\infty^2)z^4\)
and \(-z^3\), respectively, both nonzero.
If \(\alpha=0\), the roots are \(O(z^{3/2})\), as follows by
multiplying \(K_x=0\) by \(z\) and using the quadratic formula.
Both \(K_z\) values then have leading term \(-z^3\), since every
other term is \(O(z^{5/2})\) or smaller.
In both cases \(N\) has a pole at infinity. It has a zero
away from zero and infinity, a contradiction.

### Nonzero constant term

If \(\lambda=0,\ \mu\ne0\), (4) has the form
\[
 K=U(z)x^2+V(z)x+W(z),\qquad
 U=e_2-\mu/z,\quad V=e_1-\mu d_1/z,\quad
 W=e_0-\mu d_0/z-\mu^2/(4z^2).
\]
The function \(J=W-V^2/(4U)\) is rational, with a pole of order
two at zero. A finite nonzero zero of \(J'\) where \(U\ne0\)
gives a critical point at \(x=-V/(2U)\).

Consider a zero \(z_0\ne0\) of \(U\). If \(V(z_0)\ne0\),
then \(J\) has a pole there. If \(V(z_0)=0\), the derivative
\(K_x\) vanishes for all \(x\) there. For a simple zero
of \(U\), \(K_z=U'x^2+V'x+W'\) has a root in \(x\).
For a double zero, it has a root if \(V'(z_0)\ne0\) or
\(W'(z_0)=0\). In the remaining case \(V\) vanishes to
order at least two and \(W'(z_0)\ne0\), so \(J'(z_0)=W'(z_0)\ne0\).
The zero of \(U\) has multiplicity at most two because
\(zU=ze_2-\mu\) is quadratic or linear. Therefore in the
absence of a critical point of \(R\), \(J'\) has no finite
zero: any such zero would be nonzero because zero is a pole,
and the cases just considered either give a critical point
or do not give a zero of \(J'\).

If \(J\) has a pole at infinity, \(J'\) does not vanish there.
Since \(J'\) has a pole at zero, it then has a finite zero.
If \(J\) is finite at infinity, the leading \(-z^4/4\) in \(W\)
forces \(e_2\) to be a nonzero constant and \(V\) of degree two
at infinity. Indeed a linear \(e_2\) gives \(V^2/U=O(z^3)\);
if \(e_2=0\), this quotient either has degree five at infinity
or is \(O(z^3)\). With \(e_2\) constant nonzero, \(U\) has a
simple nonzero root. Either \(V\) vanishes there and gives a
critical point, or \(J\) has at least two finite poles.
In the latter case, if \(J\) has degree \(n\) as a rational
function and \(r\ge2\) finite poles, the finite poles of
\(J'\) have total order \(n+r\). At infinity \(J'\) has
a zero of order at most \(n+1\), because
\(J(z)-J(\infty)\) has order at most \(n\) there.
Thus \(J'\) has a finite zero as well. This contradiction
excludes \(\mu\ne0\).

## The fibre residue after the square substitution

We are left with
\[
             R=(zy+D(x,z))^2+E(x,z),\qquad \deg E(x,0)\le1 .
\]
Put \(T=zy+D\) on \(z\ne0\), and write \(E_0(x)=E(x,0)\)
and \(D_0(x)=D(x,0)\). The induced model
\(T^2+E(x,z)=t\) has the residue two-form, up to an overall sign,
\[
                      \omega_t=\frac{dx\wedge dz}{2zT}
\tag{5}
\]
near \(z=0,\ T\ne0\). The regular primitive of the actual
fibre pulls over to a meromorphic one-form \(\gamma_t\) on
this model by \(y=(T-D)/z\). Its derivative is (5).

Locally use \(x,z\) as coordinates, choosing the branch of
\(T\) through the point in question. By subtracting exact
meromorphic differentials remove the terms of the \(dz\)
coefficient of \(\gamma_t\) of order below \(-1\) in \(z\).
Indeed \(a_j(x)z^j dz\) with \(j<-1\) is removed by
\(d(a_j(x) z^{j+1}/(j+1))\). The remaining \(dx\) coefficient
cannot have a negative power of \(z\): differentiation of any
such term would give a pole of order at least two in (5).
Let \(h_t\) be the coefficient of \(dz/z\) on \(z=0\).
It is single valued on \(T^2+E_0(x)=t,\ T\ne0\);
the local coordinates \(x,z\) and the choice of \(T\) through
each point determine the coefficient consistently.
By (5) it satisfies
\[
                              dh_t=\frac{dx}{2T}
\]
with a consistent orientation, and its value at a point equals
the residue of \(\gamma_t\) along any transverse analytic disc
through the point. The exact terms removed along such discs
have zero residue.

Take a general \(t\), so the quartic
\(D_0(x)^2+E_0(x)-t\) has four distinct roots, none with \(D_0=0\).
At each root \(x=x_0\), hold \(x=x_0\) fixed and solve
\(T(z)^2+E(x_0,z)=t\) holomorphically with \(T(0)=D_0(x_0)\).
Since \(T(0)\ne0\) this gives a transverse disc. Its lift
\(y=(T(z)-D(x_0,z))/z\) is holomorphic at zero, so the primitive
is regular on the actual disc. Thus \(h_t=0\) at all these points.

If \(E_0\) is constant, the locus \(z=0\) in the model has
two lines \(T=\pm\sqrt{t-E_0}\). On each, \(h_t\) has a
nonzero constant derivative in \(x\), but would vanish at
the two distinct roots of \(D_0(x)=T\) on that line.
If \(E_0=a x+b\) with \(a\ne0\), the locus is a parabola
with parameter \(T\), and \(dx/(2T)=-dT/a\).
Thus \(h_t\) is affine linear with nonzero slope in \(T\)
on the connected open \(T\ne0\). The four points above have
distinct \(T\)-values, since \(T\) determines \(x\).
Both cases are impossible. This finishes the contact-four
exclusion. The rank-two consequence combines it with the
[necessary contact-four pencil](QUINTIC_KELLER_RANK_TWO_QUARTIC_PENCIL.md#two-reduced-conics).

## Reduced leading forms in top rank two

**Corollary QuinticRankTwoReducedQuarticInversion.** A degree-five
three-variable Keller map of top rank two is invertible if its
quartic target combination cancelling degree five has geometrically
reduced leading form.

In a noninvertible map, the quartic cannot be irreducible by the
[rank-two pencil restriction](QUINTIC_KELLER_RANK_TWO_QUARTIC_PENCIL.md).
It cannot have an irreducible cubic factor by the
[cubic-pencil closure](QUINTIC_KELLER_RANK_TWO_CUBIC_PENCIL_CLOSURE.md).
For a reduced leading quartic with all factors of degree at most
two, the [low-factor boundary theorem](QUINTIC_KELLER_REDUCED_BOUNDARY_PENCIL_RESIDUAL.md)
forces a union of two irreducible conics. The rank-two pencil forces
contact four, just excluded. This exhausts the reduced possibilities.
