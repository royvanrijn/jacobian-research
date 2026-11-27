# Further restrictions on cubic source-direction mates of the squarefree quartic

Let \(v=xy-1\), \(a=xv\) and \(R=az+b(x,y)\). The source degree is degree in \(z\).

**Theorem KellerSquarefreeQuarticUnitAtZeroLeadingExclusion.**
Over any characteristic-zero field, with arbitrary \(b\in k[x,y]\), \(R\)
has no polynomial Keller completion with a mate of positive source degree
whose leading coefficient is \(\alpha a+\beta\), where
\(\alpha\beta\ne0\). There is no ordinary total degree bound.

**Theorem SexticKellerSquarefreeQuarticReducedCubicCases.**
Take the [squarefree normal form](QUARTIC_AFFINE_SOURCE_COMPONENTS_THROUGH_DEGREE_SEVEN.md)
with
\[
 b=\lambda_0 y+\mu_0+xy\bigl((\lambda_1-\lambda_0)y+\mu_1-\mu_0\bigr),
 \qquad \lambda_0\lambda_1\ne0 .
\]
If \(P,Q,R\) are the components of a characteristic-zero Keller map of
ordinary degree at most six and one mate has source degree three, a
constant linear change fixing \(R\) gives source degrees \((3,2)\).
The respective leading coefficients, up to rescaling the mates, are
either \(a,a\), or \(x^3,x^2(\gamma v+\delta)\) for scalars
\(\gamma,\delta\) not both zero.

## A leading coefficient nonvanishing on both components

Work over \(\mathbb C\) by embedding a finitely generated coefficient
field. Let \(P=p(a)z^m+\cdots\), \(p(a)=\alpha a+\beta\), and let
\(\kappa\ne0\) be the constant determinant, with sign chosen to agree with
the order of the brackets below. Write \(S\) for the ring obtained from
\(\mathbb C[x,y]\) by inverting \(p(a)\) and \(mp(a)-ap'(a)\) and adjoining
\(s\) with \(s^m=p(a)\). This is finite étale over the localized ring.
It is a domain: \(p(a)\) has a simple root in the rational base
\(\mathbb C(a)\), and \(\mathbb C(x,y)=\mathbb C(a)(x)\). Its fraction
field is \(\mathbb C(s)(x)\).

Use the formal Laurent coordinate \(w=P^{1/m}=sz+O(1)\) at \(z=\infty\).
The binomial expansion and the inverse Laurent substitution have
coefficients in \(S\) since \(s\) is a unit. Put dots over the polynomials
written in \(x,y,w\). At fixed \(w\), using the two-variable Jacobian
bracket \(\{f,g\}=f_xg_y-f_yg_x\), the chain rule gives
\[
 \{\dot R,\dot Q\}=\frac{\kappa}{\dot P_z}
     =\frac{\kappa}{ms}w^{1-m}+O(w^{-m}),\qquad
 \dot R/w=c+O(w^{-1}),\quad c=a/s.
\]
The partial derivatives on the finite étale extension are the unique
extensions of the partial derivatives in \(x,y\).

We can rectify \(\dot R/w\) to \(c\) by a formal automorphism over
\(\mathbb C[[w^{-1}]]\) of \(S[[w^{-1}]]\) with constant term the identity.
Indeed,
\[
 D=(-1-2xy)\partial_x+4y^2\partial_y
 \quad\text{has}\quad D(a)=1,
 \qquad D(c)=c'(a)=\frac{mp(a)-ap'(a)}{mp(a)s}.
\]
Hence \(D_0=D/D(c)\) is regular on \(S\) and sends \(c\) to \(1\).
At order \(w^{-j}\), if \(f_j\) is the coefficient of the current error,
apply \(\exp(-w^{-j} f_j D_0)\). It cancels that coefficient and leaves
all earlier coefficients unchanged. The composition converges in the
\(w^{-1}\)-adic topology, and so does its inverse. It extends to the ring
of Laurent series with finitely many positive powers of \(w\).

Under this substitution the Jacobian factor for the bracket is a formal
unit with constant term \(1\), by the chain rule in the étale coordinates
\(x,y\). Thus, writing \(\widetilde Q\) for the transformed series, we have
\[
 \{c,\widetilde Q\}=\frac{\kappa}{ms}w^{-m}+O(w^{-m-1}).
\]
The coefficient \(h\in S\) of \(w^{-m}\) satisfies
\(\{a,h\}=\kappa/(ms c'(a))\). At fixed \(a,s\),
\(\{a,\cdot\}=-x^2\partial_x\), since \(y=(a+x)/x^2\). Consequently
\[
 h=\frac{\kappa}{ms c'(a)}\,\frac{1}{x}+f(s),
 \qquad f\in\mathbb C(s).
\]
Choose any \(s_0\) with \(s_0^m=\beta\). Both divisors \(x=0\) and
\(v=0\) have lifts at \(s=s_0\), and \(a\) vanishes simply on both.
Indeed \(p(0)\), \(mp(0)\) and \(s_0\) are nonzero, so the lifts
are étale there and \(a\) and \(s-s_0\) differ by a local unit.
The first summand is regular at the generic point of the lift of
\(v=0\), so regularity of \(h\) there forces \(f\) to be regular
at \(s_0\). The first summand has a simple pole at the generic point
of the lift of \(x=0\) because its scalar factor is nonzero at \(s_0\).
This contradicts \(h\in S\).

## The mate with linear source coefficient

We use the [cubic leading classification](SEXTIC_KELLER_SQUAREFREE_MATE_DEGREES.md).
If the cubic mate has leading coefficient \(\alpha a+\beta\), that
classification gives a degree-two mate with leading coefficient
\(\gamma a+\delta\), with \(\alpha\gamma\ne0\). The preceding
obstruction applied to each mate gives \(\beta=\delta=0\).
If the cubic leading coefficient is proportional to \(x^3\), the
classification gives the stated quadratic case or a source-linear
mate proportional, after scaling, to \(xz+d(x,y)\).

To exclude the latter case, let \(Q=xz+d(x,y)\). Its differential
cannot vanish. On \(x=0\), the \(z\)-derivative is zero and the
\(x\)-derivative can be made zero by choosing \(z=-d_x(0,y)\).
Thus \(d_y(0,y)\) never vanishes, and \(d(0,y)=c_0y+c_1\) with
\(c_0\ne0\). Write \(D_1=(d-c_0y-c_1)/x\), which has degree
at most five. The polynomial coordinate system
\[
 (x,Z,Q)=\bigl(x,z+D_1(x,y),c_0y+c_1+x(z+D_1(x,y))\bigr)
\]
has inverse \(y=(Q-c_1-xZ)/c_0,\ z=Z-D_1(x,y)\) of degree at
most ten. The [bounded coordinate completion theorem](LOW_DEGREE_COORDINATE_KELLER_MAPS_IN_THREE_VARIABLES.md)
applies since \(6\cdot 10\le100\) and would make \(R\) a
polynomial coordinate. But on the fibre \(R=\mu_1\),
\(1-xy\) is a nonconstant unit. It has no zeros there since
\(xy=1\) would imply \(\lambda_1 y=0\), and it is nonconstant
as seen on the graph over \(xv\ne0\). This fibre cannot be an
affine plane. The same argument over \(\mathbb C\) rules out a
completion over any characteristic-zero coefficient field.

These restrictions leave both quadratic cases stated in the theorem
and the cases in which both mate source degrees are at most two.
