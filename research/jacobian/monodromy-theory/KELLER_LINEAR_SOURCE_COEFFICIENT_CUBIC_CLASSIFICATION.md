# Quartic components linear in a source variable

Let \(F=(P,Q,R)\) be a polynomial Keller map of ordinary degree at most
five in three variables over a characteristic-zero field. Suppose a
quartic target component is
\[
                    R=a(x,y)z+b(x,y),\qquad a\ne0.
\]
Thus \(\deg a\le3\). The classifications below use affine source changes
over the algebraic closure and polynomial changes \(z\mapsto z+c(x,y)\)
when displaying scalar normal forms. The degree bounds for the mates are
always taken in the original affine coordinates, before these polynomial
changes.

**Proposition KellerRepeatedParallelCoefficientResidue.** In any
characteristic-zero three-variable Keller triple with
\(R=a(x)z+b(x,y)\), assume \(\xi\) is a root of \(a\).
Then \(b(\xi,y)\) is constant or affine linear with nonzero slope.
The constant case can occur only at a repeated root, and in that case
\(b_x(\xi,y)\) is a nonzero constant. This proposition places no degree
restriction on the triple or on \(a\).

**Theorem QuinticKellerParallelCoefficientCompletion.** For a quartic
component \(R=a(x)z+b(x,y)\) as above, if \(b(\xi,y)\) is nonconstant
for every root \(\xi\) of \(a\), the map \(F\) is invertible.
The inverse parametrization on every \(R\)-fibre has total degree at most
thirteen, including when \(a\) has repeated roots.

**Theorem QuinticKellerCubicSourceCoefficientResidual.** In a
noninvertible \(F\) as above, the following list exhausts the possibilities
for \(a\) and the indicated restrictions of \(b\), up to affine changes
and scaling.

1. \(a\) is univariate of degree two or three, has a repeated root
   \(\xi\), and \(b(\xi,y)\) is constant. At that root \(b_x(\xi,y)\)
   is a nonzero constant. In a frame with \(\xi=0\),
   \(R-\mu=\lambda x+x^2 H(x,y,z)\), where the quadratic part of
   \(H\) in \(y,z\) has rank at most one.
2. \(a=x^2y+\delta x-1\) and
   \(b\bmod a=\lambda xy+\mu\) with \(\lambda\ne0\).
3. \(a=x(xy-1)\). On \(x=0\) and \(xy=1\), respectively, \(b\)
   restricts to \(\lambda_0 y+\mu_0\) and \(\lambda_1 y+\mu_1\),
   with \(\lambda_0\lambda_1\ne0\).

The last two cases have cubic coefficient. A reduced cubic coefficient
whose zero set is a graph, or a union of parallel lines, forces
invertibility.

## Normal residues along a repeated line

We work over \(\mathbb C\) by embedding a finitely generated coefficient
field. Put \(\kappa=\det JF\ne0\). On the open \(a(x)\ne0\) of the
fibre \(R=t\), use \(x,y\) as coordinates. Its two-form
\[
              \omega=dx\wedge dy/a(x)
\]
has the regular primitive \(\gamma_t=\kappa^{-1}P\,dQ\) on the full
fibre. After substituting \(z=(t-b)/a\), its coefficients on the open
lie in \(\mathbb C[x,a^{-1},y,t]\). Let \(\eta=-y\,dx/a(x)\), so
\(d\eta=\omega\). The point of using the Keller primitive is that
\(t\) enters polynomially before taking residues.

The relative closed form \(\gamma_t-\eta\) can be written as
\[
       d_{x,y} G(x,y,t)+\sum_{a(\zeta)=0} c_\zeta(t)
                                      \frac{dx}{x-\zeta}
\]
with \(G\) in the same ring and \(c_\zeta\in\mathbb C[t]\).
Indeed, integrate its \(dy\) coefficient in \(y\), then use closedness
to see that the remaining \(dx\) coefficient is independent of \(y\).
Partial fraction integration in \(x\) leaves exactly the displayed
simple poles. This is the normal-slice argument in
[Wang–Yan, Lemma 3.14](https://arxiv.org/pdf/2610.36495v1) with all
multiplicities and the polynomial dependence on \(t\) retained.

Write \(a(x)=(x-\xi)^m u(x)\) with \(u(\xi)\ne0\),
and let \(p(y)=b(\xi,y)\). For \(p\) nonconstant and \(p'(y_0)\ne0\),
take the analytic disc \(z=0,\ b(x,y(x))=p(y_0)\), with
\(y(\xi)=y_0\). Pullback of the full primitive there is regular,
and pullback of an exact rational differential has zero residue.
Consequently
\[
       \operatorname{Res}_{x=\xi}\frac{y(x)}{a(x)}dx
          = c_\xi(p(y_0)).                                      \tag{1}
\]
In multiplicity one, any zero of \(p'\) would yield a critical point
of \(R\) by choosing \(z\) to kill \(R_x\). In particular \(p\)
cannot be constant in this case. In multiplicity \(m\ge2\),
a zero \(y_*\) of \(p'\) requires \(b_x(\xi,y_*)\ne0\)
by the same submersion condition, since now \(a'(\xi)=0\).

We show that at such a \(y_*\), the left side of (1) has a pole
as \(y_0\to y_*\). Choose an analytic coordinate \(v(y)\) at \(y_*\)
such that \(p(y)=p(y_*)+v(y)^e\), where \(e\ge2\) and \(v(y_*)=0\).
Set \(\epsilon=v(y_0)\), \(x-\xi=\epsilon^e s\), and
\(v(y(x))=\epsilon V(s,\epsilon)\).
Expanding \(b\) in \(x-\xi\) and dividing the branch equation by
\(\epsilon^e\) gives an analytic equation for \(V\). At \(\epsilon=0\)
it reduces to
\[
               V^e+s\,b_x(\xi,y_*)=1.
\]
The derivative in \(V\) at \(s=\epsilon=0,\ V=1\) is \(e\), so \(V\)
is holomorphic there with
\(V(s,0)=(1-s b_x(\xi,y_*))^{1/e}\).
Writing \(y\) as the analytic inverse of \(v\), the coefficient of
\((x-\xi)^j\) in \(y(x)\), for \(j\ge1\), therefore has leading
order \(\epsilon^{1-ej}\) with nonzero coefficient. The coefficient
is nonzero because the linear coefficient of \(v^{-1}\) is nonzero,
\(b_x(\xi,y_*)\ne0\), and \(\binom{1/e}{j}\ne0\).
The residue in (1) is the coefficient of \((x-\xi)^{m-1}\)
in \(y(x)/u(x)\). Its highest pole comes uniquely from
\(j=m-1\) and \(1/u(\xi)\). The right side of (1) is holomorphic
in \(y_0\), a contradiction. Thus \(p'\) has no root and \(p\)
is affine linear with nonzero slope. If \(p\) is constant and \(m\ge2\),
submersion along \(x=\xi\) requires \(b_x(\xi,y)\) to be nowhere
zero, hence constant and nonzero. This proves the proposition.

## Coordinates on the parallel-support fibres

Let \(d=\deg a\). The constant \(a\ne0\) case is a graph completion.
If \(d>0\) and every \(p(y)=b(\xi,y)\) is affine linear with
nonzero slope, then \(b(x,y)\) is a one-variable coordinate over
\(\mathbb C[x]/(a)\). To see this, work in each factor
\(\mathbb C[x]/((x-\xi)^m)\). Modulo its nilpotent maximal ideal
the substitution \(y\mapsto b(x,y)\) is invertible and affine.
An inverse lifts by successive correction modulo powers of this ideal,
using the invertible slope at each step. Combine the inverses by the
Chinese remainder theorem and choose a representative \(\psi(x,t)\)
of \(x\)-degree less than \(d\).

On \(R=t\), take the coordinates \(x,M=(y-\psi(x,R))/a(x)\).
The latter expression is polynomial in the source variables by the
modulo-\(a\) inverse identity. Its inverse parametrization on the fibre is
\[
    y=\psi(x,t)+a(x)M,\qquad
    z=\frac{t-b(x,\psi(x,t)+a(x)M)}{a(x)} .
\]
The second quotient is polynomial by the other inverse identity.
At fixed \(t\) the first has degree at most \(d+1\) in \(x,M\),
and the second degree at most \(4(d+1)-d\le13\). These bounds
also cover the constant case directly. The formulas give an
automorphism with coordinates \(R,x,M\), hence the restrictions of
\(P,Q\) to all the fibres are plane Keller maps. Each has degree
at most sixty-five. The planar degree theorem and fibrewise argument
from the [bounded coordinate-completion proof](LOW_DEGREE_COORDINATE_KELLER_MAPS_IN_THREE_VARIABLES.md)
give invertibility of \(F\).

## The squarefree cubic classification

Use the [linear-source residue restriction](KELLER_QUARTIC_TRIPLE_BOUNDARY_AND_LINEAR_COMPONENTS.md#a-residue-restriction-for-components-linear-in-a-source-variable)
on each multiplicity-one component of \(a=0\). Its normal residue
of \(dx\wedge dy/a\) is the derivative of a polynomial in the
restriction of \(b\), so has no pole on the affine normalization.
In particular distinct reduced components of a squarefree \(a\)
cannot meet in the affine plane. The residue on either smooth
branch would have a pole there. The same reasoning excludes a cusp
of a reduced cubic in the affine plane. If \(a\) is not squarefree
and has degree at most three, its repeated factor is a line.
A distinct simple line must be disjoint from it by the residue
restriction on that simple line. Thus every non-squarefree
possibility is a product of parallel lines, counted with
multiplicities, and the previous sections give item 1 of the list.

The squarefree cases of degree at most two are covered by the
[reduced quadratic completion](KELLER_QUARTIC_TRIPLE_BOUNDARY_AND_LINEAR_COMPONENTS.md#invertibility-with-a-reduced-quadratic-coefficient).
Three different affine lines of a squarefree cubic must be
parallel. For a line plus an irreducible conic, affine disjointness
forces the conic to be a hyperbola and the line to become \(x=0\)
in the frame \(xy=1\). In fact a nonconstant affine linear equation
restricts to a nonconstant polynomial on the standard parabola,
and to \(A x+B/x+C\) on the hyperbola. The latter is a unit on
the hyperbola only for an axis or an everywhere constant equation.
The line hence gives \(a=x(xy-1)\). The restriction of \(b\) to
\(x=0\) is affine linear with nonzero slope by submersion.
On \(xy=1\), the residue is a nonzero constant times \(dx/x^2\).
Its primitive has a simple pole at \(x=0\) and no other pole.
As it is polynomial in the restriction of \(b\), it must be an
affine linear polynomial in that restriction. Thus \(b\) on
\(xy=1\) is affine linear in \(1/x=y\), with nonzero slope.
This is item 3.

Finally suppose \(a\) is an irreducible cubic. The affine residue
is the adjunction differential of its projective cubic, with no poles at
smooth points of the projective closure. For a smooth projective cubic
it is holomorphic and nonzero, so has no rational primitive.
For a node, the differential on the normalization has a simple pole
at either branch and also has no rational primitive. These assertions
follow by homogenizing \(dx\wedge dy/a\) as the standard projective
residue and using the local reduced normal-crossing form at a node.
The only possibility is therefore a cusp at infinity.
In projective coordinates its normal form is \(Y^2W=X^3\),
with cusp \([0:0:1]\).

If the line at infinity is tangent to the cusp, it is \(Y=0\).
The affine cubic is a graph \(a=y-x^3\). Submersion forces
\(b=\lambda x+\mu+ac(x,y)\), \(\lambda\ne0\), since the restriction
of \(b\) to this graph has nowhere vanishing derivative.
Here \(\deg c\le1\). Taking \(Z=z+c\), the completion \((R,a,Z)\)
has polynomial inverse
\[
  x=(R-\mu-aZ)/\lambda,\qquad y=a+x^3,\qquad z=Z-c(x,y)
\]
of degree at most six, and again gives inversion of \(F\).

Otherwise the line at infinity is \(X+\rho Y=0\). Dehomogenizing
there and absorbing the quadratic and cubic terms in the expansion
of \((1-\rho Y)^3\) by an affine change of \(W\) gives
\(a=x^2y+\delta x-1\). This affine curve is isomorphic to
\(\mathbb C^\times\) by \(x\), and its residue is a nonzero constant
times \(dx/x^2\). The argument for the hyperbola now gives
\(b|_{a=0}=\alpha/x+\beta\) with \(\alpha\ne0\).
Since \(1/x=xy+\delta\) on the curve, this is item 2.
All forms and conclusions transfer to characteristic zero by the
same coefficient-field embedding and polynomial inverse descent.

## A bilinear suspension obstruction

**Theorem KellerBilinearSuspensionConstantCoefficient.** If a component
of a characteristic-zero three-variable Keller map is
\(R=A(x)+C(x)yz\), with \(C\ne0\), then \(A\) is affine linear with
nonzero slope and \(C\) is constant. In particular a quartic
component of the form \(\mu+\lambda x+x^2H(x,y,z)\), where
\(\lambda\ne0\), \(\deg H\le2\), cannot have a rank-two quadratic
part in \(y,z\).

At \(y=z=0\), the gradient of \(R\) equals \(A'(x)\,dx\).
It has no zero, so \(A(x)=\lambda x+\mu\) with \(\lambda\ne0\).
On the open \(C(x)y\ne0\) of \(R=t\), the exact fibre form is
\[
                  \omega=\frac{dx}{C(x)}\wedge\frac{dy}{y},
\]
with a regular primitive \(\gamma_t\) polynomial in \(t\) on this open.
For each root \(\xi\) of \(C\), integrate over the product of a small
circle about \(\xi\) in \(x\) and a circle about zero in \(y\).
Exactness forces \(\operatorname{Res}_{x=\xi} dx/C(x)=0\).
Partial fractions now give \(dx/C(x)=dH_0(x)\) for a rational
\(H_0\) with poles only at the roots of \(C\).

Subtract \(H_0(x)\,dy/y\) from \(\gamma_t\). The difference is
closed. Its coefficient of \(y^{-1}dy\) after removing an
exact differential is independent of \(x\), by closedness,
and is a polynomial \(s(t)\). More explicitly, integrate each Laurent
monomial of its \(dy\) coefficient whose \(y\)-exponent differs from
\(-1\) and subtract that exact differential. Closedness then makes the
remaining \(y^{-1}dy\) coefficient independent of \(x\), since an
ordinary Laurent derivative in \(y\) has zero \(y^{-1}\) term.
The other coefficient is a rational function of \(x\) alone with
polynomial dependence on \(t\).

For generic \(t\) the curve \(x=(t-\mu)/\lambda,\ z=0\) is a line in
the full fibre. On it \(\gamma_t\) is regular at \(y=0\).
Taking the \(dy/y\) residue of the difference expression gives
\[
                         H_0((t-\mu)/\lambda)+s(t)=0.
\]
Thus \(H_0(x)\) is polynomial. Since its derivative is \(1/C(x)\),
\(C\) is constant. This argument uses the regular polynomially
supplied primitive for all \(t\). A choice of primitive separately
on every generic fibre could allow the rational function of \(t\)
that the proof excludes.

For the quartic application, complete the square in the two transverse
variables of the rank-two part of \(H\). The critical point of \(H\)
in \(y,z\) is affine linear in \(x\), so this change is an affine source
change. After a linear change in \(y,z\) the scalar is
\(\mu+\lambda x+x^2B(x)+c x^2 yz\), \(c\ne0\).
The first assertion forces \(c x^2\) to be constant, a contradiction.
