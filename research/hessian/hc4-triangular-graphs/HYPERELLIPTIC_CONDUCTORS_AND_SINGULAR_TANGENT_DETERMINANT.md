# Excluding the degree-two target coordinate in directional conductors

Let \(F=(P,Q)\) be a noninvertible complex polynomial plane Keller map
in a source frame \(x,y\) with \(P\) monic of positive \(x\)-degree.
Let \(K=\mathbb C(y)\), \(A=K[P,Q]\subset K[x]\), and let \(\kappa\)
generate the exact directional conductor furnished by the
[conductor theorem](GENERIC_COORDINATE_ORDER_CONDUCTOR_AND_CLOSED_COLLISIONS.md).
In particular \(A\) has normalization \(K[x]\). Each finite line
\(y=y_0\) has the exact nonconstant specialized conductor
\(\kappa(x,y_0)\), and its actual coordinate order has colength
\(\deg_x\kappa(x,y_0)/2\).

**Theorem HyperellipticTargetCoordinateExclusion.**
There is no \(K\)-polynomial target automorphism that gives a
coordinate of \(x\)-degree two on the generic line. In particular
no reduced characteristic sequence \((2d+1,2)\) occurs for any
integer \(d\ge1\). The statement includes all drops in the conductor
degree at finite parameters. The [reflection theorem](DEGREE_TWO_ELEMENTS_AND_KELLER_REFLECTION.md)
excludes every degree-two element of the generic order, including
elements supplied by approximate roots. The argument below establishes
the target-coordinate case by controlling the action scale.

## The canonical order over the parameter line

Assume the indicated target-coordinate change exists. Affine
normalization of its quadratic coordinate gives \(T=u^2\) with
\(u=x-h(y)\), \(h\in K\). Subtract the even part of the other
coordinate as a polynomial in \(T\). The resulting generators
are \(T\) and \(u\) times a polynomial in \(T\). In an order of
the form
\[
                  K[T]\oplus u g(T)K[T]
\]
the conductor is \(g(u^2)K[u]\): for an element of \(K[u]\) to
multiply all of \(K[u]\) into the order, both it and its product
by \(u\) must belong to the order. We can therefore choose
\[
              g(u^2)=\kappa(u+h,y),\qquad V=u g(T),
\]
using a \(K^\times\)-scaling of the second target coordinate.
The original pair is obtained from \(T,V\) by a polynomial
target automorphism \(\Phi_y\) over \(K\). Its nonzero target
Jacobian determinant is a scalar \(\Lambda(y)\in K^\times\).
Write \(g(T)=\sum_{j=0}^d a_j T^j\), \(d\ge1\) and \(a_d\ne0\).

For any finite \(y_0\), a finite root of the nonconstant
\(\kappa(x,y_0)\) is the limit of a generic root \(x_+(y)\)
after a local parameter extension. Its reflected mate
\(x_-(y)=2h(y)-x_+(y)\) has the same \(P\)-value since both
have the same \(T,V\). Monicity of \(P\) bounds this mate:
it is a root of \(P(X,y)-P(x_+(y),y)\) with bounded coefficients.
Consequently \(h\) has no finite pole. It is polynomial, and
so are the coefficients \(a_j\).

The coefficient \(a_0\) is nonzero. If it were identically
zero, both generators \(T,V\) would have zero longitudinal
derivative at \(u=0\), contradicting the Keller condition at
generic parameters. If \(a_0(y_0)=0\) at a finite parameter,
a generic root of \(g(T)\) would tend to zero and the distinct
reflected mates would both tend to \(u=0\). The two original
divided differences vanish along those mates, so they would
specialize to \(P_u=Q_u=0\), again impossible. Hence \(a_0\)
is a nonzero constant.

The polynomial \(g\) is generically squarefree. If a multiple
root persisted, choose its nonzero square root \(r\) over a
finite parameter extension. At \(u=\pm r\) the \(u\)-derivative
of \(V\) vanishes, and the \(u\)-derivatives of \(T\) are
opposite. The two longitudinal derivatives of the original
pair are likewise opposite. Differentiating the pair
equalities with respect to \(y\) then gives equal parameter
derivatives at the two mates. The two Jacobian determinants
would be opposite, although each equals the same nonzero
constant \(J=\operatorname{Jac}_{x,y}(P,Q)\). The moving
translation \(u=x-h(y)\) itself has determinant one.

## The determinant at a finite singular point

Write both original coordinates as polynomials in \(u,y\).
Their even parts belong to \(\mathbb C[y][T]\). Their odd
parts are \(u\,g(T)\) times polynomials in \(K[T]\). The latter
quotients lie in \(\mathbb C[y][T]\): division by \(g(T)\)
in the power series ring \(\mathbb C[y][[T]]\) is valid
because \(a_0\) is a nonzero constant. Its coefficients agree
with those of the known polynomial quotient over \(K\).
This gives lifts
\[
   P=A_P(T,y)+V B_P(T,y),\qquad
   Q=A_Q(T,y)+V B_Q(T,y)
\]
regular in the whole parameter line. Denote the pair of
polynomial lifts by \(\Psi_y(T,V)\). Neither its ambient
invertibility nor its ambient determinant is assumed.

The relation of the canonical plane embedding is
\[
                 H_y(T,V)=V^2-Tg(T)^2=0.
\]
Thus the components of \(\Phi_y-\Psi_y\) are multiples of
\(H_y\) over \(K[T,V]\), by division by its monic quadratic
in \(V\). At a generic point \((\alpha,0)\) with
\(g(\alpha)=0\), both \(H_y\) and its two target derivatives
vanish. Therefore
\[
              \Lambda(y)=\det D\Psi_y(\alpha,0).
\tag{1}
\]

Now fix any finite \(y_0\). The polynomial \(g_0=g|_{y_0}\)
is nonconstant by the exact conductor theorem, with nonzero
constant coefficient. Choose a root \(\alpha_0\ne0\) and a
generic root \(\alpha\) specializing to it. The specialized
actual pair lies in
\(\mathbb C[u^2,u g_0(u^2)]\), by the regular lifts. Its order
and this canonical order both have colength \(\deg g_0\), so
they coincide. Consequently \(\Psi_{y_0}\) restricts to an
isomorphism from the curve \(H_{y_0}=0\) onto the specialized
image curve of the actual pair. The point \((\alpha_0,0)\)
has two-dimensional Zariski tangent space, since both
derivatives of \(H_{y_0}\) vanish. Its image has
two-dimensional tangent space as well. The restricted
isomorphism induces an isomorphism of these tangent spaces,
so the ordinary two-by-two derivative of the ambient lift
\(\Psi_{y_0}\) is invertible there. Equation (1) now shows that
\(\Lambda(y)\) has neither a pole nor a zero at \(y_0\).
This argument allows repeated roots of \(g_0\), a drop of
degree, and ramified parameter extensions. It follows that
\(\Lambda\) is a nonzero constant.

## Constant-action-scale contradiction

For a generic reflected pair \(u=-r,r\), where \(g(r^2)=0\),
the action for the original pair is \(\int_{-r}^{r}P\,dQ\).
Its \(y\)-derivative is \(-2Jr\). Differentiation gives this
by integration by parts and cancellation of the endpoint
terms, since the endpoints have equal values of both
coordinates. The target automorphism multiplies the
canonical action by \(\Lambda\): at a fixed parameter the
difference of the two target primitives is closed on the
affine target plane and hence exact. Its contribution
between two equal target endpoints is zero.

The canonical action is
\(-\int_{-r}^r 2u^2g(u^2)\,du\), by integration by parts
and \(V(\pm r)=0\). On differentiating, the boundary terms
again vanish. Dividing the resulting identity by \(-4r\)
shows that the polynomial
\[
  \sum_{j=0}^{d}\frac{\Lambda a_j'}{2j+3}\,T^{j+1}-\frac J2
\]
vanishes at all roots of \(g(T)\), where prime denotes \(d/dy\).
Those roots are simple, so the polynomial equals
\((\alpha_1 T+\beta_1)g(T)\) over \(K\). Its constant term
gives \(\beta_1=-J/(2a_0)\ne0\), while the \(T\)-coefficient
gives \(\alpha_1 a_0+\beta_1 a_1=0\). Thus \(\alpha_1\) is
polynomial. From the coefficient of \(T^{d+1}\) one has
\(\alpha_1=\Lambda a_d'/((2d+3)a_d)\). A nonzero polynomial
with polynomial logarithmic derivative is constant, so
\(a_d\) is constant and \(\alpha_1=0\). Successive
coefficients now give \(a_1=\cdots=a_d=0\), a contradiction.

This extends the [individual hyperelliptic exclusions](CONDUCTOR_TWELVE_HYPERELLIPTIC_AND_CUBIC_EXCLUSIONS.md#the-degree-two-branch)
using the tangent determinant instead of a classification
of specialized generator pairs. In conductor degree fourteen
it excludes the reduced sequence \((15,2)\).
The [degree-two approximate-root proof](CONDUCTOR_FOURTEEN_HYPERELLIPTIC_EMBEDDINGS.md)
covers the two other embeddings of the same degree semigroup at conductor fourteen.
