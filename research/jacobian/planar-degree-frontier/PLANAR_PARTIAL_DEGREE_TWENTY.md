# Excluding the planar partial-degree pair (12,18)

**Theorem PlanarCommonPartialDegreeTwenty.**
Let \(P,Q\) be polynomials over a characteristic-zero field with
nonzero constant Jacobian. If both partial degrees in one common
source variable are at most twenty, then \((P,Q)\) is invertible.

The [standard-corner reduction](SMALL_PARTIAL_GCD_STANDARD_CORNERS.md#planar-pairs)
leaves only the nondividing pair of partial degrees \((12,18)\).
Over the algebraic closure, after the shear preserving the two
leading rows, translate the common high root to zero. In
the original axes \((x,y)\), the two leading forms in weight
\((-1,\rho)\), for \(\rho=3\) or \(6\), are proportional to
\(R^2,R^3\), where
\[
 R=x^{-\rho}r(z),\qquad z=x^\rho y,\qquad r(z)=zL(z),
 \qquad \deg L=5,\quad L(0)\ne0 .
\]
Indeed the only two primitive corners in the swapped axes are
\((6,15)\) and \((6,30)\), both with other endpoint \((1,0)\)
and direction \((\rho,-1)\). Proportionality to powers of a common
form also follows from the zero highest Jacobian and the
positive weights in ratio \(2:3\). Their common monic leading
factor \(h\) in \(x\) has multiplicity \(5\rho\) at zero.
The proof below excludes both cases. It adapts the polynomial
resultant and the successor argument of the
[Weyl \((6,30)\) exclusion](../../dixmier/dixmier-formal-lifting/POLYNOMIAL_WEYL_ORDER_TWELVE_HIGH_ROOT_FRACTIONAL_EXCLUSION.md#the-order-eleven-obstruction-and-its-generalization-at-the-successors)
to a commutative pair and treats the \((6,15)\) corner with
the same argument.

## A second polynomial Jacobian identity

Scale \(P,Q\) to have leading coefficients \(h^2,h^3\) in \(y\)
and write \(J(P,Q)=\kappa\). Over the finite function-field
extension obtained by adjoining \(s=h^{1/6}\), the formal variable
\(T=P^{1/12}\) can be chosen with leading \(y\)-term \(s y\).
Formal inversion at infinity gives the same Laurent series
field in \(y^{-1}\) and \(T^{-1}\). Write \(Q=\sum c_j T^j\).
Then \(c_{18}=1\), and
\[
 \left.\frac{\partial Q}{\partial x}\right|_T
 =-\frac{\kappa}{P_y}
\]
has \(T\)-degree \(-11\). The constant field of the finite
function-field extension over the algebraically closed ground
field is the ground field. Therefore \(c_j\) is constant for
\(j>-11\). Subtract the constant and \(c_{12}P\) from \(Q\).
This preserves both initial forms because the subtracted terms
have smaller \(y\)-degrees and smaller weight \((-1,\rho)\).

Put \(F=T^6(P+c_6)\), so \(F^2=P(P+c_6)^2\) is polynomial.
If there were a nonnegative exponent of \(Q-F\) with
nonzero coefficient, take the highest such exponent \(j\).
It is not divisible by six. The polynomial \(Q^2-F^2\) has
\(y\)-degree \(18+j\) and \(y\)-leading coefficient
\(2c_j s^{18+j}\), with vanishing order
\(5\rho(18+j)/6\) at \(x=0\). Its weight \((-1,\rho)\)
is therefore at least \(\rho(18+j)/6\).
The identity \(J(P,Q^2-F^2)=2\kappa Q\) forces its highest
form to commute with \(R\) whenever that weight is greater
than \(2\rho-1\). In fact, for a form \(x^{-w}u(z)\), the
equation of commutation with \(R\) is
\(\rho r u'= w r' u\).
Since \(r\) has a simple root at zero, this requires
\(w/\rho\) to be a nonnegative integer and
\(\deg u=6w/\rho\). The lower weight bound and the
\(y\)-degree \(18+j\) force equality, contradicting \(6\nmid j\).

Now set \(F_1=F+c_{-6}T^{-6}\) and take the polynomial
\[
 W=Q^2-P(P+c_6)^2-2c_{-6}(P+c_6).
\]
We still have \(J(P,W)=2\kappa Q\).
The degree \(j\) of \(Q-F_1\) is at least \(-11\), by
the displayed derivative identity, and at most \(-1\).
Also \(j\ne-6\). The correction \(F_1^2-F^2-2c_{-6}(P+c_6)\)
is \(c_{-6}^2 P^{-1}\), of degree at most \(-12\), so
\(W\) has \(y\)-degree \(18+j\). If \(18+j>6(2\rho-1)/\rho\),
then \(j>-11\). The same leading-coefficient and
commutation argument again forces \(6\mid j\), impossible.
Thus
\[
 \deg_y W\le
 \begin{cases}10&(\rho=3),\\ 11&(\rho=6).\end{cases}
\]
Its weight is exactly \(2\rho-1\). A smaller weight
cannot give the bracket of weight \(3\rho\) in
\(J(P,W)=2\kappa Q\). A larger weight would force a
commuting initial form of weight at least \(2\rho\)
and \(y\)-degree at least twelve.

Write its initial form as \(x^{-(2\rho-1)}w(z)\). The
Jacobian identity at the top weight gives, up to orientation,
\[
 \rho r w'-(2\rho-1)r'w=\kappa_0 r^2,\qquad \kappa_0\ne0.
\]
For a root of multiplicity \(l\) in \(r\), an order \(v\le l\)
in \(w\) would have nonzero leading multiplier
\(\rho v-(2\rho-1)l\) and give order below \(2l\).
An order \(v>l+1\) would give order above \(2l\).
Thus \(w\) has order \(l+1\). In particular \(w/(z r)\)
is polynomial of degree at most four and vanishes at every
root of \(L\). There is a root \(\lambda_0\) of \(L\) of
multiplicity \(m\in\{2,3,4,5\}\).

## Increasing slopes after the root shift

Pass to \((X,Y)=(y,x)\), reversing the sign of both
Jacobian identities. In the fractional polynomial ring
\(k[X^{1/\rho},X^{-1/\rho}][Y]\) consider weights \((1,-s)\).
At \(s_0=1/\rho\) let \(\eta=X^{1/\rho}Y\).
The common factor on the first face is \(X L(\eta^\rho)\)
and the face of \(W\) is proportional to
\[
 X^{2-1/\rho}\eta^{-(2\rho-1)}w(\eta^\rho).
\]
Translate \(Y\) by \(\lambda X^{-1/\rho}\) for a root
\(\lambda\) of \(\lambda^\rho=\lambda_0\).
This shift preserves both Jacobian identities. On these
faces the monomials of smallest \(Y\)-degree now have
the following exponent expressions:
\[
 X^2\eta^{2m},\qquad X^3\eta^{3m},\qquad
 X^{2-1/\rho}\eta^{m+1}.
\]
All three coefficients are nonzero. For slightly larger
\(s\) the first two monomials are uniquely highest.
Their weights, if they still remain highest, are
\(a=2(1-m\Delta)\) and \(b=3a/2\), where \(\Delta=s-1/\rho\).
The third monomial gives a lower bound
\(\omega_0=2-1/\rho-(m+1)\Delta\) for the weight of \(W\).
For \(s>s_0\), no term of \(Y\)-degree greater than \(m+1\)
can match the weight of that third monomial, because it was
at most that weight at \(s_0\) and loses weight faster.

There must be a first slope \(s>s_0\) at which the first
or second monomial shares its highest weight with another
term of its polynomial. It occurs strictly before
\(1/\rho+1/m<1\). Otherwise at that latter slope both
weights would be zero and the Jacobian would have negative
weight, since a bracket costs \(1-s>0\). The supports are
finite, which ensures a first tie. At this slope both
weights are positive. Either \(a+b=1-s\), or the top
forms commute, since the Jacobian is scalar.

If \(a+b=1-s\), then
\[
                \rho s=\frac{5m+4\rho}{5m-1}.
\]
For \(\rho=3\) the cancelled gcd divides thirteen, and
for \(\rho=6\) it divides twenty-five. For \(m=2,3,4,5\)
in either case the reduced denominator is \(5m-1>3m\).
Equal weight monomials on the \(1/\rho\)-lattice have
\(Y\)-degrees differing by a multiple of this denominator.
At the first tie every other highest term has smaller
\(Y\)-degree than the previously highest monomial. An
equal-degree term would have the same exponent of \(X\),
and a larger degree would have exceeded the original
highest weight immediately before the tie. Since the
anchored \(Y\)-degrees are \(2m,3m\), a tie is impossible
with that denominator.

In the commuting case set \(\xi=X^sY\). Coprimality
of the weight ratio \(2:3\) and vanishing of the bracket
give the two highest forms as scalar multiples of
\(X^a f(\xi)^2,X^b f(\xi)^3\) with \(\deg f=m\).
Write the top of \(W\) as \(X^\omega U(\xi)\),
with \(\deg U\le m+1\). The identity \(J(P,W)=2\kappa Q\)
and \(\omega_0+a-(1-s)=b\) force this top to commute
with \(X^a f^2\) if \(\omega>\omega_0\). Comparing the
highest \(\xi\)-degree in that equation gives
\[
 \omega=\frac{a\,\deg U}{2m}\le\frac{a(m+1)}{2m}
 =\omega_0-\left(1-\frac{1}{\rho}-\frac{1}{m}\right)
 <\omega_0.
\]
Hence \(\omega=\omega_0\) and the nonzero bracket gives
\[
 U'-\frac{2\omega_0}{a}\frac{f'}{f}U=c f,\qquad c\ne0.
\]
Since \(2\omega_0/a>(m+1)/m>1\), at every root of
\(f\) of multiplicity \(l\), an order \(v\le l\) in
\(U\) gives a nonzero term of degree \(v-1<l\) on the
left. Thus \(U\) vanishes at each root to order at least
\(l+1\). The degree bound \(m+1\) forces \(f\) to have
only one root, with multiplicity \(m\), and \(U\) to
have order \(m+1\) there. That root is nonzero since
a tie occurred.

The consecutive \(Y\)-degrees of the resulting binomial
power in \(P\) require \(\rho s\in\mathbb Z\).
Translate by the resulting root times \(X^{-s}\)
in the same fractional ring. The three top forms
become monomials at the very same exponent points
as the three previously highest monomials. Thus
the next first tie exists, is at a larger slope,
and obeys the same argument. The bound on \(\deg U\)
persists because terms of larger \(Y\)-degree were
below the monomial at the old slope, and lose
weight faster. Since \(\rho\) is fixed and these
slopes are strictly increasing in \((0,1)\) with
\(\rho s\) integral, this is impossible. This rules
out both primitive corners. Faithful flatness
descends the inverse from the algebraic closure.

## Consequences

**Corollary SexticUnrestrictedCubicGenericFibreChart.**
Let a characteristic-zero sextic Keller triple \((P,Q,R)\)
have a polynomial parametrization \(Z(u,v)\) of the generic
fibre \(R=t\) over \(k(t)\) with geometric rank two everywhere.
If all three chart coordinates have degree at most three
in \(v\), then the triple is invertible. The coefficient
degrees in \(u\) are unrestricted.

The [generic-fibre descent](PARTIAL_DEGREE_SLICES_AND_HESSIAN_TRANSFER.md#a-three-variable-generic-fibre-criterion)
gives the planar Keller pair \(P(Z),Q(Z)\) of common
partial degree at most eighteen in \(v\).
Apply the theorem and use the same descent of the inverse.
Existence of the polynomial chart is a separate hypothesis.

**Corollary HessianDirectionalDegreeTwenty.**
Assume Tang's retained complex structural alternative.
Let \(H\) be a four-variable complex polynomial with nonzero
constant Hessian determinant. If three independent constant
directions \(w\) satisfy \(D_w^{21}H\in\mathbb C\), then
the gradient of \(H\) is a polynomial automorphism.

In the cotangent case, one of the three directions
projects to a nonzero direction in the source plane.
The coefficients of the fibre variables in \(D_w^{21}H\)
bound the partial degrees of the two planar Keller
components in that direction by twenty. The theorem
and the [cotangent inversion formula](PARTIAL_DEGREE_SLICES_AND_HESSIAN_TRANSFER.md#four-variable-hessian-potentials-with-small-directional-degrees)
give inversion.
