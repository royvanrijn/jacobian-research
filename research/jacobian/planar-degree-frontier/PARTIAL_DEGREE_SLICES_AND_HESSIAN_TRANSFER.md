# Partial degree eleven from leading-factor slices

The common partial-degree-eleven theorem has a source-normalization
step at its Guccione–Guccione–Valqui interface. Their published
standardization supplies this step. The leading-factor slice equation
gives another proof with the same normalization, and also applies to
scalar-commutator Weyl pairs.

## One translation for the two leading coefficients

**Proposition AxisLeadingFactorRootTranslation.**
Let \(k\) be algebraically closed of characteristic zero. Take either
a polynomial Keller pair \(P,Q\) in \(k[x,y]\) or a scalar-commutator
pair of polynomial differential operators in \(k[x][D;\partial_x]\).
Write their positive partial \(y\)-degrees or differential orders as
\(au,bu\), where \(a,b\) are coprime. In both settings their leading
coefficients are \(\alpha h^a,\beta h^b\) with a monic \(h\in k[x]\)
and \(\alpha,\beta\ne0\). If \(v=\deg h>u\), there is a common
translation of \(x\) such that the leading coefficients have orders
\(a\mu,b\mu\) at zero for some \(\mu>u\), and their degrees
remain \(av,bv\). In particular, after interchanging the two axes,
both endpoints of each right face have negative \((1,-1)\)-weight.
For Weyl operators use the Fourier interchange, with the sign
that preserves the commutator.

For the polynomial pair, the coefficient of \(y^{(a+b)u-1}\)
in the Jacobian gives \(b p' q-a p q'=0\) for the leading
coefficients \(p,q\). Thus \(p^b/q^a\) is constant, giving the
common factorization by coprimality and unique factorization.
The first polynomial has leading form \(p y^{au}\) at weight
\((0,1)\). By the [plane leading-symbol slice theorem](../../dixmier/dixmier-formal-lifting/WEIGHTED_HAMILTONIAN_SLICE_FULL_BOUNDARY.md),
there is a rational \(G\) with \(J(p y^{au},G)=1\).
Let \(g\in k(x)\) be the coefficient of \(y^{1-au}\) in the Laurent
expansion of \(G\) at \(y=0\). Comparing the coefficient of \(y^0\)
in the bracket gives
\[
                     B'-\frac{h'}{u h}B=1,\qquad B=-au\,p g.
\]
For Weyl pairs the common factorization and this same rational
equation are given by the
[leading-factor residue theorem](../../dixmier/dixmier-formal-lifting/POLYNOMIAL_WEYL_LEADING_FACTOR_RESIDUE.md).

The equation permits no finite pole of \(B\): at a root of \(h\)
of order \(m\ge0\), a Laurent term of initial exponent \(j<0\)
has nonzero coefficient multiplied by \(j-m/u\) at exponent
\(j-1\). Thus \(B\) is a polynomial. At every root of \(h\) the
equation forces a simple zero of \(B\), and its constant term
there is \(B'(1-m/u)=1\), excluding multiplicity \(u\).
If all multiplicities were less than \(u\), then the number
\(s\) of roots would satisfy \(s>v/u>1\). We would have
\(\deg B\ge s\). At infinity the leading term on the left
would then have degree \(\deg B-1>0\) and nonzero multiplier
\(\deg B-v/u\), a contradiction. Hence there is a root
of multiplicity \(\mu>u\).

Translate that root to zero and interchange \(x,y\).
The right-face endpoints are then \((au,a\mu),(au,av)\)
and \((bu,b\mu),(bu,bv)\). In the Weyl case the Fourier map
reverses the axes of the PBW support. Reordering a monomial
after the interchange only adds terms lowering both exponents,
which cannot contribute to the right face. This proves the
same endpoint statement there. This comparison uses actual
leading rows and does not import the planar endpoint gcd
obstruction into the Weyl algebra.

## The endpoint interface and the planar theorem

**Theorem MaxElevenStandardEndpointInterface.**
Over an algebraically closed field of characteristic zero,
the hypothesis named \(\texttt{PlaneKellerStandardEndpointGCDObstruction}\)
in the pinned Lean max-11 theorem holds. With the published
Nagata prime-gcd input, the pinned max-11 theorem implies
polynomial inversion for every planar Keller pair with both
partial degrees at most eleven in a single source coordinate.
The mathematical inversion conclusion descends to any field of
characteristic zero.

For the endpoint interface, let \(a,b>1\) be coprime and
retain the hypotheses of the
[exact Lean definition](https://github.com/dcposch/jc2-lean/blob/6ee763cd623649a5af0517dd12d7202ce2e2ea3c/max11-partial-y/EndpointGCDRoute.lean).
Thus the partial degrees are \(au,bu\), the degrees of the
top coefficients are \(av,bv\), the total degrees are \(aS,bS\),
and \(u<v\). If \(u=0\) both polynomials depend only on \(x\),
contrary to the Keller condition. After swapping \(x,y\),
the right endpoint of the first polynomial is \((au,av)\),
of negative \((1,-1)\)-weight. The two degree ratios in
Definition 4.3 of Guccione–Guccione–Valqui are \(a:b\), so
this is an \((a,b)\)-pair in their polynomial algebra \(L\).
Their Proposition 5.20 makes it standard without changing
the right endpoint or ordinary total degrees, and the
automorphism stays in \(L\). Corollary 7.9 now gives
\(\gcd(u,v)>2\) from the normalized endpoint \((u,v)\).
This is the first conclusion of Corollary 7.9, which
applies to every standard pair. Its second conclusion
has the separate global minimality interpretation for \(B\).
[Primary text, Definition 4.3, Proposition 5.20 and Corollary 7.9](https://arxiv.org/pdf/1401.1784v3)

Alternatively, apply the preceding root-translation proposition
before the swap. It preserves both total degree and the common
coefficient degrees, and makes the start of the right face
have negative \((1,-1)\)-weight as well as its end. This gives
the same standard pair required by Corollary 7.9.

[Nagata's prime total-degree gcd theorem and the Lean headline](EXTERNAL_JC2_AND_HC4_INPUTS.md)
are the two remaining proof inputs for the partial-degree conclusion.
Nagata's constant-one convention is met by rescaling a target
component. The Lean declaration still takes both classical
results as parameters. This argument supplies the mathematical
instances and does not add a formalization of the classical
proofs. For an arbitrary characteristic-zero field, apply the
result after algebraic closure. The map on coordinate rings is
then an isomorphism by faithful flatness of that field extension.

## Four-variable Hessian potentials with small directional degrees

**Corollary HessianFourThreeBoundedDirections.**
Assume the retained complex structural alternative of Tang.
Let \(H\in\mathbb C[x_1,x_2,x_3,x_4]\) have nonzero constant
Hessian determinant. Suppose there are three linearly
independent constant directions \(v\) for which \(D_v^{12}H=0\).
Then \(\nabla H\) is a polynomial automorphism. In particular,
in a linear source frame, three individual partial degrees
may be bounded by eleven with no bound in the fourth direction.

In the cotangent alternative of
[Tang's structural theorem](../../hessian/tang-hc4-audit/STRUCTURAL_REVIEW.md),
write \(H=R(x,y)+zA(x,y)+wB(x,y)\) in linear coordinates.
The two fibre directions span a plane, so one of the three
directions projects to a nonzero vector \(\bar v\) in the
\((x,y)\)-plane. The coefficients of \(z,w\) in \(D_v^{12}H\)
are \(D_{\bar v}^{12}A,D_{\bar v}^{12}B\). Both vanish.
Use \(\bar v\) as the second direction in a linear plane
frame to apply the max-11 theorem to \((A,B)\). The last
two gradient outputs recover \(x,y\) polynomially. The first
two then recover \(z,w\) by the invertible linear system
with determinant \(J(A,B)\in\mathbb C^\times\).
Invertible linear changes of source preserve polynomial
invertibility of the gradient.

## A three-variable generic-fibre criterion

**Corollary KellerThreePolynomialGenericFibreChart.**
Take a characteristic-zero Keller triple \((P,Q,R)\) in
\(k[x_1,x_2,x_3]\). Set \(K=k(t)\), with \(t\) identified
with \(R\) on the generic fibre. Suppose there are three
polynomials \(Z_i(u,v)\in K[u,v]\) giving a map into this
fibre, so \(R(Z_1,Z_2,Z_3)=t\), and the parametrization has
rank two everywhere over \(\overline K\).
If both \(P(Z)\) and \(Q(Z)\) have \(v\)-degree at most eleven,
then the triple is invertible. In particular, when
\(P,Q\) have ordinary total degrees at most eleven, it
suffices that the three \(Z_i\) be affine in \(v\) over
\(K[u]\). This includes the ordinary sextic case.

The columns of \(dZ\) span the tangent plane of \(R=t\)
at every point of the parametrization. Since the full
Jacobian of \((P,Q,R)\) is a nonzero constant, the
differentials of \(P(Z),Q(Z)\) are independent over
\(\overline K\) at every point. Their planar Jacobian is
therefore a nonzero constant, by the Nullstellensatz.
Apply the planar theorem over \(K\).

The generic fibre coordinate ring is the localization of
\(k[x_1,x_2,x_3]\) at \(k[R]\setminus\{0\}\). It is an
integral domain of transcendence degree two over \(K\).
Its homomorphism into \(K[u,v]\) induced by \(Z\) is
injective: a nonzero prime kernel would lower transcendence
degree, while the images of \(P,Q\) are algebraically
independent. Those images generate \(K[u,v]\) by max-11.
The generic fibre ring is therefore \(K[P,Q]\).

To descend generation, clear denominators in \(t=R\) in
the expression for each source coordinate as a polynomial
in \(P,Q\). If an irreducible \(\pi(t)\in k[t]\) remains in
the denominator of a reduced expression, its numerator
gives a nonzero polynomial in \(P,Q\) over \(k[t]/(\pi)\)
vanishing on the fibre \(R=t\bmod\pi\). This is impossible.
After algebraic closure that fibre is nonempty, because
\(R\) is nonconstant, and is smooth of dimension two,
because \(dR\) never vanishes. On each component the
restrictions of \(dP,dQ\) are independent by the Keller
identity. The restrictions of \(P,Q\) are therefore
algebraically independent by the differential criterion.
Hence no irreducible denominator survives, and every
source coordinate belongs to \(k[P,Q,R]\).

A polynomial coordinate completion for \(R\) is one way of
supplying the chart. The theorem allows a chart over the
generic parameter alone and requires its full polynomial
admission and rank control. The rational image fields of
the [axial sextic tuple](../monodromy-theory/SEXTIC_AXIAL_THREE_IMAGE_FIELD.md)
by themselves supply neither requirement. There is no
claim here that the unrestricted residual sextic family
has such a chart.

A useful test of the admission condition is a binary affine
component \(R=f(x)+g(x)y\) with \(g\ne0\).
If \(g\) is constant, the chart
\[
                 (x,y,z)=(u,(t-f(u))/g,v)
\]
is polynomial of rank two and affine in \(v\). This applies
to any \(P,Q\) with both original partial \(z\)-degrees at
most eleven. If \(g\) is nonconstant, \(g(x)\) and \(t-f(x)\)
are coprime in \(k(t)[x]\), since \(t\) is transcendental
over the field of definition of every root of \(g\).
On the generic fibre \(t-f(x)=g(x)y\). Substituting into a
Bézout relation makes \(g(x)\) a unit on that fibre.
Every polynomial parametrization from \(\mathbb A^2_{k(t)}\)
therefore maps \(g(x)\) to a nonzero constant. It maps
\(x\) to a constant since \(g\) is nonconstant, and the fibre
equation then also makes \(y\) constant. Such a chart has
rank at most one. The rational formula for \(y\) in this
case cannot supply the required polynomial chart.
