# Excluding degree-two elements of directional Keller orders

**Theorem KellerDirectionalDegreeTwoElementExclusion.**
Let \((P,Q):\mathbb C^2\to\mathbb C^2\) be a polynomial map with nonzero
constant Jacobian, with \(P\) monic of positive \(x\)-degree over
\(\mathbb C[y]\). If \(\overline{\mathbb C(y)}[P,Q]\) contains an element
of \(x\)-degree two, the map is a polynomial automorphism. In particular,
the generic directional order of a noninvertible map cannot have degree
semigroup \(\langle2,N\rangle\) for any odd \(N\ge3\).

The constant Jacobian lets us reduce the generator degrees using target
changes over \(\mathbb C\). The first approximate root then makes the
reduced even part a polynomial shear in the other generator. The odd part
of the Jacobian forces the shear coefficients to be constant. The resulting
pair is injective on the fixed line of the source reflection.

## The reflection and its polynomial centre

Assume the map is noninvertible and write \(K=\mathbb C(y)\). A quadratic
element over \(\overline K\) is a finite linear combination of monomials in
\(P,Q\). The conditions on this same finite list to cancel all coefficients
above degree two and prescribe leading coefficient one are linear over
\(K\). Consistency over \(\overline K\) gives consistency over \(K\).
Completing the square in this monic quadratic gives \(T=u^2\) in the order
\(A=K[P,Q]\), where \(u=x-h(y)\) and \(h\in K\).

The even part of each element of \(A\) is in \(K[T]\subset A\), so its
odd part belongs to \(A\) as well. The odd parts form \(u g(T)K[T]\) for a
nonzero polynomial \(g\) because their coefficients after division by
\(u\) form an ideal of \(K[T]\), and the normalization is \(K[u]\).
The conductor in \(K[u]\) is \(g(u^2)K[u]\). In fact, a polynomial
\(a(T)+u b(T)\) lies in the conductor precisely when it and its product
with \(u\) both lie in \(A\), which requires \(g\mid a,b\).
Thus the degree semigroup is \(\langle2,N\rangle\), where
\(N=2\deg g+1\) is the smallest odd degree. The
[exact directional conductor theorem](GENERIC_COORDINATE_ORDER_CONDUCTOR_AND_CLOSED_COLLISIONS.md)
gives a polynomial generator \(\kappa(x,y)\) of the conductor that is
nonzero and nonconstant on every finite line \(y=y_0\). In particular
\(N\ge3\). It also follows from the even and odd parts that every
element of degree less than \(N\) is even in \(u\).

To prove \(h\in\mathbb C[y]\), fix any finite \(y_0\) and a finite root of
\(\kappa(x,y_0)\). A small disk around the root with no zeros on its
boundary contains a root \(x_+\) of \(\kappa(x,y)\) for all sufficiently
nearby \(y\), by continuity of the coefficients and Rouché's theorem.
At generic nearby parameters the root \(x_+\) and its reflection
\(x_-=2h(y)-x_+\) have equal \(P,Q\) values, by the description of \(A\).
The polynomial \(P(X,y)-P(x_+,y)\) in \(X\) is monic with bounded
coefficients there, so its root \(x_-\) stays bounded. Hence \(h\) has no
pole at \(y_0\). The substitution \(x=u+h(y)\) is now a polynomial
source automorphism of determinant one.

## Constant target degree reduction

Suppose a polynomial Keller pair has longitudinal degrees \(n,m>0\)
and leading coefficients \(a(y),b(y)\). The coefficient of
\(u^{n+m-1}\) in its Jacobian is \(n a b'-m a'b\). This vanishes, so
if \(a\) is constant, then \(b\) is constant. Since the original \(P\)
is monic, the starting pair has both leading coefficients constant.
When one degree divides the other, subtract a constant multiple of a
power of the lower-degree polynomial to cancel the higher leading
term. The unchanged polynomial still has constant leading
coefficient, so the cancellation can be repeated. A degree-zero
generator would make \(A\) a one-generator polynomial order with
fraction field \(K(u)\), forcing the other generator to have degree
one. A degree-one element would give \(A=K[u]\), contrary to the
nonconstant conductor. We arrive by constant target automorphisms at
\((U,W)\) of degrees \(n>m>1\) with \(m\nmid n\).

## The degree sequence and the first approximate root

Apply the one-place-at-infinity theorem of Assi and García-Sánchez,
[Section 1, Proposition 2 and Section 6](sources/assi-garcia-sanchez-1407.0490v1.pdf).
Over \(\overline K\) the polynomial normalization has exactly one
place at infinity. The defining equation is monic of degree \(n\) in
the variable for \(W\), by \(K(U,W)=K(u)\) and \([K(u):K(U)]=n\).
Its coefficient at index \(i\) has degree at most \(\lfloor im/n\rfloor\)
in \(U\). This follows by taking elementary symmetric functions of
the conjugates of \(W\) and using their pole order \(m/n\) when \(U\)
has pole order one. The first coefficient is constant in \(U\) and
can be removed by a translation in the equation over \(K\). The
approximate roots are defined over \(K\) by uniqueness and can be
expressed in the original variables. The norm of the translated
\(W\) has degree \(m\) in \(U\), by the same unique place at infinity
and its leading pole. Section 6 identifies the intersection
multiplicities with the degrees of the evaluated polynomials for
this polynomial normalization.

Use the characteristic sequence \(r_0=n,r_1=m,\dots,r_s\) and
\(d_i=\gcd(r_0,\dots,r_{i-1})\), \(e_i=d_i/d_{i+1}\), with strictly
decreasing gcds ending at \(d_{s+1}=1\). The entries generate
\(\langle2,N\rangle\), so both \(2\) and \(N\) occur. If \(m=2\),
then \(n\) is odd and the two entries generate the semigroup, giving
\(n=N\). In this case \(W\) is affine in \(T\) over \(K\), and the
even part of \(U\) belongs to \(K[W]\).

Suppose \(m>2\). Define
\[
C_i=\sum_{j<i,\ j\ge1}(e_j-1)r_j-r_0+d_i.
\]
The conductor formula gives \(C_{s+1}=N-1\). Also
\(C_2=d_2(n/d_2-1)(m/d_2-1)\ge2d_2\), and
\[
C_{i+1}-C_i=(d_i-d_{i+1})(r_i/d_{i+1}-1)\ge0.
\]
If the gcd just before entry \(2\) is even, that entry changes it to
\(2\). The first odd entry must then be \(N\) and must be last,
since it changes the gcd to one. Its increment \(N-1\) exceeds the
remaining budget \(N-1-C_2\). Thus \(2\) is last and \(d_s\) is odd,
with \(C_s=N-d_s\). If \(r_i=N\) with \(2\le i<s\), then
\[
C_{i+1}\ge 2d_2+N-d_i+d_{i+1}>N,
\]
also impossible. Hence \(N\) is one of \(n,m\). The odd integer
\(N/d_2\) is at least three. If the other degree divided by \(d_2\)
were at least three, we would have \(C_2\ge2(N-d_2)>N\). Thus
\(n=N\) and \(m=2d_2\).

For every intermediate index \(2\le i<s\), freeness gives
\(e_i r_i\in\langle r_0,\ldots,r_{i-1}\rangle\). The smallest
positive member of this prefix is at least \(2d_i\) by induction,
starting with \(m=2d_2\), so \(r_i\ge2d_{i+1}\). The sum of all
intermediate increments is \(C_s-C_2=d_2-d_s\), which is exactly
the sum of their lower bounds \(d_i-d_{i+1}\). Consequently
\(r_i=2d_{i+1}\) at each intermediate step. In particular \(r_2<N\)
whether or not there is an intermediate step.

The first nontrivial approximate root is monic of degree \(e_1=N/d_2\)
in \(W\) with lower coefficients in \(K[U]\). Its value on the
parametrization has degree \(r_2<N\), so it is even in \(u\). The
same holds for \(W\). The longitudinal degrees of the nonzero
terms \(c_j(U)W^j\) with \(0\le j<e_1\) are distinct modulo \(N\).
Since the term \(W^{e_1}\) has degree \(2N\), no lower term can
have larger degree and the only term that can cancel it at degree
\(2N\) is a nonzero constant times \(U^2\). Thus the approximate
root has the form \(\gamma U^2+b(W)U+c(W)\), with \(\gamma\in K^\times\).
Subtract its reflected value at \(-u\). The difference is the
nonzero difference \(U(u)-U(-u)\) times
\(\gamma(U(u)+U(-u))+b(W)\) and must vanish. Therefore the even
part of \(U\) is \(-b(W)/(2\gamma)\) in this case as well.

## The fixed line

Write \(U=V+L(W,y)\), where \(V\) is odd in \(u\) and \(L\) is a
polynomial in its first argument over \(K\). The odd part \(V\)
is already in \(\mathbb C[u,y]\), as it is the odd part of the
polynomial \(U\). The identity
\[
\operatorname{Jac}_{u,y}(W,U)=
\operatorname{Jac}_{u,y}(W,V)+W_u(\partial_y L)(W,y)
\]
has a constant left side. Its first term on the right is even, and
its second is odd. The latter vanishes, so \(\partial_y L=0\)
by \(W_u\ne0\) and the transcendence of \(W\) over \(K\).
The shear is therefore a target automorphism over \(\mathbb C\).

Write \(J_*\) for the nonzero constant Jacobian of \((W,V)\).
On \(u=0\) we have \(V=0\) and
\[
J_*=-W_y(0,y)V_u(0,y).
\]
Both factors are polynomials in \(y\), so each is a nonzero
constant. Thus \(W(0,y)\) is affine of degree one. The restriction
of \((W,V)\) to the line \(u=0\) is injective, and
[Gwoździewicz's Theorem 1.1](sources/gwozdziewicz-alg-geom-9305008v1.pdf)
makes it an automorphism. The source and target changes used on the
Keller map were polynomial automorphisms, giving the contradiction.

Any degree-two target coordinate belongs to the directional order, so
the theorem gives the [target-coordinate exclusion](HYPERELLIPTIC_CONDUCTORS_AND_SINGULAR_TANGENT_DETERMINANT.md)
as a corollary. In the [two conductor-fourteen embeddings](CONDUCTOR_FOURTEEN_HYPERELLIPTIC_EMBEDDINGS.md)
\((15,6,2)\) and \((15,10,2)\), the last approximate root supplies a
degree-two element by Proposition 2. The same deduction applies to
every characteristic sequence containing an entry \(2\), with no bound
on the conductor.
