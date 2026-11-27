# The alternate embeddings of the conductor-fourteen cubic order

A polynomial order in one variable with pole semigroup
\(\langle3,8\rangle\) has numerical conductor fourteen. The
[degree-fourteen characteristic-sequence classification](CONDUCTOR_FOURTEEN_CHARACTERISTIC_SEQUENCES.md)
lists four reduced embeddings with this semigroup. The
[direct cubic–octic exclusion](CONDUCTOR_FOURTEEN_CUBIC_OCTIC_EXCLUSION.md)
treats the embedding of reduced type \((8,3)\). Here we exclude
\((8,6,3)\) for rational polynomial orders and \((12,8,3)\) for
the directional order of a polynomial plane Keller map. The
[three-stage exclusion](CONDUCTOR_FOURTEEN_THREE_STAGE_EMBEDDING.md)
treats \((12,8,6,3)\). The proofs below use the degree-three approximate root supplied by
[Assi–García-Sánchez, Propositions 2 and 4](sources/assi-garcia-sanchez-1407.0490v1.pdf).

## A sextic lower generator and an octic upper generator

**Theorem ConductorFourteenCubicOrderEightSixThreeExclusion.** A
two-generator order \(K[W,V]\) in \(K[z]\) with fraction field
\(K(z)\), over a characteristic-zero field \(K\), cannot have
reduced one-place sequence \((8,6,3)\).

In this sequence \(\deg W=6\), \(\deg V=8\), and the pole
semigroup is \(\langle3,8\rangle\). We may work over an
algebraic closure of the coefficient field for the one-place
theorem. The approximate root is monic of degree four in \(W\)
and has source degree three. An element of degree three may
be made monic and depressed by scaling and translating the
source. In degrees at most six the order has basis consisting
of the constant, this monic degree-three element and its square,
by the semigroup degree filtration. Completing the square in
the expression for \(W\), then affinely transforming \(W\),
gives
\[
 W=U^2,\qquad U=z^3+az+h.
\]
The constant \(h\) is needed since the center of the square is
fixed by the given sextic generator. The degree-three approximate
root is affine in \(U\), with nonzero slope.

Take \(V\) monic. The powers \(W^j\) for \(0\le j<4\) have
degrees distinct modulo eight, so no terms of weight above
24 can cancel in the approximate root. Here the weights of
\(W,V\) are six and eight. The leading terms are \(W^4\) and
\(-V^3\), since they are the only terms at weight 24. The
coefficient of \(V^2\) is affine in \(W\). Depress the cubic
in \(V\) by adding an affine function of \(W\) to \(V\), a
triangular target automorphism that preserves its monic leading
term. The approximate root identity now has the form
\[
 V^3+A(U)V+B(U)=0,\qquad A(-U)=A(U),\qquad
 B(U)-B(-U)=2\gamma U,\qquad\gamma\ne0.
 \tag{1}
\]
Indeed the remaining coefficients in the approximate root are
polynomials in \(W=U^2\), whereas the approximate root itself
is affine in \(U\) with nonzero slope. The polynomial in (1)
is minimal over \(K(U)\): the field \(K(z)\) has degree three
over \(K(U)\), and \(V\) of degree eight cannot belong to
\(K(U)\). In fact the pole order at infinity of any rational
function of \(U\) is a multiple of three.

Use the \(K[U]\)-basis \(1,z,z^2\) of \(K[z]\). The trace
of \(z\) over \(K(U)\) is zero and the trace of \(z^2\) is
\(-2a\), from \(z^3+az+h-U=0\). Equation (1) has trace zero
for \(V\). Thus its source presentation in this basis is
\[
 V=q(U)z+p(U)(z^2+2a/3),\qquad
 p=U^2+p_1U+p_0,\qquad q=q_2U^2+q_1U+q_0.
\]
The degree bounds follow by distinct degree residues modulo
three, and the leading coefficient of \(p\) follows from
monicity of \(V\).

Taking the norm of \(Y-V\) along the cubic equation for \(z\)
gives the coefficients
\[
\begin{aligned}
 A&=a q^2-3(U-h)pq-a^2p^2/3,\\
 B&=-(U-h)^2p^3+(U-h)ap^2q-(U-h)q^3
      -2a^3p^3/27-2a^2pq^2/3 .
\end{aligned}
\]
Write \(A_i=[U^i]A\) and \(B_i=[U^i]B\). In order,
\(A_5=B_7=A_3=0\) imply
\[
 q_2=0,\qquad p_1=2h/3,\qquad q_0=h(-4a^2+9q_1)/27 .
\]
If \(h=0\), then \(p\) is even and \(q\) odd in \(U\), and
the displayed formula gives even \(B\), contrary to (1).
If \(h\ne0\), then \(B_5=0\) gives
\(p_0=(4a^3-9a q_1-5h^2)/27\).
With this substitution \(A_1=8h^3(q_1-a^2)/27\), so
\(q_1=a^2\). Now \(B_3=64h^5/729\ne0\), again contrary
to (1). This rules out the rational polynomial embedding.

## An octic lower generator and a degree-twelve upper generator

**Theorem ConductorFourteenCubicOrderTwelveEightThreeExclusion.**
For a complex polynomial plane Keller map, in a longitudinal
frame with \(P\) monic, the exact directional order cannot
have reduced sequence \((12,8,3)\).

Apply Theorem DirectionalThreeTwoLaurentExclusion below with
\(m=4\), \(r=3\).

## The general three-to-two Laurent obstruction

**Theorem DirectionalThreeTwoLaurentExclusion.** Let \(F=(P,Q)\)
be a complex polynomial plane Keller map with \(P\) monic in
the longitudinal source coordinate \(x\). Its reduced generic
directional order cannot have one-place characteristic sequence
\((3m,2m,r)\) with \(2\le r\le m\) and \(\gcd(m,r)=1\).

Set \(K=\mathbb C(y)\) and let the reduced generators \(W,V\)
have source degrees \(3m,2m\), respectively. The one-place
approximate root for the gcd \(m\) is monic of degree three in
\(V\) and has source degree \(r\). For its terms in \(V^j\)
with \(j=0,1,2\), the weights \(2mj\) are distinct modulo
\(3m\), so no terms of weight above \(6m\) can cancel. After
normalizing \(V,W\) to be source monic, the terms of weight
\(6m\) are \(V^3-W^2\). The coefficient of \(V^2\) is
constant and that of \(V\) is at most affine in \(W\).
Complete the square in \(W\) by a triangular shear in \(V\),
then depress the cubic in \(V\) by a translation. The
resulting generators satisfy
\[
 W^2=V^3+kV+C+\eta U,\qquad \eta\ne0,\qquad \deg_x U=r ,
 \tag{2}
\]
where \(U\) is a scalar normalization of the approximate root.

Order the target pair as \((V,W)\), and write \(\lambda\)
for the determinant of the ambient target change from \(P,Q\)
to this pair. Put \(x=\rho z\) with \(\rho^{5m}=\lambda\) and
multiply the two generators by \(\rho^{-2m},\rho^{-3m}\).
The resulting target determinant is one over the finite
parameter extension. Keep \(V,W\) for these normalized
generators. Equation (2) retains its form with rescaled
coefficients and \(U\). The
[parameter identity](CONDUCTOR_TWELVE_HYPERELLIPTIC_AND_CUBIC_EXCLUSIONS.md#the-cubicseptic-branch)
gives
\[
 V_z W_y-V_y W_z=J\rho+\partial_z H(V,W)
 \tag{3}
\]
for a polynomial \(H\) over that extension. Here \(J\) is the
original nonzero constant Jacobian. The parameter derivative
field of the determinant-one target change has divergence zero,
so its contribution in (3) is Hamiltonian.

At \(z=\infty\), the monic Laurent branch of (2) gives
\[
 W=V^{3/2}+\tfrac{k}{2}V^{-1/2}+O(z^{r-3m}).
\]
The contribution of \(C+\eta U\) to this expansion has degree
at most \(r-3m\). The quadratic and later binomial corrections
have degree at most \(-5m\). Differentiation in \(y\) preserves
the degree bound on the tail and differentiation in \(z\)
lowers it. Since \(V_z,V_y\) both have degree at most \(2m-1\),
the tail makes no contribution to the polynomial part of
the Jacobian. The chain rule cancels the term \(V^{3/2}\)
and gives the exact polynomial identity
\[
 V_z W_y-V_y W_z=k'\partial_z\bigl([\sqrt V]_{\ge0}\bigr),
 \tag{4}
\]
where \([\sqrt V]_{\ge0}\) is monic of degree \(m\).

From (3)–(4), \(H(V,W)\) has degree at most \(m\). Its
possible degrees in that range are multiples of \(r\) by the
semigroup \(\langle3m,2m,r\rangle\). The condition
\(\gcd(m,r)=1\) with \(r\ge2\) excludes degree \(m\).
Thus \(\partial_z H\) has degree at most \(m-2\). Comparing
degree \(m-1\) in (3)–(4) forces \(k'=0\). This makes
\(\partial_z H=-J\rho\) a nonzero constant, contrary to the
absence of degree one in the order.

In particular the theorem excludes the conductor-18 sequence
\((15,10,3)\) by taking \(m=5\), \(r=3\).

## Exact replay

The [checker](verify_conductor_fourteen_cubic_order_alternate_embeddings.py)
reconstructs the trace-zero norm and every coefficient elimination
in the first proof. It also checks the Laurent-to-polynomial
chain-rule identity used in the second proof on formal
coordinates. The source degree bounds justify the omitted
Laurent terms in the written argument.

~~~sh
.venv/bin/python research/hessian/hc4-triangular-graphs/verify_conductor_fourteen_cubic_order_alternate_embeddings.py
~~~
