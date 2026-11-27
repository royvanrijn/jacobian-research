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

Here \(W\) has source degree twelve and \(V\) degree eight.
The monic approximate root has degree three in \(V\) and
source degree three. The terms in \(V^j\) for \(j=0,1,2\)
have degrees distinct modulo twelve, so terms of weight
above 24 cannot cancel, where \(W,V\) have weights twelve
and eight. Normalize both generators to be source monic.
The leading terms of the approximate root are \(V^3-W^2\).
Its coefficient of \(V^2\) is constant and that of \(V\)
is at most affine in \(W\). Complete the square in \(W\)
by a triangular shear in \(V\), then depress the cubic in
\(V\) by a translation. The resulting generators satisfy
\[
 W^2=V^3+kV+C+\eta U,\qquad \eta\ne0,\qquad\deg U=3.
 \tag{2}
\]
The scalars and polynomials are over the parameter field.

Order the target pair as \((V,W)\), and let \(\lambda\)
be the determinant of the ambient target change from the
original \(P,Q\) to this pair. Put \(x=\rho z\) with
\(\rho^{20}=\lambda\) and multiply the two generators
by \(\rho^{-8},\rho^{-12}\), respectively. Over the finite
parameter extension the resulting target determinant is one.
Keep \(V,W\) for the normalized generators. Equation (2)
keeps its form, with rescaled coefficients and \(U\).
The [parameter identity](CONDUCTOR_TWELVE_HYPERELLIPTIC_AND_CUBIC_EXCLUSIONS.md#the-cubicseptic-branch)
therefore gives
\[
 V_z W_y-V_y W_z=J\rho+\partial_z H(V,W)
 \tag{3}
\]
for a polynomial \(H\) over the extended parameter field.
The divergence of the polynomial parameter field of a
target automorphism with determinant one is zero, so this
parameter field is Hamiltonian as in the cubic–septic proof.

In the monic Laurent branch at \(z=\infty\), (2) gives
\[
 W=V^{3/2}+\tfrac{k}{2}V^{-1/2}+O(z^{-9}).
\]
Indeed the contribution of \(C+\eta U\) has order at most
\(3-12=-9\), and the higher powers of the binomial correction
have still smaller \(z\)-degree. The contribution of \(V^{3/2}\)
to the left side of (3) cancels by the chain rule. The
\(O(z^{-9})\) term has no contribution to the polynomial
part of that Jacobian, because \(V_z\) and \(V_y\) have degrees
at most seven. Thus the exact polynomial Jacobian is
\[
 V_z W_y-V_y W_z=k' \partial_z \bigl([\sqrt V]_{\ge0}\bigr).
 \tag{4}
\]
Here \([\sqrt V]_{\ge0}\) is the monic polynomial part of the
Laurent square root, of degree four. This follows by applying
the chain rule to \((k/2)V^{-1/2}\) and then taking the
polynomial part.

It follows from (3)–(4) that \(H(V,W)\) has degree at most
four. The pole semigroup \(\langle3,8\rangle\) has only
nongaps zero and three in this range, so
\(H(V,W)=\beta U+\text{constant}\). Equating the degree-three
coefficients of (3)–(4) forces \(k'=0\). The degree-two
coefficient then forces \(\beta=0\), contradicting the
nonzero \(J\rho\) at degree zero. This excludes the embedding
without a specialized-embedding comparison.

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
