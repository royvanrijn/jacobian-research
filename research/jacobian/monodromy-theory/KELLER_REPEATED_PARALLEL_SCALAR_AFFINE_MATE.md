# Affine mates of a polynomial scalar linear in a source direction

**Theorem KellerRepeatedParallelScalarAffineMateExclusion.**
Let \(g,A\) be univariate polynomials over a field of
characteristic zero, with \(A\) nonconstant. If the polynomial
scalar \(r=g(x)+A(x)v\) has Keller mates \(P,Q\) in the three
variables \(x,v,y\), then \(\deg_y P,\deg_y Q\ge2\).
There is no ordinary degree bound on \(g,A,P,Q\).

For \(A=x^2,\ g=x\) this covers the exceptional cubic scalar.
It also applies after a polynomial source shear with constant
nonzero Jacobian, using the free variable in the sheared frame.

## Localization by the affine mate

Work over \(\mathbb C\) by embedding the finitely generated
coefficient field. The Keller condition makes \(r\) a submersion.
Consequently at every root \(a\) of \(A\) we have \(A'(a)=0\)
and \(g'(a)\ne0\). Otherwise the two partial derivatives
\(r_v=A(x)\) and \(r_x=g'(x)+A'(x)v\) would simultaneously
vanish for some \(v\) at \(x=a\).

Put \(B=\mathbb C[x,v]\) and
\[
 D=-A(x)\partial_x+(g'(x)+A'(x)v)\partial_v .
\]
The determinant equation, up to orientation, is
\[
                 D(P)Q_y-P_yD(Q)=\kappa\ne0.
\]
The derivation fixes \(r\) and acts in \(\mathbb C(r,x)\) as
\(-A(x)\partial_x\) at fixed \(r\). Its rational constants
are \(\mathbb C(r)\).

Suppose \(Q=\alpha y+b\) with \(\alpha,b\in B\) and \(\alpha\ne0\).
Over \(B[\alpha^{-1}]\) take \(Q\) as the free coordinate
instead of \(y\). The determinant equation becomes
\(\alpha D(P)|_Q=\kappa\). Thus all positive coefficients of
\(P\) in \(Q\) are invariant. We obtain
\[
              P=\Psi(Q,r)+T,\qquad
 \Psi\in\mathbb C(r)[Q],\quad T\in B[\alpha^{-1}],
 \quad \alpha DT=\kappa.                                  \tag{1}
\]
Evaluate at \(y=0\). Then
\(T=P|_{y=0}-\Psi(b,r)\) belongs to \(B[H(r)^{-1}]\)
for a nonzero polynomial \(H\) clearing the denominators
of the coefficients of \(\Psi\). This ring is preserved by
\(D\), so (1) makes \(\alpha\) a unit in \(B[H(r)^{-1}]\).

## The components of a scalar fibre

For each \(c\in\mathbb C\) put
\[
 L_c(x)=\prod_{\substack{A(a)=0\\g(a)=c}}(x-a),
 \qquad J_c=(r-c)/L_c(x).
\]
At a root of \(L_c\) the polynomial \(g-c\) has a simple
zero. Thus \(J_c\) is primitive and linear in \(v\), and is
irreducible in \(B\). The distinct \(x-a\) for \(a\) in
this product and \(J_c\) are the full factorization of \(r-c\).
In particular \(r-c\) vanishes to order one on each component.

Unique factorization of \(H(r)\) now gives \(\alpha\) as a
nonzero scalar times a product of powers of these fibre
components. Substituting \(J_c=(r-c)/L_c(x)\) shows
\(\alpha=a_r(r)a_x(x)\) with nonzero rational factors.
Thus (1) has the form
\(\partial_x T=R(r)W(x)\) with nonzero rational factors.
There is a rational primitive \(V(x)\in\mathbb C(x)\) of \(W\).
Indeed one may specialize \(T/R\) at any value of \(r\)
for which its rational expression in \(x\) is defined.
Absorbing the remaining invariant into \(\Psi\), we get
\[
                         P=\Psi(Q,r)+R(r)V(x).             \tag{2}
\]
The function \(V\) is nonconstant.

Choose a root \(a\) of \(A\) and put \(c=g(a)\).
Let \(i\) be the multiplicity of \(x-a\) in \(\alpha\),
\(j\) the multiplicity of \(J_c\), and \(m\ge2\) the
multiplicity of \(a\) in \(A\). In (2), \(R\) has pole
order \(j\) at \(c\), and \(W\) has order \(j-i-m\) at
\(x=a\), since \(\partial_x T=-\kappa/(A\alpha)\).
If \(j=0\), then \(W\) has a pole of order at least two,
so \(V\) has a pole at \(x=a\). Whenever \(V\) has a pole
at \(a\), the two orders of \(R(r)V(x)\) differ:
the pole at \(x=a\) has order greater than \(j\), while on
\(J_c=0\) the order of pole is just \(j\).
If \(V\) is regular at \(a\), then \(j>0\). Replace \(V\)
by \(V-V(a)\) and absorb the constant multiple of \(R\)
into \(\Psi\). Now \(R(r)V(x)\) has pole order \(j\) on
\(J_c=0\) and strictly smaller pole order on \(x=a\).
In either case the orders differ and at least one is positive.

The restriction of \(Q\) to either component is nonconstant.
If \(\alpha\) does not vanish identically there, use its
\(y\)-coefficient. Otherwise the determinant equation gives
\(-\overline{P_y}\,D(\bar b)=\kappa\) on that component,
forcing \(\bar b\) nonconstant. This reduction is valid
because \(D(x-a)=-A(x)\) is divisible by \(x-a\) and
\(D(J_c)=A L_c'J_c/L_c\) is divisible by \(J_c\).

If \(\Psi\) has a coefficient pole at \(r=c\), let \(k\) be
the largest such order. Its leading polar coefficient
is a nonzero complex polynomial in \(Q\). It remains
nonzero on either component because \(Q\) restricts
nonconstantly to each integral component. Therefore
\(\Psi(Q,r)\) has the same pole order \(k\) along both.
If the coefficients are regular at \(c\), it has poles
on neither. This contradicts regularity of \(P\) in (2).

Finally, a mate independent of \(y\) would reduce the
determinant equation to a product of its \(D\)-derivative
and the other mate's \(y\)-derivative. Both factors would
be units, so the other mate would be affine in \(y\),
which has already been excluded.
