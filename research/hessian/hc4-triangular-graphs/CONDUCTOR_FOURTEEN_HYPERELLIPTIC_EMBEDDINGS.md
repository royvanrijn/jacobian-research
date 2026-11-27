# Excluding the degree-two approximate roots at conductor fourteen

**Theorem ConductorFourteenHyperellipticEmbeddingExclusion.** Let
\(F=(P,Q)\) be a complex polynomial plane Keller map in a source frame
with \(P\) monic in the longitudinal variable of positive degree.
The exact generic order cannot have reduced characteristic sequence
\((15,6,2)\) or \((15,10,2)\). Together with the
[degree-two target-coordinate exclusion](HYPERELLIPTIC_CONDUCTORS_AND_SINGULAR_TANGENT_DETERMINANT.md),
this exhausts the embeddings with degree semigroup
\(\langle2,15\rangle\) in
[conductor degree fourteen](CONDUCTOR_FOURTEEN_CHARACTERISTIC_SEQUENCES.md).

Both exclusions follow over the algebraic closure of the generic
parameter field. The necessary approximate-root equation has no
solution with a nonzero degree-two term.

## Normalizing the approximate-root equation

Write \(K=\overline{\mathbb C(y)}\). Take reduced generators \(U,W\) of
degrees \(15,m\), where \(m=6\) or \(10\), and normalize them to be
monic. The degree-two approximate root gives an element
\(T=u^2\) after centering the longitudinal parameter and making an
affine change of values. Even and odd parts of any polynomial in the
order belong to the order because its even part is a polynomial in
\(T\). The odd parts are \(u\) times an ideal in \(K[T]\). Its monic
generator \(G(T)\) has degree seven, since the smallest odd member
of \(\langle2,15\rangle\) is \(15\). Thus \(W=q(T)\), with \(q\)
monic of degree \(k=m/2\), and the odd part of \(U\) is \(uG(T)\).

Use [Assi–García-Sánchez, Proposition 2 and its approximate-root setup](sources/assi-garcia-sanchez-1407.0490v1.pdf).
Its applicability to the actual reduced order is checked in the
[conductor-fourteen classification](CONDUCTOR_FOURTEEN_CHARACTERISTIC_SEQUENCES.md).
The degree-two approximate root, which is affine in \(T\), is monic of degree \(d=15/k\) in \(W\),
with lower coefficients in \(K[U]\). The degrees \(jm\) for
\(0\le j<d\) have distinct residues modulo \(15\), so terms
above longitudinal degree \(30\) cannot cancel. At degree \(30\)
the only possible cancellation is between \(W^d\) and \(U^2\).
Consequently the approximate-root identity is quadratic in
\(U\) with constant nonzero coefficient of \(U^2\).
Complete the square in \(U\) by a polynomial in \(W\). The
square of the resulting generator is even in \(u\). Its odd part
is still \(uG(T)\), so its even part is zero. A translation in
\(W\) depresses the degree-\(d\) polynomial on the right.
The resulting necessary identity is
\[
 T G(T)^2
   =q(T)^d+\sum_{i=1}^{d-2} c_i q(T)^i+c_0+\eta T,
 \qquad \eta\ne0 ,
 \qquad (k,d)=(3,5)\ \text{or}\ (5,3).
\tag{1}
\]
The nonzero coefficient \(\eta\) is required by the longitudinal
degree two of the approximate root.

## Exact coefficient certificates

Write \(q(T)=T^k+\sum_{j=1}^k p_jT^{k-j}\) and
\(F=q^d+\sum_{i=1}^{d-2}c_iq^i\). Successive cancellation of the
coefficients of degrees \(14,13,\ldots,8\) determines \(G\)
uniquely. Starting at \(G=T^7\), the correction at step \(j\) is
\(\frac12[T^{15-j}](F-TG^2)T^{7-j}\).
After these seven steps the residual \(R=TG^2-F\) has degree at
most seven. Equation (1) is equivalent to
\[
 R_2=\cdots=R_7=0,\qquad R_1\ne0,
 \qquad R_i=[T^i]R.
\tag{2}
\]

Scaling \(T\) by \(s\in K^\times\) and renormalizing the leading
coefficients sends \(p_j\) to \(s^{-j}p_j\) and \(\eta\) to
\(s^{-14}\eta\). If the first nonzero \(p_j\) is at index \(j\),
we may therefore set it to \(1\) and the preceding coefficients
to zero. We also retain the chart with every \(p_j=0\).
For each of the four charts when \(k=3\) and the six charts
when \(k=5\), the [exact replay](verify_conductor_fourteen_hyperelliptic_embeddings.py)
generates the polynomials of (2) from the coefficient recursion
over \(\mathbb Q\). With the auxiliary variable \(z\), Singular
constructs coefficients expressing \(1\) as a polynomial linear
combination of
\[
 R_2,\ R_3,\ R_4,\ R_5,\ R_6,\ R_7,\ zR_1-1 .
\]
The replay multiplies the returned coefficients by these
generators and verifies the sum equals \(1\) in the polynomial
ring in each chart. Thus (2) has no solution in any of the charts
over a characteristic-zero field, as asserted.
