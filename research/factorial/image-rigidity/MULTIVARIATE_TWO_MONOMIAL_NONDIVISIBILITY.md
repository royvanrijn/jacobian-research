# Consecutive nondivisibility for two monomials in several variables

For any number of variables and any positive Gamma shapes, the moment
polynomial of a two-monomial support at order \(m\) cannot divide its
successor at order \(m+1\). In particular, every multivariate ordinary
two-monomial polynomial satisfies Strong Factorial at terminals 2 and 3
over every characteristic-zero field.

## An inequality using the endpoint coefficients

Take distinct exponent vectors \(\alpha,\beta\in\mathbb Z_{\ge0}^s\)
and positive shapes \(\gamma_i\). Set
\[
 C_{m,k}=\prod_{i=1}^{s}(\gamma_i)_{(m-k)\alpha_i+k\beta_i},
 \qquad Q_m(t)=\sum_{k=0}^m \binom m k C_{m,k}t^k .
\]
At ordinary shapes \(\gamma_i=1\), these are the actual factorial moment
polynomials for the two monomials, after dividing an order-\(m\)
moment of \(A X^\alpha+B X^\beta\) by \(A^m\), with \(t=B/A\).

For \(m\ge1\), the ratio
\[
 E_m=\frac{C_{m,1}C_{m,m-1}}{C_{m,0}C_{m,m}}
\]
strictly decreases with \(m\), starting at \(E_1=1\).
Each coordinate's contribution is unchanged by swapping its two
exponents. For exponents \(a\le b\) in that coordinate and \(d=b-a\),
the contribution equals
\[
 \frac{(\gamma+ma)_d}{(\gamma+mb-d)_d}.
\]
If \(d=0\) this is one. If \(d>0\), each quotient
\((\gamma+ma+j)/(\gamma+mb-d+j)\), \(0\le j<d\),
is positive and strictly decreasing for \(m\ge1\).
Its derivative has numerator \(-d(\gamma+a+j)\).
The vectors are distinct, so at least one coordinate decreases strictly.

Suppose \(Q_m\) divides \(Q_{m+1}\). Write \(C_k=C_{m,k}\) and
\(C'_k=C_{m+1,k}\). Their quotient must be \(A+Bt\) with
\(A=C'_0/C_0>0\) and \(B=C'_{m+1}/C_m>0\).
Comparing the coefficients of \(t\) and \(t^m\) gives
\[
 (m+1)C'_1=m A C_1+B C_0,\qquad
 (m+1)C'_m=A C_m+m B C_{m-1}.
\]
The product of these equalities and the arithmetic-geometric mean
inequality give
\[
 (m+1)^2 C'_1 C'_m
 \ge AB\bigl(m\sqrt{C_1 C_{m-1}}+\sqrt{C_0 C_m}\bigr)^2
 \ge AB(m+1)^2 C_1 C_{m-1}.
\]
The second inequality uses \(E_m\le1\). Dividing by
\((m+1)^2 C'_0 C'_{m+1}\) gives \(E_{m+1}\ge E_m\),
a contradiction. The argument also applies when the two compared
coefficients coincide at \(m=1\).

At terminal 2 the lower polynomial is linear, so a common root would
imply the excluded divisibility. At terminal 3 its quadratic
discriminant is \(4(C_{2,1}^2-C_{2,0}C_{2,2})<0\) by \(E_2<1\).
Any common root of \(Q_2,Q_3\) would therefore give a common conjugate
root and make \(Q_2\) divide \(Q_3\).
At ordinary shapes the coefficients are rational, so rational Bézout
transfers both conclusions to every characteristic-zero field.
A nonzero monomial has nonzero factorial moments at every positive order.

More generally, irreducibility over the coefficient field of either
polynomial of a consecutive pair suffices for coprimality. An
irreducible lower polynomial would have to divide the upper one,
and an irreducible upper polynomial has degree too large to be
a divisor of the lower one.

## All orders when both private exponent sums are at least two

At ordinary shapes write
\[
 D_\alpha=\sum_{\beta_i=0}\alpha_i,\qquad
 D_\beta=\sum_{\alpha_i=0}\beta_i .
\]
If \(D_\alpha,D_\beta\ge2\), all consecutive pairs are coprime,
even when the pure factorial moments of the two monomials differ.

For a private exponent \(e\ge1\), the sequence
\(J_e(j)=(ej)!/(j!)^e\) is log-convex, since its consecutive ratio is
\(\prod_{\ell=1}^e(e-(e-\ell)/(j+1))\), which increases with \(j\).
As \(J_e(0)=1\), this gives
\[
 \frac{(ek)!}{((em)!)^{k/m}}
 \le \left(\frac{k!}{(m!)^{k/m}}\right)^e .
\]
A coordinate shared by the two monomials contributes at most one
after normalization by its geometrically interpolated endpoint
factorials, by log-convexity of \(j!\). Thus after the positive
rescaling \(z=t(C_{m,m}/C_{m,0})^{1/m}\), the polynomial
\(Q_m(t)/C_{m,0}\) has the form
\(1+z^m+\sum_{k=1}^{m-1}c_{m,k}z^k\) with
\(0<c_{m,k}\le\binom{m}{k}^{-1}\).
Indeed the private coordinates contribute at most
\(\bigl(k!(m-k)!/m!\bigr)^2\) to the normalized factorial ratio.

For \(m\ge3\), let \(\epsilon_m=\sum_{k=1}^{m-1}\binom m k^{-1}\).
Then \(\epsilon_m\le4(m-2)/(m(m-1))\le2/3\), by using the
first and last terms and bounding each other term by
\(2/(m(m-1))\).
For a root with \(|z|\le1\), we have \(|1+z^m|\le\epsilon_m\).
For \(|z|>1\) apply the same argument to the reciprocal polynomial.
In either case its argument \(\theta\), chosen in \([0,\pi]\)
by conjugation, satisfies
\[
 |m\theta-(2j_m+1)\pi|\le e_m,\qquad e_m=\arcsin\epsilon_m,
\]
for some integer \(j_m\ge0\). This is the angular width of the disk
centered at \(-1\) with radius \(\epsilon_m\).

If a root belongs to both orders \(m,m+1\), its argument
is the same after the two positive scalings. Subtracting the
two angular inequalities gives
\(|\theta-2(j_{m+1}-j_m)\pi|\le e_m+e_{m+1}<\pi\).
Thus \(j_{m+1}=j_m\) and \(\theta\le e_m+e_{m+1}\).
For \(0\le x\le2/3\), convexity of \(\arcsin x\) gives
\(\arcsin x\le1.1x\), since
\(\sin(11/15)\ge 11/15-(11/15)^3/6>2/3\).
Consequently
\[
 (2j_m+1)\pi\le(m+1)e_m+m e_{m+1}<8.8<3\pi.
\]
For the strict bound, each summand is below \(4.4\)
by \(\epsilon_n\le4(n-2)/(n(n-1))\). Therefore \(j_m=0\).

A root of \(Q_m\) at an angle \(0<\theta\le\pi/m\) is impossible,
since each term has nonnegative imaginary part and the linear term
has positive imaginary part. Angle zero is excluded by positivity.
On the other hand \(e_n<\pi/(n-1)\) for every \(n\ge4\).
For \(n=4\) use \(\epsilon_4=2/3\). For \(n\ge5\), bound the first
two and last two terms exactly and bound the rest by
\(\binom n3^{-1}\) to obtain
\[
 (n-1)\epsilon_n
 \le 2+\frac{8n-34}{n(n-2)}\le\frac83.
\]
The last comparison reduces to \(2(n-7)^2+4>0\).
Hence \((n-1)e_n\le44/15<\pi\).
With \(n=m+1\) and \(j_m=j_{m+1}=0\), the angular inequality
forces \(\theta\le(\pi+e_{m+1})/(m+1)<\pi/m\), a contradiction.
Terminals 2 and 3 were covered above, and rational Bézout gives
the assertion over any characteristic-zero field.

In particular, all disjoint nonconstant supports are coprime
at all orders. The private exponent proof applies when both total
degrees are at least two. If one is degree one, orient it as the
first monomial. Independence of its variable from the second
monomial gives
\[
 Q_m(t)/m!=\sum_{k=0}^{m}\frac{(k\beta)!}{k!}\,t^k,
 \qquad (k\beta)!:=\prod_{i=1}^s(k\beta_i)!.
\]
The successive polynomials on the right differ by a single
nonzero term of top degree and have constant term one, so they
are coprime.



## Varying one Gamma shape

Choose a coordinate on which the exponents differ. Swap the two
monomials, if needed, so that \(d_i=\beta_i-\alpha_i>0\) in that
coordinate. Hold the other positive shifts in a field \(K\) and
replace \(\gamma_i\) by an indeterminate \(z\). In the monic reciprocal
of \(Q_m(t)/C_{m,0}\), every nonleading coefficient is a nonzero
constant in \(K\) times \((z+\alpha_i m)_{d_i k}\), where \(1\le k\le m\).
It has valuation precisely one at \(z+\alpha_i m\).
The reciprocal is therefore Eisenstein over \(K[z]\) at this prime
and irreducible over \(K(z)\). This proves consecutive coprimality
when that shift is transcendental over \(K\), including for supports
with differences of both signs. The ordinary specialization has all
shifts equal to one. The remaining multivariate assertion is coprimality
at every terminal \(N\ge4\) for arbitrary ordinary two-monomial supports.

