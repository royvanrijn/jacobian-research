# Consecutive moment polynomials for two monomials in several variables

For any number of variables and any positive real Gamma shapes, the
moment polynomial of a two-monomial support at order \(m\) cannot
divide its successor at order \(m+1\). At terminals 2, 3, and 4 the
consecutive polynomials are coprime over \(\mathbb C\). In particular,
ordinary multivariate binomials satisfy Strong Factorial at these three
terminals over every characteristic-zero field.

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

## Terminal 4 for arbitrary supports

**Theorem.** Let \(\alpha,\beta\) be distinct nonnegative integer
exponent vectors, and let every Gamma shape \(\gamma_i\) be a positive
real number. Then \(Q_3\) and \(Q_4\) have no common complex root.

We use three coefficient invariants, which are unchanged by scaling
the variable by a positive real number or multiplying each polynomial
by a nonzero scalar. They are also unchanged by reversing both
polynomials.
\[
 c=\sqrt{\frac{C_{3,1}C_{3,2}}{C_{3,0}C_{3,3}}},
 \qquad
 u=\sqrt{\frac{C_{4,1}C_{4,3}}{C_{4,0}C_{4,4}}},
 \qquad
 v=\frac{C_{4,2}}{\sqrt{C_{4,0}C_{4,4}}}.
\]
They satisfy
\[
 0<c<1,\qquad c^{3/2}\le u<c,\qquad
 \frac{v^2}{u^2}\le\frac{2+u^2}{3}.                         \tag{1}
\]
Indeed \(c^2=E_3\) and \(u^2=E_4\), so the strict inequalities
follow from the endpoint ratios above. In an individual row with
ordered exponents \(a\le b\) and gap \(d=b-a>0\), the negative logarithm
of each factor of \(E_m\) is
\[
 \log\!\left(1+\frac{d(m-1)}{\gamma+ma+j}\right),\qquad 0\le j<d .
\]
The inner fraction is increasing and concave in \(m\ge1\) since its
derivative is \(d(\gamma+a+j)/(\gamma+ma+j)^2\). The logarithm preserves
concavity here. Interpolating between \(m=1\), where the loss is zero,
and \(m=4\) gives \(-\log E_4\le-\tfrac32\log E_3\).
For the last inequality in (1), in this row put \(A=\gamma+4a\).
The factors in \(v^2/u^2\) and \(u^2\) compare as
\[
 \frac{A+d+j}{A+2d+j}
 \le \frac{2+(A+j)/(A+3d+j)}{3}.
\]
Multiply over the rows and \(j\), using
\(\prod_i\frac{2+x_i}{3}\le\frac{2+\prod_i x_i}{3}\) for \(0<x_i\le1\).
This product inequality follows by induction from
\((2+x)(2+y)\le 3(2+xy)\). Equal-exponent rows contribute one.

We also need a bound for the imbalance between the two interior
coefficient ratios of the cubic:
\[
 \left|\log\frac{C_{3,2}^3 C_{3,0}}{C_{3,1}^3 C_{3,3}}\right|
 \le 2\log\frac{1+c^2}{2c}.                                \tag{2}
\]
For one row with increasing exponents, set \(x_j=d/(\gamma+3a+j)\).
Its factors in the ratio on the left are \((1+x_j)^2/(1+2x_j)\),
while its factors in \(c^2\) are \(1/(1+2x_j)\).
Reversing one row inverts the first ratio and leaves the second
unchanged. Thus the left side of (2) is at most the sum of the logarithms
of \((1+x_j)^2/(1+2x_j)\) over the sorted rows. Writing
\(z_j=1+2x_j\ge1\), this bound follows from
\[
 \prod_j\frac{1+z_j}{2}\le\frac{1+\prod_j z_j}{2},
 \qquad \prod_j z_j=c^{-2}.
\]
The two-factor inequality is \((1+z)(1+w)\le 2(1+zw)\) for
\(z,w\ge1\), so it again iterates over any number of factors.

Suppose there is a common root. The Gamma integral makes \(Q_4(t)>0\)
for real \(t\), since the fourth power of the corresponding real
binomial is nonnegative and is positive on a nonempty open subset
of the positive orthant. Therefore the common root is nonreal.
Scale the variable by its modulus and make both polynomials monic.
Their common conjugate factor and complementary factors have the form
\[
 (t^2+h t+1)(t+s),\qquad
 (t^2+h t+1)(t^2+p t+q),
 \qquad -2<h<2,\quad s,q>0,\quad p^2<4q .
\]
Here the last inequality holds because the quartic has no real roots.
By reversing the polynomials if necessary, arrange \(s\ge1\). Put
\(S=(s+s^{-1})/2\), \(Z=(\sqrt q+1/\sqrt q)/2\), and \(P=p/\sqrt q\).
Coefficient comparison now gives
\[
 9c^2=1+h^2+2hS,\qquad
 16u^2=P^2+h^2+2hPZ,\qquad 6v=2Z+hP.                     \tag{3}
\]

If \(h\le0\), positivity of the cubic coefficients \(s+h\) and
\(1+hs\) gives \(h>-1\). Since \(S\ge1\), the cubic equation in (3)
gives \(3c\le1+h\). Also \(Z\ge1\) and \(P<2\), so
\(6v\ge2+2h\ge6c\). This contradicts \(v\le u<c\) from (1).

Assume \(h>0\). Set \(h_*=3c-1\), \(x=1-c\), and
\(\Delta=h_*-h\). The cubic equation in (3) gives \(h\le h_*\), so
\(0<x<2/3\), and
\[
 \Delta=\frac{h(s+s^{-1}-2)}{2+h_*+h}.                     \tag{4}
\]
We claim that \(\Delta\le 4x^2/27\). Log-convexity of the cubic
normalized coefficients gives
\((s+h)^2\le3(1+hs)\), or \(s^2-hs+h^2\le3\).
Since the left side at \(s=2\) equals \(3+(h-1)^2\) and increases
for \(s\ge2\), we have \(s\le2\).
Write \(w=s^{1/3}\). Substituting the cubic coefficients in (2)
and applying the tangent bound for the concave power \(2/3\) gives
\[
 \frac{w(w^3+h)}{1+hw^3}-1
 \le \left(\frac{1+c^2}{2c}\right)^{2/3}-1
 \le \frac{x^2}{3c}.                                     \tag{5}
\]
The left side equals
\((w-1)(w+1)w(w+w^{-1}-h)/(1+hw^3)\).
We have \((w^2+w+1)/(w+1)<2\) for \(1\le w\le2^{1/3}\).
On the same interval,
\[
 \frac{1+hw^3}{(1+h)w^{5/2}}
 \le \frac{1+2w^3}{3w^{5/2}}\le1 .
\]
The first inequality uses \(h\le2\) and \(w^3\ge1\).
For the second, the derivative of \(w^{-5/2}+2w^{1/2}\) is
\(w^{-7/2}(w^3-5/2)<0\). Since
\(\sqrt{s+s^{-1}-2}=(w-1)(w^2+w+1)/w^{3/2}\), (5) yields
\[
 \sqrt{s+s^{-1}-2}\le \frac{2(1+h)x^2}{3c(2-h)}.
\]
Using (4), \(h\le h_*\), and \(2-h\ge3x\) proves the claimed bound:
\[
 \Delta\le\frac{2h(1+h)x^2}{81c^2}
 \le\frac{2(3-1/c)x^2}{27}\le\frac{4x^2}{27}.             \tag{6}
\]

The quartic equations in (3) give a lower bound on \(v\).
Its positive \(t\) and \(t^3\) coefficients imply
\(P>-h\min(\sqrt q,1/\sqrt q)\).
If \(P\le0\), then \(P(P+2hZ)\le0\) and \(4u\le h\).
But \(u\ge c^{3/2}>h/4\), using \(c>1/3\) and \(h\le3c-1\le2c\).
Consequently \(P>0\) and
\[
 P=-hZ+\sqrt{16u^2+h^2(Z^2-1)} .
\]
At fixed \(h,u\), the expression \(2Z+hP\) increases in \(Z\ge1\).
Its derivative is
\[
 2-h^2+\frac{h^3 Z}{\sqrt{16u^2+h^2(Z^2-1)}}
 \ge 2-h^2+\frac{h^3}{4u}\ge 2-h^2+\frac{h^3}{4}>0 .
\]
The first bound uses \(4u>h\), and the last polynomial decreases
to zero as \(h\) increases from \(0\) to \(2\).
We obtain \(6v\ge2+h(4u-h)\) by evaluating at \(Z=1\).

Finally put \(y=1-u\). By (1), \(x\le y\le3x/2<1\), where
the upper bound follows from the tangent inequality for \(c^{3/2}\).
Using \(\sqrt{1-z}\le1-z/2\) on the last bound in (1) gives
\[
 6v\le6(1-y)\sqrt{1-(2y-y^2)/3}
      \le6-8y+3y^2-y^3.
\]
The difference of the lower and upper bounds for \(6v\), using
\(h=2-3x-\Delta\), is
\[
 -(3x+\Delta)^2+4y(3x+\Delta)-3y^2+y^3 .
\]
It is concave in \(y\le1\). At \(y=x\) it equals
\(x^3-2x\Delta-\Delta^2>0\) by (6). At \(y=3x/2\) it equals
\(9x^2/4+27x^3/8-\Delta^2>0\), also by (6).
This contradicts (3) and proves the theorem.

For ordinary shapes the two polynomials have rational coefficients,
so Bézout over \(\mathbb Q\) extends the conclusion to every
characteristic-zero field. This argument places no bound on the
exponents or on the number of variables.

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
at every terminal \(N\ge5\) for arbitrary ordinary two-monomial supports.

