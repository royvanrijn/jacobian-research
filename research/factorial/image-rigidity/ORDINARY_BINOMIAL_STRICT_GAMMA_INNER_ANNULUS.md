# A strict inner annulus for Gamma-power moment polynomials

For \(\gamma>0\) and integers \(m,d\ge1\), put
\[
T_{m,\gamma,d}(t)=\sum_{k=0}^{m}\binom mk(\gamma)_{dk}t^k,
\qquad
L(\gamma,H)=\gamma+\frac{4\gamma H}{33\gamma+5H},
\]
where \((x)_k=x(x+1)\cdots(x+k-1)\). Every zero \(t\) of
\(T_{m,\gamma,d}\) with \(m\ge2\) satisfies
\[
(L(\gamma,d(m-1)))_d < |1/t| < (\gamma+d(m-1))_d .
\tag{1}
\]
The upper bound is the retained Gamma-product annulus. The improvement
in the lower bound comes from combining the differential equation at a
minimum-modulus root with a second-moment estimate for all roots.

## The individual Gamma factor

Write \(h=m-1>0\) and
\[
G_{m,\gamma}(z)=\sum_{k=0}^{m}\binom mk(\gamma)_kz^{m-k}.
\]
Use (A1)–(A4) in the [Gamma annulus proof](GAMMA_PRODUCT_SHARP_LOWER_ANNULUS_AND_SHAPE_SEPARATION.md#1-the-classical-factor-bound-with-an-elementary-proof)
and the real-part estimate in the
[Gamma matrix disk lemma](GAMMA_PRODUCT_ROOT_MODULUS.tex).
The zeros \(z_j\) are simple and satisfy
\(|z_j|<\gamma+h\) and \(\Re z_j\le -\gamma+h\).
Choose \(z_i\) of minimum modulus \(r\). Taking the real part of the
differential equation at \(z_i\) gives
\[
 r-\gamma\ \ge\ -\Re z_i-\gamma
   =\sum_{j\ne i}\frac{|z_j|^2-r^2}{|z_i-z_j|^2}\ \ge\ 0.
\tag{2}
\]
In particular \(\Re z_i\le-\gamma\). The first two coefficients
after the leading coefficient are \(m\gamma\) and
\(\binom m2\gamma(\gamma+1)\). By Vieta,
\[
 \sum_j z_j=-m\gamma,\qquad
 \sum_j z_j^2=m^2\gamma^2-mh\gamma(\gamma+1)
              =m\gamma^2-mh\gamma.
\]
Using Cauchy–Schwarz on the real parts,
\[
 \sum_j|z_j|^2
 =2\sum_j(\Re z_j)^2-\Re\sum_j z_j^2
 \ge m\gamma^2+mh\gamma .
\tag{3}
\]

Set \(s=h/\gamma\), \(\tau=4/(33+5s)\), and \(Q=2\tau+s\tau^2\).
Suppose \(r\le\gamma+h\tau\), so
\(r^2-\gamma^2\le\gamma h Q\). For each \(j\),
\[
 |z_j+\gamma|^2<(4+s)\gamma h,\qquad
 |z_i+\gamma|^2\le r^2-\gamma^2\le Q\gamma h .
\]
The first inequality follows by expanding the square and using the modulus
and real-part bounds. The second uses \(\Re z_i\le-\gamma\).
The inequality \((u+v)^2\le(5/4)u^2+5v^2\) gives
\[
 |z_j-z_i|^2 < \gamma h\!\left[\tfrac54(4+s)+5Q\right].
\tag{4}
\]
For \(\ell=33+5s\), completing a square gives
\[
 \frac{Q}{\tau}(1+5\tau)
 =\frac{14+148/\ell-2640/\ell^2}{5}
 \le\frac{14+148^2/(4\cdot2640)}{5}
 =\frac{10609}{3300}<\frac{13}{4}.
\]
Since \(1/\tau-\tfrac54(4+s)=13/4\), this is equivalent to
\(1-Q>\tau[\tfrac54(4+s)+5Q]\). In particular \(1-Q>0\).
The sum of the numerators in (2) is at least \(m\gamma h(1-Q)\)
by (3). Combining (2) and (4) gives
\[
 r-\gamma > \frac{m(1-Q)}{\tfrac54(4+s)+5Q}
             >m\tau>h\tau ,
\]
a contradiction. Thus every root of \(G_{m,\gamma}\) has modulus
strictly greater than \(L(\gamma,m-1)\).

## The product bound and its ordinary consequences

The rising-factorial identity
\((\gamma)_{dk}=d^{dk}\prod_{j=0}^{d-1}((\gamma+j)/d)_k\)
expresses the reciprocal \(x^mT_{m,\gamma,d}(1/x)\) as the scaled
binomial convolution of the \(d\) Gamma factors with shapes
\((\gamma+j)/d\). The reciprocal form of
[closed-disc Grace composition](GAMMA_PRODUCT_SHARP_LOWER_ANNULUS_AND_SHAPE_SEPARATION.md#2-every-positive-shape-a-sharper-product-annulus)
multiplies their strict lower modulus bounds. For fixed \(h\),
the shift \(L(\alpha,h)-\alpha\) increases in \(\alpha\).
Replacing the shift of each factor by its value at \(\gamma/d\),
then scaling by \(d^d\), yields (1).

For the ordinary actual pair \(P_{N-1,a,d}=T_{N-1,a(N-1)+1,d}\)
and \(P_{N,a,d}=T_{N,aN+1,d}\), with integer \(a\ge0\), \(d\ge1\),
\(N\ge2\), (1) gives coprimality whenever
\[
 a+\frac{4(aN+1)d(N-1)}{33(aN+1)+5d(N-1)}\ge d(N-2).
\tag{5}
\]
Indeed the upper-order lower radius is then at least the lower-order
upper radius, with strict separation of the root sets. This includes
equality in (5). For \(N=2\) the lower-order polynomial has degree one
and the upper bound in (1) is an equality at its root. A rational
Bézout identity gives coprimality over every characteristic-zero
coefficient field.

In particular, for \(N\ge25\), the ordinary pair is coprime at every
point of the [ordinary containing region](ORDINARY_BINOMIAL_RESIDUAL_REGION.md)
satisfying
\[
 25a\ge (22N-47)d .
\tag{6}
\]
Indeed (6) implies \(a\ge15d\), so \(5d(N-1)/(aN+1)\le1/3\).
The shift in (5) is then at least \(3d(N-1)/25\), and this plus
the lower bound on \(a\) in (6) gives (5).

For any nonconstant monic rational common factor \(F(x)\) of the
actual reciprocal pair with degree \(k\) and \(N\ge3\), (1) improves
the Archimedean norm restriction to
\[
 \bigl((L(aN+1,d(N-1)))_d\bigr)^k < |F(0)|
   <\bigl((a(N-1)+1+d(N-2))_d\bigr)^k .
\tag{7}
\]
Take the product of the moduli of the roots of \(F\), using the lower
bound for \(P_{N,a,d}\) and the upper bound for \(P_{N-1,a,d}\).
Both bounds in (7) are rational and can be compared to integer norms
by clearing denominators.

## A stronger bound at larger shape-to-order ratios

For the individual factor \(G_{m,\gamma}\), the bound improves to
\(r>\gamma+h/6\) when \(\gamma\ge15h\). To prove it, suppose instead
that \(r-\gamma\le h/6\), and set
\[
 \eta=\frac{r^2-\gamma^2}{\gamma h},\qquad
 s=\frac h\gamma,\qquad q_j=|z_j|^2-r^2.
\]
Then \(0\le\eta\le 1/3+s/36\le181/540<1\).
Writing \(u_j=z_j+\gamma\), the identities used in (3) give
\(\sum_j u_j=0\) and \(\sum_j u_j^2=-mh\gamma\).
The retained modulus and real-part bounds also give
\(|\Re u_j|\le h\). Therefore
\[
 \sum_j|z_j-z_i|^2
 =\sum_j|u_j|^2+m|u_i|^2
 \le m\gamma h(1+2s+\eta).
\]
By (3), \(\sum_j q_j\ge m\gamma h(1-\eta)\), and the retained
upper modulus bound gives \(q_j\le\gamma h(2+s-\eta)\).
Use weighted Cauchy–Schwarz in (2) to obtain
\[
 r-\gamma
 \ge\frac{(\sum_j q_j)^2}{\sum_j q_j|z_j-z_i|^2}
 \ge m\,\frac{(1-\eta)^2}{(2+s-\eta)(1+2s+\eta)} .
\]
The last rational function decreases as either \(s\) or \(\eta\)
increases for \(s\ge0\), \(0\le\eta<1\). In \(\eta\), this
follows by logarithmic differentiation since
\(-2/(1-\eta)+1/(2+s-\eta)<0\), and the remaining term is negative.
At \(s=1/15,\eta=181/540\) the function is
\(128881/741455>1/6\). This contradicts \(r-\gamma\le h/6\).

Thus the same multiplicative Grace argument gives
\[
 (\gamma+d(m-1)/6)_d < |1/t|
 \qquad\text{if }\gamma\ge15d(m-1),\quad m\ge2 .
\tag{8}
\]
In this range (8) can replace the lower radius in (1) and (7).
In particular, at \(N\ge25\), the ordinary pair is coprime on the
larger strip
\[
 6a\ge(5N-11)d .
\tag{9}
\]
Indeed (9) gives \(a\ge15d\) and
\(a+d(N-1)/6\ge d(N-2)\), so the strict upper-order lower
radius in (8) is larger than every lower-order root modulus.
