# Matched tanh transport controlled by excess outer loss

Let \(n\ge2\), \(N=n+1\), and \(R=n^2/(n^2-1)\). Use an ordinary
two-monomial profile of distinct exponent columns in any finite
dimension. On unequal row \(i\) let \(\ell_i\ge0\) be the smaller
exponent, \(d_i\ge1\) the gap, and \(\epsilon_i\in\{-1,1\}\) its
orientation. The row contributions to the lower and upper
reflected outer losses are
\[
 U_i=\frac12\sum_{h=1}^{d_i}
       \log\frac{n\ell_i+h+(n-1)d_i}{n\ell_i+h},
 \qquad
 W_i=\frac12\sum_{h=1}^{d_i}
       \log\frac{N\ell_i+h+nd_i}{N\ell_i+h}.
\]
Write \(U=\sum_iU_i,\ W=\sum_iW_i\), and \(\chi=1-1/n-1/N>0\).
With the [matched tanh coordinates](MULTIVARIATE_MATCHED_TANH_MOMENTS.md),
let \(M_p=\sum_{i,h}\epsilon_i T_{i,h}^p\) and
\(M_p'=\sum_{i,h}\epsilon_i(T_{i,h}')^p\). For a uniform
constant put
\[
 \begin{aligned}
 L_i&=\ell_i+d_i,&
 b_i&=\frac{(n-1)d_i}{2n\ell_i+(n-1)d_i},\\
 B&=\max_i b_i,&
 E&=\max_i\frac{L_i b_i}{d_i}.
 \end{aligned}
\]
Both \(B\) and \(E\) are at most one.

**Theorem.** For every real \(p\ge1\),
\[
 \left|M'_p-R^pM_p\right|
 \le R^p B^{p-1}
     \max\left\{\frac{4pN}{n-1}E^2,\
               \frac{2p(p+1)(2n+1)}{\chi(n-1)^2}E^3\right\}
     (W-RU).
 \tag{1}
\]
The same bound applies to the unsigned moments by taking all
\(\epsilon_i=1\). This uses the integer block of Gauss shifts on
each row. It ties the transport error to the surplus in the outer
loss ratio shared by the two orders. Private rows are included.

For instance, for \(p=2\) and \(n\ge6\) the first entry
of the maximum in (1) dominates, so the signed skew
\(P=M_2/2\) obeys
\[
 \left|P'-R^2P\right|
       \le R^2\frac{4N}{n-1}\,B E^2(W-RU).
 \tag{2}
\]
The inequality between the two entries follows from
\(2(n^2-n-1)(n-1)\ge 3n(2n+1)\) for \(n\ge6\). After
subtracting the right side the cubic is \(2n^3-10n^2-3n+2\),
whose value and derivative are positive at \(6\), with positive
second derivative thereafter.
For fully shared rows with lower reciprocal cap
\(d_i/(n\ell_i+1)\le t\), admissible choices are
\[
 B=\frac{(n^2-1)t}{2n+(n^2-1)t},
 \qquad
 E=\frac{(n-1)(1+Nt)}{2n+(n^2-1)t}.
 \tag{3}
\]
Indeed \(d_i/\ell_i\le Nt\), and each expression for a row is
increasing in this ratio. Either bound can also use the actual
largest factor reached in the contraction below.

## Two contributions to the row gain

Fix a row and omit its index. With \(y=1/j\), write
\[
 -\log u_j=\frac{d(1-y)}2 J(y),\qquad
 J(y)=\sum_{h=1}^d\int_0^1\frac{ds}{D(h,s)},\qquad
 D(h,s)=\ell+(1-y)ds+yh .
\]
The [outer-ratio proof](MULTIVARIATE_MATCHED_OUTER_RATIO_MONOTONICITY.md)
shows how to pair the summands in \(-J'(y)\) after writing
\(s=(r-b)/d\), \(1\le r\le d,\ 0\le b\le1\).
For \(h>r\) set \(a=h-r\) and denote the paired denominators by
\(D_1,D_2\), with \(D_2-D_1=(1-2y)a\). On
\(1/N\le y\le1/n\), both are positive and at most \(L=\ell+d\).
The paired sum is
\[
 a(D_1^{-2}-D_2^{-2})+b(D_1^{-2}+D_2^{-2}).
\]
The \(b\)-terms, including the diagonal, integrate to at least
\(d/(2L^2)\), since \(D_1,D_2\le L\) and \(\int_0^1 b\,db=1/2\).
The \(a\)-term is at least \(2(1-2y)a^2/L^3\) by the mean
value theorem. The integration weight of each \((h,r)\) term in
\(-J'\) is \(db/d\), and
\(\sum_{1\le r<h\le d}(h-r)^2=d^2(d^2-1)/12\).
It follows that
\[
 -J'(y)\ge \frac{d}{2L^2}
               +\frac{(1-2y)d(d^2-1)}{6L^3}.
\]
Integrating from \(1/N\) to \(1/n\), whose length is \(1/(nN)\),
and using \(W-RU=\tfrac{dn}{2N}[J(1/N)-J(1/n)]\) on the
single row, gives
\[
 W-RU\ge G_++G_-,
 \quad G_+=\frac{d^2}{4N^2L^2},
 \quad G_-=\frac{\chi d^2(d^2-1)}{12N^2L^3}.
 \tag{4}
\]
Here the average of \(1-2y\) on the integration interval
is \(\chi\).

## Signed transport on the integer grid

Put \(A=(n-1)/2,\ c=n/N\), and \(f(s)=A/(n\ell/d+A+s)\).
Then \(T_h=f(h/d)\), and direct substitution of the upper
coordinate gives
\[
 \frac{T_h'}R=f\!\left(\tfrac12+c(h/d-\tfrac12)\right).
\]
For \(\lambda\in[c,1]\) write \(f_{\lambda,h}=f(1/2+
\lambda(h/d-1/2))\). All these values are at most the row
bound \(b_i=A/(n\ell/d+A)\) given above. The derivative
\(f'=-f^2/A\) gives
\[
 R^{-p} \sum_h(T_h')^p-\sum_h T_h^p
    =\frac pA\int_c^1
         \sum_h(h/d-\tfrac12) f_{\lambda,h}^{p+1}\,d\lambda.
\]
Decompose \(h/d-1/2=m+v_h\), with \(m=1/(2d)\) and
\(v_h=(2h-d-1)/(2d)\). The \(m\)-contribution is nonnegative
and bounded above by \(p b_i^{p+1}/(n^2-1)\). The \(v_h\)
contribution is nonpositive, since its grid is symmetric and
\(f\) decreases. Pairing \(v_h\) and \(-v_h\) and bounding the
derivative of \(f^{p+1}\) shows
\[
 0\le-\sum_h v_h f_{\lambda,h}^{p+1}
 \le \frac{(p+1)\lambda}{A}\, b_i^{p+2}
                                 \sum_h v_h^2.
\]
Now \(\sum_h v_h^2=(d^2-1)/(12d)\). Hence the magnitude
of the \(v_h\)-contribution is at most
\[
 \frac{p(p+1)(2n+1)}{6(n^2-1)^2}
                  \frac{d^2-1}{d}\,b_i^{p+2}.
\]
The first majorant equals
\(\frac{4pN}{n-1} b_i^{p-1}(L_i b_i/d_i)^2 G_+\).
The second equals
\(\frac{2p(p+1)(2n+1)}{\chi(n-1)^2}
 b_i^{p-1}(L_i b_i/d_i)^3 G_-\), including its
vanishing at \(d=1\).
Thus the absolute row difference is at most the
larger uniform coefficient in (1), before the factor \(R^p\),
times \(G_++G_-\). Sum over rows with either orientation and
use (4).

An additional consequence of (4) is the bound
\[
 W-RU\ge \frac{\chi}{12N^2}\sum_{i,h} T_{i,h}^3 .
\]
In fact \(T_{i,h}\le b_i\le d_i/L_i\), so the \(G_-\) term
controls the fraction \(1-d_i^{-2}\) of the row's cubic moment
with the displayed constant, and \(G_+\) controls the rest
since \(L_i\ge d_i\ge1\). This complements the
[bound in terms of the squared signed second moment](MULTIVARIATE_MATCHED_OUTER_RATIO_MONOTONICITY.md#coupling-to-the-lower-order-tanh-moments).
The inequalities keep the shape information needed in paired
coefficient enclosures. They do not by themselves compare root
arguments.

