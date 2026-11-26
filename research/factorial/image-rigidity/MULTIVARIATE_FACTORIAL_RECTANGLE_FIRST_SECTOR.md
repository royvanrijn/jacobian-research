# First-sector angular separation from the exact factorial rectangle

For distinct nonnegative integer exponent columns write
\[
 P_m(z)=\sum_{k=0}^m\binom mk
           \prod_i((m-k)a_i+kb_i)!\,z^k.
\]
Remove equal rows, whose factors are positive scalars. The argument
below compares the two consecutive coefficient arrays by exchanging
the order of two integer products. It gives a first-sector test in
terms of the actual stationary curvature and a separate test with
exactly one private unit row. Both tests compare root arguments at
arbitrary relative radii.

## The exact adjacent coefficient ratio

Take \(n\ge2\), \(N=n+1\), \(L=\log(N/n)\), and compare
\(A=P_n\) with \(B(z)=P_N(\kappa z)\) for \(\kappa>0\).
Write \(R_k=B_k/A_k\) for the ratio of the actual coefficients
at \(0\le k\le n\). For now assume all rows are shared.
Put \(d_i=b_i-a_i\), \(D_i=|d_i|\),
\(\ell_i=\min(a_i,b_i)\), and \(\Delta=\sum_i d_i\).
Let
\[
 {\cal J}_i=
 \begin{cases}
 \{1,\ldots,D_i\}, & d_i>0,\\
 \{0,\ldots,D_i-1\}, & d_i<0 .
 \end{cases}
\]
For \(0\le k<n\) the exact identity is
\[
 \log\frac{R_{k+1}}{R_k}
 =\log\frac{N-k}{n-k}+\log\kappa+\Delta L-\mathcal A_k,\qquad
 \mathcal A_k=\sum_i D_i
  \sum_{j\in{\cal J}_i}
    \int_{(k+j/D_i)/N}^{(k+j/D_i)/n}
                   \frac{dt}{a_i+d_i t}.                         \tag{1}
\]
In fact \(R_k=N\kappa^k g(k)/(N-k)\), with
\(g(k)=\prod_i(na_i+d_i k+1)_{a_i}\).
For \(d_i>0\), putting \(x_i=na_i+d_i k\) gives the
factorial rectangle identity
\[
 \frac{(x_i+d_i+1)_{a_i}}{(x_i+1)_{a_i}}
       =\prod_{j=1}^{d_i}\frac{x_i+a_i+j}{x_i+j}.
\]
For \(d_i<0\) the corresponding product is
\(\prod_{j=0}^{D_i-1}(x_i-j)/(x_i+a_i-j)\).
For a factor in either product use \(s=k+j/D_i\).
Its logarithm is
\(\operatorname{sgn}(d_i)\int_n^N a_i/(m a_i+d_i s)\,dm\).
Subtract its contribution \(\operatorname{sgn}(d_i)L\);
the remaining integral is
\(-D_i\int_{s/N}^{s/n}(a_i+d_i t)^{-1}\,dt\), by \(t=s/m\).
This proves (1).

We use the following first-sector fact. If \(R_k\) strictly
increases, then \(A,B\) have no common root at
\(0<|\arg z|\le2\pi/N\). Indeed at any root
\(z=r e^{i\theta}\) of \(A\) in the positive sector,
\(\sum_{k=0}^n A_k r^k\sin(k\theta)=0\).
The sines change sign just once, from positive to negative,
and both signs occur. Subtract a value of \(R\) between the
last positive and first negative index from the weighted sum
for \(\Im B(z)\). The resulting sum is strictly negative,
and the extra term \(B_Nr^N\sin(N\theta)\) is nonpositive.
The same conclusion holds if only the first ratio increment
can be zero, since both nonzero sine signs occur at indices
at least one. Conjugation treats the negative sector.
[Direct first-sector coefficient-ratio test](MULTIVARIATE_DIRECT_RATIO_FIRST_SECTOR.md)

## Fully shared rows at subcritical stationary curvature

Write
\[
 \begin{split}
 F(t)&=\sum_i\frac{d_i^2}{a_i+d_i t},\qquad
 c=\max_{0\le t\le1}t(1-t)F(t),\\
 V&=\sum_i d_i^2/\ell_i,\qquad
 \eta=\max_i D_i/\ell_i,\qquad
 \widehat\eta=\max_i (D_i-1)/\ell_i .
 \end{split}
\]

**Theorem.** If \(n>\widehat\eta\) and
\[
                 c\left[\frac{N}{n-\widehat\eta}
                                    +\frac{2+\eta}{N}\right]\le1 ,
                                                               \tag{2}
\]
the consecutive moments have disjoint root arguments in
\(0<|\arg z|\le2\pi/N\), for arbitrary signed degree
difference. There are no restrictions on the shared dimension.

For a specified relative radius \(\kappa\), a more general
sufficient test is
\[
 B_0+A_0/2
 \ \le\ 1+|\log\kappa+\Delta L|\,(B_0/A_0+1/2),
 \qquad A_0=V/N,\quad B_0=cN/(n-\widehat\eta).                    \tag{3}
\]
To prove both, take \(q=n-k\) and \(y=(k+1)/N\) in (1).
Each integration interval has length at most \(y/n\).
For a positive gap only \(t<y\) can increase
\((a_i+d_i t)^{-1}\) above its value at \(y\), and then
\(y-t\le(1-1/D_i)/N\). For a negative gap only \(t>y\)
can do so, and then \(t-y\le(1-1/D_i)/n\).
In both cases the denominator is at least
\((1-\widehat\eta/n)(a_i+d_i y)\). Thus
\[
 \mathcal A_k\le
 \min\left\{\frac{V}{N},\
            \frac{yF(y)}{n-\widehat\eta}\right\}
 \le \min\left\{A_0,\frac{B_0}{q}\right\}.                       \tag{4}
\]
The first bound uses the sharper maximum interval length \(1/N\)
and \(a_i+d_i t\ge\ell_i\). The second uses \(1-y=q/N\).

Reverse the columns and conjugate if necessary. This preserves the
angle of a hypothetical contact and sends \(\kappa\) to
\(\kappa^{-1}\) and \(\Delta\) to \(-\Delta\).
The quantities \(A_0,B_0\) are unchanged, since \(F\) is reflected.
We may therefore assume that \(d=\log\kappa+\Delta L\ge0\).
By strict convexity of \(1/x\),
\(\log(1+1/q)>1/(q+1/2)\). The quotient
\[
 \frac{\min\{A_0(q+1/2),\,B_0(1+1/(2q))\}}{1+d(q+1/2)}
\]
is maximized at \(q=B_0/A_0\) among positive real \(q\).
The first expression in the minimum divided by the denominator
increases up to that point, and the second decreases after it.
Condition (3) makes this maximum at most one, so all increments
in (1) are strictly positive. For (2), use
\(V\le4(1+\eta/2)c\), which follows by evaluating \(t(1-t)F(t)\)
at \(t=1/2\). Condition (2) gives \(B_0+A_0/2\le1\)
and hence (3) for every \(\kappa>0\).

In particular \(\widehat\eta=0\) for ordinary unit gaps.
Then \(c\le n/(n+4)\) suffices in (2), since \(\eta\le1\).
For a fixed strict subcritical band \(c\le c_*<1\), (2)
holds from a computable order independent of the profile.
Indeed the individual maxima of
\(t(1-t)d_i^2/(a_i+d_i t)\) are \((\sqrt{a_i}-\sqrt{b_i})^2\).
Their bounds by \(c_*\), together with \(\ell_i\ge1\), give
\(\eta\le c_*+2\sqrt{c_*}\).

## One private unit on one side

**Theorem.** Suppose there is exactly one private row, of
exponent one on its nonzero side. All other unequal rows
are shared. With \(V=\sum_{\rm shared} d_i^2/\ell_i\), the
consecutive moments have disjoint root arguments in
\(0<|\arg z|\le2\pi/N\) whenever \(V\le n\log(N/n)\).
In particular \(V\le4/5\) suffices at every \(n\ge2\).

Orient the private row as \((0,1)\) and put
\(d=\log\kappa+\Delta_{\rm sh}L\), where \(\Delta_{\rm sh}\)
is the signed gap sum over shared rows. The private row
contributes one to \(g(k)\), so (1) still holds with \(d\)
replacing \(\log\kappa+\Delta L\) and with only shared rows
in \(\mathcal A_k\). Its loss is at most
\((k+1)V/(nN)\) by the exact interval lengths. Thus the
ratio increment is at least
\[
                d+\log(1+1/q) -(N-q)V/(nN).
\]
Here \(V\le nL<1\). The last two terms decrease as \(q\)
increases on \(1\le q\le n\), since their derivative is
\(-1/[q(q+1)]+V/(nN)<0\).
Their minimum is \(L-V/(nN)\), and is attained only at \(k=0\).

In reverse orientation the private row is \((1,0)\), with
\(g(k)=N-k\) for that row. It exactly cancels the binomial
factor \(N/(N-k)\) in the coefficient ratio, up to the
constant \(N\). The reversed ratio increments are consequently
at least \(-d-V/N\) by (1) for the reversed shared rows.
If \(d\ge -L+V/(nN)\), the original ratios increase, with
strict increments at every \(k\ge1\). Otherwise
\(-d>L-V/(nN)\ge V/N\), so all reversed increments are
strictly positive. This proves angular separation at every
relative radius. Finally \(nL>2n/(2n+1)\ge4/5\) for \(n\ge2\),
by strict convexity of \(1/x\) on \([n,N]\).

The two theorems concern the first sector. Angles beyond
\(2\pi/N\) may have multiple sine sign changes, and (1) with
the same monotonicity test gives no comparison there.
