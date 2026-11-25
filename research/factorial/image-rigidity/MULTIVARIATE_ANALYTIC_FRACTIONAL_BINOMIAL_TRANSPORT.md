# From matched fractional transport to alternating binomial sums

Let \(a,b\) be distinct ordinary exponent columns with all unequal
rows shared. Put \(n\ge2\), \(N=n+1\), \(d_i=|b_i-a_i|\), and
\(\ell_i=\min(a_i,b_i)>0\) on the unequal rows. Use
\[
 A_j(s)=\frac1j\sum_i\left[\log\Gamma(1+j((1-s)a_i+sb_i))
            -(1-s)\log\Gamma(1+ja_i)-s\log\Gamma(1+jb_i)\right].
\]
Then \(e^{jA_j(k/j)}\) is the normalized factorial multiplier
with the binomial coefficient removed. The definitions at
\(s=0,1\) give \(A_j(0)=A_j(1)=0\).

Let \(0<r<\min_i\ell_i/d_i\). Set
\[
 B_r=\sum_i\frac{d_i^2}{\ell_i-rd_i},\qquad
 D_r=\frac1{2nN}\sum_i\frac{d_i^2}{(\ell_i-rd_i)^2},\qquad
 C_r=\frac{(1+r)^2}{2},\qquad
 M_r=N C_r D_r\,e^{N C_rB_r}.
\]
The following comparison concerns the actual order-\(N\) array
and the degree-\(N\) lift of the lower Gamma chord:
\[
 Q(z)=\sum_{k=0}^N\binom Nk
       \bigl(e^{N A_N(k/N)}-e^{N A_n(k/N)}\bigr)z^k .
\]
For all complex \(z\), with \(a_z=\min(1,|z|)\), one has
\[
 |Q(z)|
 \le M_r\sum_{m=1}^N\binom Nm a_z^m|1+z|^{N-m}
                        \frac{m!}{(Nr)^m}.                    \tag{1}
\]
In particular
\[
                    |Q(-1)|\le M_r\frac{N!}{(Nr)^N}.           \tag{2}
\]
If \(U,W\) are the matched outer losses and \(R=n^2/(n^2-1)\),
the [rowwise excess-loss bound](MULTIVARIATE_TANH_EXCESS_LOSS_TRANSPORT.md#two-contributions-to-the-row-gain)
also gives
\[
  D_r\le \frac{2N}{n}
       \max_i\left(\frac{\ell_i+d_i}{\ell_i-rd_i}\right)^2(W-RU).
                                                               \tag{3}
\]

## Analytic Gamma chords

The chords \(A_j\) are analytic on a neighborhood of the closed
complex \(r\)-neighborhood of \([0,1]\). Throughout that neighborhood,
\(t_i(s)=(1-s)a_i+sb_i\) has positive real part at least
\(\ell_i-rd_i\). The [fractional-index Gauss
integral](MULTIVARIATE_FRACTIONAL_INDEX_TRANSPORT.md#the-comparison-kernel)
applies there by analyticity. In particular
\[
  A_j''(s)=\sum_i d_i^2 \int_0^\infty
       e^{-t_i(s)u}\,g(u/j)\,du,\qquad g(x)=\frac{x}{e^x-1}.
\]
On \(x\ge0\) one has \(0<g(x)\le1\) and
\(0\le-g'(x)\le1/2\). For the derivative upper bound,
clearing denominators reduces the inequality to
\(2x\le e^x-e^{-x}\). Thus \(\|A_j''\|_\infty\le B_r\)
and \(\|A_N''-A_n''\|_\infty\le D_r\).
The latter uses
\(|g(u/N)-g(u/n)|\le u/(2nN)\) before integrating.
Formula (3) follows by summing the row lower bounds
\(d_i^2/[4N^2(\ell_i+d_i)^2]\) for the excess loss.

If \(H(0)=H(1)=0\) and \(\|H''\|_\infty\le D\) on this
neighborhood, either real endpoint has derivative bounded
by \(D/2\) by integrating along \([0,1]\). Each complex \(s\)
is within \(\tfrac12+r\) of at least one endpoint. Taylor's
integral formula from that endpoint, along a segment in the
neighborhood, yields \(\|H\|_\infty\le C_rD\).
Apply this to each chord and their difference. For
\(f=e^{NA_N}-e^{NA_n}\), integrate the derivative of
\(e^{N((1-t)A_n+tA_N)}\) over \(0\le t\le1\).
This gives \(\|f\|_\infty\le M_r\), with both endpoints zero.

## Newton differences at the root scale

Write \(\Delta_h f(x)=f(x+h)-f(x)\). Newton interpolation at
\(h=1/N\) and the binomial formula give
\[
  Q(z)=\sum_{m=1}^N\binom Nm z^m(1+z)^{N-m}\Delta_h^m f(0).
\]
The missing term of order zero is \(f(0)=0\). The \(m\)-fold
integral for the forward difference uses \(f^{(m)}\) at points
in \([0,m/N]\). Each such point has a complex disk of radius
\(r\) inside the domain. Cauchy's estimate gives
\[
            |\Delta_h^m f(0)|\le\frac{m!M_r}{(Nr)^m}.
\]
Apply the same proof to \(f(1-s)\) and \(1/z\), then multiply by
\(z^N\), to get the bound with \(1\) in place of \(|z|\) in the
\(m\)-th power. Taking the smaller bound proves (1).
Setting \(z=-1\) retains only \(m=N\), proving (2).

This is a bound after exponentiation and alternating summation
of the actual Gamma chords at their respective fractional knots.
An angular theorem still requires separation of the lower
polynomial from its degree-\(N\) lift and a comparison of (1) to
the local modulus on contours around the relevant roots.
