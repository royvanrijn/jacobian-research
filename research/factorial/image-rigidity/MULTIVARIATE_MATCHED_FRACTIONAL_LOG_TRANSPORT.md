# Matched transport of assembled log coefficients on the fractional grid

Let \(n\ge2\), \(N=n+1\), \(R=n^2/(n^2-1)\),
\(\chi=1-1/n-1/N\), and \(\Delta=W-RU\).
Use the matched tanh factors, outer losses, and path bounds \(B,E\)
from the [full-kernel transport theorem](MULTIVARIATE_MATCHED_TANH_KERNEL_TRANSPORT.md).
Set
\[
 \lambda=\frac{n}{n-1},\qquad b=\lambda B,\qquad d=1-b
\]
and assume \(b<1\). Define \(V_3=\sum_{i,h}T_{i,h}^3\) at the
lower order. Let \(\ell_{j,k}\) be the logarithm of the
endpoint-normalized factorial coefficient after removal of
\(\binom jk\), as in the
[Cayley expansion](MULTIVARIATE_CAYLEY_OPERATOR_RAY_TEST.md#ordinary-factorial-moments).

**Theorem.** There are explicit continuous functions \(A_n,A_N\)
on \([0,1]\) with \(A_j(k/j)=\ell_{j,k}/j\) at their own knots,
such that
\[
 \|A_N-A_n\|_\infty
 \le
 \Delta\left(\frac{K_0}{2}+\frac{K_1}{4N}
                         +\frac{K_2}{120N^3}\right)
            +\frac{19\lambda^3 V_3}{30n^3d^3}.             \tag{1}
\]
Here
\[
 C_+=\frac{4N E^2}{n-1},\qquad
 C_-=\frac{2(2n+1)E^3}{\chi(n-1)^2},\qquad
 K_m=\max(C_+a_m,C_-b_m)
\]
and the norm entries are
\[
 \begin{array}{c|cc}
 m & a_m/\lambda & b_m/\lambda\\ \hline
 0 & 1/d       & (2-b)/d^2\\
 1 & 2/d^2     & 4/d^3\\
 2 & 8b/d^3    & 24b/d^4 .
 \end{array}
\]
Using \(V_3\le 12N^2\Delta/\chi\) from the
[excess-loss theorem](MULTIVARIATE_TANH_EXCESS_LOSS_TRANSPORT.md),
the last term in (1) is bounded by
\(38\lambda^3N^2\Delta/(5\chi n^3d^3)\).
The estimate compares assembled, normalized log coefficients
at a common fractional index through the specified extensions.
The two adjacent grids have no common interior knot.

## Centered midpoint sums

For a function \(g\) on \([0,1]\) write
\[
 {\cal C}g(s)=\int_0^s g(x)\,dx-s\int_0^1 g(x)\,dx .
\]
In a cell \([\ell/j,(\ell+1)/j]\) put
\(m=(\ell+1/2)/j\), \(v=x-m\), and define
\[
 L_jg(x)=g(m)+v g'(m),\qquad
 Q_j(x)=\frac{v}{24j^2}-\frac{v^3}{6}.
\]
The cell integral of \(L_jg\) is \(g(m)/j\).
The primitive \(Q_j\) vanishes at all knots and satisfies
\[
 Q_j'(x)=\frac1{24j^2}-\frac{v^2}{2},\qquad
 \|Q_j\|_\infty=\frac1{72\sqrt3\,j^3}<\frac1{120j^3}.
\]
For \(g\in C^3([0,1])\) define the corrected interpolation
\[
 {\cal A}_j g(s)={\cal C}(L_jg)(s)-Q_j(s)g''(s).             \tag{2}
\]
It agrees at the knots with the normalized centered midpoint sums.
It also agrees with \({\cal C}\) on all quadratics, as the
expansion below shows.

Taylor's integral remainder about the midpoint gives
\[
 L_jg(x)=g(x)-\frac{v^2}{2}g''(x)+r_j(x),\qquad
 |r_j(x)|\le\frac{|v|^3}{3}\|g'''\|_\infty .
\]
Indeed the error from replacing \(g''(m+t v)\) by \(g''(x)\)
in the integral remainder is bounded by
\(|v|^3\|g'''\|_\infty\int_0^1(1-t)^2dt\).
Therefore \(\int_0^1|r_j|\le\|g'''\|_\infty/(96j^3)\).
Using the expression for \(Q_j'\) and integrating by parts in (2),
\[
 {\cal A}_j g
  ={\cal C}g-\frac1{24j^2}{\cal C}(g'')
             -{\cal C}(Q_jg''')+{\cal C}r_j.              \tag{3}
\]
The two remainder terms have total norm at most
\(7\|g'''\|_\infty/(480j^3)\).
This uses \(\|{\cal C}f\|_\infty\le\|f\|_\infty/2\)
for the first term and \(\|{\cal C}f\|_\infty\le\int_0^1|f|\)
for the second. Also \({\cal C}(g'')(s)\) is the error of
linear interpolation of \(g'\), so its absolute value is
at most \(s(1-s)\|g'''\|_\infty/2\).
Thus (3) gives
\[
 \|({\cal A}_N-{\cal A}_n)g\|_\infty
       \le \frac{19}{480n^3}\|g'''\|_\infty.               \tag{4}
\]
For a changing function \(h\), formula (2) also gives directly
\[
 \|{\cal A}_N h\|_\infty
 \le \tfrac12\|h\|_\infty+\frac{\|h'\|_\infty}{4N}
                             +\frac{\|h''\|_\infty}{120N^3}.       \tag{5}
\]

## Applying the midpoint comparison to factorial logs

The identity
\[
 1-\frac{2\ell}{j-1}
   =\frac{j}{j-1}\left(1-2\frac{\ell+1/2}{j}\right),
 \qquad \frac{N}{n} T'_{i,h}
              =\lambda (T'_{i,h}/R)
\]
puts both orders on the same continuous kernel scale.
Write \(z(s)=1-2s\), \(t_{n,i,h}=T_{i,h}\), and
\(t_{N,i,h}=T'_{i,h}/R\). At each order use
\[
 g_j(s)=\sum_{i,h}\epsilon_i
                   \log(1-\epsilon_i\lambda t_{j,i,h} z(s)),
 \qquad A_j={\cal A}_j g_j .
\]
For an increasing factor the sampled logarithm at the midpoint
is \(\log(1-(1-2\ell/(j-1))T_{j,i,h})\). Its difference from
the original Gauss logarithm is constant in \(\ell\) and disappears
under centering. For a decreasing factor the negative reflection
\(g(s)\mapsto-g(1-s)\) gives the centered sum at \(j-k\).
Indeed its first \(k\) sampled terms are minus the last \(k\)
of the increasing factor, and the full sum has reversed sign.
This proves \(A_j(k/j)=\ell_{j,k}/j\) by (2).

Put \(h=g_N-g_n\). For \(m=0,1,2\) apply kernel transport to
\[
 \phi_{m,\epsilon,s}(t)=
       \partial_s^m\{\epsilon\log(1-\epsilon\lambda t z(s))\}.
\]
This kernel can differ by row. The proof of kernel transport
bounds each row separately, so the shared derivative bounds \(a_m,b_m\)
give \(\|h^{(m)}\|_\infty\le K_m\Delta\).
For clarity, with \(u=\epsilon\lambda t z(s)\) the expressions
needed for the derivative bounds are
\[
\begin{array}{c|cc}
 m & \partial_t\phi_{m,\epsilon,s}
       & 2\partial_t\phi_{m,\epsilon,s}
                            +t\partial_t^2\phi_{m,\epsilon,s}\\\hline
 0 &-\lambda z/(1-u)&-\lambda z(2-u)/(1-u)^2\\
 1 &2\lambda/(1-u)^2&4\lambda/(1-u)^3\\
 2 &-8\epsilon\lambda^2t/(1-u)^3
                         &-24\epsilon\lambda^2t/(1-u)^4 .
\end{array}
\]
Here \(|u|\le b\), and \((2-u)/(1-u)^2\) increases for \(u<1\).
Finally \(|\partial_s^3\{\epsilon\log(1-u)\}|
\le16\lambda^3t^3/d^3\), so
\(\|g_n'''\|_\infty\le16\lambda^3V_3/d^3\).
Decomposing \(A_N-A_n=({\cal A}_N-{\cal A}_n)g_n+{\cal A}_Nh\)
and using (4)–(5) gives (1).

For the Gaussian-subtracted log correction, add
\(jU_j s(1-s)/(j-1)\) to \(A_j(s)\).
The added terms differ by \((N/n)\Delta s(1-s)\).
Also (2) preserves the exact centered linear and quadratic kernels,
so the signed quadratic direction is
\(-\lambda^2(\sum_{i,h}\epsilon_i t_{j,i,h}^2)
\,s(1-s)(1-2s)/3\).
This assembles the logarithmic coefficient transport. Estimating
the alternating sums after exponentiation at the angular gap scale
is a further step.

