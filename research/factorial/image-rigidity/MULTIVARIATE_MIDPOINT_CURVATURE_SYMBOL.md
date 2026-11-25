# Signed midpoint curvature in the normalized factorial symbol

For distinct fully shared ordinary exponent columns, omit the equal
rows and set
\[
 d_i=b_i-a_i,\qquad M_i=\frac{a_i+b_i}{2},\qquad
 S_*=\sum_i\frac{d_i^2}{M_i},\qquad
 \delta=\max_i\frac{|d_i|}{M_i}.
\]
The signed cubic sum and the two nonnegative correction sums are
\[
 K=\sum_i\frac{d_i^3}{M_i^2},\qquad
 B_4=\sum_i\frac{d_i^4}{M_i^3},\qquad
 \Lambda=\sum_i\frac{d_i^2}{M_i^2}.
\]
Write \(x_i(v)=a_i+d_i v\), and define the factorial log symbol
relative to the Gaussian of curvature \(S_*\) by
\[
 \begin{split}
 E_j(v)=\sum_i\big[&
   \log\Gamma(jx_i(v)+1)
  -(1-v)\log\Gamma(ja_i+1)-v\log\Gamma(jb_i+1)\big]\\
  &+\frac{jS_*}{2}v(1-v).
 \end{split}
\]
Thus the exact normalized coefficients are
\(\binom jk e^{-S_* k(j-k)/(2j)} e^{E_j(k/j)}\).
The Gaussian splitting in the
[ordinary diffuse lift](MULTIVARIATE_DIFFUSE_POLYNOMIAL_ANGLE_LIFT.md)
uses this type of analytic symbol.

**Lemma.** Suppose \(\delta\le1/60\), and put
\(L=(11/2)|K|+34B_4\). With suprema on \(|v|\le5\), for every
\(n\ge1\) and \(N=n+1\),
\[
 \|E_j''\|\le jL+2\Lambda\quad(j=n,N),\qquad
 \|E_N''-E_n''\|\le L+\frac{2\Lambda}{n^2}.
\]
Consequently, for \(\varepsilon_n=L+2\Lambda/n\),
\[
 \sup_{|v|\le5}|E_j(v)|\le15j\varepsilon_n,\qquad
 \sup_{|v|\le5}|E_N(v)-E_n(v)|\le15\varepsilon_n .
\]
When \(K=0\) and \(n\ge1/\delta\), one may take
\(\varepsilon_n\le36 S_*\delta^2\). The identity \(K=0\) holds
for an exactly mirrored set of rows and also permits signed
cancellation without exact mirroring.

Let \(t=v-1/2\). On this disk
\(\Re x_i\ge (109/120)M_i\) and \(|t|\le11/2\), so the Gamma
logarithms use the branch analytic in the right half-plane.
Expanding the reciprocal to its linear term gives the exact identity
\[
 \sum_i\frac{d_i^2}{x_i(v)}-S_*=
    -tK+t^2\sum_i\frac{d_i^4}{M_i^3(1+td_i/M_i)}.
\]
Its modulus is at most \(L\), since
\((121/4)/(109/120)<34\).

For real \(j\ge1\), write \(h_j(x)=j^2\psi'(jx+1)\).
The [trigamma series](https://dlmf.nist.gov/5.15.E1)
and sum–integral comparison give
\[
 |h_j(x)-j/x|\le\frac1{(\Re x)^2}.
\]
Comparing the series at the right endpoints to the integral of
\((jx+t)^{-2}\) on successive unit intervals bounds their
difference by the integral of \(2/(j\Re x+t)^3\). Multiplication
by \(j^2\) yields the claim.

The derivative in \(j\) has a sharper sum–integral error because
its integrand vanishes at the left endpoint. In fact
\[
 \partial_j h_j(x)=\sum_{k\ge1}\frac{2jk}{(jx+k)^3},
 \qquad \int_0^\infty\frac{2jt}{(jx+t)^3}\,dt=\frac1x.
\]
Differentiation of the series is uniform on compact subsets of
\(\Re x>0\) and positive \(j\). For
\(f(t)=2jt/(jx+t)^3\), \(f(0)=0\) and \(f(t)\to0\). Thus the
sum of the composite trapezoid values on the unit intervals
equals the displayed series. The integral form of the
trapezoid error on a unit interval has kernel
\((t-k)(k+1-t)/2\), bounded by \(1/8\). Writing \(r=\Re x\), we
obtain
\[
 \left|\partial_jh_j(x)-\frac1x\right|
 \le\frac18\int_0^\infty
 2j\left(\frac6{(jr+t)^4}+\frac{12t}{(jr+t)^5}\right)dt
 =\frac3{4j^2r^3}.
\]

Differentiate \(E_j\) twice in \(v\) to get
\(E_j''=\sum_i d_i^2h_j(x_i)-jS_*\). The reciprocal identity and
the bound for \(h_j\) give the first inequality. Apply the bound
for \(\partial_jh_j\), then integrate from \(n\) to \(N\), for
the second. Here \(M_i\ge1\) gives
\(\sum_i d_i^2/M_i^3\le\Lambda\).
Both \(E_j\) and its adjacent difference vanish at \(v=0,1\).
For such an analytic function \(F\), integrating over \([0,1]\)
gives \(|F'(0)|\le\|F''\|/2\), then integration along \([0,v]\)
gives \(\|F\|\le15\|F''\|\) on the disk. This proves the log
bounds. Finally \(B_4\le S_*\delta^2\) and
\(\Lambda\le S_*\delta\), since the nonzero integer \(d_i\) have
\(|d_i|\ge1\). These inequalities give the signed-cancellation
specialization.
