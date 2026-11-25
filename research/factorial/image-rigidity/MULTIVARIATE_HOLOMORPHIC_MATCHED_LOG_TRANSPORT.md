# Holomorphic transport of consecutive normalized factorial logs

Let \(n\ge2\), \(N=n+1\), and let \(a,b\) be distinct ordinary
nonnegative exponent vectors. Omit equal rows and put
\[
 d_i=b_i-a_i,\qquad L_i=\max(a_i,b_i),\qquad
 \eta=\max_i |d_i|/L_i,\qquad x_i(s)=a_i+d_is .
\]
Let \(U,W\) be the lower and upper reflected outer losses in the
[matched tanh transport theorem](MULTIVARIATE_TANH_EXCESS_LOSS_TRANSPORT.md).
Write \(R=n^2/(n^2-1)\), \(\Delta=W-RU\), and \(\nu=1/(3N)\).
On the convex complex domain \(\Re x_i(s)>-\nu\) for every row, define
\[
 A_j(s)=\frac1j\sum_i\left[
   \log\Gamma(1+jx_i(s))-(1-s)\log\Gamma(1+ja_i)
                       -s\log\Gamma(1+jb_i)\right],
 \qquad j=n,N .
\]
Use the Gamma logarithm real on the positive axis. At each own knot
\(s=k/j\), \(jA_j(s)\) is the exact log of the endpoint-normalized
factorial multiplier after removal of \(\binom jk\).

Let \(D=A_N-A_n\) and \(B(s)=-D(s)/(s(1-s))\), with removable
endpoint values. For \(\tau>0\) assume
\(\Re x_i(v)+\nu\ge\tau L_i\) at \(v=0,1,s\) in every row.
For each integer \(m\ge0\),
\[
 \big|B^{(m)}(s)\big|
 \le \frac{2N}{n}\frac{m!}{m+2}
                 \frac{\eta^m}{\tau^{m+2}}\Delta.                 \tag{1}
\]
The estimate also holds along the convex hull of \(0,1,s\).
It covers private rows, with a complex neighborhood of width
limited by \(\nu/|d_i|\) at a private endpoint. At real points of
the domain \(B\) is positive and strictly log convex. Its even
derivatives are strictly positive.

If all rows are shared and \(\Re x_i(v)\ge\tau_0 L_i\) at
\(v=0,1,s\), with \(\tau_0>0\), the same positivity and strict
log-convexity statements hold for \(B_0-B\) on this positive
real domain, where
\[
 B_0(s)=\frac{\sum_i\big[\log x_i(s)-(1-s)\log a_i-s\log b_i\big]}
                   {2nN\,s(1-s)} .
\]
The analytic difference satisfies
\[
 \big|B_0^{(m)}(s)-B^{(m)}(s)\big|
 \le \frac{m!}{12}\left(\frac1{n^2}-\frac1{N^2}\right)
       \sum_i\frac{|d_i|^{m+2}}{(\tau_0 L_i)^{m+3}} .                \tag{2}
\]
Thus the excess outer loss controls a holomorphic symbol whose samples
come from the actual factorial logs. The common scale in (1) includes
every derivative, with extra factors of the relative row gap.

## The exact curvature kernel

Write \(h(t)=t/(e^t-1)\), with \(h(0)=1\). The trigamma integral,
or termwise integration of the trigamma series, gives
\[
 j\,\psi'(1+jx)=\int_0^\infty e^{-xu}h(u/j)\,du
 \qquad (\Re x>-1/j).
\]
Indeed the integral before the change of variable \(u=jt\) is
\(\int_0^\infty e^{-jxt}t/(e^t-1)\,dt
=\sum_{q\ge1}(jx+q)^{-2}\), with absolute convergence.
Consequently
\[
 D''(s)=\sum_i d_i^2\int_0^\infty e^{-x_i(s)u}q_n(u)\,du,
 \qquad q_n(u)=h(u/N)-h(u/n).                                    \tag{3}
\]
For \(t>0\) direct differentiation gives \(h'(t)<0\), since
\(e^t(1-t)\) strictly decreases from one. Moreover
\[
 0<h'(t)+\tfrac12\le t/6 .
\]
To check the latter inequalities, put \(v=t/2\) and differentiate
\(h(t)+t/2=v\coth v\). The derivative is
\((\sinh(2v)-2v)/(4\sinh^2 v)>0\). The upper bound reduces to
\(3\sinh t\le t\cosh t+2t\), whose power series difference has
zero linear and cubic terms and positive higher odd terms.
It follows that, for \(u>0\),
\[
 0<q_n(u)<\frac{u}{2nN},\qquad
 0<\frac{u}{2nN}-q_n(u)
              \le \frac{u^2}{12}(n^{-2}-N^{-2}).                 \tag{4}
\]
This uses integration on \([u/N,u/n]\). For a bound with
exponential decay, the
[fractional-index kernel](MULTIVARIATE_FRACTIONAL_INDEX_TRANSPORT.md#the-comparison-kernel)
gives \(-h'(t)\le e^{-t/3}/2\), and consequently
\(0<q_n(u)\le e^{-\nu u}u/(2nN)\).

Let \(\mathcal T=\{(\alpha,\beta):\alpha,\beta\ge0,\
\alpha+\beta\le1\}\). Integrating first in \(\beta\) and then
in \(\alpha\) gives the chord identity
\[
 B(s)=\int_{\mathcal T} D''(\alpha s+\beta)\,d\alpha\,d\beta .
\]
It holds off the endpoints and then by analytic continuation at them.
The formula for \(B_0\) is the same with \(q_n(u)\) in (3) replaced
by \(u/(2nN)\), since the corresponding second derivative is
\(\sum_i d_i^2/(2nN x_i^2)\).

Differentiate \(m\) times under the integrals, which is valid on
compact subsets of the domain by exponential domination. The
resulting factors are \((-\alpha d_i u)^m\). Since
\(\int_{\mathcal T}\alpha^m\,d\alpha\,d\beta
=1/((m+1)(m+2))\), the exponential bound for \(q_n\) gives
\[
 |B^{(m)}(s)|
 \le \frac{m!}{2nN(m+2)}
                   \sum_i\frac{|d_i|^{m+2}}{(\tau L_i)^{m+2}} .
\]
The rowwise gain \(G_{+,i}\) in the
[integer-grid excess bound](MULTIVARIATE_TANH_EXCESS_LOSS_TRANSPORT.md#two-contributions-to-the-row-gain)
gives \(\sum_i d_i^2/L_i^2\le4N^2\Delta\), proving (1).
The second bound in (4) gives (2) by integrating \(u^{m+2}\).
For real \(s\), both \(B\) and \(B_0-B\) are positive mixtures of
real exponentials in \(s\). Their even derivatives are positive
by the integral formulas. For either function, its second
logarithmic derivative is the variance of \(-\alpha d_i u\) under
the corresponding integral weights normalized by the function.
This variance is strictly positive because at least one row
has nonzero \(d_i\) and the associated kernel is positive on
\(\mathcal T^\circ\times(0,\infty)\).

## A paired multiplier on the same degree

The normalized moment polynomial is
\[
 F_j(z)=\sum_{k=0}^j \binom jk e^{jA_j(k/j)}z^k .
\]
Define the two degree-\(n\) polynomials
\(L_0=F_N-zF_N'/N\) and \(L_1=F_N'/N\), so
\(F_N=L_0+zL_1\). They are obtained from the exact \(F_n\) by
multiplying its coefficient at \(k\) by \(e^{g_\epsilon(k/n)}\),
where
\[
 g_\epsilon(s)=N A_N((ns+\epsilon)/N)-nA_n(s),\qquad \epsilon=0,1 .
\]
For \(\epsilon=0\) this follows from
\((1-k/N)\binom Nk=\binom nk\), and for \(\epsilon=1\) it follows
from \((k+1)\binom N{k+1}/N=\binom nk\).

Set \(r_\epsilon=\epsilon-s\). Taylor's integral identity gives
\[
 \begin{split}
 g_\epsilon(s)={}& A_n(s)+r_\epsilon A_n'(s)
           -N s(1-s)B(s)+r_\epsilon D'(s)\\
 &+\frac{r_\epsilon^2}{N}
       \int_0^1 (1-t)A_N''(s+t r_\epsilon/N)\,dt .
 \end{split}                                                     \tag{5}
\]
The shift stays in the convex hull used in (1). Equation (3) with
\(h(u/N)\) in place of \(q_n(u)\) bounds \(|A_N''|\) there by
\(\tau^{-1}\sum_i d_i^2/L_i\), since \(h(t)=e^{-t/2}t/(2\sinh(t/2))
\le e^{-t/2}\) and \(1/(2N)>\nu\). Equation (1) controls \(D'\) through
\(D=-s(1-s)B\). The principal terms of the two multipliers at the same
fractional position therefore share \(A_n(s)+(1/2-s)A_n'(s)-Ns(1-s)B(s)\),
and differ by the two signs of \(A_n'(s)/2\).

Equations (1) and (5) supply holomorphic bounds for exact grouped
multipliers at both orders. The [alternating binomial estimate](MULTIVARIATE_ANALYTIC_FRACTIONAL_BINOMIAL_TRANSPORT.md)
uses the exact analytic chords to preserve cancellation through
exponentiation and finite differences. For the angular comparison
the remaining estimate must be relative to the reference root
derivative or contour modulus. The [matched angular models](MULTIVARIATE_MATCHED_SECOND_ORDER_ANGLES.md)
give one quantitative target for that paired remainder in the
weak-curvature window.

## A gap-dependent angular derivative for the exact difference

Let \(G(z)=\sum_{k=0}^j g_k z^k\), \(j\ge2\), be a real reciprocal
reference with simple roots \(e^{i\theta_\ell}\) on the unit circle.
Write \(g>0\) for its minimum circular angular spacing. Fix one root at
angle \(\theta\) and define
\[
 c_\ell=\cot((\theta-\theta_\ell)/2)\quad(\theta_\ell\ne\theta),\qquad
 K=\tfrac12\sum_\ell c_\ell,\quad Q=\sum_\ell |c_\ell|^3,\quad
 H_i=\frac{a_i+b_i}{2}+\nu-\left(\frac12+\frac{\pi^2}{2jg}\right)|d_i|.
\]
Suppose \(H_i>0\) for every unequal row. For the analytic root
branch of
\[
                G_t(z)=\sum_{k=0}^j g_k e^{tjD(k/j)}z^k
\]
starting at \(e^{i\theta}\), the derivative of its argument satisfies
\[
 |\dot\theta(0)|\le
       \frac{|K|}{2nNj}\sum_i\frac{d_i^2}{H_i^2}
       +\frac{Q}{8nNj^3}\sum_i\frac{|d_i|^4}{H_i^4}.               \tag{6}
\]
In particular \(j=N\) compares the lifted lower chord and exact upper
chord at the level of their first variation whenever the reference
satisfies the circle and spacing hypotheses.

Let \(f(v)=e^{-ijv/2}G(e^{iv})\), which is real on the real axis.
The unit-circle factorization gives, for real \(x\),
\[
 \frac{f(\theta+ix)}{f'(\theta)}
  =2i\sinh(x/2)\cosh(x/2)^{j-1}
       \prod_{\ell}\bigl(1+i c_\ell\tanh(x/2)\bigr).
\]
In either direction from \(\theta\), the \(k\)-th root within
circular distance \(\pi\) has distance at least \(kg\). Assign
a possible antipodal root to just one direction. Since
\(\sin(\delta/2)\ge\delta/\pi\) for \(0\le\delta\le\pi\),
its \(|c_\ell|\) is at most \(\pi/(kg)\).
Set \(p=\pi|x|/(2g)\). Using \(|\tanh(x/2)|\le|x|/2\) and
bounding both lists of distances, the modulus of the product
is at most \(\prod_{k\ge1}(1+p^2/k^2)\le e^{\pi p}\).
For the last bound sum the decreasing function
\(\log(1+p^2/u^2)\) at the positive integers and compare with its
integral from zero to infinity. The integral is \(\pi p\) by
integration by parts. The other factors have combined modulus
at most \(|x|e^{j|x|/2}\).

The real part of this ratio is bounded in modulus by its modulus
times the absolute sine of the sum of the product arguments.
With \(w=\tanh(x/2)\), the inequality
\(|\arctan z-z|\le |z|^3/3\) for real \(z\) gives
\[
 \left|\sum_\ell\arctan(c_\ell w)\right|
            \le 2|K|\,|w|+Q|w|^3/3 .
\]
Thus
\[
 \left|\Re\frac{f(\theta+ix)}{f'(\theta)}\right|
 \le e^{(j/2+\pi^2/(2g))|x|}
               \bigl(|K|x^2+Q|x|^4/24\bigr).                    \tag{7}
\]

Write \(s=1/2+v\). Integrate the kernel (3) twice to express
\(D(s)-D(1/2)-vD'(1/2)\) as
\[
 \sum_i\int_0^\infty
     e^{-(a_i+b_i)u/2}q_n(u)\,
                   \frac{e^{-d_iuv}-1+d_iuv}{u^2}\,du .
\]
On the Fourier spectrum of \(f\), \(v=(k-j/2)/j\).
Hence the exponential term shifts \(\theta\) to
\(\theta+i d_i u/j\). The constant term vanishes at the root,
and the linear term contributes a purely imaginary ratio to
\(f'(\theta)\). Implicit differentiation gives the derivative
of the root argument as the real part of minus the perturbed
centered value divided by \(f'(\theta)\). Combining (7) with
\(q_n(u)\le e^{-\nu u}u/(2nN)\) and integrating \(u e^{-H_i u}\)
and \(u^3 e^{-H_i u}\) proves (6).

This comparison preserves a common root-product denominator.
To use it for a finite angular displacement, the spacing and
location estimates must hold along the relevant deformation.
Separation from the lower degree also requires a comparison with
the upper-grid lift before transporting its multiplier.
