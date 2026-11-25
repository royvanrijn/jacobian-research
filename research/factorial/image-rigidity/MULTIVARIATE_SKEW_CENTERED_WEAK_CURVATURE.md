# Skew-centered weak-curvature angular separation

Let \(a,b\) be distinct nonnegative integer exponent vectors and let \(X_i\)
be independent unit-rate Gamma variables of positive shapes \(\gamma_i\).
Set \(U=X^a\), \(V=X^b\), \(d_i=b_i-a_i\), and
\(P_j(t)=\mathbb E(U+tV)^j\). For \(N=n+1\ge3\) put
\[
 \lambda_i(j)=\gamma_i+\tfrac j2(a_i+b_i),\qquad
 V_j=\sum_i\frac{d_i^2}{\lambda_i(j)},\qquad \sigma=\sqrt{V_n}.
\]

**Theorem.** If \(\sigma\le 10^{-5}N^{-5/2}\), the roots of \(P_n\) and
\(P_N\) have disjoint arguments, and both polynomials have simple roots.
In particular, for ordinary factorial moments the two polynomials
are coprime over every characteristic-zero field when
\[
       \sum_{i:d_i\ne0}\frac{d_i^2}{\max(a_i,b_i)}
                      \le\frac{1}{3\cdot 10^{10}N^4}.                 \tag{1}
\]
At \(n=1\) consecutive coprimality follows directly from strict
Cauchy–Schwarz for the two distinct monomials.

The proof uses the Hermite circles and midpoint tilt from the
[weak-curvature theorem](MULTIVARIATE_WEAK_CURVATURE_ALL_TERMINALS.md).
Absorbing the term linear in the cubic cumulant into a real displacement
of each Hermite center preserves the angular intervals of those circles.

## Midpoint representation

Under expectation \(\mathbb E_j\), use independent Gamma shapes
\(\lambda_i(j)\) and set
\[
 L_j=\sum_i d_i\log X_i,\quad
 \mu_j=\mathbb E_j L_j,\quad W_j=(L_j-\mu_j)/\sigma,\quad
 K_j(z)=\log\mathbb E_j e^{zW_j},\quad \rho=V_j/\sigma^2 .
\]
Here \(\rho\in[n/N,1]\subset[2/3,1]\). For
\(t=-\exp(-\mu_j+\sigma x)\) the zeros of \(P_j(t)\) are represented
by the zeros of
\[
 F_j(x)=\mathbb E_j\left[(x+W_j)^j
        \left(\frac{\sinh(\sigma(x+W_j)/2)}{\sigma(x+W_j)/2}\right)^j
                \right].
\]
In fact \((U+tV)^j\) factors into a nonvanishing scalar exponential,
the Gamma tilt \((UV)^{j/2}\), and the displayed integrand. Introduce
\(h=\sigma\sum_{\ell=1}^jT_\ell\), where the \(T_\ell\) are independent
uniform variables on \([-1/2,1/2]\), independent of \(W_j\). The uniform
moment generating function gives
\[
 F_j(x)=j![z^j]\mathbb E_h\exp\big((z+h)x+K_j(z+h)\big),
 \qquad \mathbb E h=0,\quad |h|\le\sigma N/2.                          \tag{2}
\]

Fix the order \(j\). For \(K_j\), take the analytic branch at zero. Its
quadratic coefficient differs from \(\rho/2\) by at most \(\sigma^2/2\),
and its Taylor coefficients of order \(k\ge3\) have absolute value at
most \((2/k)\sigma^{k-2}\). Indeed, the Gamma polygamma series gives
\(1/\lambda\le\psi'(\lambda)\le1/\lambda+1/\lambda^2\). Also
\(\lambda_i(j)^{-1}\le\sigma^2\) and
\(|d_i|/(\lambda_i(j)\sigma)\le\sigma\) for each \(d_i\ne0\).
For \(k\ge3\), differentiating the trigamma series gives the bound
\[
 |[z^k]K_j(z)|\le\frac{\sigma^{k-2}}{k}
              \sum_i\frac{d_i^2}{\sigma^2}\psi'(\lambda_i(j))
 \le\frac{1+\sigma^2}{k}\sigma^{k-2}.
\]
Writing \(c=[z^3]K_j(z)\), these estimates give
\[
 |c|\le 2\sigma/3,\qquad
 K_j(z)=\rho z^2/2+c z^3+E(z),\qquad
 |E(z)|\le\tfrac12\sigma^2|z|^2+\sigma^2|z|^4
 \quad(\sigma|z|\le1/2).                                               \tag{3}
\]
The last bound sums the geometric majorant with \(2/k\le1/2\) for
\(k\ge4\). The disk \(\sigma|z|\le1/2\) lies inside the nonvanishing
Gamma product's shape strip. [DLMF, polygamma series](https://dlmf.nist.gov/5.15.E1)

## Cancellation on each Hermite disk

Let \(G_j(x)=j![z^j]\exp(xz+\rho z^2/2)\), whose roots are imaginary
rotations of the scaled probabilists' Hermite roots. The two
orders have anchors separated by at least \(1/(200N^{5/2})\).
Set \(r=1/(800N^{5/2})\) and consider \(|w-z_0|=r\) about an
anchor \(z_0\) of \(G_j\). The Hermite estimates give
\[
 |z_0|\le2\sqrt N,\qquad |w|\le3\sqrt N,\qquad
 \left|\frac{G_j^{(k)}(w)}{G_j(w)}\right|
                 \le 150N^2(6\sqrt N)^k\quad (1\le k\le j),
 \quad
 \left|\frac{G'_j(w)}{G_j(w)}\right|\le\frac2r.                         \tag{4}
\]
The last bound uses the root-circle logarithmic derivative estimate
\(1/r+80\sqrt N(1+\log N)\le2/r\).
All bounds follow from the
[Hermite comparison and circle estimates](MULTIVARIATE_WEAK_CURVATURE_ALL_TERMINALS.md#quantitative-comparison-of-the-hermite-pairs).

The equation \(\rho G_j''+wG_j'-jG_j=0\) and its derivative give
\[
 G_j'''(w)=A(w)G_j'(w)+B(w)G_j(w),\qquad
 A(w)=\frac{j-1}{\rho}+\frac{w^2}{\rho^2},\qquad
 B(w)=-\frac{jw}{\rho^2}.
\]
Choose \(s=-c A(z_0)\). Since \(c\) is real and \(z_0\) is imaginary,
\(s\) is real. From (3), (4), and \(\rho\ge2/3\),
\(|s|\le(2\sigma/3)(3N/2+9N)=7\sigma N\).

Set \(R=12\sqrt N\). At every point \(|z|\le R\) in (2) with
\(x=w+s\), write the exponent as \(wz+\rho z^2/2+T\), and put
\[
 T_0=sz+h(w+\rho z)+c z^3.
\]
Then \(\lvert T_0\rvert\le\sigma R^3\). To see this, bound the three
terms by \(7\sigma N R\), \(\sigma N R\), and \(2\sigma R^3/3\).
Since \(|h|\le R\) and \(\sigma |z+h|\le2\sigma R\le1/2\), (3) gives
\[
\begin{aligned}
 |T-T_0|
 &\le |hs|+\tfrac12 \rho h^2
          +|c|\,|(z+h)^3-z^3|+|E(z+h)|\\
 &\le \sigma^2\left(\tfrac{29}{8}N^2+\tfrac73 NR^2
                                 +2R^2+16R^4\right)
 \le 20\sigma^2R^4 .
\end{aligned}
\]
Here the cubic difference is at most \(7R^2|h|\) by factoring
the difference of cubes. The hypothesis gives
\(\sigma R^3\le1/4\) and \(20\sigma R\le1\). Thus
\(|T|\le2\sigma R^3\le1/2\). From the exponential Taylor remainder,
\[
 |e^T-1-T_0|
 \le 20\sigma^2R^4+\tfrac12e^{1/2}(2\sigma R^3)^2
 \le 5\sigma^2R^6,
\]
using \(e^{1/2}<2\) and \(R^2\ge432\).

Averaging gives \(\mathbb E T_0=sz+cz^3\). The remainder is analytic in
\(z\) on a neighborhood of \(|z|\le R\), uniformly for the bounded set of
\(h\), so Cauchy's coefficient estimate applies to its expectation. On
extracting the degree-\(j\) coefficient in (2) and using (4), this yields
\[
 \left|\frac{F_j(w+s)-G_j(w)-sG_j'(w)-cG_j'''(w)}{G_j(w)}\right|
 \le 5\sigma^2 R^6(1+150N^2)<\frac14 .                                 \tag{5}
\]
Indeed each coefficient of degree \(k\) in the remainder is bounded by
\(5\sigma^2R^6/R^k\), and \((6\sqrt N)/R=1/2\). For the last strict
inequality use \(5\cdot12^6\cdot151<10^{10}/4\).

The linear correction remaining in (5) has ratio to \(G_j(w)\) equal to
\[
       c\left((A(w)-A(z_0))\,\frac{G_j'(w)}{G_j(w)}+B(w)\right).
\]
The two terms in parentheses are bounded by
\((45/2)\sqrt N\) and \((27/4)N^{3/2}\), respectively. The first bound
uses \(|w-z_0|=r\) and the last estimate in (4). Multiplication
by \(|c|\le2\sigma/3\) gives a bound of \(10\sigma N^{3/2}<1/4\) for
the sum. Rouché's theorem now places exactly one zero, counted
with multiplicity, of \(F_j(w+s)\) on the open disk \(|w-z_0|<r\).

In \(x=w+s\) the disk has center \(z_0+s\) and its imaginary interval
is still \((\Im z_0-r,\Im z_0+r)\). All these intervals are disjoint
within and across the two orders because the Hermite anchors are
separated. They give all the roots of \(P_j\) through the exponential
substitution, which has nonzero derivative and is injective on the
strip \(|\Im x|\le3\sqrt N\) since \(6\sigma\sqrt N<2\pi\). Both orders
use the same \(\sigma\). Their root arguments are therefore disjoint.

For ordinary factorial moments, put
\(S=\sum_{i:d_i\ne0}d_i^2/\max(a_i,b_i)\). Since
\(\lambda_i(n)\ge(n/2)\max(a_i,b_i)\) and \(n\ge2N/3\),
\(\sigma^2\le2S/n\le3S/N\). Condition (1) implies the hypothesis.
The resulting rational Bézout identity gives coprimality over
every characteristic-zero field. At \(n=1\), the only root of
\(P_1\) is \(-\mathbb E U/\mathbb E V\), while
\(P_2(t)=\mathbb E(U+tV)^2>0\) on the real line because the two
monomials are not proportional under Gamma measure. \(\square\)

