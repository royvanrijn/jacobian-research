# A sharper all-terminal weak-curvature band

Use the notation and arbitrary positive Gamma shapes of the
[explicit weak-curvature proof](MULTIVARIATE_WEAK_CURVATURE_ALL_TERMINALS.md).
For \(N=n+1\ge3\), put
\[
   \sigma^2=\sum_i\frac{(b_i-a_i)^2}{\gamma_i+n(a_i+b_i)/2}.
\]
If \(\sigma\le(3\cdot10^7 N^{5/2})^{-1}\), the two moment
polynomials have simple roots and disjoint sets of root arguments.
In particular the ordinary pair is coprime whenever
\[
                     S\le 10^{-16} N^{-4}.                  \tag{1}
\]

## The Hermite gap

For \(\rho\in[n/N,1]\), compare the roots of \(\mathrm{He}_n\)
with those of \(\mathrm{He}_N\) scaled by \(\sqrt\rho\).
Their distance is at least
\[
                          \delta=\frac1{200 N^{3/2}}.       \tag{2}
\]
The [Hermite comparison in the original proof](MULTIVARIATE_WEAK_CURVATURE_ALL_TERMINALS.md#quantitative-comparison-of-the-hermite-pairs)
gives the same-set spacing \(b_0=1/(20\sqrt N)\), an unscaled
distance \(1/(10\sqrt N)\) from a lower-order root to either
adjacent larger-order root, and the following scaled comparison.
For each positive lower-order root \(x\), let \(y\) be its
unscaled positive outer neighbor. Then
\(\sqrt\rho\,y-x\ge x/(3N^2)\).
The inner positive neighbor only moves farther below \(x\) under
scaling, so its unscaled bound suffices. The zero and negative
neighbors have separation at least \(b_0/2\), as in the original
proof, and reflection covers the negative lower roots.

For the outer neighbor, if \(y\le\sqrt N/20\), then
\[
 (1-\sqrt\rho)y\le y/N\le1/(20\sqrt N).
\]
Subtract this from the unscaled gap to get \(1/(20\sqrt N)\).
If \(y>\sqrt N/20\) and \(x\ge\sqrt N/40\), the scaled comparison
gives \(1/(120N^{3/2})\). Finally if
\(y>\sqrt N/20\) and \(x<\sqrt N/40\), then \(\sqrt\rho>3/4\)
gives separation greater than \(\sqrt N/80\).
These three cases give (2).

## The relative error on the larger Hermite circles

Keep the same midpoint tilt, centered entire functions \(F_j\),
and reference Gaussian polynomials \(G_j\) of orders \(j=n,N\).
Take circles of radius \(\delta/4\) about each imaginary Hermite
root. The reciprocal-distance sum of the other Hermite roots is
still bounded by \(80\sqrt N(1+\log N)\), by same-set spacing.
Thus on those circles
\[
 \left|G_j'/G_j\right|\le
   800N^{3/2}+80\sqrt N(1+\log N)\le880N^{3/2}.
\]
With \(D_0=6\sqrt N\), this improves the derivative bound in
the original relative-error argument to
\[
            \left|G_j^{(k)}/G_j\right|
                  \le150N D_0^k\qquad(0\le k\le j).
\]
The \(k=0,1\) bounds are immediate. The induction uses the same
differentiated Hermite equation with normalized coefficients
bounded by \(3/4\) and \(1/24\).

Set \(R=2D_0=12\sqrt N\) as before. The original cumulant
bounds hold for the new hypothesis
\(\sigma\le(3\cdot10^7N^{5/2})^{-1}\): in particular
\(|h|\le\sigma N/2\) implies \(|z+h|\le2R\) and
\(\sigma|z+h|\le1/2\) for \(|z|\le R\), and
\(24\sigma R^3<1/2\). The same coefficient argument now gives
\[
  \left|F_j/G_j-1\right|
       \le48\sigma R^3(1+150N)
       \le48\cdot 12^3\cdot151\,\sigma N^{5/2}<1/2 .
\]
Rouché therefore gives one zero on the interior of each
circle's disk. The common midpoint scale \(\sigma\) and the
disjoint imaginary-coordinate ranges of the disks in (2) give
disjoint physical root arguments under the exponential change
of variable. The angle coordinates lie within one period,
as \(6\sigma\sqrt N<2\pi\).
For ordinary shapes \(\sigma^2\le W/n\le2S/n\).
Under (1) this is at most \(3\cdot10^{-16}N^{-5}\), which
implies the stated hypothesis on \(\sigma\).
