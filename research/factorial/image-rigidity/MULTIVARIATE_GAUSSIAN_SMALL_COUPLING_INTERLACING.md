# An explicit small-coupling Gaussian reference comparison

For \(a>0\) put
\[
 B_j(z;2a)=\sum_{k=0}^j\binom jk e^{-(a/j) k(j-k)}z^k .
\]
Let \(n\ge2\) and \(N=n+1\). These consecutive polynomials have
simple and strictly interlacing unit-circle roots if
\[
                         0<a\le 10^{-4}N^{-3/2}.             \tag{1}
\]
The two heat times are \(a/n\) and \(a/N\), the slowest physical
upper time from the
[Gaussian reference comparison](MULTIVARIATE_GAUSSIAN_REFERENCE_ANGULAR_INTERLACING.md).
That comparison's outward monotonicity and inner inequalities
extend the conclusion to every upper time in \([a/N,a/n]\)
with the lower time fixed at \(a/n\).

The proof uses an exact Gaussian comparison for the sum of
uniform variables in the sine representation. Its variance
does not change with \(j\). Absorbing this common variance
gives a common rescaling of the Hermite anchors and an error
quadratic in \(a\).

## An auxiliary variable with degree-independent variance

For \(j=n,N\) put \(\lambda_j=\sqrt{2a/j}\).
In the signed angular coordinate \(z=-e^{ix}\), the real centered
function for \(B_j(z;2a)\) is
\[
 f_j(x)=e^{-(a/j)\partial_x^2}\sin^j(x/2).
\]
This follows by applying the Fourier multiplier term by term
to the finite trigonometric expansion. If \(Z\) is standard
normal, that multiplier is also given by expectation of the
shift \(x+i\lambda_j Z\). Thus, with \(v=x/\lambda_j\) and
\(D_{\sin}(w)=\sin(w/2)/(w/2)\),
\[
 F_j(v):=(2/\lambda_j)^j f_j(\lambda_j v)
     =\mathbb E\left[(v+iZ)^j D_{\sin}(\lambda_j(v+iZ))^j\right].
\]
Let \(T_1,\ldots,T_j\) be independent uniform variables on
\([-1/2,1/2]\), independent of \(Z\), and let
\(h=\lambda_j\sum_{\ell=1}^j T_\ell\).
The uniform characteristic function gives the displayed sine
factor. Tilt the Gaussian \(Z\) by \(e^{-hZ}\) and extract a
coefficient at power \(j\). This gives
\[
 F_j(v)=j![z^j]e^{vz-z^2/2}\,
                    \mathbb E_h e^{ih(v-z)+h^2/2}.            \tag{2}
\]
The extraction variable in (2) is independent of the original binomial
variable. In particular the shifted Hermite in (2) comes from
\(\mathbb E_Z e^{-hZ}(v+iZ)^j
=e^{h^2/2}\mathrm{He}_j(v-ih)\).

Put \(c=a/6\) and
\[
 D(w)=\frac{\sinh(w/2)}{w/2},\qquad
 M_h(\zeta)=D(\lambda_j\zeta)^j,\qquad
 A(\zeta)=M_h(\zeta)\exp(-c\zeta^2/2).
\]
The variance of \(h\) is \(\lambda_j^2 j/12=c\).
For complex \(\zeta\) the uniform-variable bound gives
\[
 |M_h(\zeta)|\le M_h(\Re\zeta)\le
                     \exp(c(\Re\zeta)^2/2).
\]
Indeed \(D(w)\le e^{w^2/24}\) on the real line by comparing
coefficients with \((2k+1)!\ge6^k k!\).
Consequently \(|A(\zeta)|\le\exp(c(\Im\zeta)^2/2)\).

In (2) write \(w=i(v-z)\). An independent standard normal gives
\(e^{h^2/2}=\mathbb E_G e^{Gh}\). Complete the square in
\(\mathbb E_G M_h(G+w)\), and move the contour to
\(G=r+cw/(1-c)\) with \(r\in\mathbb R\). The integrand is entire,
and the above subGaussian bound gives Gaussian decay uniformly
on the bounded-height contour pieces at infinity. The result is
\[
 \mathbb E_h e^{hw+h^2/2}
  = (1-c)^{-1/2} e^{c w^2/(2(1-c))}\,(1+\eta_j(w)),\qquad
 1+\eta_j(w)=\mathbb E_r A\!\left(r+\frac{w}{1-c}\right),     \tag{3}
\]
where \(r\) is real normal of variance \(1/(1-c)\).

Set \(v=\sqrt{1-c}\,y\), \(z=\sqrt{1-c}\,\xi\) in the
coefficient extraction. Combining the three quadratic terms
in (2)–(3) gives the exact expression
\[
 F_j(\sqrt{1-c}\,y)=
   (1-c)^{-(j+1)/2}e^{-cy^2/2}\,
   j![\xi^j]e^{y\xi-\xi^2/2}
          \bigl(1+\eta_j(i\sqrt{1-c}(y-\xi))\bigr).             \tag{4}
\]
With \(\eta_j=0\), the nonvanishing prefactor multiplies
\(\mathrm{He}_j(y)\).

## Uniform error on the Hermite circles

Use \(|y|\le3\sqrt N\) and \(|\xi|\le12\sqrt N\).
For the argument of \(A\) in (3) this gives
\[
 \zeta=r+i(y-\xi)/\sqrt{1-c},\qquad
 |\zeta|\le28\sqrt N\quad\text{if } |r|\le12\sqrt N.
\]
Here \(c\le1/16\). On the latter interval
\(|\lambda_j\zeta|^2\le2352a<1/4\) by (1) and \(N/j\le3/2\).
For \(|u|\le1/2\) the series for \(D\) gives
\[
                  |\log D(u)-u^2/24|\le |u|^4/100.
\]
For example, \(|D(u)-1-u^2/24|\le |u|^4/1800\),
\(|D(u)-1|\le |u|^2/20\), and the logarithm's higher terms
are bounded by \(|D(u)-1|^2/[2(1-|D(u)-1|)]\).
It follows that
\[
 |j\log D(\lambda_j\zeta)-c\zeta^2/2|
 \le 4a^2|\zeta|^4/(100j)<40000a^2N.
\]
The exponential remainder then gives \(|A(\zeta)-1|\le80000a^2N\)
on this interval. For \(|r|>12\sqrt N\), the bound
\(|\Im\zeta|\le16\sqrt N\) and \(aN\le10^{-4}\) give
\(|A(\zeta)|<2\). A Gaussian exponential-moment bound gives
\(\Pr(|r|>12\sqrt N)\le2e^{-50N}\).
Together these give
\[
 \sup_{\lvert\xi\rvert\le12\sqrt N}
   |\eta_j(i\sqrt{1-c}(y-\xi))|
            \le 80000a^2N+6e^{-50N}.                         \tag{5}
\]

On each circle of radius \(r_0=1/(800N^{5/2})\) around a
root of \(\mathrm{He}_j\), its derivatives satisfy
\[
  \left|\mathrm{He}_j^{(k)}(y)/\mathrm{He}_j(y)\right|
          \le150N^2(6\sqrt N)^k\qquad(1\le k\le j).
\]
These are the unit-variance
[Hermite derivative estimates](MULTIVARIATE_WEAK_CURVATURE_ALL_TERMINALS.md#a-relative-error-bound-on-hermite-circles),
rotated to the real line. Applying Cauchy's coefficient bound
to (5), the ratio of the correction in (4) to its Hermite term
has absolute value at most
\[
 (1+150N^2)(80000a^2N+6e^{-50N})<\tfrac12 .
\]
The last inequality uses (1) and \(\sum_{k\ge1}2^{-k}=1\)
for the Cauchy derivative weights. Rouché gives exactly one
simple zero in each circle. The function has real symmetry
in \(y\), so these zeros are real.

The common normalization \(x/(\lambda_n\sqrt{1-c})\)
takes the anchors to the roots of \(\mathrm{He}_n\) and
the roots of \(\mathrm{He}_N\) scaled by \(\sqrt{n/N}\).
Their consecutive gaps are at least \(1/(200N^{5/2})\)
by the
[quantitative Hermite comparison](MULTIVARIATE_WEAK_CURVATURE_ALL_TERMINALS.md#quantitative-comparison-of-the-hermite-pairs).
The rescaled disks have radius at most \(r_0\).
Thus both the relative angular order and the simplicity
persist for the exact roots. The corresponding \(x\) all
lie in \((-\pi,\pi)\), since
\(|x|\le3\lambda_j\sqrt{1-c}\sqrt N\le3\sqrt{3a}<\pi\).
The map \(x\mapsto-e^{ix}\) therefore produces all the
distinct roots of the stated binomial polynomials.

