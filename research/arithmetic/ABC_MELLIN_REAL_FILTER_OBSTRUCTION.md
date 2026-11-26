# A sharp Mellin cost for filtering close sums

The real filter for modular target popularity restricts \(a+b-c\) to a short interval. On any fixed comparable range of \(a/c\) and \(b/c\), imposing that filter by an absolutely summable superposition of separated multiplicative characters has a sharp power cost. A dual curved measure proves the lower bound for any such filter that separates the exact relation from a nearby displaced curve.

[Real-filtered popularity tests](ABC_SHORT_INTERVAL_RADICAL_FILTER.md#consequences-for-real-filtered-popularity-tests) · [Signed modular modes](ABC_MODULAR_POPULARITY_NONPRINCIPAL_MODES.md#unit-tests-at-the-smaller-real-filtered-prime-scale)

## The coefficient norm

Let \(I\) be a compact interval of positive length inside \((-\infty,0)\), let \(g(x)=\log(1-\exp x)\), and fix \(0<\delta<1/4\). Suppose \(\sigma\) is a finite complex Borel measure on \(\mathbb R^2\), and set
\[
F_\sigma(x,y)=\int_{\mathbb R^2}
                e^{2\pi i(\xi x+\eta y)}\,d\sigma(\xi,\eta).
\]
If
\[
\operatorname{Re}F_\sigma(x,g(x))\ge 1,\qquad
\operatorname{Re}F_\sigma(x,g(x)+2\delta)\le \tfrac12
\quad (x\in I),
\]
then
\[
                  \|\sigma\|_{\mathrm{TV}}\gg_I\delta^{-1/2}.
                                                        \tag{MF1}
\]
The exponent is sharp, even for a smooth cutoff which equals one on the full original arc and vanishes whenever \(|y-g(x)|\ge\delta\) in a neighborhood of that arc.

Choose a fixed nonnegative smooth bump \(\chi\) supported in the interior of \(I\) with integral one. Push \(\chi(x)\,dx\) forward to the graph of \(g\) to obtain \(\mu\). Its Fourier transform satisfies
\[
|\widehat\mu(\xi,\eta)|\ll_I (1+|\xi|+|\eta|)^{-1/2}.
\]
Indeed, the first frequency range is \(|\xi|\le 2(1+\sup_{\operatorname{supp}\chi}|g'|)|\eta|\). Here \(g''=-\exp(x)/(1-\exp(x))^2\) is bounded away from zero on the support. For \(|\eta|>1\), the set where the phase derivative is at most \(|\eta|^{1/2}\) in magnitude has length \(O_I(|\eta|^{-1/2})\). On its complement, integrate by parts and use monotonicity of the derivative to bound the variation of its reciprocal by \(O(|\eta|^{-1/2})\). In the complementary range the phase derivative is bounded below in magnitude by \(|\xi|/2\). One integration by parts gives \(O_I((1+|\xi|)^{-1})\), using the bounded second derivative and \(|\eta|\ll_I |\xi|\). Near the origin use the mass bound. The constants can absorb the fixed choice of bump.

Subtract from \(\mu\) its translate by \((0,2\delta)\) to obtain \(\nu\). Fourier transformation gives
\[
|\widehat\nu(\xi,\eta)|
 \ll_I \min(1,\delta|\eta|)(1+|\xi|+|\eta|)^{-1/2}
 \ll_I \delta^{1/2}.
\]
For the last inequality split at \(1+|\xi|+|\eta|=\delta^{-1}\). On the other hand \(\operatorname{Re}\int F_\sigma\,d\nu\ge1/2\). Fubini gives \(\int F_\sigma\,d\nu=\int\widehat\nu(-\xi,-\eta)\,d\sigma(\xi,\eta)\). The sup norm bound proves (MF1), including when the frequencies are discrete.

For the upper bound, extend \(g\) smoothly outside a larger compact interval inside \((-\infty,0)\). Let \(w\) be smooth with compact support there and equal to one near \(I\). Take \(\rho\) smooth, \(0\le\rho\le1\), supported on \((-1,1)\) and equal to one near zero. Then
\[
                  K_\delta(x,y)=w(x)\rho((y-g(x))/\delta)
\]
satisfies both requirements. The support has area \(O_I(\delta)\), and the \(L^2\) norm of each derivative of order \(s\) is \(O_{I,s}(\delta^{1/2-s})\) by the chain rule. Plancherel and Cauchy–Schwarz give \(\int_{|\zeta|\le\delta^{-1}}|\widehat K_\delta(\zeta)|\,d\zeta\ll_I\delta^{-1/2}\). For the annulus \(R\le|\zeta|\le2R\) with \(R\ge\delta^{-1}\), the same inequalities applied with derivatives of order two give
\[
\int_{R\le|\zeta|\le2R}|\widehat K_\delta(\zeta)|\,d\zeta
                      \ll_I R^{-1}\delta^{-3/2}.
\]
Dyadic summation gives \(\|\widehat K_\delta\|_1\ll_I\delta^{-1/2}\), so Fourier inversion supplies a measure \(\sigma\) of the matching norm.

This square-root scale is consistent with the classical growth of absolute Fourier norms for a smooth curved phase. The dual witness here uses two curves to make a statement for any filter with the specified contrast. [Lebedev, *Quantitative estimates in Beurling–Helson type theorems*, introduction](https://arxiv.org/pdf/1112.5677)

## Consequence for the modular estimate

Use \(x=\log(a/c)\), \(y=\log(b/c)\), and \(X/2<c\le X\). On the fixed comparable range, a displacement of order \(\delta\) in \(y\) changes \(a+b-c\) by order \(X\delta\), since \(\exp(g(x))=1-\exp x\) is bounded below on \(I\). For a filter of absolute width \(X^\theta\), choose \(\delta\) a sufficiently large fixed multiple of \(X^{-1+\theta}\) so the displaced arc lies outside the pass band. In the scale \(N=X^L\), the unavoidable coefficient norm then has order at least
\[
                   X^{(1-\theta)/2}=N^{(1-\theta)/(2L)}
\]
up to a fixed constant. Each separated frequency does preserve multiplicative weights since its phase factors into a power of \(a\), a power of \(b\), and a power of \(c\). An argument applying a uniform saving \(N^{-\kappa}\) at each such frequency and then integrating absolute values therefore has a coefficient loss larger than its saving whenever \(L\) is sufficiently small and fixed with \(\kappa\) fixed. The height limit is taken at that budget. The lower bound tests a continuum of comparable ratios. Control on only the actual radical-restricted sample, or signed cancellation across the frequencies, would require a different estimate. Multilinear finite field estimates for product dictionaries address the unfiltered frequencies. [Bourgain, *On Exponential Sums in Finite Fields*, Theorem A](https://www.math.ias.edu/files/avi/Bourgain_expon2.pdf)

The additive smoothing \(\rho((a+b-c)/X^\theta)\) itself has an additive Fourier representation with bounded total variation by rescaling \(\widehat\rho\). For a factored letter, its additive phase depends on the product of the slots and thus does not provide the separated slot weights used by the multiplicative group moment bound. [Logarithmic group moments](ABC_REAL_STRIP_GROUP_DICTIONARIES.md#logarithmic-moments)
