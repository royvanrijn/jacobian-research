# Exterior rigidity for the Mehler equation

Fix \(0<r<1\), let \(T_r\) be the one-dimensional Gaussian
noise operator, and let \(q\) be a bounded measurable function.
If \(\lambda\) is real with \(|\lambda|<1\) and
\[
                  T_rq(x)=\lambda q(x)
              \quad\text{for almost every }x>S,
                                                        \tag{EL1}
\]
then \(q=0\) almost everywhere on the whole line. The conclusion
allows arbitrary bounded data on the other side of \(S\).

## Iteration at expanding observation points

Put \(f=T_rq-\lambda q\), so \(f\) is bounded and supported
in \((-\infty,S]\). Denote convolution by the normal density of
standard deviation \(a\) by \(P_a\), with \(P_0\) the identity.
The semigroup relation \(T_uT_v=T_{uv}\) gives
\[
  P_{\sqrt{1-r^{2n}}}q(x)
     =\lambda^n q(x/r^n)
      +\sum_{k=1}^n \lambda^{k-1}
         P_{\sqrt{1-r^{2(n-k)}}} f(x/r^k).
\]
For almost every \(x\), the first term tends to zero. Each
fixed summand converges and the boundedness of \(f\) gives
a summable dominating sequence, with the convention
\(\lambda^0=1\) if \(\lambda=0\). Thus
\[
         P_1 q(x)=\sum_{k=1}^{\infty}
                      \lambda^{k-1}P_1f(x/r^k).          \tag{EL2}
\]
Both sides are continuous, so this holds everywhere. For
\(x/r>S\) and \(x>0\), the support of \(f\) gives
\[
  |P_1 q(x)|
       \le \frac{\|f\|_\infty}{1-|\lambda|}
                    \exp\bigl(-(x/r-S)^2/2\bigr).
                                                        \tag{EL3}
\]
In particular the decay on the positive half-line is
strictly faster than the variance-one Gaussian rate.

## One-sided Gaussian deconvolution

The one-sided Hardy theorem in the critical case says that
if \(g\) is tempered, \(g(x)e^{x^2/2}\) on the positive
half-line is tempered, and \(\widehat g(\xi)e^{\xi^2/2}\)
is tempered under the unitary Fourier convention, then
\(g=P_1 \mu\) for a tempered distribution \(\mu\) supported
in \((-\infty,0]\). In Demange's Theorem 1.6.8 take
\(u\mapsto g(\sqrt{2\pi}\,u)\), using the Fourier convention
with kernel \(e^{-2\pi i u\eta}\). Both critical weights
become \(e^{\pi u^2}\) and \(e^{\pi\eta^2}\), respectively.
The Gaussian \(e^{-\pi u^2}\) is the rescaled variance-one
density, and the dilation preserves the supporting half-line.
[Bruno Demange, *Uncertainty principles associated to non-degenerate quadratic forms*](https://www.numdam.org/item/10.24033/msmf.431.pdf)

For every real \(a\), apply the theorem to \(g(x)=P_1q(x+a)\).
The real-axis hypothesis follows from (EL3), including
after translation. The Fourier-side hypothesis holds
since multiplication by the reciprocal Gaussian recovers
the tempered Fourier transform of \(q(\,\cdot+a)\).
Injectivity of Gaussian convolution on tempered
distributions identifies \(\mu=q(\,\cdot+a)\). Indeed the
Fourier multiplier is nonzero and can be inverted
against every compactly supported smooth test function.
Consequently \(\operatorname{supp}q\subseteq(-\infty,a]\)
for every \(a\), which proves (EL1).

The special case \(\lambda=r\) is the fixed-pressure
linearization of an isolated two-label graph at a
horizontal interface. For a vertical height perturbation
\(h(x)=b+\epsilon q(x)\), the derivative of the stationary
contrast there is a positive multiple of
\(T_rq(x)-r q(x)\). The two terms come from the height
of the source interface and the shifted observation,
respectively. The proof therefore excludes bounded
solutions of the linear equation with residual supported
on the left, even when the height perturbation changes
sign infinitely often. The exact graph equation has a
nonlinear residual on the right. Shifted-contour rate
improvement controls that residual for bounded isolated
graphs and gives exact straightness.
[Isolated interface ends](THREE_CELL_ISOLATED_END_RIGIDITY.md)
[Nonlinear straightening by shifted-contour interpolation](THREE_CELL_TAIL_FLATNESS.md#improving-the-gaussian-rate-by-shifting-the-inversion-contour)

## Exterior angular regions in any dimension

Let \(q\) be bounded on \(\mathbb R^d\), \(d\ge 1\), and write
\(f=T_rq-\lambda q\) for the \(d\)-dimensional Mehler operator.
Suppose there is a unit direction \(e\) with
\[
 \liminf_{t\to\infty}
       \frac{\operatorname{dist}(te,\operatorname{ess\,supp}f)}{t}
       > r.                                                \tag{EL4}
\]
Then \(q=0\). For example, (EL4) holds if the equation holds
on a half-space, or on a translated circular cone of half-angle
\(\theta\) with \(\arcsin r<\theta\le\pi/2\).

The derivation of (EL2) works with isotropic \(d\)-dimensional
Gaussian convolution \(P_1^{(d)}\). Let \(y\perp e\) be fixed
and choose \(r<c<1\) strictly below the limit inferior.
For every \(k\ge1\) and sufficiently large \(t\), the distance
of \((te+y)/r^k\) from the support of \(f\) is at least
\(c t/r^k\). This follows from the distance bound on the
central ray and the triangle inequality, allowing a margin
between \(c\) and the limit inferior to absorb \(|y|/t\).
For standard \(d\)-dimensional \(Z\) and \(0<\eta<1\),
\(\Pr(\|Z\|\ge u)\le \eta^{-d/2}e^{-(1-\eta)u^2/2}\)
by the exponential moment of \(\|Z\|^2\). Choose \(\eta\)
so that \((1-\eta)c^2>r^2\). The series in (EL2) bounds
\(P_1^{(d)}q(te+y)\) by a constant times
\(\exp(- (1-\eta)c^2 t^2/(2r^2))\).

For each fixed \(y\) the function \(t\mapsto P_1^{(d)}q(te+y)\)
is the one-dimensional Gaussian convolution of the bounded
longitudinal function obtained by smoothing \(q\) only in
the transverse coordinates at \(y\). The argument using
translated one-sided Hardy bounds makes this slice zero.
Thus \(P_1^{(d)}q=0\), and Gaussian injectivity gives \(q=0\).
