# Matched nonlinear transport for all fully shared mirrored profiles

For every fixed fully shared ordinary mirrored profile, place its
exact order-\(n+1\) log-factorial symbol on the order-\(n\)
fractional grid. The roots of that polynomial can be paired with
the roots of the actual order-\(n\) polynomial with angular
displacement \(O(n^{-2})\). Both polynomials have simple unit-circle
roots after positive rescaling, for all sufficiently large \(n\).
The constants and a sufficient terminal threshold are explicit
below.

The [positive Laplace increment](MULTIVARIATE_MATCHED_LAPLACE_NONLINEAR_TRANSFER.md#positive-curvature-increment-for-ordinary-shared-rows)
and its fixed-grid root argument apply. A sharper real-shift
bound uses the full circular mesh, and the second Poisson
moment keeps the required Laplace rate low. The transport here
holds at one fixed degree. An across-degree comparison at a
common symbol remains a separate step.

## The exponential rate from the full circular mesh

Let the reference have degree \(m\ge2\), all roots on the unit
circle, and circular angular mesh at least \(\gamma/m\), with
\(0<\gamma\le2\pi\). Put \(M=\pi/\gamma\), \(D=2+4M^2\), and choose
\[
           0<a_0\le\min\{1/32,1/(16M)\}.
\]
Fix a reference root and use \(F\), \(q\), \(b\), and
\(\widetilde F=e^{ibw}F\) from the
[nonlinear fixed-degree estimate](MULTIVARIATE_MATCHED_LAPLACE_NONLINEAR_TRANSFER.md#a-positive-shift-estimate-at-fixed-degree).
The same harmonic bound gives \(|b|\le V_m=M(1+\log m)\).
For the \(j\)-th circular neighbor of the fixed root, the absolute
sine of the half-gap is at least
\(\sin(\gamma\min(j,m-j)/(2m))\). Consequently for real \(h\) and
\(v=|h|/(2m)\),
\[
 \log\left|\prod_{\zeta\ne z}A_{q_\zeta}(h/(2m))\right|
  \le m I(v),\qquad
 I(v)=\int_0^{1/2}
        \log\left(1+\frac{\sinh^2v}{\sin^2(\gamma x/2)}\right)\,dx .
\]
Indeed the logarithm in the integral decreases with \(x\).
There are at most two roots at each circular neighbor index,
so their contributions to the half-log modulus are bounded
by the right-endpoint Riemann sum.

For \(v>0\), differentiation under the integral and the
substitution \(t=\tan(\gamma x/2)\) give
\[
          I'(v)=\frac4\gamma
                 \arctan\left(\coth v\,\tan(\gamma/4)\right)
                   \le \frac{2\pi}{\gamma}.
\]
At \(\gamma=2\pi\) the formula is interpreted by the limit.
The integral has \(I(0)=0\) by dominated convergence, since its
integrand near \(x=0\) is bounded on a compact \(v\)-range by
\(\log(1+\mathrm{const}/x^2)\). It follows that the product has
modulus at most \(e^{M|h|}\).

Define
\[
       K=M+\frac1{2m}+\frac{a_0D}{2},\qquad
       B=\exp(4a_0^2D+a_0),\qquad b_0=\tfrac12\exp(-4a_0^2D).
\]
The same complex product estimate after the gauge gives
\[
 |\widetilde F(h+y)|\le B(|h|+2a_0)e^{K|h|}\quad (|y|\le2a_0),
 \qquad
 |\widetilde F(y)|\ge b_0|y|\quad (|y|\le2a_0).
\]
In the upper bound, \(M\) is the real product rate, \(1/(2m)\)
comes from the sinh prefactor, and \(a_0D/2\) comes from the
complex perturbation estimate. The sign and size assumptions
on \(a_0\) give the same logarithm-series bounds as before.
Indeed \(Mm\ge1\) and \(M a_0\le1/16\).
Cauchy's estimate followed by integration along the real segment
gives, for \(|y|\le a_0\),
\[
 |\widetilde F(h+y)-\widetilde F(y)|
        \le (B/a_0)|h|(|h|+2a_0)e^{K|h|}.
\]
The polynomial prefactor is handled by the second moment of the
positive jump law.

## Weighted jump moments and the effective threshold

For the positive measure \(\nu\) in the nonlinear transport, set
\[
        \beta=m\int (|h|+h^2)e^{K|h|}\,\nu(dh),\qquad
        \rho=\frac{2B}{a_0 b_0}\beta e^{(K+V_m+1)\beta}.
\]
Then \(\rho\le a_0\) suffices for the same angular bound \(\rho/m\).
In fact, truncate the jumps as in the nonlinear proof and write
\(\beta_1=m\int |h|e^{K|h|}\,d\nu\) and
\(\beta_2=m\int h^2e^{K|h|}\,d\nu\) on the truncation.
The sum \(J\) of the absolute jumps has
\[
 \mathbb E[(J^2+2a_0J)e^{KJ}]
  =e^{m\int(e^{K|h|}-1)\,d\nu}
          \left(\beta_2+\beta_1^2+2a_0\beta_1\right)
  \le \beta e^{(K+1)\beta}.
\]
The equality follows by two derivatives of the exponential
moment. The bound uses \(2a_0\le1\),
\(\beta_1+\beta_2\le\beta\), and \(1+\beta\le e^\beta\).
The Fourier factor used to remove the gauge has modulus at least
\(e^{-V_m\beta}\). Comparing on \(|y|=\rho\) and taking the
truncation limit proves the asserted Rouché pairing.

Now let \(L=\min_i\ell_i/d_i>0\) over the unequal shared rows
of a mirrored ordinary profile. The
[Gegenbauer mesh bound](MULTIVARIATE_MATCHED_LAPLACE_NONLINEAR_TRANSFER.md#actual-mirrored-reference-profiles)
and Schur–Szegő composition give the mesh lower bound
\[
        \frac{\gamma}{m},\qquad
        \gamma=\frac{2\pi}{1+3L/2},
        \qquad m\ge \max\{2,4/L\}.
\]
Use one Gegenbauer factor with ratio \(L\). Its parameter satisfies
\(\alpha\le mL+1\) and \(\alpha\ge1\), and the mesh lower bound
is \(2\pi/(m+\alpha)\) before composition. Choose \(a_0\) as above
with the additional condition \(a_0\le L/(8D)\). These choices give
\[
       K\le K_*=\tfrac12+\tfrac{15}{16}L
             < a_i=\ell_i/d_i+\tfrac12
               \quad\hbox{for every row}.
\]
Indeed \(M=\tfrac12+3L/4\),
\(1/(2m)\le L/8\), and \(a_0D/2\le L/16\).
The positive row-density bound of the nonlinear proof now gives
\[
       \beta\le\frac{m}{2n(n+1)} R,\qquad
       R=\sum_i\left(\frac1{a_i-K_*}+\frac1{(a_i-K_*)^2}\right).
\]
Take \(m=n\) and the actual lower-order polynomial as reference.
For example a sufficient explicit check is
\[
 \frac{2B}{a_0 b_0}\frac{R}{2(n+1)}
       \exp\!\left(\frac{(K_*+V_n+1)R}{2(n+1)}\right)\le a_0 .
\]
The left side tends to zero and the root shift obtained by dividing
it by \(n\) is \(O(n^{-2})\). The mirrored symbol difference keeps
the coefficients real and reciprocal. Reciprocal conjugation
preserves each log disk, so the unique paired root in that disk
lies on the unit circle.

## Why the unadjusted stable compression has a larger defect

For a degree-\(N\) circle-rooted polynomial \(P_N\), the
Cayley derivative from the cut at \(1\) gives
\[
                   J_N=(1-z)P_N'+NP_N .
\]
The Cayley substitution \(z=(x-i)/(x+i)\) turns this into the
derivative of \((x+i)^N P_N((x-i)/(x+i))\) after removal of the
nonzero factor \((x+i)^{N-1}\). Thus \(J_N\) has degree
\(n=N-1\) and strictly interlaces the roots of \(P_N\) when they are simple
and none lies at the cut.

This compression does not match the common-symbol fractional
grid up to the \(O(n^{-1})\) scale in full log coefficients
afforded by the positive jump increment. For a fixed \(C^2\)
symbol \(L\), put \(b_{N,k}=e^{N L(k/N)}\). The normalized
coefficients of \(J_N\), apart from a constant, are
\(b_{N,k}+b_{N,k+1}\). At \(s=k/n\) the limiting log defect
from \(e^{nL(s)}\) is uniform on \([0,1]\) and equals
\[
       D_L(s)= L(s)-sL'(s)+\log(1+e^{L'(s)}),
\]
by Taylor expansion at \(s\). Consequently
\[
    D_L'(s)=\left(\frac{e^{L'(s)}}{1+e^{L'(s)}}-s\right)L''(s).
\]
For a fully shared log-factorial symbol, \(L''>0\) at the
endpoints and the slopes there are finite. The derivative
of the defect is therefore positive at zero and negative at
one. Its nonaffine part cannot be removed by a scalar and
positive rescaling of the polynomial. The same obstruction
holds for any fixed positive relative weight on the second
summand. It leaves open whether the compression defect has
a directional angular inequality suitable for interlacing.
