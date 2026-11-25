# Eventual diffuse coprimality at the critical coupling

Let
\[
 B_n(z;s)=\sum_{k=0}^n\binom nk e^{-s k(n-k)/(2n)}z^k .
\]
At \(s=4\), the polynomials \(B_n,B_{n+1}\) are coprime for all
sufficiently large \(n\). This completes the exactly critical
parameter in the [diffuse Curie–Weiss comparison](MULTIVARIATE_DIFFUSE_LOW_TEMPERATURE_INTERLACING.md).
The argument below keeps \(s=4\) fixed.

The roots lie on the unit circle, so it suffices to consider their
principal absolute arguments \(0\le\theta\le\pi\).
The regular-arc estimate already excludes common roots when
\(\theta\) is bounded away from zero. We treat the bounded
\(n^{3/4}\theta\) window and the remaining small arguments.
[Unitary Hermite circle theorem](https://doi.org/10.5802/ahl.227)

## The bounded quartic window

At \(a=s/2=2\), the exact Gaussian integral after removing
\(e^{in\theta/2}\) has phase
\[
 F(t,\theta)=-\frac12-\frac{t^2}{8}
             +\log\!\left(2\cosh\frac{t+i\theta}{2}\right).
\]
Use the logarithm with the determination real at \(\theta=0\)
on the real line. It exists wherever
\(|\Im(t+i\theta)|<\pi\), since the hyperbolic cosine then has
positive real part. For real \(\theta\) the phase satisfies
\(F(-\bar t,\theta)=\overline{F(t,\theta)}\).

Put \(\lambda=n^{-1/4}\) and \(\theta=\lambda^3 v\). Allow
\(\lambda\) to vary over positive reals and use the real integral
\[
 J(v,\lambda)=\frac1\lambda\int_{\mathbb R}
    \exp\!\left(\lambda^{-4}
      [F(t,\lambda^3v)-F(0,\lambda^3v)]\right)\,dt .
\]
At \(\lambda=n^{-1/4}\) this has the same zeros in \(v\) as the
normalized polynomial. For real \(v\) it is real by reflection.

At \(\theta=0\) the real phase has a unique maximum at \(t=0\).
One can see this from the derivative of the horizontal log modulus
in the regular-arc proof, with \(P_1(r)=r\coth(r/2)>2\) for \(r>0\).
Near zero its difference from the maximum is
\(-t^4/192+O(t^6)\). The retained tails \(|t|\ge r_0\) for any
sufficiently small fixed \(r_0>0\) therefore have strict loss for
small \(\theta\), with Gaussian decay at infinity. In \(|t|\le r_0\),
a mixed derivative bound gives
\[
 \Re[F(t,\theta)-F(0,\theta)]
   \le -c t^4+C|\theta||t|
\]
with \(c>0\). Under \(t=\lambda x\), the exponent is consequently
bounded in real part by \(-cx^4+C|v||x|\).

Taylor expansion of the scaled exponent gives the limit
\[
 \frac{ivx}{4}-\frac{x^4}{192}.
\]
For the derivative bounds, the coefficients of \(t,t^2,t^3,t^4\)
in \(F(t,\theta)-F(0,\theta)\) are respectively
\[
 \frac{i\theta}{4}+O(\theta^3),\quad
 O(\theta^2),\quad O(\theta),\quad
 -\frac1{192}+O(\theta^2).
\]
Their remainders are analytic at zero. The remaining terms have
the form \(t^5 R(t,\theta)\) with \(R\) and its first derivatives
bounded on a fixed neighborhood. Substitution
\(t=\lambda x,\theta=\lambda^3v\) and division by \(\lambda^4\)
show that the \(v\)- and \(\lambda\)-derivatives of the scaled
exponent are bounded by \(C(1+|x|^5)\) for \(|\lambda x|\le r_0\)
and bounded \(v\). In the derivative of \(\lambda x^5 R\),
use \(|\lambda x|\le r_0\) for the term differentiating \(R\)
in its first argument. Endpoint terms at \(|x|=r_0/\lambda\)
are exponentially small. The omitted tails and their derivatives
are also exponentially small, since differentiated factors grow
at most polynomially in \(|t|\) and \(\lambda^{-1}\).
Dominated convergence gives uniform convergence of \(J,\partial_vJ\)
on bounded \(v\)-intervals, with \(\partial_\lambda J\) uniformly
bounded there. The limiting integral is
\[
 J(v,0)=\int_{\mathbb R}e^{-x^4/192+i v x/4}\,dx .
\]

Every real zero of this limit is simple. Rescale it to
\(Y(\omega)=\int_{\mathbb R}e^{-x^4+i\omega x}\,dx\).
Integration by parts gives \(Y'''=(\omega/4)Y\).
The function is real and even for real \(\omega\), and
\(Y(0)>0\). If \(Y\) and \(Y'\) vanished at \(\omega_0>0\),
then \(Y''(\omega_0)\ne0\) by uniqueness for this ODE.
After changing the sign if needed, \(Y''(\omega_0)>0\).
Immediately to the right \(Y>0\), so \(Y'''=(\omega/4)Y>0\).
The second derivative remains at least its initial positive value
as long as \(Y>0\), which then persists and forces at least quadratic
growth. This contradicts \(Y(\omega)\to0\), which also follows by
integration by parts.

On any bounded \(v\)-interval there are only finitely many limit
zeros. The positive ones have disjoint neighborhoods with
\(|\partial_v J|\) bounded below for small \(\lambda\).
Each neighborhood has one real root \(v_k(\lambda)>0\) of
\(J\), and away from the neighborhoods \(J\) has no zeros
on the bounded interval for small \(\lambda\). The implicit
function theorem gives \(\partial_\lambda v_k=O(1)\).
Consequently
\[
 \frac{d}{d\lambda}\bigl(\lambda^3 v_k(\lambda)\bigr)
  =3\lambda^2v_k(\lambda)+\lambda^3v_k'(\lambda)>0
\]
for sufficiently small \(\lambda\). The angles at \(n,n+1\)
cannot agree on a branch. Different branches cannot agree either:
their scaled \(v\)-coordinates are separated by a positive margin,
whereas a common angle would make those coordinates differ by
a factor \((1+1/n)^{3/4}\). This proves the exclusion in each
bounded \(n^{3/4}\theta\) window.

## Intermediate small positive arguments

Use the left saddle \(-r+iy\) from the regular-arc construction.
Here \(a=2\), \(C=\cos u=2\sinh r/r-\cosh r\), and
\(y=r\sin u/\sinh r\). Taylor expansion gives
\[
 u=\frac r{\sqrt3}+O(r^3),\qquad
 y=\frac r{\sqrt3}+O(r^3),\qquad
 \theta=u-y=\frac{2}{9\sqrt3}r^3+O(r^5).
\]
Center the reflection-normalized phase at \(iy\), writing
\[
 F_r(\zeta)=-\frac12-\frac{(iy+\zeta)^2}{8}
              +\log\!\left(2\cosh\frac{\zeta+iu}{2}\right),
 \qquad F_s=F_r(-r).
\]
The limiting polynomial in the scaled coordinate \(w=\zeta/r\)
is
\[
 \frac{192}{r^4}\bigl(F_r(rw)-F_r(0)\bigr)
   \longrightarrow -w^4+2w^2+\frac{4i}{\sqrt3}(3w-w^3)
\]
locally uniformly in \(w\). Indeed the first four Taylor
coefficients in \(\zeta\) are
\[
 \frac{ir^3}{16\sqrt3}+O(r^5),\quad
 \frac{r^2}{96}+O(r^4),\quad
 -\frac{ir}{48\sqrt3}+O(r^3),\quad
 -\frac1{192}+O(r^2),
\]
respectively. In particular
\(\Im F_s=-r^4/(24\sqrt3)+o(r^4)\), since \(F_r(0)\) is real.

There is a uniform loss estimate on the original horizontal
left half-line through the saddle. Its log modulus, apart
from a constant, is the function
\[
 \mathcal F_r(Y)=-Y/8+\tfrac12\log(\cosh\sqrt Y+C),
 \qquad Y=(\Re\zeta)^2 .
\]
The function is analytic in \(Y\) near zero and its second
derivative at \(r=Y=0\) is \(-1/96\).
At \(Y=r^2\) the first derivative is zero by stationarity.
On a small fixed \(0\le Y\le r_0^2\) the second derivative
is bounded above by a negative constant for small \(r\).
Thus, for real \(\zeta\) with \(|\zeta|\le r_0\),
\[
 \Re(F_r(\zeta)-F_s)\le -c(\zeta^2-r^2)^2.
\]
Off this fixed interval the original horizontal has uniform
strict loss, including at \(r=0\), by the critical horizontal
maximum and continuity. The far tail has Gaussian decay.

Let \(m=nr^4\) and
\(\Psi_r(w)=(F_r(rw)-F_s)/r^4\).
Shift the original contour to the horizontal through \(iy\)
and split at \(\zeta=0\). The normalized polynomial is twice
the real part of the left packet. This representation also
gives its real-order interpolation in the logarithm strip.
The left packet is
\[
 I_n=\sqrt{\frac n{8\pi}}\,r e^{nF_s} K(m),\qquad
 K(m)=\int_{-\infty}^0 e^{m\Psi_r(w)}\,dw .
\]
The scaled phase has a nondegenerate saddle at \(w=-1\).
Indeed \(A_r=-F_r''(-r)/r^2\) is bounded and has real part
bounded below by a positive constant, by the uniform loss
and the scaled Taylor coefficients. The third derivative of
\(\Psi_r\) is uniformly bounded near \(w=-1\). Its real loss
is at least \(c(w^2-1)^2\) when \(|rw|\le r_0\). In that region
the same Taylor coefficients give \(|\Psi_r(w)|\le C(1+|w|^4)\).
The endpoint \(w=0\) has strict loss in the scaled phase.

Gaussian comparison on a fixed neighborhood of \(w=-1\)
and exponential loss elsewhere give
\[
 \sqrt m K(m)=C_r+O(m^{-1/2}),\qquad C_r=\sqrt{2\pi/A_r}.
\]
The amplitudes \(C_r\) are bounded away from zero. Differentiation
in \(m\) is uniform:
\[
 \frac{d}{dm}(\sqrt m K(m))
   =\frac{K(m)}{2\sqrt m}
      +\sqrt m\int_{-\infty}^0\Psi_r(w)e^{m\Psi_r(w)}\,dw
   =O(m^{-1}).
\]
In the saddle neighborhood, the differentiated integrand
has the additional factor \(\Psi_r(w)=O((w+1)^2)\), giving
the stated order after normalization. On the other bounded
and intermediate pieces the loss dominates the polynomial
bound for \(\Psi_r\). The far-left part has exponential
loss in \(m/r^4\), absorbing the rescaled lengths and all
powers of \(r^{-1}\) from differentiation.

Thus \(\sqrt m K(m)=C_r(1+\eta_r(m))\) with
\(\eta_r=O(m^{-1/2})\) and \(\eta_r'=O(m^{-1})\).
The change in argument of the nonzero left packet from
\(n\) to \(n+1\) is
\[
 \Im F_s+O(r^4/m)
   =-\frac{r^4}{24\sqrt3}+o(r^4)+O(r^4/m).
\]
The real prefactors and \(C_r\) contribute no change in
argument. For small \(r\) and sufficiently large \(m\),
the change lies strictly between \(-\pi\) and zero, so the
left packets cannot both be purely imaginary. This excludes
a common root. Since \(m=nr^4\) is comparable to
\((n^{3/4}\theta)^{4/3}\), the bounded quartic window covers
bounded \(m\), and the regular arcs handle fixed positive
distance from zero. This gives eventual coprimality at \(s=4\).
