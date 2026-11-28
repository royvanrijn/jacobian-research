# The limiting angle bound in the open-semicircle geometry

Under the exact-three bounded wedge-end hypotheses of the
[acute-end projection analysis](THREE_CELL_ACUTE_END_PROJECTION.md),
suppose all three outward directions of an actual balanced
planar maximizer lie in an open semicircle. The smallest acute
angle \(\theta_0\) satisfies
\[
                   \tan^2\theta_0+\tan^4\theta_0<r^2.
\]
The inequality with a weak sign is proved by the
[open-semicircle angle bound](THREE_CELL_ACUTE_END_PROJECTION.md#an-angle-bound-in-the-open-semicircle-geometry).
We exclude equality by identifying the two additional critical
terms in the transverse projection equation.

Keep the notation of (AU1)–(AU5). Equality gives
\(\lambda_0+\kappa_0=\Lambda\).
If the successive gaps are unequal, use the extreme arm adjacent
to the larger gap. It has no neighbor at the globally smallest
angle, so its height rate exceeds \(\lambda_0\) by (AU1).
Its direct leakage has strictly larger rate than \(\lambda_0\)
before longitudinal convolution and strictly larger rate than
\(\kappa_0\) after it. The residual for this arm in (AU3) is
therefore faster than \(\Lambda\). In fact the observation
square terms are faster by \(2\lambda_0>\Lambda\), the observation
shift times the convolved direct error is faster by the two
strict rate gains, and the direct affine replacement in (AU2)
is faster by the same square bound. Both steps of critical
Hardy from the weak angle bound apply to this arm. The exact
projection reconstructed from its affine models has limit zero
or sign opposite the incident spectator at \(-\infty\), whereas
the geometric limit has the spectator sign.

## Critical reconstruction for equal gaps

Now suppose both successive gaps equal \(\theta_0\).
Put \(m=\tan\theta_0>0\) and \(\alpha=r\sin\theta_0>0\).
Equality gives \(m^2+\alpha^2=r^2\), since
\(m^2+m^4=r^2\). In particular \(\theta_0<\pi/4\).
The two extreme arms have mutual angle \(2\theta_0\) and
its direct source rate exceeds \(\Lambda\), since
\[
         \tan^2(2\theta_0)=\frac{4m^2}{(1-m^2)^2}
                          >m^2+m^4.
\]
Thus, for either extreme arm, only the nearest neighbor
contributes an inverse affine model at or below the critical
rate. Work on one extreme arm in its lower-minus-upper
orientation. Let \(\eta=\pm1\) be the sign of its nearest step
contribution, and use its shift \(\tau\) from (AU1).
Set \(c=\tau\cos\theta_0\), \(c_1=\psi'(d)\), and
\[
  X(t)=\Phi(-(m t+\tau)/\sigma),\qquad
  Y(t)=\Phi(-(\alpha t+c)/\sigma).
\]
Mehler convolution gives \(T_r X=Y\). By (AU1), \(u=h-b\)
equals \(-\eta X/c_1\) up to a faster quadratic Gaussian.
Consequently \(T_r\widetilde u=-\eta Y/c_1\) up to a faster
Gaussian, where \(\widetilde u\) is \(u\) cut off on a late
right half-line. The faster rate here is measured against
\(\kappa_0\). The omission on the left has rate \(\Lambda\)
with a linear shift, and \(\Lambda>\kappa_0\).
At observation height \(b\), the vertical derivative of the
nearest signed convolved step equals
\[
         \frac{r\cos\theta_0}{\sigma}\,
                     \phi((\alpha t+c)/\sigma)
\]
up to a faster Gaussian with rate greater than \(\kappa_0\).
The upper step has positive contrast and its probability
increases with observation height. The lower step has negative
contrast and its probability decreases, giving the same
positive sign. The affine replacement error and its derivative
are faster by the thin-strip estimate of (WG3).
[Affine step replacement, WG3](THREE_CELL_WEDGE_END_ASYMPTOTICS.md#leading-interaction-between-the-affine-directions)

Expand the exact identity \(T_r[C(\cdot,b+u(t))](t)=\beta\)
through first order in \(u(t)\). The binary part of the
smoothed vertical derivative is
\(-r c_1-r\psi''(d)T_r\widetilde u+
O(T_r\widetilde u^2)\), up to the bounded-left error.
Using \(V=Q+D\) and \(Q=c_1u+O(u^2)\) shows that the
right residual in (AU3) is
\[
 \begin{split}
 F_{\rm err}(t)=&
       \frac{r\psi''(d)}{c_1^2}\,X(t)Y(t)
      +\frac{\eta r\cos\theta_0}{\sigma c_1}\,
                           X(t)\phi((\alpha t+c)/\sigma)\\
     &+O(e^{-(\Lambda+\delta)t^2})                        \tag{AB1}
 \end{split}
\]
for some \(\delta>0\). The terms displayed arise respectively
from \(r\psi''(d)u\,T_r\widetilde u\) and the negative of \(u\)
times the smoothed step derivative. All other terms are faster.
In detail, (AU2) is faster since \(2\lambda_0>\Lambda\),
as is the observation square. The strip derivative remainder
is quadratic in the source displacement, so its longitudinal
smoothing has rate strictly above \(\kappa_0\), and it is
multiplied by \(u\). The far arm's direct term is kept in
\(C_{\rm aff}\) and has rate above \(\Lambda\), while its
smoothed derivative is faster than that of the nearest arm.

Each product in (AB1) is an exact Mehler convolution of a
bounded source profile supported on a left half-line.
For \(XY\), take independent normal source coordinates
\((P,Q)\) with means \((mt,\alpha t)\) and variance \(\sigma^2\)
each. The event is \(P\le-\tau,\ Q\le-c\).
Rotate so the first coordinate is along \((m,\alpha)/r\), which
makes its mean \(rt\) and the other mean zero. At fixed
first coordinate \(s\), the event is impossible for
\(s>(-m\tau-\alpha c)/r\). Its conditional probability tends to
one as \(s\to-\infty\) because both components of \((m,\alpha)\)
are positive. This conditional probability is the required
source profile for \(XY\).

For \(X\phi((\alpha t+c)/\sigma)\), fix \(Q=-c\) and integrate
the joint conditional density over \(P\le-\tau\), with the
prefactor \(\sigma\). After the same rotation the first
coordinate \(s\) has Jacobian \(ds/dP=m/r\). The profile
is a constant multiple of a transverse Gaussian density
evaluated at an affine function of \(s\), on the same
left half-line. This affine function has nonzero slope
because \(\alpha>0\). The profile is bounded and tends to zero
as \(s\to-\infty\).

At equality the first inverse of the nearest affine
CDF has rate
\[
 \frac{\lambda_0}{R-2\sigma^2\lambda_0}
       =\frac{1}{2\sigma^2m^2}>\Lambda ,
 \qquad R=r^2.
\]
Thus \(J\) has just this one inverse term, its terminal
error is faster than \(\Lambda\), and \(J\) tends to
\(-r\eta\) at \(-\infty\).
The first Hardy step for \(H=V-J\) applies to (AB1) and
(AU4) at rate \(\Lambda\), with the linear shifts
absorbed by translation. Hence \(H\) is supported on the
left. On the forward tail, subtract the two exact
convolutions in (AB1) from \(T_rH\). The remaining
convolution is supercritical, so the second Hardy step
identifies \(H\) globally with the corresponding source
profiles. Their limits, together with the affine model
\(J\), give
\[
     \eta-\beta=-r\eta-\frac{r z_\beta}{2\phi(z_\beta)},
       \qquad z_\beta=\Phi^{-1}((1+\beta)/2).              \tag{AB2}
\]
The left side is the reverse geometric limit of \(V\).
The last term is \(r\psi''(d)/c_1^2\).
The second critical product makes no contribution to this
limit.

Since \(\eta-\beta\) has sign \(\eta\), (AB2) forces \(\beta\)
to have sign \(-\eta\).
Putting \(p=|\beta|>0\) gives
\(r z_p/(2\phi(z_p))-p=1+r\).
Its left side vanishes at \(p=0\) and has strictly increasing
derivative \(r(1+z_p^2)/(4\phi(z_p)^2)-1\) on \(0<p<1\).
It therefore takes the given positive value at most once.
The two extreme pressures are opposite with equal magnitude,
and the middle pressure is zero by the cyclic pressure sum.

For the middle arm, the two direct affine contributions
cancel exactly. Their angles are opposite, the middle
pressure offset is zero, and the two outer offsets are
opposite. Equation (AU1) now gives the middle height a rate
strictly larger than \(\lambda_0\). In (AU3) its affine
forcing is zero and its residual is faster than \(\Lambda\),
using the improved middle rate in both observation terms
and (AU2) for the direct replacement. The first Hardy
step gives \(V\) supported on the left, and \(T_r V\)
then has rate strictly greater than \(\Lambda\).
Hardy at variance \(\sigma^2\) forces \(V=0\) globally.
The middle pair score is consequently equal to its
pressure on its entire supporting line. Its transverse
derivative stays bounded away from zero in a fixed tube
about this line on the tail, by the wedge comparison.
The actual graph must coincide with the line. A straight
wedge end cannot have acute neighbors in the exact-three
maximizer geometry, yielding the contradiction.
[The straight-end exclusion](THREE_CELL_WEDGE_END_ASYMPTOTICS.md#geometry-when-at-least-one-wedge-end-is-straight)
