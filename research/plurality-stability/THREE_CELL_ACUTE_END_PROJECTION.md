# Unsmoothed asymptotics at the closest acute ends

Consider an actual balanced planar three-cell Gaussian
noise-stability maximizer at \(0<r<1\). Suppose that outside
a compact set its interface consists of exactly three ends,
one for each pair, in distinct outward directions. Each
end is a bounded transverse graph about its outward direction
in the sense of the wedge-end analysis. Assume at least one
pair of outward directions makes an acute angle. Write
\[
 R=r^2,\qquad \sigma^2=1-R,\qquad
 \Lambda=\frac{R}{2\sigma^2},\qquad
 a(\xi)=\frac{R\xi}{1+2\sigma^2\xi}.
\]
Let \(\theta_0\) be the smallest angle between distinct
acute directions, and set \(\lambda_0=\tan^2\theta_0/(2\sigma^2)\),
\(\kappa_0=a(\lambda_0)=\Lambda\sin^2\theta_0\).

In the coordinates of any reference end, the actual height
displacement has, at rate \(\lambda_0\), the negative of the
direct transverse Gaussian leakage of its closest acute
neighbors. More precisely, use the reference lower-minus-upper
contrast, set \(d=(1-r)b\), and let \(\theta_j\) be the relative
signed direction of another arm \(j\), with upward unit normal
\((-\sin\theta_j,\cos\theta_j)\) whenever \(\cos\theta_j>0\).
Write \(b_j\) for that arm's affine offset in its own normal
coordinate. Then for some \(\delta>0\),
\[
 \begin{split}
 h(t)-b
   &=-\frac1{\psi'(d)}
      \sum_{|\theta_j|=\theta_0}\operatorname{sgn}(\theta_j)
         \Phi\!\left(-\frac{|\tan\theta_j|\,t+\tau_j}{\sigma}\right)
      +O(e^{-(\lambda_0+\delta)t^2}),\\
 \tau_j&=\operatorname{sgn}(\theta_j)
                 (b_j/\cos\theta_j-rb),\qquad
 \psi(z)=2\Phi(z/\sigma)-1 .                              \tag{AU1}
 \end{split}
\]
An empty sum gives the strictly faster decay.
If \(\lambda_0\ge\Lambda\), the exact conditional transverse
projection of each pair contrast on its pressure line equals
the pressure on a sufficiently late source half-line.
[Wedge flattening and acute-arm geometry](THREE_CELL_WEDGE_END_ASYMPTOTICS.md)

At a reference end with a unique nearest arm making the
global minimum angle, the leading term in (AU1) has fixed
sign. Its first two derivatives also give the leading slope
and curvature. The fixed-partition correlation derivative
along this arm is nonconstant for every pressure and either
orientation of its acute neighbor. Consequently in the
exact-three-end geometry, a nonsimplex higher-order first tie
is excluded for every angular configuration.

## The exact transverse projection

For a reference pair let \(f=1_{\rm lower}-1_{\rm upper}\),
\(\beta=\psi(d)\), and define the conditional contrast at a
fixed source longitudinal coordinate by
\[
 C(s,z)=\int f(s,w)p_\sigma(w-rz)\,dw,\qquad
 V(s)=C(s,b)-\beta .
\]
Use \(u=h-b\) and \(Q(s)=\psi(d+u(s))-\psi(d)\) on the
reference tail. Extend \(Q\) by zero on the left.
For sufficiently large \(s\), the reference binary
step and the entering arms give
\(V(s)=Q(s)+D(s)\), where \(D(s)\) is the sum of direct
vertical conditional masses of their far-side steps,
with sign \(+1\) above and \(-1\) below. This uses the
ordered label transitions on vertical slices in (WG3).
There is no wedge error in this definition of \(V\).

Let \(C_{\rm aff}(s)\) be the sum of these signed tail
probabilities when every entering arm is replaced by
its affine line. For a single entering arm the direct
unsigned term is
\(\Phi(-(|\tan\theta_j|s+\tau_j)/\sigma)\).
Each has quadratic rate at least \(\lambda_0\), with
a possible linear shift. The same rate bound, up to
arbitrarily small loss, holds for the actual direct
step probabilities and their first derivatives with
respect to the transverse observation \(z\) near \(b\).
Indeed the actual vertical intersection height is
\(s\tan\theta_j+b_j/\cos\theta_j+O(1)\).

If all three height displacements have a common Gaussian
rate \(\lambda>0\), then, for every \(\epsilon>0\),
\[
      |D(s)-C_{\rm aff}(s)|
           \le C_\epsilon e^{-(\lambda_0+\lambda-\epsilon)s^2}
                          \qquad(s\to+\infty).             \tag{AU2}
\]
At a given vertical slice \(s\), the actual-arm
intersection is at longitudinal parameter
\(\ell=s/\cos\theta_j+O(1)\) in its own coordinates.
Its vertical displacement from the affine intersection
equals \((h_j(\ell)-b_j)/\cos\theta_j\).
Integrating the conditional vertical density over this
small interval yields the sum of rates
\(\tan^2\theta_j/(2\sigma^2)+\lambda/\cos^2\theta_j\),
up to linear shifts and arbitrarily small rate losses.
This proves (AU2).

Taylor expansion in the observation shift gives a residual
whose direct forcing is \(D\), rather than the slower
longitudinal smoothing of \(D\).
Set \(k=\psi''(d)/\psi'(d)\).
For \(s\) on the forward tail,
\[
 C_z(s,b)=-r\psi'(d)-rk V(s)+O(Q(s)^2+E_1(s)),
\]
where \(E_1\) bounds the absolute direct step
probabilities and their first vertical derivatives.
This follows from the binary derivative
\(-r\psi'(d+u(s))\) and
\(\psi'(d+u)=\psi'(d)+kQ+O(Q^2)\).
On the left use an additional bounded term supported
on a left half-line. The \(z\) second derivative
of \(C\) is uniformly bounded by conditional Gaussian
derivative moments.

On the graph, \(T_r[C(\cdot,b+u(t))](t)=\beta\),
where \(T_r\) in this formula acts in the longitudinal
source coordinate. Consequently
\[
 (1-rku(t))T_rV(t)
  =r\psi'(d)u(t)+O\bigl(u(t)^2+
                  |u(t)|T_r(Q^2+E_1)(t)\bigr).
\]
Here \(E_1\) includes the bounded-left term.
The prefactor stays away from zero. Subtract \(rV\),
using \(V=Q+D\) and \(Q=\psi'(d)u+O(u^2)\).
By (AU2) this gives
\[
 \begin{split}
 T_rV-rV&=-r C_{\rm aff}+F_{\rm err},\\
 |F_{\rm err}(t)|&\le C_\epsilon e^{-d_\lambda t^2},\qquad
 d_\lambda=\min\{\mu(\lambda),\lambda+\kappa_0-\epsilon,
                                    \lambda_0+\lambda-\epsilon\},
                                                               \tag{AU3}\\
 \mu(\lambda)&=\min\{2\lambda,\lambda+a(2\lambda)\}.
 \end{split}
\]
Indeed Gaussian convolution turns the positive-half-line
rate \(\lambda_0-\epsilon'\) of \(E_1\) into a rate
arbitrarily close to \(\kappa_0=a(\lambda_0)\).
The bounded-left part has quadratic rate \(\Lambda\)
with a possible linear shift and causes no loss here.
The pure square term and the term \(u\,T_r Q^2\)
give \(\mu\).
All quantities in (AU3) are bounded on the left.

## Removing the known direct forcing

Each single affine term in \(C_{\rm aff}\) has the form
\(C^{(0)}(s)=\eta\Phi(-(\alpha s+\tau)/\omega)\), where
\(\alpha>0\), \(\omega>0\), and \(\eta=\pm1\).
Its quadratic rate is \(\gamma=\alpha^2/(2\omega^2)\).
If \(\gamma<\Lambda\), it has an explicit bounded
Mehler inverse \(C^{(1)}\). Use the CDF with parameters
\(\alpha/r,\ \tau,\ \sqrt{\omega^2-\sigma^2\alpha^2/R}\)
in place of \(\alpha,\tau,\omega\). Independent Gaussian
variances give \(T_r C^{(1)}=C^{(0)}\).
The new rate is \(\gamma/(R-2\sigma^2\gamma)>\gamma/R\).
Repeat until the current rate is at least \(\Lambda\).
This requires finitely many steps, since below \(\Lambda\)
the rates increase by at least the factor \(1/R\).
When the initial rate is already at least \(\Lambda\),
take no inverse steps.

Define \(J\) by summing \(-\sum_{n=1}^N r^n C^{(n)}\)
over the affine terms, with its own stopping index for
each term. Direct telescoping gives
\[
       T_rJ-rJ=-r C_{\rm aff}+\mathcal E,\qquad
       |\mathcal E(t)|\le C e^{-\Lambda t^2+C't}
                         \quad(t\to+\infty).             \tag{AU4}
\]
All terms of \(J\) have rates strictly larger than
\(\lambda_0\). If \(\lambda_0\ge\Lambda\), then \(J=0\).
Set \(H=V-J\); its residual \(T_rH-rH\) has the two
bounds in (AU3) and (AU4).

## The simultaneous bootstrap

Initially all three ends have a positive common exponent
by (WG1). For any \(0<\lambda<\lambda_0\) already reached
for the heights, \(H\) has a Gaussian right bound
\(C e^{-\lambda t^2}\) and is bounded on the left.
If \(d_\lambda<\Lambda\), apply the shifted-contour step
(TF8)–(TF14) to \(H\) with this \(d_\lambda\).
The positive-real variance-one smoothing bound follows
from the Mehler telescoping series. In that series the
bounded-left residual has rate strictly greater than
\(1/2\), and the critical term in (AU4) has rate at least
\(1/2\) with a possible linear shift.
The complex bound on positive real parts still holds.
The right Gaussian bound gives
\(|P_1H(x+iy)|\le Ce^{-\lambda x^2/(1+2\lambda)+y^2/2}\)
there, since the bounded-left part has the faster
rate \(1/2\) with a shift. The sector and the vertical
contour shift only use positive real parts.
[Shifted Gaussian inversion, TF8–TF14](THREE_CELL_TAIL_FLATNESS.md#improving-the-gaussian-rate-by-shifting-the-inversion-contour)

With \(\epsilon=0\), the gain inequality
\(d_\lambda(1+2\sigma^2\lambda)>(1+R)\lambda\)
is strict for \(0<\lambda\le\lambda_0\) on each of the
three branches in (AU3). For \(\mu\) this is (TF14).
For \(\lambda+\kappa_0\), the difference is
\[
     (R-2\sigma^2\kappa_0)(\lambda_0-\lambda)
                          +2\sigma^2\lambda^2>0 .
\]
For \(\lambda_0+\lambda\), use \(\lambda_0+\lambda\ge2\lambda\).
After fixing a positive starting exponent and taking
\(\epsilon\) small enough, the gain remains strict
throughout the compact interval from that exponent to
\(\lambda_0\). Thus the shifted-contour step improves
the right rate of \(H\). Recover a common height
exponent from \(Q=V-D=H+J-D\) up to the direct rate
\(\lambda_0\). The linear shifts in \(D\) can be
absorbed by taking any slightly smaller exponent.

If \(d_\lambda\ge\Lambda\), then the positive-real
telescoping series gives \(|P_1H(t)|\le C e^{-t^2/2+C't}\).
The critical one-sided Gaussian deconvolution theorem,
after translation, implies that \(H\) is supported
on a left half-line. This conclusion also holds with
equality at the rate \(\Lambda\). Otherwise the
continuous rate update increases every exponent below
\(\lambda_0\) and extends with a strict gain at
\(\lambda_0\), if \(d_{\lambda_0}<\Lambda\).
If the iteration approaches the boundary \(d_\lambda=\Lambda\)
from below, the contour update tends to \(\Lambda\).
An infinite sequence trapped at a limit \(\lambda\le\lambda_0\)
there would have \(\lambda=\Lambda=d_\lambda\),
contrary to the strict gain inequality.
These observations allow finitely many updates until
\(H\) either vanishes on a tail or has a rate strictly
larger than \(\lambda_0\).
[Critical one-sided deconvolution](GAUSSIAN_EXTERIOR_LINEAR_RIGIDITY.md#one-sided-gaussian-deconvolution)

We have therefore obtained \(V=J+H\) with a rate
strictly larger than \(\lambda_0\).
When \(\lambda_0\ge\Lambda\), the Hardy conclusion
must occur before passing \(\Lambda\) and \(J=0\), so
\(V\) vanishes exactly on a tail. Equation (AU2) and
\(Q=V-D\) now identify the leading term of \(Q\) as
\(-C_{\rm aff}\) at rate \(\lambda_0\).
The quadratic error in inverting \(Q\) to get \(u\)
is faster, using any achieved height rate above
\(\lambda_0/2\). This proves (AU1).
All fixed derivatives of \(u\) are bounded by the
implicit score equation. The Taylor interpolation argument
in (WG7) therefore preserves a strict rate gap when
differentiating the error any fixed number of times.

## Correlation contact

At an end with a unique globally closest acute neighbor,
apply the score identity (WG7) with upper-minus-lower
orientation, where \(E\) denotes the signed outside
contribution to the upper-minus-lower score.
The left side is faster than \(\kappa_0\), by (AU1)
and \(\lambda_0>\kappa_0\).
The strip term in (WG7) is also faster. Indeed it is bounded by \(C T_r u^2(t)\) and
\(a(2\lambda_0)>\kappa_0\).
The affine replacement errors and all other arms have
strictly faster rates by (WG3). So (WG7) gives
\[
             r(g(t)-g_\infty)=
                 \mathcal B_b E(t,b)
                   +O(e^{-(\kappa_0+\delta)t^2})          \tag{AU5}
\]
for some \(\delta>0\).
Here the operator term on the right has a single
nonzero signed Gaussian leading term of rate
\(\kappa_0\), from the second vertical derivative of
the nearest affine tail. Its coefficient is nonzero
because \(\cos\theta_0>0\).
Thus \(g\) is nonconstant on the reference tail.
[The score equation and signed leading second derivative, WG7](THREE_CELL_WEDGE_END_ASYMPTOTICS.md#a-signed-derivative-obstruction-at-a-wedge-end)

For the open-semicircle case, one of the extreme
directions belongs to a globally closest acute pair
and has a unique acute neighbor at that angle.
Equation (AU5) excludes constant \(g\) on that pair.
Hence the mixed contact inequality excludes an
endpoint where the deficit's first and second
correlation derivatives both vanish. The positively
spanning and antipodal cases already reduce to the
centered simplex or to the global reflection result.
The argument does not preclude an open-semicircle
maximizer with a transverse or ordinary quadratic
correlation contact.
[Geometry and earlier contact reductions](THREE_CELL_WEDGE_END_ASYMPTOTICS.md#geometry-when-at-least-one-wedge-end-is-straight)
[Mixed contact inequality](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#a-first-crossing-obstruction-at-the-level-of-hermite-weights)

## An angle bound in the open-semicircle geometry

Suppose the three directions lie in an open semicircle, and
put \(x=\tan^2\theta_0\) for the smallest acute angle.
An actual maximizer in this geometry must have \(x+x^2\le r^2\).

Indeed, suppose \(x+x^2>r^2\), which is equivalent to
\(\lambda_0+\kappa_0>\Lambda\). Taking the common height
rate sufficiently close to \(\lambda_0\) in (AU3), the
rate of \(F_{\rm err}\) is strictly larger than \(\Lambda\).
All branches in (AU3) exceed \(\Lambda\) at the limiting
rate, since \(a(2\lambda_0)>\kappa_0\) and
\(\lambda_0>\kappa_0\).
The first Hardy step for \(H=V-J\) therefore makes \(H\)
supported on a left half-line. Hence on the forward tail
\(T_rH=F_{\rm err}-\mathcal E\).
All terminal terms of \(\mathcal E\) of rate strictly
above \(\Lambda\) now belong to a faster Gaussian error.
A terminal term at rate \(\Lambda\) is
\(r^{N_j+1}\eta_j T_r1_{(-\infty,a_j]}\) for some \(a_j\),
because its CDF slope divided by its standard deviation
is \(r/\sigma\). It follows that the bounded function
\[
       H+\sum_{\text{critical terminals }j}
                    r^{N_j+1}\eta_j1_{(-\infty,a_j]}
\]
has variance-\(\sigma^2\) convolution of rate strictly
above \(1/(2\sigma^2)\) on the positive half-line.
Applying critical one-sided Hardy after every real shift
forces that function to be zero.

For either extreme arm of the open semicircle, all entering
neighbors lie on the same side, so their direct signs
\(\eta_j\) are a common \(\eta\). The formula for \(V=J+H\)
then gives a limit at \(-\infty\) which is zero or of sign
\(-\eta\). But the reverse direction of the extreme arm
lies strictly in the major angular sector of its incident
spectator label. This is the lower reference label when
\(\eta=+1\) and the upper reference label when \(\eta=-1\).
The exact conditional source projection therefore satisfies
\(V(s)\to\eta-\beta\) as \(s\to-\infty\). The transverse
Gaussian has fixed mean \(rb\), so it lies in the reverse
pure-label sector with probability tending to one. Since
\(|\beta|<1\), this limit has strict sign \(\eta\), a
contradiction.
[One-sided Gaussian deconvolution](GAUSSIAN_EXTERIOR_LINEAR_RIGIDITY.md#one-sided-gaussian-deconvolution)
