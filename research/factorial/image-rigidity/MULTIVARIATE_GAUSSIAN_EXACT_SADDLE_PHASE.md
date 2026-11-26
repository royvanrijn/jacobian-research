# Exact saddle phase and microscopic Gaussian root bounds

For \(t>0\) let
\[
 Q_j(z;t)=\sum_{k=0}^j\binom jk e^{-tk(1-k/j)}z^k .
\]
This is the Gaussian polynomial at raw heat time \(t/j\).
The roots are simple on the unit circle by the
[Gaussian heat theorem](MULTIVARIATE_GAUSSIAN_REFERENCE_ANGULAR_INTERLACING.md).
Write \(0<\theta_{j,1}<\cdots<\theta_{j,\lfloor j/2\rfloor}<\pi\)
for their positive arguments, excluding the structural root at
\(\pi\) in odd degree. Define
\[
 \theta_*=
 \begin{cases}
  2\arccos\sqrt{t/2}-\sqrt{t(2-t)},&t<2,\\
  0,&t\ge2.
 \end{cases}
\]
Every \(\theta_{j,i}\) exceeds \(\theta_*\), by the
[finite-degree zero-free sector](MULTIVARIATE_GAUSSIAN_SHARP_CLOCK_BULK_CONTACT.md#a-zero-free-sector-at-every-subcritical-degree)
below criticality and coefficient positivity otherwise.

**Theorem.** There is a real increasing saddle phase
\(\phi_t(\theta)\) from \(0\) to \(\pi/2\) on
\((\theta_*,\pi)\), independent of the degree, for which
\[
       (i-\tfrac12)\pi <j\phi_t(\theta_{j,i})<i\pi.                  \tag{1}
\]
For every \(t\) there are positive explicit constants
\(c_t,C_t\) such that
\[
 c_t\left(\frac{i-1/2}{j}\right)^{\!\beta}
    <\theta_{j,i}-\theta_*
    <C_t\left(\frac{i}{j}\right)^{\!\beta},\qquad
 \beta=\begin{cases}2/3,&t<2,\\3/4,&t=2,\\1,&t>2.\end{cases}        \tag{2}
\]
These statements hold for every degree and index. The proof also
gives uniform bounds through the critical crossover. In addition, for
\(0<t\le6\) the consecutive root angles strictly interlace at each
lower root with \(\phi_t(\theta_{j,i})\ge\pi/4\). This last statement
concerns the sharp physical pair \(t/j,t/(j+1)\).

## A global contour through the right saddle

Fix \(0<\theta<\pi\), and abbreviate
\[
 G(z)=e^{-z^2/(4t)}\cosh((z+i\theta)/2),\quad
 I(x,y)=\arctan\!\left(\tanh(x/2)\tan((\theta+y)/2)\right)
                 -xy/(2t).
\]
For \(x>0\) and \(-\theta<y<\pi-\theta\), \(I\) is the imaginary
part of a branch of \(\log G(x+iy)\). Direct differentiation gives
\[
 I_x=-\frac y{2t}+\frac{\sin(\theta+y)}{2(\cosh x+\cos(\theta+y))},
 \qquad
 I_y=-\frac x{2t}+\frac{\sinh x}{2(\cosh x+\cos(\theta+y))}.
\]
Thus \(I_{yy}>0\) and \(I_{xx}=-I_{yy}<0\) in the open strip.

Let \(x_c>0\) solve \(x_c\tanh(x_c/2)=t\).
If \(t>2\) let \(x_0>0\) solve \(x_0\coth(x_0/2)=t\).
Otherwise set \(x_0=0\).
The minimum of \(I(x,\cdot)\) on
\([-\theta,\pi-\theta]\) is at the lower endpoint for \(0<x\le x_0\)
and at the upper endpoint for \(x\ge x_c\). In the remaining
interval it is at \(y=u(x)-\theta\), where
\[
 u(x)=\arccos\left(\frac{t\sinh x}{x}-\cosh x\right),\qquad
 \Theta(x)=u(x)-\frac{x}{\sinh x}\sin u(x).                        \tag{3}
\]
These assertions follow from the sign of \(I_y\) at the two
endpoints and its strict increase in \(y\).

The argument of \(\arccos\) in (3) is strictly decreasing there.
In fact \(t\le x\coth(x/2)\) on this interval, and its derivative
is at most
\[
 \frac{\cosh x+1}{x\sinh x}\big(x-\sinh x\big)<0 .
\]
Consequently \(\Theta\) strictly increases, since \(x/\sinh x\)
strictly decreases and \(u\) strictly increases.
Its lower endpoint is \(\theta_*\) and its upper endpoint is \(\pi\).
The minimum phase \(m_\theta(x)\) has derivative
\((\theta-\Theta(x))/(2t)\) in the interior case.
At the two endpoint regimes it has derivative
\(\theta/(2t)\) and \(-(\pi-\theta)/(2t)\), respectively.

When \(\theta>\theta_*\), the unique maximum of \(m_\theta\)
is at \(x_*\in(x_0,x_c)\) with \(\Theta(x_*)=\theta\).
Let \(y_*=u(x_*)-\theta>0\) and
\(\phi=\phi_t(\theta)=m_\theta(x_*)\).
The two partial derivatives of \(I\) vanish at this saddle.
Also
\[
 0<\phi<\theta/2,\qquad
                      \phi_t'(\theta)=\frac{x_*}{2t}.             \tag{4}
\]
The first inequality uses \(m_\theta(0)=0\) by continuity, and the
upper bound uses \(m_\theta(x_*)\le I(x_*,0)<\theta/2\).
For the derivative, the formula for \(I(x,u(x)-\theta)\)
makes its partial \(\theta\)-derivative equal to \(x/(2t)\) and
the \(x\)-derivative vanishes at the maximum. At the upper endpoint
\(x_*\to x_c\), \(y_*\to0\), so \(\phi_t(\pi-)=\pi/2\).
At the lower endpoint the phase tends to zero.

On \(0<x<x_*\), take the solution of \(I(x,y)=\phi\) above
the minimizing \(y\). It exists because the upper endpoint value
\(\pi/2-x(\pi-\theta)/(2t)\) decreases in \(x\) and exceeds
\(\phi\) at \(x_*\). On \(x>x_*\), use the solution below
the minimum. It exists because the lower endpoint value
\(x\theta/(2t)\) increases and exceeds \(\phi\) at \(x_*\).
The two branches join smoothly. To check this, write near the saddle
\[
 I(x,y)=m_\theta(x)+(y-u(x)+\theta)^2 B(x,y),\qquad
 \phi-m_\theta(x)=(x-x_*)^2 C(x)
\]
with positive analytic \(B,C\). The selected branch solves
\(y-u(x)+\theta=-(x-x_*)\sqrt{C(x)/B(x,y)}\), which has a smooth
solution by the implicit-function theorem.

This contour \(z=x+i y(x)\) runs from \(i(\pi-\theta)\) to
the positive real direction at infinity, and \(y\) strictly
decreases. Indeed, near zero the level equation gives
\[
 \frac{\pi-\theta-y(x)}{x}\longrightarrow\cot\phi,\qquad I_x>0 .
\]
For the first limit use the arctangent equation for \(I=\phi\)
and \(\tanh(x/2)\sim x/2\). Substituting the same limit in the
formula for \(I_x\) gives the claimed sign.
At infinity \(y(x)=t(\theta-2\phi+O(e^{-x}))/(x-t)\), so \(I_x<0\).
There are no other zeros of \(I_x\) on the selected contour
besides the saddle. For a fixed \(y\), such a zero maximizes
\(I(\cdot,y)\) by strict concavity. If \(y\ne y_*\), this
maximum exceeds \(\phi\), since \(I(x_*,y)>I(x_*,y_*)\) by strict
convexity in \(y\). If \(y=y_*\), the unique maximum is at \(x_*\).
Since \(y'=-I_x/I_y\), with \(I_y>0\) on the left branch and
\(I_y<0\) on the right, \(y'<0\) on both branches.
In particular \(0<y(x)<\pi-\theta\).

## The finite-degree phase and its labels

By the centered Gaussian integral, the polynomial times the factor
\(e^{-ij\theta/2}\) is proportional by a positive constant to
\(\int_{\mathbb R} G(z)^j\,dz\).
Deform this line to the right contour and its reflection across the
imaginary axis. They join at a zero of \(G\), and Gaussian decay
makes the vertical closing sides vanish. The integrand is entire,
so the deformation is valid for every integer \(j\).
On the right contour \(G(z)^j=e^{ij\phi}|G(z)|^j\).
By conjugation, the full integral becomes
\[
 2\Re\big(e^{ij\phi}(a_j-i b_j)\big),\qquad
 a_j=\int_0^\infty |G(x+iy)|^j\,dx,\quad
 b_j=\int_0^\infty |G(x+iy)|^j(-y')\,dx.                           \tag{5}
\]
Both integrals are strictly positive and finite.

Put \(\chi_j=\arg(a_j-i b_j)\in(-\pi/2,0)\).
The continuous phase \(\psi_j=j\phi+\chi_j\) determines the zeros
through its half-integer multiples of \(\pi\). For continuity
one can equivalently express \(e^{ij\phi}(a_j-i b_j)\)
as the integral on the positive real axis minus the integral on
the vertical segment from zero to \(i(\pi-\theta)\).
Near \(\theta=\pi\), \(b_j/a_j\to0\), since the total vertical
variation is \(\pi-\theta\) and \(a_j\) is bounded below by its
integral over \(1\le x\le2\). At the lower endpoint \(\phi\to0\),
and \(\psi_j>-\pi/2\) throughout. Thus continuity forces the
phase through all levels \((i-\tfrac12)\pi\) for
\(1\le i\le\lfloor j/2\rfloor\).
The exact root count assigns precisely one root to each such
level. Since \(\chi_j\in(-\pi/2,0)\), this proves (1).

### Explicit constants for the edge scales

For \(t<2\), put \(u_c=\arccos(t-1)\), \(M=\cosh x_c\), and
\[
 L=\tfrac12(1-t/2)\min\{1,2-t\},\qquad
 K=M\left(\frac{\pi}{2\min\{\sin u_c,\cos(u_c/2)\}}+\frac16\right).
\]
Then \(L x^2\le\Theta(x)-\theta_*\le Kx^2\).
For the lower bound, differentiating \(u\) in (3) gives
\(u'(x)\ge(1-t/2)x\).
This uses \(x\cosh x-\sinh x\le x^2\sinh x/2\), obtained by
integrating \(v\sinh v\). Also
\(1-(x/\sinh x)\cos u\ge\min\{1,2-t\}\).
For the upper bound,
\(\cos u_c-\cos u\le M x^2/2\).
The sine factorization of this cosine difference gives
\[
 u-u_c\le\frac{\pi(\cos u_c-\cos u)}
                       {2\min\{\sin u_c,\cos(u_c/2)\}}.
\]
Then \((u-\sin u)-(u_c-\sin u_c)\le2(u-u_c)\) and
\(1-x/\sinh x\le Mx^2/6\) give the upper estimate.

At \(t=2\), \(u_c=0\) and \(Lx^3\le\Theta(x)\le Kx^3\)
holds with
\[
 l=1/\sqrt3,\quad U=\pi\sqrt{\cosh x_c}/2,\qquad
 L=\frac{2l^3}{3\pi^2},\quad K=\frac{U^3+U\cosh x_c}{6}.
\]
Indeed \(x^2/6\le1-\cos u\le x^2\cosh(x_c)/2\) by the
positive Taylor series of \(1+\cosh x-2\sinh(x)/x\).
This yields \(l x\le u(x)\le Ux\).
Now use \(2u^3/(3\pi^2)\le u-\sin u\le u^3/6\) in (3)
and the same bound on \(1-x/\sinh x\).

In either case, if \(d=2\) or \(3\) is the respective power and
\(\Delta=\theta-\theta_*\), equation (4) gives
\[
 \frac{d\Delta^{1+1/d}}{2t(d+1)K^{1/d}}
       \le\phi_t(\theta)
       \le\frac{d\Delta^{1+1/d}}{2t(d+1)L^{1/d}}.
\]
For \(t>2\), integrate
\(x_0/(2t)<\phi_t'(\theta)<x_c/(2t)\).
These bounds and (1) give (2).

### Uniform critical crossover

For \(3/2\le t\le5/2\) the root bounds can be made uniform as
\(t\) approaches \(2\). With \(q=i/j\), the symbol \(\asymp\)
below denotes two-sided bounds by absolute positive constants
for every degree and index:
\[
 \begin{array}{ll}
 t\le2:&
 \displaystyle \theta_{j,i}-\theta_*
        \asymp \max\!\{(2-t)^{1/6}q^{2/3},\,q^{3/4}\},\\[3pt]
 t\ge2:&
 \displaystyle \theta_{j,i}
        \asymp \min\!\{q/\sqrt{t-2},\,q^{3/4}\}.
 \end{array}                                                     \tag{7}
\]
At \(t=2\) omit the infinite entry in the minimum.

To verify this, put \(C(x)=t\sinh(x)/x-\cosh x\) and
\(C(0)=t-1\). Its derivative has the exact expansion
\[
 -C'(x)=\sum_{k\ge1}\left(1-\frac{t}{2k+1}\right)
                       \frac{x^{2k-1}}{(2k-1)!}.
\]
On \(0<x<x_c\), this is between \(x/6\) and \(x\cosh4\),
as \(x_c<4\) in the stated time interval.
For \(t\le2\), it gives
\(C(0)-C(x)\asymp x^2\), and hence
\(u\asymp\sqrt{2-t+x^2}\) by \(1-\cos u\asymp u^2\)
on \([0,\pi]\).
Since \(u_c\le\pi/3\), the sine factorization as above gives
\(u-u_c\asymp x^2/u\).
Thus
\((u-\sin u)-(u_c-\sin u_c)\asymp u^3-u_c^3\asymp x^2u\).
The remaining nonnegative term \((1-x/\sinh x)\sin u\) is
at most a constant times \(x^2u\).
We obtain
\(\Theta(x)-\theta_*\asymp x^2\sqrt{2-t+x^2}\).

For \(t>2\), the same derivative bound gives
\(x_0^2\asymp t-2\) and \(1-C(x)\asymp x^2-x_0^2\) above \(x_0\).
Consequently \(u\asymp\sqrt{x^2-x_0^2}\).
Now \(\Theta(x)\asymp x^2u\).
For the lower bound, if \(u\le\pi/2\) then
\((1-x/\sinh x)\sin u\) is bounded below by a
positive constant times \(x^2u\).
If \(u\ge\pi/2\), use \(u-\sin u\asymp u^3\) and \(x\le4\).
The upper bound follows from \(u\le\text{const}\cdot x\).

In the subcritical case invert the estimate for \(\Theta\) to get
\(x_*\asymp\min\{(2-t)^{-1/4}\Delta^{1/2},\Delta^{1/3}\}\)
for \(\Delta=\theta-\theta_*\).
Integrating (4), using the last half of the interval for the
lower bound, gives
\(\phi\asymp\min\{(2-t)^{-1/4}\Delta^{3/2},\Delta^{4/3}\}\).
In the supercritical case,
\(\theta\asymp x_*^2\sqrt{x_*^2-x_0^2}\) gives
\(x_*\asymp x_0+\theta^{1/3}\), by splitting at \(x_*=2x_0\).
Then \(\phi\asymp \sqrt{t-2}\,\theta+\theta^{4/3}\).
Inverting these bounds and using (1) proves (7).

## A strict adjacent phase increment for \(\phi\ge\pi/4\)

Assume \(0<t\le6\). Put \(A=\pi-\theta\) and write a point on the
right contour as \(z=w+iA\), \(w=x-is\), \(s=A-y(x)>0\).
In polar coordinates \(w=R e^{-i\alpha}\),
\(0<\alpha=\arctan(s/x)<\pi/2-\phi\).
To see the latter inequality, use
\(\cosh((z+i\theta)/2)=i\sinh(w/2)\) and the
[hyperbolic sine product](https://dlmf.nist.gov/4.36.E1).
Along the contour,
\[
 \phi =\frac\pi2-\alpha-\frac{x(A-s)}{2t}
       -\sum_{k\ge1}
          \arctan\frac{R^2\sin(2\alpha)}
                         {4\pi^2 k^2+R^2\cos(2\alpha)}.            \tag{6}
\]
The product denominators have positive real part since \(s<A<\pi\).

For \(\phi\ge\pi/4\) all the contour rays satisfy \(\alpha\le\pi/4\).
At fixed such \(\alpha>0\), the imaginary phase as a function of
\(R\) is strictly convex. In fact the quadratic contribution in
(6) has second derivative \(\sin(2\alpha)/(2t)\).
For \(K_k=4\pi^2 k^2\), the second derivative of the \(k\)-th
arctangent is
\[
 \frac{2K_k\sin(2\alpha)
       [K_k^2-2K_k R^2\cos(2\alpha)-3R^4]}
      {[K_k^2+2K_k R^2\cos(2\alpha)+R^4]^2}
 < \frac{2\sin(2\alpha)}{K_k}
 \quad(R>0).
\]
The series differentiates locally uniformly. Its sum is less than
\(\sin(2\alpha)/12\), giving strict convexity also at \(t=6\).

Let \(\alpha_*\) be the saddle's ray. For every \(\alpha\ne\alpha_*\)
the ray meets the horizontal line \(y=y_*\) at a point whose
phase is less than \(\phi\), by strict concavity of \(I(\cdot,y_*)\).
Thus on that ray a stationary point of the radial phase
cannot lie at level \(\phi\). For \(\alpha=\alpha_*\) the saddle
is the unique point at level \(\phi\) by radial strict convexity.
Therefore the polar angle \(\alpha(x)\) has no stationary point
on the contour except possibly the saddle, since
\[
 \frac{d\alpha}{dx}
    =\frac{x I_x - s I_y}{(x^2+s^2) I_y}
\]
away from that point and its numerator is the radial phase
derivative times \(R\). It strictly decreases on both sides:
near zero (6) gives
\(\alpha=\pi/2-\phi-Ax/(2t)+O(x^2)\), and at infinity
\(\alpha\) tends to zero through positive values.
It follows that \(0<-y'=s'\le s/x<\cot\phi\) along the graph,
including at the saddle by continuity. From (5) we obtain
\[
 \phi-\pi/2<\chi_j<0,\qquad
 2\phi-\pi/2<\psi_{j+1}-\psi_j<\pi/2 .
\]
The increment is strictly positive. Since each phase level has
a unique root, the corresponding angles are strictly interlaced.
For \(\phi<\pi/4\) the exact phase representation and the root
bounds still apply, but this estimate gives no strict lower sign
for the adjacent increment.



There is a quantitative gap away from the rank cutoff.
When \(\phi_t(\theta_{j,i})\ge\pi/4\), the improved sector for
\(\chi_j\) gives
\(\phi_t(\theta_{j,i})<i\pi/(j+1)\).
Let \(N=n+1\), \(0<\varepsilon<1/4\),
\(n\ge1/\varepsilon\), and \(i\ge(1/4+\varepsilon)n\).
The lower bound in (1) makes both the degree-\(n\) and
degree-\(N\) phases at index \(i\) at least \(\pi/4\).
Using the improved upper bound and then
\(\phi_t'(\theta)\le x_c/(2t)\) gives
\[
 \theta_{n,i}-\theta_{N,i}
    >\frac{2t}{x_c}
        \left(\frac{(i-1/2)\pi}{n}-\frac{i\pi}{n+2}\right)
    \ge\frac{2\pi t\varepsilon}{x_c(n+2)}.
\]
The same upper bound on the degree-\(n\) phase at index \(i\)
and the lower bound in (1) at the next upper root give
\[
 \theta_{N,i+1}-\theta_{n,i}>\frac{\pi t}{x_c(n+1)}
\]
whenever that next positive root exists.

For \(0<t\le2\), the decrease of \(\alpha(x)\) and the same
slope bound \(0<-y'<\cot\phi\) hold at every positive saddle
phase. To verify the radial convexity in the remaining rays
\(\alpha>\pi/4\), put \(H(w)=\sinh(w/2)/(w/2)\).
Its reciprocal square has the exact positive integral
\[
  H(w)^{-2}=\mathbb E e^{-V w^2},\qquad
  V=\sum_{k\ge1}V_k,\quad V_k\sim{\rm Gamma}(2,4\pi^2 k^2),
\]
with independent variables and the second parameter denoting rate.
This follows from the same hyperbolic sine product and the gamma
Laplace transform. The sum has finite expectation, and the integral
is absolutely convergent for \(|\Im w|\le\pi\).
For \(w=x-is\) with \(x>0\) and \(0<s<\pi\), it gives
\[
 \Im H(w)^{-2}
   \le 2xs\,\mathbb E[V e^{V(s^2-x^2)}]
   <2xs\,\mathbb E[V e^{V\pi^2}]
   =xs/2.
\]
The last expectation is \(1/4\), by differentiating
\(H(-is)^{-2}=[(s/2)/\sin(s/2)]^2\) at \(s=\pi\) and dividing
by \(2\pi\).
The radial phase therefore satisfies
\[
 I_{RR}
   =\frac{\sin(2\alpha)}{2t}
        - \frac{\Im H(w)^{-2}}{R^2}
   >\sin(2\alpha)\left(\frac1{2t}-\frac14\right)\ge0 .
\]
The same strict-convexity argument then applies to the whole right
contour. Its resulting sector still has width \(\pi/2-\phi\),
which by itself leaves the smaller-phase adjacent sign open.
