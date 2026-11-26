# Polynomial angular gaps at the Gaussian physical clocks

Write
\[
 R_j(z;\tau)=\sum_{k=0}^j\binom jk e^{-\tau k(j-k)}z^k,\qquad
 Q_j(z;t)=R_j(z;t/j).
\]
For every \(n\ge1\), \(N=n+1\), \(t>0\), and
\(\tau'\in[t/N,t/n]\), each root argument of \(R_n(\cdot;t/n)\)
has circular distance from every root argument of \(R_N(\cdot;\tau')\)
greater than
\[
                       \frac{t}{36 nN(1+t)} .                    \tag{1}
\]
For the sharp outer comparison we use a derivative of the exact
contour phase. For the equal-time inner comparison at \(t\ge2\)
we use the root motion equation and the individual angular mesh.

## Sharp outer gap

In the [exact saddle phase](MULTIVARIATE_GAUSSIAN_EXACT_SADDLE_PHASE.md)
notation, fix \(\theta\) between the zero-free edge and \(\pi\).
Put \(A=\pi-\theta\) and
\(f(z)=-z^2/(4t)+\log\cosh((z+i\theta)/2)\).
For an integer \(m\ge1\) the complex half integral may be taken
along \(z=x+iA\), \(0<x<\infty\):
\[
 B_m(\theta)
 =i^m\int_0^\infty e^{-m(x+iA)^2/(4t)}\sinh^m(x/2)\,dx .
\]
Along this line the hyperbolic factor is independent of \(\theta\).
Differentiating under the integral and deforming both integrals back
to the right saddle contour yields
\[
 B_m'(\theta)=\frac{im}{2t}\int_{\mathrm{right}}z e^{mf(z)}\,dz .
 \tag{2}
\]
The Gaussian factor controls the horizontal ends uniformly on
compact angle intervals, and \(e^{mf}\) is entire. The phase
\(\psi_m\) is a continuous argument of \(B_m\).
On the right contour \(z=x+iy(x)\), the phase of \(e^{mf}\)
is constant and \(dz=(1+iy')\,dx\) has positive real and
nonpositive imaginary parts. In particular
\[
              \int_{\mathrm{right}} |e^{mf(z)}\,dz|
                        \le\sqrt2\,|B_m| .
 \tag{3}
\]

Choose \(a=t+2\). Then \(a>x_*\), since \(x_*<x_c\) where
\(x_c\tanh(x_c/2)=t\), and \(x_c-t=2x_c/(e^{x_c}+1)<1\).
For \(x\ge a\),
\[
 2I_y=-x/t+\frac{\sinh x}{\cosh x+\cos(\theta+y)}
      \le-x/t+\coth(x/2)\le-1/t .
\]
The last bound uses \(e^x-1\ge2t\) and \(x\ge t+2\).
Thus \(U'=I_y(1+y'^2)\le-1/(2t)\) there, where \(U=\Re f\).
The [right-branch contour bound](MULTIVARIATE_GAUSSIAN_SUPERCRITICAL_PHASE_RESERVE.md#dyadic-width-and-the-phase-reserve)
also gives \(|y'|\le1\) on this part.
Since \(e^{mU(x+q)}\le e^{-mq/(2t)}e^{mU(x)}\) for \(x\ge a\)
and \(q\ge0\), integration gives
\[
 \int_{\mathrm{right}} (x-a)_+ |e^{mf(z)}\,dz|
 \le \sqrt2\int_a^\infty (x-a)e^{mU(x)}\,dx
 \le \frac{2\sqrt2 t}{m}\int_a^\infty e^{mU(x)}\,dx
 \le \frac{2\sqrt2 t}{m}|B_m| .
\]
In the second step write the excess \(x-a\) as a tail integral
and apply the preceding exponential bound. In the third use that
the horizontal amplitude is the real part of \(B_m\) after rotation
by \(-m\phi\). Together with \(|z|\le a+\pi+(x-a)_+\), (2) and
(3) prove
\[
 |\psi'_m(\theta)|\le\frac{\sqrt2 m}{2t}(t+2+\pi)+\sqrt2
                    \le 6m(1+1/t).
 \tag{4}
\]
This controls the argument of the normalized half integral
without estimating the modulus of the centered polynomial.

At the positive degree-\(n\) root angle \(\theta_{n,i}\), the
[all-coupling phase reserve](MULTIVARIATE_GAUSSIAN_SUPERCRITICAL_PHASE_RESERVE.md)
gives
\(\psi_N(\theta_{n,i})-\psi_n(\theta_{n,i})
>(\pi/2-\log4)/n>1/(6n)\).
The positive sharp upper angle \(\theta_{N,i}\) lies below
\(\theta_{n,i}\) and obeys
\(\psi_N(\theta_{N,i})=\psi_n(\theta_{n,i})\) by the root labels.
By (4), \(\theta_{n,i}-\theta_{N,i}>t/[36nN(1+t)]\).

## Equal-raw inner gap

Use positive heights from the negative axis in decreasing order,
\(x_i(\tau)\) in degree \(n\) and \(y_i(\tau)\) in degree \(N\), both
at raw time \(\tau=t/n\). For \(0<t\le2\), the
[equal-raw Poisson estimate](MULTIVARIATE_GAUSSIAN_QUANTITATIVE_PHYSICAL_GAPS.md#a-polynomial-inner-gap-at-the-equal-raw-clock)
gives every cross-angle gap at least
\(2\sqrt{t/(18 n^3N)}\) when \(n\ge2\). This exceeds (1).

For \(t\ge2\), consider an inner gap between positive heights
\(\delta=x_i-y_{i+1}>0\). Let \(F_N(s)\) be the raw heat velocity
of the upper height \(y_{i+1}\) when that height is placed at \(s\)
and its reflected height is placed at \(-s\), with all other
upper heights fixed. At \(s=x_i\), the same-index cancellations
in the
[equal-raw contact argument](MULTIVARIATE_GAUSSIAN_REFERENCE_ANGULAR_INTERLACING.md#first-contact-across-adjacent-degrees)
give
\[
 \dot x_i-F_N(x_i)
 \ \ge\ \cot(x_i/2)+2\tan(x_i/2)\ >2 .
 \tag{5}
\]
In detail, the paired force \(F(x,a)=2\sin x/(\cos a-\cos x)\)
increases with \(a\) on either side of \(x\). The pair
\(y_{k+1}\le x_k\) is on the same side of \(x_i\) by interlacing.
The self-reflected terms are equal. For odd \(n\) the lost center
contributes \(\cot(x_i/2)\) to the lower-minus-upper force.
For even \(n\) the unmatched lower smallest positive height
contributes more than \(2\cot(x_i/2)\), and the upper center
subtracts only \(\cot(x_i/2)\). In either parity the extra largest
height \(y_1>x_i\) contributes negatively to the upper force,
with absolute value at least \(2\tan(x_i/2)\).

The upper angular mesh is at least \(\mu=\pi/N\) for \(t\ge2\)
by the [individual Gaussian mesh bound](MULTIVARIATE_GAUSSIAN_LAPLACE_ANGLE_RESPONSE.md#angular-gaps-and-velocities-of-the-reference),
since its normalized time is \(N\tau\ge2\).
Suppose \(\delta\le\mu/4\). While \(s\) moves from \(y_{i+1}\)
to \(x_i\), any upper argument at cyclic rank distance \(k\) from
\(y_{i+1}\) has circular separation from \(s\) at least \(k\mu/2\).
Indeed, its original separation is at least \(k\mu\), and the
motion reduces it by at most \(\delta\), or by at most \(2\delta\)
for the moving reflected root. Each ordinary cotangent force
derivative has absolute value at most
\(\frac12(2\pi/(k\mu))^2\). There are at most two arguments per
cyclic rank distance. For the reflected term \(\cot s\), the
argument changes twice as fast and adds at most \(2\pi^2/\mu^2\)
to the same bound. Summing with \(\sum_{k\ge1} k^{-2}\le2\) gives
\[
 |F_N'(s)|\le10\pi^2/\mu^2=10N^2 .
\]
Equations (5) and \(\dot y_{i+1}=F_N(y_{i+1})\) therefore imply
\(\dot\delta>2-10N^2\delta\) while \(\delta\le\pi/(4N)\).
At \(t=2\) the Poisson bound exceeds \(1/(10N^2)\). This value
is a strict barrier at every later raw time, by the derivative
inequality at a first downward crossing. Hence the equal-raw
positive inner gaps exceed \(1/(10N^2)>t/[36nN(1+t)]\).

## Full physical interval

The structural root at \(-1\) in the odd degree is at least half
the circular mesh of the even root set away from every even root,
by reflection. That even mesh is at least
\(\min\{\pi,\sqrt{32t/\pi}\}/N\) throughout the physical interval.
This proves (1) for a structural root as well, including \(n=1\).
As the upper raw time grows from \(t/N\) to \(t/n\), its positive
heights move strictly outward. The outer distance can only grow
from its sharp value, while the inner distance is bounded below
by its equal-raw value. These neighboring gaps and the structural
ones give (1) on each semicircle. For two nonstructural roots on
opposite semicircles the circular distance is at least the distance
between their reflected arguments on the same semicircle, so the
same bound applies.
