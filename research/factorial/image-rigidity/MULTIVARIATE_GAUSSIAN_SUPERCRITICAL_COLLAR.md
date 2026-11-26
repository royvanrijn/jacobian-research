# A verified supercritical collar for sharp Gaussian clocks

For
\[
 Q_j(z;t)=\sum_{k=0}^j\binom jk e^{-t k(1-k/j)}z^k
\]
and \(N=n+1\), \(n\ge2\), the sharp physical pair
\(Q_n(\cdot;t),Q_N(\cdot;t)\) strictly interlaces at every
\[
                     2\le t\le 2+\frac{8^{-N}}{4000nN}.          \tag{1}
\]
There is also one positive collar valid at every degree.
If \(n_0\ge2\) is any threshold in the
[eventual Gaussian theorem](MULTIVARIATE_DIFFUSE_CRITICAL_TRANSITION.md)
with \(s_-=4\), it suffices that
\[
                     2\le t\le
             2+\frac{8^{-(n_0+1)}}{4000n_0(n_0+1)}.             \tag{2}
\]
The threshold \(n_0\) from that theorem is unevaluated.
The finite-degree estimate (1) requires no threshold and no
checks in individual degrees.

## A phase reserve at the critical time

Use the notation of the
[exact saddle phase](MULTIVARIATE_GAUSSIAN_EXACT_SADDLE_PHASE.md)
and the
[subcritical and critical interlacing proof](MULTIVARIATE_GAUSSIAN_SUBCRITICAL_SHARP_CLOCK_INTERLACING.md).
In the latter proof the left-width estimate with constant \(25/9\)
holds also for \(\phi\le\pi/4+1/1000\). Only the saddle-angle
table changes. Here are its rational upper checks. Put
\(\xi=1/(2q_{\min})+1/1000\), \(p=\xi/(1-\xi^2/2)\).
For each indicated interval, the bound on \(u'\) is
\[
 \frac12\left[(q_{\max}+1)\frac{1+p}{1-p}
                             +\frac76 q_{\max}\right].
\]
Its exact values are
\[
 \begin{array}{c|c}
 q\text{ interval}&\text{upper bound on }u'\\
 [1,5/4]&174698911/35855952\\
 [5/4,3/2]&33672383/8297592\\
 [3/2,2]&141969928/32927973\\
 [2,5/2]&16962547/3826664
 \end{array}
\]
and each is below \(49/10\). As in the original argument, this
gives \(\sin\beta_*>1/10\) and hence the same left-width bound.
There are no other changes to the contour inequalities.

At a positive degree-\(n\) root angle, in this phase range the
lower phase increment is therefore
\[
       \psi_N-\psi_n>
          \frac{\pi/2-25/18}{n}>\frac1{10n}.
\]
The upper bound is \(\phi+\pi/2\le3\pi/4+1/1000<7\pi/8\).
Above this phase range, the exact adjacent-sector estimate gives
\(\psi_N-\psi_n>2\phi-\pi/2>1/500\) and
\(\psi_N-\psi_n<\pi/2\).
Thus every positive lower root satisfies
\[
          \frac1{500n}<\psi_N-\psi_n<\frac{7\pi}{8}.
\]
Since its own phase is an odd multiple of \(\pi/2\), the absolute
cosine of \(\psi_N\) there is greater than \(1/(1000n)\).

## A finite-degree angular gap at \(t=2\)

The exact right contour has
\(
 G(z)=e^{-z^2/8}\cosh((z+i\theta)/2).
\)
On the segment with \(1\le\Re z\le2\) its modulus exceeds \(1/4\),
using
\[
 e^{-\Re(z^2)/8}\ge e^{-1/2},\qquad
 |\cosh((z+i\theta)/2)|\ge\sinh(1/2)>1/2.
\]
Therefore \(a_N>4^{-N}\).
The real centered integral for \(Q_N\) has positive common factor
\(2^N e^{-N/2}\sqrt{N/(8\pi)}>1/4\) multiplying
\(2\Re[e^{iN\phi}(a_N-i b_N)]\).
At any positive lower root angle the absolute centered degree-\(N\)
polynomial is greater than \(4^{-N}/(2000n)\).
On the full real line, the derivative of the centered degree-\(N\)
polynomial is at most \((N/2)2^N\) in modulus, term by term in its
finite trigonometric expansion.
The mean value theorem then gives separation from every upper root
by more than
\[
                             g_n=\frac{8^{-N}}{1000nN}.           \tag{3}
\]
Conjugation covers negative lower angles.
When the lower degree is odd its structural root is \(-1\).
The even upper roots nearest to it are a conjugate pair whose
circular gap is at least \(\pi/N\) by the
[Gaussian mesh bound](MULTIVARIATE_GAUSSIAN_LAPLACE_ANGLE_RESPONSE.md#angular-gaps-and-velocities-of-the-reference)
at \(t=2\). Their distances from \(-1\) in angle are therefore
at least \(\pi/(2N)>g_n\).

The [normalized Gaussian angular speed bound](MULTIVARIATE_GAUSSIAN_LAPLACE_ANGLE_RESPONSE.md#angular-gaps-and-velocities-of-the-reference)
keeps each root's speed at most \(\pi/4\) on \(t\ge2\).
For a pair of roots the change of their angular gap is therefore
at most \((\pi/2)|t-2|\). Equation (3) makes this smaller than
the initial gap throughout (1), so the interlacing persists.
For \(n\ge n_0\), eventual coprimality and the unit-circle location
exclude every contact along \(t\ge2\). Interlacing established at
\(t=2\) therefore holds there at those degrees.
For \(n<n_0\), (2) lies within (1). This proves the uniform collar.
