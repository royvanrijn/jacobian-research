# Quantitative angular gaps for the physical Gaussian clocks

Let
\[
 R_j(z;\tau)=\sum_{k=0}^j\binom jk e^{-\tau k(j-k)}z^k,\qquad
 Q_j(z;t)=R_j(z;t/j).
\]
For \(N=n+1\), \(n\ge2\), and \(0<t\le2\), compare the lower raw
time \(t/n\) with any upper raw time
\[
                         t/N\le\tau'\le t/n .
\]
The consecutive roots strictly interlace in the symmetric angular
ordering. Every two root arguments from opposite degrees have
circular distance at least
\[
                  g_n(t)=\frac{(\sqrt t/12)^N}{1000nN}.            \tag{1}
\]
Each individual root set has circular mesh at least
\[
        \gamma(t)/N,\qquad \gamma(t)=\min\{\pi,\sqrt{32t/\pi}\}.
\]
The mesh is from the
[Gaussian angular response estimate](MULTIVARIATE_GAUSSIAN_LAPLACE_ANGLE_RESPONSE.md#angular-gaps-and-velocities-of-the-reference)
at normalized strength at least \(t\).

In particular, the upper normalized strength \(N\tau'\) may
exceed \(2\). The subcritical and critical theorem is needed
only at the slowest upper time \(t/N\). Moving the upper clock
to \(\tau'\) increases its outer heights from the negative
axis, while its inner gaps are bounded below by the gaps at
the equal raw clock.

## The sharp outer gap

At the sharp pair \(t/n,t/N\), use the
[exact saddle phase](MULTIVARIATE_GAUSSIAN_EXACT_SADDLE_PHASE.md).
The [phase reserve calculation](MULTIVARIATE_GAUSSIAN_SUPERCRITICAL_COLLAR.md#a-phase-reserve-at-the-critical-time)
applies unchanged for every \(0<t\le2\).
In detail, the expanded saddle-angle table for
\(\phi\le\pi/4+1/1000\) uses only \(0<t\le2\) and
\(q=\sinh(x_*)/x_*<5/2\), as does the rest of the
[left-width proof](MULTIVARIATE_GAUSSIAN_SUBCRITICAL_SHARP_CLOCK_INTERLACING.md).
For larger \(\phi\), the adjacent-phase sector test applies.
Thus at every positive lower root angle,
\[
             \frac1{500n}<\psi_N-\psi_n<7\pi/8 ,
 \qquad |\cos\psi_N|>\frac1{1000n}.
\]
The right contour integrand has the one-factor modulus
\[
 |G(x+iy)|=e^{-(x^2-y^2)/(4t)}
                |\cosh((x+i(y+\theta))/2)|.
\]
On \(\sqrt t\le x\le2\sqrt t\) this exceeds
\(e^{-1}\sinh(\sqrt t/2)>\sqrt t/6\).
Consequently its positive horizontal integral \(a_N\) exceeds
\(\sqrt t\,(\sqrt t/6)^N\).
The common factor of the real centered Gaussian integral is
\(2^N e^{-Nt/4}\sqrt{N/(4\pi t)}\). Its product with \(a_N\)
is greater than \((\sqrt t/6)^N/4\).
Multiplying by \(2|\cos\psi_N|\) proves that at a positive lower
root angle the absolute centered upper polynomial is greater
than \((\sqrt t/6)^N/(2000n)\).
Its derivative with respect to angle is at most \((N/2)2^N\),
term by term. The mean value theorem gives (1) at the sharp
pair for every positive lower root, and conjugation covers
negative lower angles.

If the lower polynomial has the structural root \(-1\), the
even upper root pair nearest \(-1\) is conjugate and separated
by at least its own circular mesh \(\gamma(t)/N\).
Each of its two distances from \(-1\) is at least
\(\gamma(t)/(2N)>g_n(t)\), using \(\gamma(t)\ge\sqrt t\).

## A polynomial inner gap at the equal raw clock

At the common raw time \(\tau=t/n\), put
\(q=e^{-\tau}\) and \(u=\tanh(\tau/2)\).
An exact derivative identity is
\[
 \frac{d}{dz}R_N(z;\tau)
       =Ne^{-n\tau}R_n(ze^\tau;\tau).
\]
It follows by shifting the binomial index, since
\((k+1)(N-k-1)=k(n-k)+n-k\).
If \(z=e^{i\theta}\) is a lower root, the upper derivative vanishes
at \(qz\) inside the unit circle.
Taking the real part of the logarithmic derivative there gives
\[
   \sum_{\alpha} P_q(\theta-\alpha)=N,\qquad
   P_q(d)=\frac{1-q^2}{1-2q\cos d+q^2}
         =\frac{u}{u^2+(1-u^2)\sin^2(d/2)} ,
\]
where \(\alpha\) ranges over all upper arguments at the equal clock.
Indeed \(\Re(qz/(qz-e^{i\alpha}))=(1-P_q(\theta-\alpha))/2\).
Every Poisson summand is at least \(u\).
For any one upper argument, the others therefore give
\(P_q(\theta-\alpha)\le N-nu\), and it follows exactly that
\[
 \sin^2\!\frac{\theta-\alpha}{2}
    \ \ge\ \frac{u(1-nu)}{(N-nu)(1+u)}.                           \tag{2}
\]
Here \(u<1/n\).

For a numerical-free simplification, \(u\ge\tau/3\) and
\(1-nu\ge1/(4n^2)\). For \(0\le v\le1/2\),
\(\tanh v\ge v-v^3/3\) follows by integrating
\(1-\tanh^2 v\ge1-v^2\).
In particular \(\tanh v\ge(11/12)v\), which in the same
integral gives \(\tanh v\le v-v^3/4\).
Apply the upper bound at \(v=1/n\) and use
\(\tau/2\le1/n\). Also
\((N-nu)(1+u)\le3N/2\).
Thus each circular distance at the equal clock is at least
\[
                        2\sqrt{\frac{t}{18 n^3 N}} .              \tag{3}
\]
This is larger than (1).

Finally, a positive lower height from \(-1\) stays between its
outer and inner upper neighbors throughout the physical interval.
The outer distance can only increase from the sharp clock.
The inner distance can only decrease to its value at the equal
clock. A structural lower zero has its smallest neighboring
upper distances at the sharp clock, and structural upper zeros
are fixed. Comparisons with these neighboring upper roots
give the bound for all upper roots in circular distance.
Equations (1) and (3) therefore apply to the physical pair
without replacing its upper clock by the lower normalized
strength.
