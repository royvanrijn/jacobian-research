# Terminal-seven separation for curvature at least 73/10

For distinct ordinary exponent columns \(a,b\) in arbitrary dimension,
put \(S=\sum_{i:a_i\ne b_i}(b_i-a_i)^2/\max(a_i,b_i)\).
If \(S\ge73/10\), the factorial moments \(Q_6,Q_7\) have disjoint
root arguments after independent positive rescalings. In particular
they are coprime over every characteristic-zero field. The bound
allows any pattern of private rows and any total-degree difference.

Use the [fixed-terminal coefficient estimate](MULTIVARIATE_FIXED_TERMINAL_CURVATURE.md#large-curvature-and-the-endpoint-radii).
At order \(j\), normalize both endpoint coefficients to one. Each
interior coefficient is at most
\[
    \binom jk \exp\!\left[-\frac{73\,k(j-k)}{20j}\right].
\]
The [exact replay](verify_multivariate_terminal_seven_high_curvature.py)
replaces this exponential by the reciprocal of the positive
Taylor sum for its inverse through degree 32. Each interior
coefficient ranges from zero to this rational upper bound. The
constant and leading terms strictly dominate on the respective
circles of radii \(1/8\) and \(8\), so all roots remain between those
circles for every array in either coefficient box.

The five ray cosines in order of increasing angle are
\[
       \frac{2167}{2500},\quad \frac{2969}{5000},\quad
       \frac{14}{125},\quad-\frac{3303}{10000},\quad-\frac{379}{500}.
\]
For each ray and both orders, use the
[rational Bernstein projection test](MULTIVARIATE_TERMINAL_SIX_JOINT_RAY_BANDS.md#moving-rays-and-exact-replay)
to exclude a zero at every radius in the verified annulus. A straight
segment in each box to its anchor \(1+t^j\) preserves the counts in
the sectors bounded by the rays and their conjugates.

For the anchors, let \(x=\cos\theta\). For \(0<\theta<\pi\), the sign
of \(e^{-3i\theta}(1+e^{6i\theta})\) is that of \(8x^3-6x\).
The sign of \(e^{-7i\theta/2}(1+e^{7i\theta})\) is that of
\(8x^3-4x^2-4x+1\), after dropping the positive factor
\(2\cos(\theta/2)\). Sign changes at the rays and at the two ends
give the upper-sector counts \(7,6,7,6,7,6\), with numbers labeling
orders and one root per sector. The seventh-order anchor has one
remaining real root at \(-1\).

For the actual polynomials, the first and last sectors are counted
together with their lower reflections. The last joined sector
contains one septic root and two sextic roots. The single septic
root there is real by conjugation. The sextic moment is positive
for every real coefficient ratio, by integrating the sixth power
of the real binomial against the product exponential measure.
Thus the two sextic roots in the last joined sector are nonreal.
All root arguments are disjoint. A rational Bézout identity gives
the claimed coprimality in characteristic zero.
