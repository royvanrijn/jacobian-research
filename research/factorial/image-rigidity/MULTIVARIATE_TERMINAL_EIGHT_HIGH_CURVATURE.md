# Terminal-eight separation for curvature at least 15/2

For distinct ordinary exponent columns \(a,b\), let
\(S=\sum_{i:a_i\ne b_i}(b_i-a_i)^2/\max(a_i,b_i)\).
If \(S\ge15/2\), the order-seven and order-eight factorial
moments have disjoint complex root arguments. Hence they are
coprime over every characteristic-zero field. The result
places no restriction on the signs of the gaps, the private
degrees, or the dimension.

Use the [total-curvature coefficient bound](MULTIVARIATE_FIXED_TERMINAL_CURVATURE.md#large-curvature-and-the-endpoint-radii)
with \(j=7,8\). Normalize both endpoints to one. For each interior
coefficient the bound is \(\binom jk e^{-15 k(j-k)/(4j)}\).
Replace the exponential by the reciprocal of the positive Taylor
lower bound on its inverse. The resulting rational coefficient
box contains the whole segment from the actual array to \(1+t^j\).
The constant and leading terms strictly dominate on the respective
circles of radii \(1/8\) and \(8\) throughout these boxes.

The [exact checker](verify_multivariate_terminal_eight_high_curvature.py)
uses the following ray cosines, in increasing angular order:
\[
 \frac{4511}{5000},\quad \frac{687}{1000},\quad\frac{3037}{10000},
 \quad-\frac{841}{10000},\quad-\frac{5081}{10000},\quad-\frac{799}{1000}.
\]
It applies the
[rational Bernstein tests](MULTIVARIATE_TERMINAL_SIX_JOINT_RAY_BANDS.md#moving-rays-and-exact-replay)
and proves strict nonvanishing for both orders throughout the
annulus on all six rays in 82 radial cells.

For \(x=\cos\theta\), the sign of
\(e^{-7i\theta/2}(1+e^{7i\theta})\) for \(0<\theta<\pi\) is the
sign of \(8x^3-4x^2-4x+1\) after dividing out \(2\cos(\theta/2)\).
The sign of \(e^{-4i\theta}(1+e^{8i\theta})\) is the sign of
\(8x^4-8x^2+1\). The checked signs at the rays and the two
ends give one upper anchor root per sector, in the alternating
order \(8,7,8,7,8,7,8\). The remaining seventh-order root is
negative real. Root counts continue along the coefficient
segments, joining the first and last upper sectors to their
lower reflections. In the last joined sector the single septic
root is real by conjugacy. The two eighth-order roots are
nonreal since the eighth factorial moment is positive for
every real coefficient ratio. This proves angular separation
and coprimality by a rational Bézout identity.
