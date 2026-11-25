# Terminal-eight separation for light rows with outer mean at most .10

Let all unequal rows of distinct ordinary exponent columns be shared.
Write \(\ell_i=\min(a_i,b_i)\), \(d_i=|b_i-a_i|\) and let \(u\)
be the normalized order-seven outer geometric mean. If
\(0<u\le.10\) and \(\max_i d_i/(7\ell_i+1)\le1/8\), the order-seven
and order-eight factorial moment polynomials have disjoint root
arguments, up to independent positive rescalings. Rational Bézout
gives coprimality over every characteristic-zero field.

The [checker](verify_multivariate_terminal_eight_small_outer.py)
uses the local coefficient enclosures from the
[paired terminal-eight proof](MULTIVARIATE_TERMINAL_EIGHT_PAIRED_LIGHT.md).
It bounds each logarithmic coefficient correction independently
by summing the shared parameter intervals. This yields rational
exponent bounds \([a,b]\) for each normalized multiplier
\([u^b,u^a]\), including the Gaussian exponent in that interval.
The exponent lower endpoints are strictly positive, so at \(u=0\)
the interior multipliers all vanish.

The five parameter intervals use continuous affine ray cosines
and exact radial Bernstein projection checks in the annulus
\(1/8<|t|<8\). The checker proves annulus dominance in the whole
coefficient box and verifies the signs for both pure binomial
anchors at \(u=0\). On the unit circle the relevant cosine
quotients are \(8x^3-4x^2-4x+1\) and \(8x^4-8x^2+1\), with
\(x=\cos\theta\). Their signs give upper sectors in order
\(8,7,8,7,8,7,8\), with one extra negative real septic root.

For each \(u\), move along the straight coefficient segment to the
symmetric Gaussian arrays and then along those arrays to \(u=0\).
The verified boxes and ray exclusions preserve the sector counts.
The unique septic root in the last joined sector is real by
conjugacy; the actual octic pair is nonreal by positivity of the
even factorial moment on the real axis. At \(u=.10\) the rays
also match those of the paired proof exactly.

Replay from the repository root with

~~~sh
python3 -O research/factorial/image-rigidity/verify_multivariate_terminal_eight_small_outer.py
~~~
