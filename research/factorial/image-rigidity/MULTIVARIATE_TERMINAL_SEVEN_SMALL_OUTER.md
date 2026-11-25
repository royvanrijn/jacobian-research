# Terminal-seven separation for light shared rows with outer mean at most .10

Let \(a,b\) be distinct ordinary exponent columns whose unequal
rows are all shared. Write \(\ell_i=\min(a_i,b_i)\) and
\(d_i=|b_i-a_i|\) on those rows, and let \(u\) be the
order-six normalized outer geometric mean. If
\[
 0<u\le1/10,\qquad \max_i d_i/(6\ell_i+1)\le1/8 ,
\]
the order-six and order-seven factorial moment polynomials
have disjoint root arguments up to independent positive
rescalings. Rational Bézout gives coprimality over every
characteristic-zero field.

For the order-six Gauss factors use the cap \(t=1/8\).
The corresponding order-seven reciprocals satisfy
\(x'\le6t/(7-t)\) by the
[linked reciprocal comparison](MULTIVARIATE_TERMINAL_SEVEN_PAIRED_LIGHT.md#coefficient-enclosures).
Apply the exact independent exponent bounds of the
[one-private-linear proof](MULTIVARIATE_TERMINAL_SEVEN_PRIVATE_LINEAR.md)
to the shared array alone. In particular, if the exponent
interval is \([a,b]\), the normalized multiplier lies in
\([u^b,u^a]\). The bounds cover the symmetric Gaussian
exponents \(e_{6,k}=k(6-k)/5\) and
\((36/35)e_{7,k}=(36/35)k(7-k)/6\).
The [checker](verify_multivariate_terminal_seven_small_outer.py)
verifies the local skew convexity and the positive exponent
enclosures and uses exact integer-power bounds for all
powers of the interval endpoints. At \(u=0\) all allowed
interior multipliers are zero.

The checker uses five affine cosine curves with rational
knots, matching the paired-light rays at \(u=.10\).
It checks constant and leading-term dominance at radii
\(1/8\) and \(8\) in the full coefficient box for
\(0\le u\le.10\). The rational Bernstein projection checks
then exclude all five rays throughout that annulus and
all five knot intervals. The replay used \(275\) sextic
and \(334\) septic radial cells.

At a fixed \(u\), join the actual coefficient array to
the symmetric Gaussian array in the same coefficient
box by a straight segment. Along the Gaussian arrays
let \(u\) decrease to zero with the same moving rays.
The two limiting anchors are \(1+t^6\) and \(1+t^7\).
For \(x=\cos\theta\), their cosine quotients on the
unit circle have signs given by \(8x^3-6x\) and
\(8x^3-4x^2-4x+1\), respectively. The checker evaluates
their signs between the rays and at \(x=\pm1\) to place
the successive upper roots in sectors of orders
\(7,6,7,6,7,6\), with one additional negative real
septic root. The certified exclusions preserve these
sector counts along the path. In the last joined sector
the actual septic root is real by conjugacy, and the
sextic pair is nonreal by positivity of the sixth
factorial moment on the real line. Their arguments
are disjoint.

Replay from the repository root with

~~~sh
python3 -O research/factorial/image-rigidity/verify_multivariate_terminal_seven_small_outer.py
~~~
