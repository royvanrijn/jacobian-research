# Effective relative-gap bands for fully shared terminal-six profiles

Let \(a,b\) be ordinary exponent columns with
\(a_i,b_i>0\) on every unequal row. Put
\[
 S=\sum_i\frac{(b_i-a_i)^2}{\max(a_i,b_i)},\qquad
 \eta=\max_i\frac{|b_i-a_i|}{\max(a_i,b_i)} .
\]

**Theorem.** The order-five and order-six factorial
moments are coprime in each of the following ranges.
The bounds are uniform in the dimension and the
exponent sizes.

| Shared curvature \(S\) | Maximum relative gap \(\eta\) |
| --- | ---: |
| \([5/4,3/2]\) | \(1/500\) |
| \([3/2,2]\) | \(1/100\) |
| \([2,3]\) | \(1/50\) |
| \([3,6]\) | \(1/25\) |

For \(S\ge6\) the
[high-curvature theorem](MULTIVARIATE_TERMINAL_SIX_PRIVATE_ANGLES.md#total-curvature-at-least-six)
already applies to arbitrary relative gaps.
The new ranges make the diffuse comparison effective
on moderate curvature bands.

## Coefficient boxes at a common curvature

Write \(w_j(k)=k(j-k)/(2j)\). After normalizing the
endpoints to one by separate positive changes of variable,
the coefficients at the two orders obey
\[
 \binom jk e^{-S w_j(k)/(1-\eta)}
 \le c_{j,k}\le
 \binom jk e^{-S w_j(k)}.                                  \tag{1}
\]
In fact the second difference of a shared row of gap
\(d\), with minimum \(l\) and maximum \(r\), lies between
\(d^2/(jr)\) and \(d^2/(jl)\) on the logarithmic scale.
The factorial-ratio proof is given in the
[two-private linear case](MULTIVARIATE_TERMINAL_SIX_TWO_PRIVATE_LINEAR.md#a-two-sided-loss-estimate-in-both-orders).
Since \(l\ge(1-\eta)r\), the endpoint-zero quadratic
comparison and summation give (1).

Divide each displayed band into closed subintervals
of width \(1/200\). In each box use the left
curvature endpoint for the upper coefficient bound
and the right endpoint for the lower bound. Both
orders use the same subinterval. The verifier raises
each bound to the order \(j\), encloses the resulting
exponential by a positive Taylor sum through degree
64 with a geometric tail bound, and encloses the
\(j\)-th root with rational denominator \(10000\).

The anchor is the reciprocal polynomial with coefficients
\(\binom jk \exp(-s_* w_j(k))\), where \(s_*\) is the
subinterval midpoint. The four ray cosines are rational
linear interpolations of the entries below, evaluated
at \(s_*\).

| Curvature | First | Second | Third | Fourth |
| --- | ---: | ---: | ---: | ---: |
| \(5/4\) | \(-.181\) | \(-.429\) | \(-.760\) | \(-.861\) |
| \(3/2\) | \(-.060\) | \(-.350\) | \(-.708\) | \(-.850\) |
| \(2\) | \(.135\) | \(-.210\) | \(-.640\) | \(-.807\) |
| \(3\) | \(.402\) | \(-.020\) | \(-.502\) | \(-.768\) |
| \(4\) | \(.577\) | \(.120\) | \(-.403\) | \(-.721\) |
| \(5\) | \(.674\) | \(.210\) | \(-.330\) | \(-.680\) |
| \(6\) | \(.740\) | \(.278\) | \(-.280\) | \(-.648\) |

The ordering and all ray projections are checked exactly.

## Exact ray and anchor verification

The constant and leading terms dominate on radii \(1/8\)
and \(8\), respectively. It suffices to certify the rays
inside that radial annulus. For each ray the lower
bounds of the positive and negative projections are
polynomials in the radius with rational coefficients,
as in the
[projection formula](MULTIVARIATE_TERMINAL_SIX_PRIVATE_ANGLES.md#exact-ray-inequalities).
On a radial interval, a positive Bernstein coefficient
list proves strict positivity of that projection.
Exact bisection uses the de Casteljau averaging rule.
The [verifier](verify_multivariate_terminal_six_diffuse_bands.py)
checks that at least one projection is strictly positive
on every terminal radial interval.

At the anchor, put \(X=z+z^{-1}\). Dividing the quintic
by \((1+z)z^2\) and the sextic by \(z^3\) gives the
quadratic and cubic in
[the two-private linear proof](MULTIVARIATE_TERMINAL_SIX_TWO_PRIVATE_LINEAR.md#a-linked-curvature-band-and-its-symmetric-anchor).
The anchor uses the present binomial-weighted
coefficients in those same expressions.
The verifier checks the signs \(+,-,-,+\) for the
quadratic and \(-,-,+,+\) for the cubic at the four
doubled ray cosines. For the cubic it also checks
negative sign at \(X=-2\); the sign at \(X=2\) is
positive from the original positive coefficients.
The degrees and signs therefore locate every root
of each reduced polynomial in the required
subintervals of \((-2,2)\).

The sector counts from the anchor carry to the full
coefficient box by straight interpolation. Counting
conjugate sectors together around the negative real
axis leaves one real quintic root there. The original
even factorial moment is strictly positive for real
parameters, so this real root cannot match a sextic
root. The other sectors have disjoint root arguments.
Rational Bézout now proves the theorem over every
characteristic-zero field.
