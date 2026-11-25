# Terminal-nine retained rows above three smaller reciprocal caps

Let the unequal exponent rows be fully shared, \(x_i=d_i/(8a_i+1)\),
and \(u\) their order-eight normalized outer geometric mean.
Suppose all \(x_i\le1/8\). In each listed band, if at least one
\(x_i\) exceeds the specified cap, the order-eight and order-nine
factorial moments have disjoint root arguments.

| Outer mean \(u\) | Cap |
|---|---:|
| \([31/50,7/10]\) | \(3/40\) |
| \([7/10,19/25]\) | \(11/200\) |
| \([3/4,17/20]\) | \(7/200\) |

Retain every row above the smaller cap. The
[exact replay](verify_multivariate_terminal_nine_retained_bands.py)
keeps the same signed retained list at both orders.

## Finite enumeration and tighter tail cap

The [shared high-curvature bound](MULTIVARIATE_TERMINAL_NINE_SHARED_HIGH_CURVATURE.md)
allows \(S<53/10\). Also \(-\log u\ge7S/16\), so a band
with lower endpoint \(u_-\) has a second curvature budget
\((16/7)\log(1/u_-)\).
Search all rows above the cap within the smaller budget, with
both orientations and possible multiplicities up to reversal.
A retained list's squared outer mean is the product of
\(\prod_{h=1}^d(8a+h)/(8a+7d+h)\) over its rows.
All residual factors are at most one, so the retained square
is at least \(u_-^2\). Use this to prune the enumeration.

With retained outer mean \(u_C\) and remaining budget \(B\),
the tail mean is in
\[
 \left[\max\{u_-/u_C,(8/15)^B\},
       \min\{1,u_+/u_C\}\right]
\]
by the [terminal-nine retained estimates](MULTIVARIATE_TERMINAL_NINE_RETAINED_ROWS.md#retained-factors-and-the-tail).
On a tested interval with tail lower endpoint \(v_-\), each
residual row has own outer mean at least \(v_-\).
For a unit gap this implies
\[
                 v_-^2\le\frac{8a+1}{8a+8}.
\]
Combine the resulting exact integer lower bound on \(a\)
with the assumed reciprocal cap. For \(d\ge2\) each of the
\(d\) outer factors is at most \((8a+d)/(8a+8d)\).
Thus \(v_-\le(8a+d)/(8a+8d)\) and, for \(v_->1/8\),
\[
                 \frac{d}{8a+1}\le\frac{1-v_-}{8v_--1}.
\]
These separate unit and larger-gap estimates tighten the
cap in the joint coefficient models. The Gaussian midpoint
anchor still uses the exact same retained factors.

## Annuli and paired rays

The endpoint-normalized multipliers of each factorial row
are log convex as functions of the coefficient index, and
the Gaussian tail has the same property.
For their products, if \(b=F_{j,1}\) and \(c=F_{j,j-1}\), then
\(c\le F_{j,k-1}/F_{j,k}\le b^{-1}\).
For the interior adjacent coefficients, the binomial ratio
ranges from \(2/(j-1)\) to \((j-1)/2\).
At the two ends the full ratios are \(1/(jb)\) and \(jc\).
Thus whenever \(b,c>(j-1)/20\), all these ratios are strictly
between \(1/10\) and \(10\), since \(j\le9\).
The [annular argument](MULTIVARIATE_TERMINAL_SEVEN_RETAINED_BANDS.md#annuli-anchors-and-projections)
places the roots in that annulus. The ratio inequalities
are linear after cross multiplication and hence also hold
on the straight segment from anchor to actual.
When the stricter analogous test with \(1/8,8\) holds,
the checker first tries the independent coefficient bounds
on the narrower radial interval. Otherwise use the paired
projection on the wider annulus.

The interval tail bounds and endpoint attenuation are as
in the retained-row proof. The independent coefficient
bounds or the common log-moment model verify a strict
projection for both the exact arrays and Gaussian anchors,
with the retained coefficients present in each.
Subdivide the tail interval if necessary.
Certified midpoint root balls and rational rays then give
the alternating sector count as in the retained-row proof.
The sole ninth-order root in the sector joining across the
negative axis is real, while the eighth-order moment is
positive there. This gives disjoint root arguments and
characteristic-zero coprimality.
