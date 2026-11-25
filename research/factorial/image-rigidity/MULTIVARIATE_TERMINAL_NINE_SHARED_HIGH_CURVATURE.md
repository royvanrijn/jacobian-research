# Shared terminal nine at curvature at least 53/10

Let the unequal rows of two ordinary exponent columns be fully shared.
For each such row write \(a_i\ge1\) for its smaller exponent and
\(d_i>0\) for its gap. If
\[
                  S=\sum_i\frac{d_i^2}{a_i+d_i}\ \ge\ \frac{53}{10},
\]
the order-eight and order-nine factorial moments have disjoint root
arguments under independent positive rescalings. They are coprime over
every characteristic-zero field.

The [unrestricted terminal-nine certificate](MULTIVARIATE_TERMINAL_NINE_HIGH_CURVATURE.md)
covers \(S\ge1493/200\). For the remaining range, use the order-eight outer
geometric mean \(u\) of the multipliers after endpoint normalization and
removal of the binomial coefficients. The
[checker](verify_multivariate_terminal_nine_shared_high_curvature.py)
encloses both coefficient arrays using this one parameter.

## Candidate rows and linked envelopes

Write \(s_i=d_i^2/(a_i+d_i)\). The squared outer mean on row \(i\) is
\[
 v_i^2=\prod_{h=1}^{d_i}
                  \frac{8a_i+h}{8a_i+7d_i+h}.
\]
The [general outer-loss bound](MULTIVARIATE_TERMINAL_SEVEN_SHARED_HIGH_CURVATURE.md#a-row-restriction-from-the-total-outer-loss)
at \(n=8\) gives
\[
 v_i^2\le e^{-7s_i/8},\qquad
 u^2\le v_i^2 e^{-(7/8)\max(53/10-s_i,0)},\qquad
 u\le e^{-371/160}<99/1000.
\]
Indeed each log loss for \(v_i^2\) is at least
\(7d_i/(8a_i+7d_i+h)\ge 7d_i/(8(a_i+d_i))\).
A row whose largest order-eight reciprocal \(d_i/(8a_i+1)\) exceeds
the floor \(t=1/8\) has
\[
 a_i<d_i/(8t)-1/8,\qquad d_i<(1493/200)(1+1/(8t)).
\]
The checker enumerates this finite rectangle, imposing the exact
curvature and reciprocal cutoffs. On an interval of \(u\), it discards
a row if the squared interval minimum exceeds the squared profile
cap above. The squared cap is enclosed from above by dividing the
rational \(v_i^2\) by a positive Taylor lower bound for the inverse
exponential. Rows below the reciprocal floor have order-nine reciprocals
at most \(8t/(9-t)\), since the matched Gauss factors satisfy
\(x'=8x/(9-\gamma x)\) for \(0<\gamma\le1\).
The separate maxima over the surviving candidates, together with the
floor bounds, give reciprocal caps for each of the two orders.

Apply the [mean-loss, skew and outer-ratio formulas](MULTIVARIATE_TERMINAL_SEVEN_SHARED_HIGH_CURVATURE.md#joint-coefficient-bounds)
at orders eight and nine. The checker proves convexity of each noncentral
skew in \(s=\log(1+(j-1)x)/2\) on \(0\le x\le1\). It uses the exact
log weights in that reference, differentiates twice with the positive
scale \((1+(j-1)x)d/dx\), and checks the rational coefficients after
\(x=z/(1+z)\). The enumerated caps and the floor bounds are within
this interval. For the mean-loss ratio the lower bound is its value
at the reciprocal cap and the upper bound is \(k(j-k)/(j-1)\) at
order \(j\). The skew-to-loss ratio is bounded by its value at the cap.
The exact outer-log ratio bounds from the matched Gauss factors turn
the order-nine losses into powers of the same \(u\). Rational logarithm
enclosures give the coefficient intervals.

Intersect these intervals with the
[total curvature caps](MULTIVARIATE_FIXED_TERMINAL_CURVATURE.md#large-curvature-and-the-endpoint-radii)
and the [adjacent endpoint inequality](MULTIVARIATE_ORDINARY_ADJACENT_COEFFICIENT_ATTENUATION.md#an-all-order-endpoint-inequality),
which gives attenuation at most \((152/153)^{Sk^2}\) at distance \(k\)
from either endpoint. The checker uses the outward bound
\(\exp(-(53/10)k^2/153)\), since \(\log(153/152)\ge1/153\).
The attenuation is applied to the normalized multiplier, before
putting back the appropriate binomial coefficient.

## Exact ray replay

The symmetric comparison multipliers are
\[
                u^{k(8-k)/7},\qquad u^{(64/63)k(9-k)/8}.
\]
They lie in the log-loss envelopes. The bound \(u\le e^{-371/160}\)
gives both total curvature caps. At the same endpoint distance \(k\),
the order-nine exponent exceeds the order-eight exponent by \(k^2/63\).
Thus their adjacent ratio is at most \(\exp(-(53/10)k^2/144)\) and
satisfies the attenuation bounds from both endpoints. The comparison
arrays and the actual arrays are in the same coefficient boxes, so
their straight segments stay there. The checker verifies that every
root of any represented polynomial has modulus between \(1/8\) and \(8\).

The seven ray cosines interpolate affinely in \(u\) between the
checker’s rational knots. The checker bisects the outer-mean bands
as needed and tests each ray for both orders by
[strict rational Bernstein projections](MULTIVARIATE_TERMINAL_SIX_JOINT_RAY_BANDS.md#moving-rays-and-exact-replay).
It also verifies that the rays remain ordered. At \(u=0\) the
comparison polynomials are \(1+t^8\) and \(1+t^9\), and exact binomial
cosine tests place eight upper roots in successive sectors in the order
\(9,8,9,8,9,8,9,8\). The successful replays visit 14,738 radial cells.
Continuation along the comparison polynomials and then along straight
coefficient segments preserves the counts. The last sector joined
with its reflection contains two octic roots and one order-nine root.
The single order-nine root is real by conjugacy. The octic has
positive value at every real argument, by the even-moment
representation, so its pair is nonreal. This separates every
root argument. Rational Bézout then transfers coprimality to all
characteristic-zero fields.
