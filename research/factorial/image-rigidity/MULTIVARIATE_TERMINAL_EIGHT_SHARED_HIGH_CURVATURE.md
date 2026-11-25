# Shared terminal eight at curvature at least 19/4

Let the unequal rows of two ordinary exponent columns be
fully shared. Let \(a_i\) be the smaller exponent, \(d_i\) the
gap, and put \(S=\sum_i d_i^2/(a_i+d_i)\). If \(S\ge19/4\), the
order-seven and order-eight factorial moments have disjoint root
arguments under independent positive rescalings. In particular
they are coprime over every characteristic-zero field.

For \(S\ge15/2\) use the
[fixed-ray theorem](MULTIVARIATE_TERMINAL_EIGHT_HIGH_CURVATURE.md).
Otherwise let \(u\) be the order-seven outer geometric mean
after endpoint normalization and removal of the binomial factor.
The [checker](verify_multivariate_terminal_eight_shared_high_curvature.py)
uses this same \(u\) for both coefficient arrays.

## Ordinary row caps

Take \(n=7\) in the
[general outer-loss bound](MULTIVARIATE_TERMINAL_SEVEN_SHARED_HIGH_CURVATURE.md#a-row-restriction-from-the-total-outer-loss).
For a row of curvature \(s_i=d_i^2/(a_i+d_i)\) with own
outer mean \(v_i\), it gives
\[
 u\le v_i\exp\!\left[-\frac37\max(19/4-s_i,0)\right],
 \qquad u\le \exp(-57/28)<.132.
\]
At reciprocal floor \(t\), the integer candidates obey
\(a_i<d_i/(7t)-1/7\) and
\(d_i<(15/2)(1+1/(7t))\). Their squared row means are rational.
The first displayed inequality excludes a candidate over a given
interval of \(u\) by comparing its lower endpoint squared with
\(v_i^2\exp[-(6/7)\max(S_0-s_i,0)]\), where \(S_0=19/4\).
The checker uses a positive Taylor lower bound on the inverse
exponential. The factor \(6/7\) comes from squaring the \(u\)
bound.
Rows below the floor have order-seven reciprocal at most
\(1/7\) and order-eight reciprocal at most \(1/(8-1/7)\),
by the matched Gauss relation
\(x'=7x/(8-\gamma x)\), \(0<\gamma\le1\).
The surviving candidates above \(t=1/7\) give the remaining cap
at each order.

When \(u\ge23/250\), first separate a possible row of unordered
exponents \((1,3)\). Its squared outer mean is \(6/35\) and
its curvature is \(4/3\). If any other row exceeds reciprocal
\(1/8\), enumerate it with \(t=1/8\) and use both exact row
means. The squared outer mean of the whole profile then has
upper bound
\[
 (6/35)v_i^2\exp[-(6/7)\max(19/4-4/3-s_i,0)] .
\]
The checker verifies that this is strictly less than
\((23/250)^2\) for every candidate, also excluding a second
\((1,3)\) row. The squared residual outer mean is at least
\((23/250)^2/(6/35)\), which is strictly above
\(\exp[-(6/7)(15/4)]\) by an exact Taylor check.
The residual curvature is therefore less than \(15/4\).
This case is covered by the independent
[marked singleton certificate](MULTIVARIATE_TERMINAL_EIGHT_MARKED_SINGLETONS.md).
For the remaining profiles in \(u\ge23/250\) remove \((1,3)\)
from the row candidates. For \(u\le23/250\) keep the full list.

## Joint coefficient bounds and ray replay

Use the
[general mean and skew formulas](MULTIVARIATE_TERMINAL_SEVEN_SHARED_HIGH_CURVATURE.md#joint-coefficient-bounds)
at orders seven and eight. The formulas for the mean losses
and outer ratio hold at any degree. For the skew magnitude,
the checker tests nonnegativity of the explicit second
derivative in \(s=\log(1+(j-1)x)/2\) on \(0\le x\le1\) for
each reflected pair at the two orders. Indeed the rational
change of variable \(x=z/(1+z)\) clears denominators to
nonnegative coefficient arrays, with positive denominator
constant. All surviving reciprocals lie in that interval.
Rational logarithm intervals at the cap now enclose each
coefficient as powers of \(u\).

Intersect these bounds with the
[total curvature cap](MULTIVARIATE_FIXED_TERMINAL_CURVATURE.md#large-curvature-and-the-endpoint-radii)
and the [adjacent endpoint bound](MULTIVARIATE_ORDINARY_ADJACENT_COEFFICIENT_ATTENUATION.md#an-all-order-endpoint-inequality).
At a fixed distance \(k\) from either endpoint the latter
gives an attenuation at most \((119/120)^{(19/4)k^2}\)
before restoring the binomial factor. The checker uses the
larger \(\exp(-(19/4)k^2/120)\), by
\(\log(120/119)\ge1/120\).

As symmetric comparison arrays take normalized multipliers
\(u^{k(7-k)/6}\) and \(u^{(49/48)k(8-k)/7}\).
These exponents are in the joint envelopes.
The bound \(u\le e^{-3(19/4)/7}\) gives the total curvature
cap for both arrays. The fixed-distance exponent difference
is \(k^2/48\), so the adjacent attenuation is at most
\(\exp(-(19/4)k^2/112)\). The Gaussian arrays and the actual
arrays therefore lie in the same coefficient boxes.
The boxes are convex and the total curvature caps put all
roots at radii between \(1/8\) and \(8\).

The checker uses six rational piecewise affine ray cosines
and the strict Bernstein projections from the
[terminal-eight fixed-ray proof](MULTIVARIATE_TERMINAL_EIGHT_HIGH_CURVATURE.md).
At \(u=0\) both anchors are binomial \(1+t^j\), with
successive upper sectors \(8,7,8,7,8,7,8\). Exact evaluations of
the two binomial cosine polynomials verify these counts at
the first knot. The later ray exclusions hold in 5,434 accepted
radial cells. Continuation along the symmetric arrays and then
along straight coefficient segments preserves the counts.
Joining the last sector to its reflection gives one real
septic root and two octic roots, which are nonreal by positivity
of the even moment on real fields. Thus every argument is
separated, and rational Bézout gives coprimality in characteristic zero.
