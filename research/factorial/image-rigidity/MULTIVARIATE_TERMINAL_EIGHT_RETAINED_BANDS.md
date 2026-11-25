# Terminal eight above the middle light-row caps

For fully shared ordinary exponent columns let
\(x_i=d_i/(7a_i+1)\) be the order-seven reciprocal, where \(a_i\)
is the smaller exponent and \(d_i>0\) is the gap.
Let \(u\) be the normalized order-seven outer mean. Suppose all
\(x_i\le1/8\). In each of the following bands, if at least one
\(x_i\) exceeds the listed cap, the order-seven and order-eight
moments have disjoint root arguments.

| Outer mean \(u\) | Cap | Nonempty retained lists | Tail intervals | Radial cells |
|---|---:|---:|---:|---:|
| \([13/25,3/5]\) | \(2/25\) | 5 | 55 | 4,031 |
| \([3/5,33/50]\) | \(13/200\) | 12 | 94 | 7,323 |
| \([33/50,18/25]\) | \(9/200\) | 27 | 220 | 16,711 |

The [exact replay](verify_multivariate_terminal_eight_retained_bands.py)
keeps the entire list above the cap in both moment arrays.

## Enumeration and residual parameters

The shared high-curvature theorem covers \(S\ge19/4\). In the
remaining range use the lower outer loss
\(-\log u\ge(3/7)S\) from the
[linked order-seven bound](MULTIVARIATE_TERMINAL_EIGHT_SHARED_HIGH_CURVATURE.md#ordinary-row-caps).
On a band \([u_-,u_+]\) the checker therefore uses a budget of
\(\min(19/4,(7/3)L_+)\), where \(L_+\) is an exact upper bound for
\(-\log u_-\). Retain every row with reciprocal strictly above the
listed cap. Its outer-mean square is
\(\prod_{h=1}^{d}(7a+h)/(7a+6d+h)\), independent of orientation.
The product of these squares must be at least \(u_-^2\), since
all residual outer means are at most one. Both this product
constraint and the finite curvature budget are checked when
collecting the row lists up to column reversal and permutation.

For a retained list with outer mean \(u_{\rm ret}\) and remaining
budget \(B\), the tail mean belongs to
\[
 \left[\max\{u_-/u_{\rm ret},(7/13)^B\},
       \min\{1,u_+/u_{\rm ret}\}\right].
\]
The lower curvature bound on the tail and the fixed-distance
attenuation are as in the
[cap-\(1/8\) retained-row certificate](MULTIVARIATE_TERMINAL_EIGHT_RETAINED_ROWS.md).
The checker tightens the tail reciprocal cap on each subinterval
and uses the exact common retained factors of the two arrays.

## Anchor and strict ray tests

The comparison anchor is the retained array times the symmetric
Gaussian tail at the subinterval midpoint. Its certified root
balls must have upper arguments alternating \(8,7,8,7,8,7,8\).
The replay attempts exact independent ray projections on
the annulus \((1/8,8)\) once both outer normalized multipliers
exceed \(7/16\). Log convexity then bounds the inner consecutive
coefficient ratios between \(1/8\) and \(8\) in both orders, as
in the [terminal-seven annular test](MULTIVARIATE_TERMINAL_SEVEN_RETAINED_BANDS.md#annuli-anchors-and-projections).
At the two ends, the ratios are \(1/(jF_1)\) and \(jF_{j-1}\).
The retained list is nonempty, so its normalized multipliers are
strictly below one at both indices. These ratio inequalities
hold along the straight anchor-to-actual segment.

The fallback uses the
[paired light-row model](MULTIVARIATE_TERMINAL_EIGHT_PAIRED_LIGHT.md)
for the tail on \((1/12,12)\), with all allowed signs of the linked
lower-order skew. Multiply each coefficient contribution by its
fixed retained multiplier. Both the anchor and actual
tail arrays are represented, so a strict linear projection
also holds along their straight segment. Every normalized
multiplier lies in \([0,1]\), and \((13/12)^8<2\) gives strict
boundary dominance on this larger annulus. Bisect the tail
interval when the ray tests fail. Strict exclusion on the rays
and the relevant radii preserves the anchor's sector counts;
the last sector joins to its reflection and contains one real
septic root and two nonreal octic roots. This proves the claimed
separation and characteristic-zero coprimality.
