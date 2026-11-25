# Terminal seven across two middle shared reciprocal bands

Let the unequal rows of distinct ordinary exponent columns be
fully shared. Write \(u\) for the order-six outer geometric mean
of the endpoint-normalized factorial multipliers, and \(x_*\) for
the largest \(d_i/(6a_i+1)\) with lower exponent \(a_i\) and gap \(d_i\).

**Theorem.** If \(31/50\le u\le3/4\) and \(x_*\le1/8\), the
order-six and order-seven factorial moments have disjoint root arguments
under independent positive normalizations. In particular they are
coprime over every characteristic-zero field.

Use the [large curvature proof](MULTIVARIATE_TERMINAL_SEVEN_SHARED_HIGH_CURVATURE.md)
if \(S\ge17/4\). Otherwise split the outer means into the bands in this
table. The last column is the reciprocal cap handled by the
[paired light-row theorem](MULTIVARIATE_TERMINAL_SEVEN_PAIRED_LIGHT.md).

| Order-six outer mean | Cap |
| --- | ---: |
| \([.62,.72]\) | \(1/16\) |
| \([.72,.75]\) | \(7/200\) |

On each band retain every row above its cap. An empty list
is covered by the light-row theorem. The
[checker](verify_multivariate_terminal_seven_retained_bands.py)
covers every nonempty list, using exact rational projections and
certified complex-ball anchors. It checks nine retained lists in the
first band and 49 in the second, with 850 tail intervals and
54,242 radial cells in total.

## Bounding the retained lists and the remainder

For a band \([u_-,u_+]\) the total curvature is at most
\[
 \min\left\{\frac{17}{4},\frac{12}{5}\log(1/u_-)\right\}.
\]
This follows from the rowwise outer loss in the
[high-curvature proof](MULTIVARIATE_TERMINAL_SEVEN_SHARED_HIGH_CURVATURE.md#a-row-restriction-from-the-total-outer-loss).
The replay uses an outward rational bound on this budget.
It enumerates candidate rows above the band cap and at or
below \(1/8\), as in the
[retained-row theorem](MULTIVARIATE_TERMINAL_SEVEN_RETAINED_ROWS.md).
The squared outer mean of a row is a rational factorial product
strictly below one. A retained list whose product of these squared
means is smaller than \(u_-^2\) cannot occur and can be dropped.
Enumerate all remaining signed multiplicities up to global
column reversal. No bound on the number of rows below the cap
is needed.

If \(u_C\) is the retained outer mean and \(v\) the residual
outer mean, then \(u=u_Cv\). If \(B\) is the residual curvature
budget, the retained-row proof gives \(v\ge(7/12)^{13B/12}\).
The replay intersects this bound with
\([u_-/u_C,u_+/u_C]\) and \(v\le1\), using outward endpoints.
On any residual interval with lower endpoint \(v_-\), each
residual reciprocal satisfies \(1+5x\le v_-^{-2}\). The replay
uses the band cap or this smaller bound rounded upwards.
The retained integer factors at both orders are calculated
from their factorials. The bounded tail intervals and retained
factors are then used together, holding the same retained rows
and the same residual outer mean at both orders.

## Annuli, anchors, and projections

The first annulus test uses the retained coefficients alone
to give strict endpoint dominance on radii \(1/8\) and \(8\).
When this test gives no bound, use that both the actual and
Gaussian remainder arrays are log convex before restoring the
binomial coefficients. For the full arrays the same is true
after multiplying by the retained multipliers. If the two outer
normalized multipliers exceed \(3/8\), all consecutive full
coefficient ratios are strictly between \(1/8\) and \(8\).
Indeed if the outer multipliers are \(b=F_1,c=F_{j-1}\),
log convexity bounds \(F_{k-1}/F_k\) between \(c\) and \(1/b\).
For \(2\le k\le j-1\) the binomial coefficient adds a ratio
between \(2/(j-1)\) and \((j-1)/2\), and the two extreme
ratios are \(1/(jb)\) and \(jc\). Here \(j=6,7\) and \(b,c\le1\).
The ratio inequalities are linear and hold also along the straight
segment between the actual and Gaussian arrays. They enclose all
roots in the annulus: after scaling the larger radius to one,
multiply a polynomial with increasing positive coefficients by
\(1-t\) and compare the leading term to the remaining coefficient
differences outside the unit disk. Reversing the polynomial gives
the lower radius.

For each residual interval the two anchor arrays consist of the
retained rows times the midpoint Gaussian remainder, as in the
retained-row theorem. Complex balls certify alternating upper
sectors and isolate rational ray proposals in the intervening gaps.
The checker first attempts the independent rational coefficient
bounds with the two endpoint attenuation links. When those
bounds give no signed projection it applies the
[paired Bernstein model](MULTIVARIATE_TERMINAL_SEVEN_PAIRED_LIGHT.md#coefficient-enclosures)
to the remainder and multiplies every coefficient contribution
by the fixed retained factor. Both the Gaussian anchor and the
actual tail arrays are represented by this model. Their straight
segments stay in it, since the Gaussian displacement parameters
and the exponential remainders can all be scaled by the segment
parameter. The scalar interval for each joint contribution is
applied to its whole projected polynomial before checking signs.
This retains the reflected skew cancellations in both orders.

Each accepted tail interval has strict rational or outward
dyadic signs over the entire radial annulus. Subdivide otherwise.
The verified rays and the root bounds preserve all sector counts
along the homotopy. The joined sector at the negative real axis
contains one septic root, necessarily real, and two sextic roots
which are nonreal by positivity of the even factorial moment.
Hence the arguments are disjoint. A rational Bézout identity
gives the field-independent assertion.
