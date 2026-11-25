# Terminal nine for all 93 nonempty retained shared profiles

Suppose all unequal rows are shared and \(S<53/10\).
Let \(a_i\) be the smaller exponent and \(d_i>0\) the gap.
If some order-eight reciprocal \(d_i/(8a_i+1)\) exceeds \(1/8\),
then the eighth and ninth ordinary factorial moments have
disjoint root arguments after independent positive rescalings.
Therefore they are coprime over every characteristic-zero field.

## Retained factors and the tail

Retain all rows with reciprocal above \(1/8\) and keep their
exact factorial factors at both orders. The bounds
\(a_i<d_i-1/8\) and
\(d_i^2/(a_i+d_i)\le53/10\) give a finite candidate list.
More generally the residual budget \(B\) gives \(d_i<2B\)
and \(a_i<d_i-1/8\) for each retained row under that budget.
Mark an increasing row after reversing the columns if needed.
Enumerate residual multiplicities in both directions and
identify a list with its simultaneous reversal. The
[exact replay](verify_multivariate_terminal_nine_retained_rows.py)
finds 93 nonempty possibilities.

For a fixed retained list write \(S_C\) for its curvature and
\(v\) for the order-eight outer mean of its light residual rows.
Their reciprocal cap implies \(d_i\le a_i\), so the
[outer loss estimate](MULTIVARIATE_TERMINAL_NINE_SMALL_OUTER.md#joint-bounds-from-the-same-outer-mean)
gives \(v\ge(8/15)^{53/10-S_C}\).
Since each residual outer factor square has denominator
\(1+7x\), any residual reciprocal also satisfies
\(x\le(v^{-2}-1)/7\). On each residual interval the checker
tightens its reciprocal bound using the interval's lower endpoint
and an outward rational rounding of this last inequality.

## Two ray bounds on the same coefficient arrays

For the two orders, enclose the residual means and skews by
the [joint outer bounds](MULTIVARIATE_TERMINAL_NINE_SMALL_OUTER.md#joint-bounds-from-the-same-outer-mean)
and multiply them by the exact retained factors.
At each fixed distance \(k\) from an endpoint their adjacent
attenuation is at most \((152/153)^{S_C k^2}v^{k^2/96}\),
by applying the
[all-order endpoint inequality](MULTIVARIATE_ORDINARY_ADJACENT_COEFFICIENT_ATTENUATION.md#an-all-order-endpoint-inequality)
to the retained and residual rows separately. The symmetric
anchor uses the exact retained multipliers and Gaussian residual
powers at the midpoint \(v_0\). These meet the joint envelope
and the attenuation bounds because the Gaussian adjacent
exponent difference is \(k^2/63\).
On each fixed retained list the retained-only upper bounds
give strict annular endpoint dominance at radii \(1/8,8\).
Residual multipliers are at most one. Thus both the anchors
and the actual arrays lie in that annulus, also along the
straight coefficient path.

When independent coefficient enclosures fail, use the
[paired exponential projection](MULTIVARIATE_TERMINAL_NINE_PRIVATE_LINEAR.md)
for the residual with the exact retained coefficients in place
of the private core. This model uses a scalar interval for
each common shape contribution and also keeps the quadratic
products inside the projection. The midpoint Gaussian
coefficients belong to the same enclosures. For the paired
ray tests the radial interval \([1/10,10]\) contains the
certified annulus. In either test, a strict projection sign
passes along the straight coefficient path.

The checker isolates the midpoint anchor roots with certified
complex balls and demands the alternating upper order
\(9,8,9,8,9,8,9,8\). It proposes rational ray cosines strictly
between those exact root-angle intervals, and verifies strict
Bernstein projection signs for both orders at every radius
on the applicable interval. It bisects \(v\) whenever it
cannot certify the rays. All 93 lists pass in 280 outer-mean
cells and 25,837 radial cells, with no exceptional lists.
The unique remaining ninth-order root is negative real, and
the even moment cannot have a real root. The continued sector
counts therefore give the asserted angular separation.
