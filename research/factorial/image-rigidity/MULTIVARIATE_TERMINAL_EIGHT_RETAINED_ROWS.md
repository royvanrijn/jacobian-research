# Terminal eight for 242 retained shared-row lists

Let the unequal rows of two distinct ordinary exponent columns be
fully shared. Write \(a_i\) for the smaller exponent, \(d_i>0\) for the
gap, and \(s_i=d_i^2/(a_i+d_i)\). Suppose \(S=\sum_i s_i<19/4\) and retain
every row with \(x_i=d_i/(7a_i+1)>1/8\).

The order-seven and order-eight moments have disjoint root arguments
when the retained list is nonempty and differs, up to column reversal
and row permutation, from the three singleton lists with unordered
exponents \((1,3)\), \((2,4)\), or \((3,6)\). The
[retained-row checker](verify_multivariate_terminal_eight_retained_rows.py)
replays 242 lists, 774 tail intervals, and 69,172 radial cells.

## Enumeration and tail constraints

The inequalities \(7a+1<8d\) and \(d^2/(a+d)\le B\), where \(B\) is
the remaining curvature budget, give explicit finite bounds on \(a,d\).
Mark a retained row, orient it as increasing by column reversal, and
enumerate all multiplicities of the other admissible rows in both
orientations. Sorting and taking the smaller of the two reversed
lists gives 245 classes. Each retained row enters both orders by the
exact powers of its normalized factorial multiplier as in the
[terminal-seven enumeration](MULTIVARIATE_TERMINAL_SEVEN_RETAINED_ROWS.md#retained-rows-and-the-residual-outer-mean).

For the remainder write \(v\) for its order-seven outer mean and
\(B=19/4-\sum_{\mathrm{retained}}s_i\).
Here \(8d\le7a+1\) implies \(d/a\le1\). The same integral calculation
as at terminal seven, now with \(r=d/a\) and \(6/7\) replacing \(5/6\),
gives
\[
  -\log v_i\le s_i \frac{1+r}{2r}\log(1+6r/7)
       \le s_i\log(13/7).
\]
Consequently \(v\ge(7/13)^B\), including \(v=1\) for an empty
remainder. Each residual reciprocal satisfies \(1+6x\le v^{-2}\).
Over a subinterval \([v_-,v_+]\) use the smaller of \(1/8\) and an
outward rational bound on \((v_-^{-2}-1)/6\). The matching formula
\(x'=7x/(8-\gamma x)\), \(0<\gamma\le1\), gives the order-eight cap.
Apply the [shared joint mean and skew envelopes](MULTIVARIATE_TERMINAL_EIGHT_SHARED_HIGH_CURVATURE.md#joint-coefficient-bounds-and-ray-replay)
at these caps to the remainder and multiply by the exact retained
multipliers.

For an index at distance \(k\) from either endpoint, the retained
endpoint attenuation is at most
\((119/120)^{k^2 S_{\rm retained}}\).
Since \(S_{\rm tail}\ge-\log v/\log(13/7)\) and
\((120/119)^{75}>13/7\), the remaining factor is at most
\(v^{k^2/75}\). The checker intersects the bounds from the two endpoints
with its coefficient boxes.

## Anchors and exact rays

At a tail midpoint \(v_0\), use Gaussian residual multipliers
\(v_0^{k(7-k)/6}\) and \(v_0^{(49/48)k(8-k)/7}\).
The joint envelopes contain these exponents and their attenuation
\(v_0^{k^2/48}\) passes the tail bound. Thus the exact retained
multipliers times these arrays, and the actual moments, lie in the
same convex coefficient box.

Endpoint dominance with retained coefficients and all residual
multipliers bounded by one puts both arrays in the annulus with radii
\(1/8\) and \(8\). Certified complex root balls for each anchor give
the upper-sector order \(8,7,8,7,8,7,8\). Rational rays proposed
between these arguments are accepted only when rational Bernstein
projections certify strict signs over the full radial interval.
Bisecting tail intervals when these tests fail gives the stated finite
replay. Strict exclusion along the rays and both annulus boundaries
preserves each sector count along the coefficient segments.

The last sector joined to its reflection contains one septic root and
two octic roots. The septic root there is real by conjugation, while
the octic moment is positive at real arguments by its even integral
representation. The disjoint root arguments imply coprimality over
every characteristic-zero field by rational Bézout.
