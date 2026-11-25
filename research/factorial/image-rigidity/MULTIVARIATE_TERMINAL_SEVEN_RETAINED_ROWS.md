# Terminal seven for 635 retained shared-row lists

Let the distinct ordinary exponent columns be fully shared on
their unequal rows. At lower order six write \(a_i\) for the
smaller exponent, \(d_i>0\) for the gap, and
\(s_i=d_i^2/(a_i+d_i)\). Suppose \(S=\sum_i s_i<17/4\) and retain
all rows with \(x_i=d_i/(6a_i+1)>1/8\).

The order-six and order-seven moments have disjoint root
arguments if the resulting list is nonempty and differs from
the following four lists, with lists identified up to column
reversal and row permutation.

- One row with unordered exponents \((1,2)\).
- Two \((1,2)\) rows in the same orientation.
- Two \((1,2)\) rows in opposite orientations.
- One row with unordered exponents \((2,4)\).

The three \((1,2)\) exceptions are separately covered by the
[marked-row certificate](MULTIVARIATE_TERMINAL_SEVEN_MARKED_LIGHT_TAIL.md).
The [retained-row replay](verify_multivariate_terminal_seven_retained_rows.py)
uses rational intervals for its ray projections and certified
complex balls for the finite Gaussian anchor roots. It certifies
635 lists, 2,237 tail intervals, and 180,302 radial cells.

## Retained rows and the residual outer mean

The integer bounds in the
[shared high-curvature proof](MULTIVARIATE_TERMINAL_SEVEN_SHARED_HIGH_CURVATURE.md#a-row-restriction-from-the-total-outer-loss)
give a finite enumeration. Mark any row above the reciprocal
threshold and reverse the columns if needed so the marked row
is increasing. Enumerate arbitrary multiplicities of all other
candidate rows, in both orientations, under the remaining
curvature budget. This gives 639 nonempty classes up to reversal
and reordering.

At each order \(j\), the \(j\)-th power of the retained
normalized multiplier at index \(k\) is the product over these
oriented rows of
\[
         \frac{((ja+kd)!)^j}{((ja)!)^{\,j-k}((j(a+d))!)^{\,k}}
\]
for an increasing row, with the index reflected for a decreasing
row. This retains exactly the same integer rows in both arrays.

For the remainder let \(v\) be the order-six outer mean. The
reciprocal bound gives \(8d_i\le6a_i+1\), hence \(d_i/a_i\le6/7\).
For \(r=d_i/a_i\), the row's contribution \(U_i\) to \(-\log v\)
satisfies
\[
 \frac{U_i}{s_i}\le\frac{1+r}{2r}\log(1+5r/6)
                 \le\frac{13}{12}\log(12/7).
\]
The last expression before the inequality is increasing in \(r\)
by writing it as \(\frac5{12}\int_0^1(1+r)/(1+5rt/6)\,dt\).
If \(B=17/4-\sum_{\rm retained}s_i\), then
\(v\ge(7/12)^{13B/12}\). Empty remainders have \(v=1\).

Every residual Gauss factor satisfies \(1+5x\le v^{-2}\), since
\(v^{-2}\) is the product of factors \(1+5x\).
The replay improves the cap on a tail interval whenever its
lower endpoint gives a stronger bound. It bounds the residual
arrays by the
[joint envelopes](MULTIVARIATE_TERMINAL_SEVEN_SHARED_HIGH_CURVATURE.md#joint-coefficient-bounds)
and multiplies them by the exact retained multipliers and the
binomial coefficients.

It also uses the endpoint attenuation from each end. The
retained part contributes the factor \((90/91)^{k^2S_{\rm retained}}\)
at any fixed endpoint distance \(k\). The remaining curvature
is at least \(-\log v/[(13/12)\log(12/7)]\). Its contribution
is bounded by \(v^{k^2/60}\), using the exact inequality
\((91/90)^{720}>(12/7)^{13}\).

## Anchor and exact rays

At the midpoint \(v_0\) of a tail interval use a Gaussian
remainder with multipliers \(v_0^{k(6-k)/5}\) and
\(v_0^{(36/35)k(7-k)/6}\) at the two orders. The joint envelopes
contain these exponents. Their attenuation at fixed distance
\(k\) is \(v_0^{k^2/35}\), which is within the residual cap
\(v_0^{k^2/60}\). The exact retained arrays already satisfy
their own attenuation bound. Hence each Gaussian anchor and
actual array, as well as their straight segment, lies inside
the bounded coefficient box.

The retained coefficients alone give strict endpoint dominance
on radii \(1/8\) and \(8\), since every residual multiplier is
at most one. The replay encloses the anchor roots using the
[certified FLINT root isolation](https://flintlib.org/doc/acb_poly/index.html)
and verifies strict alternating upper arguments. It proposes
rational rays between these arguments and accepts them only
after exact ordering and strict rational Bernstein signs over
the full radial interval. Unresolved tail intervals are bisected.
On each accepted interval, continuation along the coefficient
segments preserves the sector counts. The final joined sector
has one real septic root and two sextic roots, which are nonreal
by positivity of the even factorial moment on the real axis.
Thus the actual arguments separate. Rational Bézout identities
give coprimality over every characteristic-zero field.
