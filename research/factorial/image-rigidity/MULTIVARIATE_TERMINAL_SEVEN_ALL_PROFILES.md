# Terminal-seven coprimality for every ordinary two-monomial profile

Let \(a,b\) be distinct columns of nonnegative integer exponents in
any finite dimension, and write their factorial moment polynomials as
\[
 Q_j(t)=\sum_{k=0}^j\binom jk\prod_i((j-k)a_i+kb_i)!\,t^k .
\]
Remove the equal rows, which contribute a nonzero scalar at each
order. A row is private if one exponent is zero.

**Theorem.** For every such profile, \(Q_6\) and \(Q_7\) have
disjoint complex root arguments under independent positive
changes of scale. They are coprime over every
characteristic-zero field.

## Private rows

If there is precisely one private linear row and all others
are shared, apply the [one-private-linear theorem](MULTIVARIATE_TERMINAL_SEVEN_PRIVATE_LINEAR.md).
Every other nonempty private allocation has total private degree at
least two and is covered by the
[private-degree theorem](MULTIVARIATE_TERMINAL_SEVEN_PRIVATE_ROWS.md).

## Fully shared rows

Put \(S=\sum_i d_i^2/(a_i+d_i)\), writing \(a_i>0\) for the
lower exponent in row \(i\) and \(d_i>0\) for the gap, and put
\(x_*=\max_i d_i/(6a_i+1)\).
Apply the [shared high-curvature theorem](MULTIVARIATE_TERMINAL_SEVEN_SHARED_HIGH_CURVATURE.md)
when \(S\ge17/4\). Otherwise, if \(x_*>1/8\), retain every row
above \(1/8\). The [635-list retained certificate](MULTIVARIATE_TERMINAL_SEVEN_RETAINED_ROWS.md)
covers all such nonempty lists apart from four. A single \((1,2)\)
row or a pair of them in either relative orientation is covered by
the [marked light-tail certificate](MULTIVARIATE_TERMINAL_SEVEN_MARKED_LIGHT_TAIL.md).
A single \((2,4)\) row is covered by the
[second marked certificate](MULTIVARIATE_TERMINAL_SEVEN_MARKED_24_TAIL.md).

It remains to take \(x_*\le1/8\). With \(u\) the order-six
outer geometric mean after removing binomial factors and
normalizing both endpoints to one, \(0<u\le1\). The linked
separation results cover the following closed intervals.

| Order-six outer mean | Separation proof |
| --- | --- |
| \(0<u\le1/10\) | [Small outer mean](MULTIVARIATE_TERMINAL_SEVEN_SMALL_OUTER.md) |
| \(1/10\le u\le12/25\) | [Reciprocal cap \(1/8\)](MULTIVARIATE_TERMINAL_SEVEN_LIGHT_CAP_EXTENSION.md) |
| \(12/25\le u\le31/50\) | [Cross-order skew transfer](MULTIVARIATE_TERMINAL_SEVEN_CROSS_SKEW.md) |
| \(31/50\le u\le3/4\) | [Two retained bands](MULTIVARIATE_TERMINAL_SEVEN_RETAINED_BANDS.md) |
| \(3/4\le u\le22/25\) | [Late Cayley bridge](MULTIVARIATE_TERMINAL_SEVEN_LATE_BRIDGE.md) |
| \(22/25\le u\le19/20\) | [Upper Cayley bridge](MULTIVARIATE_TERMINAL_SEVEN_UPPER_BRIDGE.md) |
| \(19/20\le u\le1\) | [Near-flat regime](MULTIVARIATE_TERMINAL_SEVEN_NEAR_FLAT.md) |

Each normalization uses a positive rescaling of the variable
and a nonzero scalar. The original \(Q_6,Q_7\) therefore have
no common complex root. Their rational Bézout identity extends
to any characteristic-zero field. The two pure monomial moments
are nonzero at each order, so the two remaining projective
coefficient choices are also covered.

Together with the [proof through terminal six](MULTIVARIATE_TERMINAL_SIX_ALL_PROFILES.md),
this proves multivariate SFC(2) at every terminal from two through seven.
