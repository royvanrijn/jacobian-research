# Terminal eight for every ordinary two-monomial profile

Let \(a,b\) be distinct nonnegative integer exponent columns of
any finite dimension and let
\[
 Q_j(t)=\sum_{k=0}^j\binom jk
             \prod_i((j-k)a_i+kb_i)!\,t^k .
\]
**Theorem.** The polynomials \(Q_7,Q_8\) have disjoint root
arguments, including after separate positive rescalings of the
variable. In particular \(\gcd(Q_7,Q_8)=1\) over every
characteristic-zero field.

Remove equal rows, since they multiply each moment by a
nonzero rational scalar. The sum of the private degrees is
either at least two or exactly one whenever there is a private
row. These alternatives are covered by the
[private-degree certificate](MULTIVARIATE_TERMINAL_EIGHT_PRIVATE_ROWS.md)
and the [single-private-linear certificate](MULTIVARIATE_TERMINAL_EIGHT_PRIVATE_LINEAR.md).

For a fully shared profile write \(S=\sum_i d_i^2/(a_i+d_i)\),
where \(a_i\) now denotes the smaller entry in each row and
\(d_i>0\) its gap, and put \(x_*=\max_i d_i/(7a_i+1)\).
The [shared high-curvature certificate](MULTIVARIATE_TERMINAL_EIGHT_SHARED_HIGH_CURVATURE.md)
covers \(S\ge19/4\).
If \(S<19/4\) and \(x_*>1/8\), retain each row with reciprocal
above \(1/8\). The [retained-list certificate](MULTIVARIATE_TERMINAL_EIGHT_RETAINED_ROWS.md)
and the [three singleton certificates](MULTIVARIATE_TERMINAL_EIGHT_MARKED_SINGLETONS.md)
cover every such nonempty list and its arbitrary light remainder.

If \(x_*\le1/8\), use the normalized order-seven outer mean
\(0<u<1\). All intervals in the table are closed at their
positive endpoints.

| Outer-mean range | Proof |
|---|---|
| \(0<u\le1/10\) | [Small outer means](MULTIVARIATE_TERMINAL_EIGHT_SMALL_OUTER.md) |
| \([1/10,12/25]\) | [Paired light rows](MULTIVARIATE_TERMINAL_EIGHT_PAIRED_LIGHT.md) |
| \([12/25,13/25]\) | [First paired bridge band](MULTIVARIATE_TERMINAL_EIGHT_PAIRED_BRIDGE.md) |
| \([13/25,18/25]\) | [Paired bridge](MULTIVARIATE_TERMINAL_EIGHT_PAIRED_BRIDGE.md), [continuation](MULTIVARIATE_TERMINAL_EIGHT_PAIRED_BRIDGE_72.md), and [above-cap retained lists](MULTIVARIATE_TERMINAL_EIGHT_RETAINED_BANDS.md) |
| \([18/25,22/25]\) | [Late shared bridge](MULTIVARIATE_TERMINAL_EIGHT_LATE_BRIDGE.md) |
| \([22/25,19/20]\) | [Upper shared bridge](MULTIVARIATE_TERMINAL_EIGHT_UPPER_BRIDGE.md) |
| \([19/20,1)\) | [Near-flat comparison](MULTIVARIATE_TERMINAL_EIGHT_NEAR_FLAT.md) |

The paired bridge and its continuation use smaller caps
\(2/25,13/200,9/200\) on their last three intervals.
The retained-list proof covers precisely the complement of each
of those caps up to \(1/8\). The last three ranges in the table
are valid for unrestricted reciprocals. Distinct columns imply
\(u<1\) by the strict factorial log convexity of a varying row.
Every case gives strict separation of the arguments.
The polynomials have integer coefficients, so their rational
Bézout identity transfers to every characteristic-zero field.
