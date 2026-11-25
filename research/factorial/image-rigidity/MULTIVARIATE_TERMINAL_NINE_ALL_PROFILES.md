# Terminal nine for every ordinary two-monomial profile

Let \(a,b\) be distinct nonnegative integer exponent columns in
any finite dimension and put
\[
 Q_j(t)=\sum_{k=0}^j\binom jk
             \prod_i((j-k)a_i+kb_i)!\,t^k .
\]
**Theorem.** The root arguments of \(Q_8,Q_9\) are disjoint,
including after independent positive rescalings. In particular
\(\gcd(Q_8,Q_9)=1\) over every characteristic-zero field.

Delete equal rows because they multiply each polynomial by a
nonzero scalar. The total private degree is the sum of the nonzero
exponents in the unequal rows that have a zero exponent. If it is
positive, use the [private degree at least two proof](MULTIVARIATE_TERMINAL_NINE_PRIVATE_ROWS.md)
or the [single private linear proof](MULTIVARIATE_TERMINAL_NINE_PRIVATE_LINEAR.md).

For a fully shared profile let \(\ell_i\) and \(d_i>0\) be the
smaller exponent and the exponent gap of row \(i\). Put
\[
 S=\sum_i\frac{d_i^2}{\ell_i+d_i},\qquad
 x_*=\max_i\frac{d_i}{8\ell_i+1}.
\]
The [shared high-curvature proof](MULTIVARIATE_TERMINAL_NINE_SHARED_HIGH_CURVATURE.md)
covers \(S\ge53/10\). For \(S<53/10\) and \(x_*>1/8\), the
[retained row proof](MULTIVARIATE_TERMINAL_NINE_RETAINED_ROWS.md)
covers every nonempty retained list together with its shared remainder.

For \(x_*\le1/8\), let \(u\) be the order-eight outer geometric
mean of the normalized coefficients after division by the binomial
factors. Distinct columns give \(0<u<1\) by strict factorial log
convexity in a varying row. The following ranges cover this interval.

| Range for \(u\) | Proof |
|---|---|
| \((0,.09]\) | [Small outer means](MULTIVARIATE_TERMINAL_NINE_SMALL_OUTER.md) |
| \([.09,.48]\) | [Matched light rows](MULTIVARIATE_TERMINAL_NINE_SHARED_MATCHED_LIGHT.md) and [their extension](MULTIVARIATE_TERMINAL_NINE_SHARED_MATCHED_EXTENSION.md) |
| \([.48,.55]\) | [Matched upper light rows](MULTIVARIATE_TERMINAL_NINE_SHARED_MATCHED_UPPER.md) at cap \(1/8\) |
| \([.55,.62]\) | [Matched upper light rows](MULTIVARIATE_TERMINAL_NINE_SHARED_MATCHED_UPPER.md) at cap \(1/10\) and [the retained complement](MULTIVARIATE_TERMINAL_NINE_RETAINED_MID_BAND.md) |
| \([.62,.70]\) | [Matched upper light rows](MULTIVARIATE_TERMINAL_NINE_SHARED_MATCHED_UPPER.md) at cap \(3/40\) and [the retained complement](MULTIVARIATE_TERMINAL_NINE_RETAINED_BANDS.md) |
| \([.70,.76]\) | [Matched upper light rows](MULTIVARIATE_TERMINAL_NINE_SHARED_MATCHED_UPPER.md) at cap \(11/200\) and [the retained complement](MULTIVARIATE_TERMINAL_NINE_RETAINED_BANDS.md) |
| \([.76,.88]\) | [Shared middle bridge](MULTIVARIATE_TERMINAL_NINE_MIDUPPER_BRIDGE.md) |
| \([.88,.95]\) | [Upper shared bridge](MULTIVARIATE_TERMINAL_NINE_UPPER_BRIDGE.md) |
| \([.95,1)\) | [Near-flat comparison](MULTIVARIATE_TERMINAL_NINE_NEAR_FLAT.md) |

The last three ranges hold without a reciprocal cap. Every case
separates the arguments of all nonzero roots. Both constant and
leading coefficients are positive, so there is no root at a
coordinate point of the homogeneous coefficient line. Complex
coprimality gives a rational Bézout identity because the original
coefficients are integers. It transfers to every characteristic-zero
field.
