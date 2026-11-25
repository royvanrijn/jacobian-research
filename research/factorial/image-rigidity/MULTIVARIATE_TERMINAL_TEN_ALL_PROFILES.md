# Terminal ten for every ordinary two-monomial profile

For two distinct nonnegative integer exponent columns \(a,b\) in
any finite dimension, put
\[
 Q_j(z)=\sum_{k=0}^j\binom jk\prod_i((j-k)a_i+kb_i)!\,z^k .
\]

**Theorem.** The ninth and tenth polynomials \(Q_9,Q_{10}\) have
disjoint root arguments, including under independent positive
rescalings. In particular \(\gcd(Q_9,Q_{10})=1\) over every
characteristic-zero field.

Delete equal rows, which contribute a nonzero scalar. If the total
private degree is at least two, use the [private-degree theorem](MULTIVARIATE_TERMINAL_TEN_PRIVATE_DEGREES.md).
If the total private degree is one, use the [retained-profile theorem](MULTIVARIATE_TERMINAL_TEN_RETAINED_PROFILES.md).
The latter also covers all fully shared profiles whose maximal
order-nine reciprocal exceeds \(1/8\).

For every remaining profile let \(u\) be the normalized order-nine
outer mean. At least one unequal row gives \(0<u<1\).
The following closed bands cover this interval, with the strict
endpoint \(u<1\) supplied by the profile itself.

| Outer mean | Proof |
|---|---|
| \(0<u\le .76\) | [Matched light ranges](MULTIVARIATE_TERMINAL_TEN_SHARED_LIGHT.md) and [their retained complements](MULTIVARIATE_TERMINAL_TEN_RETAINED_MID_BANDS.md) |
| \([.76,.88]\) | [Upper-middle retained comparison](MULTIVARIATE_TERMINAL_TEN_RETAINED_UPPER_MID.md) |
| \([.88,.93]\) | [Upper bridge](MULTIVARIATE_TERMINAL_TEN_UPPER_BRIDGE.md) |
| \([.93,1)\) | [Stable near-flat continuation](MULTIVARIATE_TERMINAL_TEN_STABLE_UPPER.md) |

In each case all root arguments are disjoint. The two constant terms
and leading coefficients are nonzero, so neither pure monomial direction
can be a projective zero. The coefficients of the original moment
polynomials are integers, so complex coprimality gives a rational
Bézout identity, which transfers to every characteristic-zero field.
