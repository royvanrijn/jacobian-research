# Terminal eleven for every ordinary two-monomial profile

For two distinct nonnegative integer exponent columns \(a,b\) in
any finite dimension, define
\[
 Q_j(z)=\sum_{k=0}^j\binom jk\prod_i((j-k)a_i+kb_i)!\,z^k .
\]

**Theorem.** The tenth and eleventh polynomials \(Q_{10},Q_{11}\)
have disjoint root arguments, including under independent positive
rescalings. In particular they are coprime over every
characteristic-zero field.

Delete equal rows, which contribute a nonzero scalar at each order.
If the total private degree is at least two, the
[private-degree proof](MULTIVARIATE_TERMINAL_ELEVEN_PRIVATE_DEGREES.md)
applies. If the total private degree is one, the
[retained-profile theorem](MULTIVARIATE_TERMINAL_ELEVEN_RETAINED_PROFILES.md)
applies. It also covers fully shared profiles with a reciprocal
\(d/(10\ell+1)>1/8\) on some unequal row of smaller exponent \(\ell\)
and gap \(d\).

For the remaining fully shared profiles let \(u\) be the normalized
order-ten outer mean. Since the columns are distinct, \(0<u<1\).
The following bands cover all possible \(u\).

| Outer mean | Proof |
| --- | --- |
| \((0,.76]\) | [Light-row bands](MULTIVARIATE_TERMINAL_ELEVEN_SHARED_LIGHT.md), [middle retained profiles](MULTIVARIATE_TERMINAL_ELEVEN_RETAINED_MID_BANDS.md), and the [three tightened retained complements](MULTIVARIATE_TERMINAL_ELEVEN_RETAINED_HIGH.md) |
| \([.76,.88]\) | [Cayley comparison with cap \(1/55\)](MULTIVARIATE_TERMINAL_ELEVEN_MIDDLE_CAYLEY.md) and the [37 upper-middle retained lists](MULTIVARIATE_TERMINAL_ELEVEN_RETAINED_HIGH.md) |
| \([.88,.94]\) | [Cayley comparison with cap \(1/70\)](MULTIVARIATE_TERMINAL_ELEVEN_UPPER_CAYLEY.md) and the [four upper retained rows](MULTIVARIATE_TERMINAL_ELEVEN_RETAINED_HIGH.md) |
| \([.94,1)\) | [Unrestricted Cayley bridge](MULTIVARIATE_TERMINAL_ELEVEN_UPPER_SEAM.md) and [stable upper comparison](MULTIVARIATE_TERMINAL_ELEVEN_STABLE_UPPER.md) |

The constant and leading coefficients of \(Q_j\) are nonzero. Thus
disjoint complex root arguments give coprimality and the integer
coefficient arrays give a rational Bézout identity. This identity
transfers the result to every characteristic-zero field.
