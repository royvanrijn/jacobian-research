# Terminal twelve for every ordinary two-monomial profile

For distinct nonnegative integer exponent columns \(a,b\) in any
finite dimension, put
\[
 Q_j(z)=\sum_{k=0}^j\binom jk\prod_i((j-k)a_i+kb_i)!\,z^k .
\]

**Theorem.** The eleventh and twelfth factorial moment polynomials
have disjoint root arguments, including after independent positive
rescalings. In particular they are coprime over every
characteristic-zero field.

Delete equal rows, which contribute a nonzero scalar at each order.
The [private-degree theorem](MULTIVARIATE_TERMINAL_TWELVE_PRIVATE_DEGREES.md)
covers total private degree at least two. The
[retained-profile theorem](MULTIVARIATE_TERMINAL_TWELVE_RETAINED_PROFILES.md)
covers total private degree one and every fully shared profile with
some reciprocal \(d/(11\ell+1)>1/8\).

In the remaining fully shared case let \(u\) be the normalized
order-eleven reflected outer mean. The columns are distinct, so
\(0<u<1\). For each band the light or Cayley proof applies under the
stated cap. The [retained-band theorem](MULTIVARIATE_TERMINAL_TWELVE_RETAINED_BANDS.md)
covers all reciprocals above the cap and at most \(1/8\).

| Outer mean \(u\) | Cap source |
| --- | --- |
| \((0,.66]\) | [Light-row bands](MULTIVARIATE_TERMINAL_TWELVE_SHARED_LIGHT.md), with their listed caps |
| \([.66,.70]\) | [Cayley bridge](MULTIVARIATE_TERMINAL_TWELVE_CAYLEY_BRIDGE.md), cap \(1/25\) |
| \([.70,.76]\) | [Conditional Cayley comparison](MULTIVARIATE_TERMINAL_TWELVE_CAYLEY_LOWER.md), cap \(3/100\) |
| \([.76,.95]\) | [Upper Cayley bands](MULTIVARIATE_TERMINAL_TWELVE_CAYLEY.md), with their listed caps |
| \([.95,1)\) | [Unrestricted stable upper range](MULTIVARIATE_TERMINAL_TWELVE_CAYLEY.md) |

For the tightened light slabs the retained coverage can be applied
in layers. On \([.36,.40]\) the caps are \(1/8,17/200,3/40,3/50\).
On \([.47,.48]\) they are \(1/8,3/40,13/200\), and on
\([.50,.56]\) they are \(1/8,3/40,11/200,1/20\).
Likewise the bridges on \([.66,.70]\), \([.70,.76]\), and
\([.855,.88]\) use the retained complements between their caps and
the larger parent caps in the retained table. This exhausts every
ordinary profile, including the common endpoints.

The constant and leading coefficients of each \(Q_j\) are nonzero.
Disjoint complex root arguments therefore give coprimality, and the
integer coefficient arrays give a rational Bézout identity.
This transfers the result to every characteristic-zero field.
