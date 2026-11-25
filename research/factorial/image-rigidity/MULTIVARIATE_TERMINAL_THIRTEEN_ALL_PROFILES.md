# Terminal thirteen for every ordinary two-monomial profile

For distinct nonnegative integer exponent columns \(a,b\) in any
finite dimension, put
\[
 Q_j(z)=\sum_{k=0}^j\binom jk\prod_i((j-k)a_i+kb_i)!\,z^k .
\]

**Theorem.** The twelfth and thirteenth factorial moment polynomials
have disjoint root arguments, including after independent positive
rescalings. In particular they are coprime over every
characteristic-zero field.

Delete equal rows, which contribute a nonzero scalar at each order.
The [retained-profile theorem](MULTIVARIATE_TERMINAL_THIRTEEN_RETAINED_PROFILES.md)
covers every profile with a private row and every fully shared profile
with some reciprocal \(d/(12\ell+1)>1/8\). Here \(\ell\) is the smaller
exponent on the row and \(d\) is the positive difference.

In the remaining fully shared case let \(u\) be the normalized
order-twelve reflected outer mean. Since the columns differ,
\(0<u<1\). The following direct tests cover all values of \(u\)
under their respective reciprocal caps.

| Outer mean \(u\) | Direct test |
| --- | --- |
| \((0,.66]\) | [Light rays](MULTIVARIATE_TERMINAL_THIRTEEN_SHARED_LIGHT.md), with the listed caps |
| \([.66,.76]\) | Cap \(1/50\) by the [light rays](MULTIVARIATE_TERMINAL_THIRTEEN_SHARED_LIGHT.md) on \([.66,.70]\) and \([.71,.72]\), and the [root disks](MULTIVARIATE_TERMINAL_THIRTEEN_SHARED_DISKS.md) on \([.69,.71]\) and \([.72,.76]\) |
| \([.76,.95]\) | [Upper Cayley bands](MULTIVARIATE_TERMINAL_THIRTEEN_CAYLEY.md), with the listed caps |
| \([.95,1)\) | [Unrestricted upper range](MULTIVARIATE_TERMINAL_THIRTEEN_CAYLEY.md) |

The [retained bands](MULTIVARIATE_TERMINAL_THIRTEEN_RETAINED_BANDS.md)
fill the intervals above the primary caps through \(1/8\).
The following refinements connect the lower direct caps and the
primary retained caps wherever the latter are larger.

- On \([.47,.50]\) and \([.52,.54]\), the retained bands fill
  \(1/22\) through \(1/20\). The
  [smaller-cap complements](MULTIVARIATE_TERMINAL_THIRTEEN_RETAINED_FALLBACKS.md)
  also fill \(1/30\) through \(1/22\) on \([.54,.55]\).
- The smaller-cap complements fill \(1/40\) through \(1/20\)
  on \([.56,.60]\), and \(1/40\) through \(1/30\)
  on \([.60,.63]\). The
  [retained bridge bands](MULTIVARIATE_TERMINAL_THIRTEEN_RETAINED_BRIDGES.md)
  continue \(1/40\) through \(1/30\) on \([.63,.66]\).
- The smaller-cap complements fill \(1/50\) through \(1/30\)
  on \([.66,.68]\). The retained bridge bands fill \(1/50\)
  through \(1/36\) on \([.68,.72]\), and the retained bands fill
  \(1/36\) through \(1/30\) there. On \([.72,.76]\) the bridge
  bands fill \(1/50\) through \(1/34\), the retained cap there.

Every cut in the outer mean lies in the adjacent closed bands,
and each direct cap includes equality. Thus these tests cover all
ordinary profiles. Both polynomials have nonzero constant and
leading coefficients. Disjoint complex root arguments yield
coprimality, and an integer coefficient array gives a rational
Bézout identity. This transfers the conclusion to every
characteristic-zero field.
