# Terminal-six coprimality for every ordinary two-monomial profile

Let \(a,b\) be distinct columns of nonnegative integer exponents in
any finite dimension. Define
\[
 Q_j(t)=\sum_{k=0}^j\binom jk
              \prod_i((j-k)a_i+kb_i)!\,t^k .
\]
A variable is private if exactly one of the entries in its row is
zero. Rows with equal entries contribute a scalar and can be removed.

**Theorem.** The polynomials \(Q_5\) and \(Q_6\) are coprime over
every characteristic-zero field. Their complex root arguments are
disjoint even under independent positive rescalings of the variable.

## Profiles with a private variable

If private variables occur on exactly one side, use the
[one-private-side theorem](MULTIVARIATE_TERMINAL_SIX_ONE_PRIVATE_SIDE.md).
Otherwise let \(P_a,P_b>0\) be the sums of the private exponents
on the two respective sides. If either sum is at least three,
use the [private-degree theorem](MULTIVARIATE_TERMINAL_SIX_PRIVATE_ANGLES.md).
The same theorem covers sums \((2,2)\), \((2,1)\) and \((1,2)\).
The [two-private linear theorem](MULTIVARIATE_TERMINAL_SIX_TWO_PRIVATE_LINEAR.md)
covers \((1,1)\). These possibilities exhaust the positive integer
private degrees.

## Fully shared profiles

Define the curvature and largest order-five Gauss reciprocal by
\[
 S=\sum_i\frac{(b_i-a_i)^2}{\max(a_i,b_i)},\qquad
 x_*=\max_i\frac{|b_i-a_i|}{5\min(a_i,b_i)+1}.
\]
Use the [combined curvature theorem](MULTIVARIATE_TERMINAL_SIX_HYBRID_CURVATURE.md)
if \(S\ge4\). For \(S<4\) and \(x_*>1/8\), retain every row whose
reciprocal exceeds \(1/8\). The
[retained-row certificate](MULTIVARIATE_TERMINAL_SIX_RETAINED_ROWS.md)
handles every nonempty retained list except a single unordered
\((1,2)\) row. The
[marked-row tail certificate](MULTIVARIATE_TERMINAL_SIX_MARKED_ROW_TAIL.md)
handles that sole exception.

It remains to take \(x_*\le1/8\). At order five remove the
binomial coefficient and normalize both endpoints to one. Let
\(u\) be the square root of the product of the two outer
normalized multipliers. Thus \(0<u\le1\). The cases are

| Order-five outer mean | Linked separation |
| --- | --- |
| \(0<u\le7/25\) | [First coupled ray band](MULTIVARIATE_TERMINAL_SIX_JOINT_RAY_BANDS.md) at cap \(1/8\) |
| \(7/25\le u\le18/25\) | [Paired Bernstein projections](MULTIVARIATE_TERMINAL_SIX_PAIRED_LIGHT_ROWS.md) for \(u\le3/5\), and the [diffuse Cayley certificate](MULTIVARIATE_TERMINAL_SIX_DIFFUSE_CAYLEY.md) for \(u\ge3/5\) |
| \(18/25\le u\le9/10\) | [Intermediate Cayley certificate](MULTIVARIATE_TERMINAL_SIX_INTERMEDIATE_CAYLEY.md) |
| \(9/10\le u\le1\) | [Near-flat extension](MULTIVARIATE_TERMINAL_SIX_NEAR_FLAT_EXTENSION.md) |

The reciprocal cap uses the largest reciprocal in each row, so all
Gauss factors of both orders are within the bounds required by the
linked certificates. These intervals include their endpoints and
cover every possible \(u\).

## Projective coefficients

Each branch gives disjoint complex root arguments for the normalized
polynomials. Normalization multiplies the variable by a positive
number and the polynomial by a nonzero scalar. Thus \(Q_5\) and
\(Q_6\) have no common complex root. Their coefficients are rational,
so a rational Bézout identity extends the conclusion to every
characteristic-zero field. Both pure monomial coefficients of each
moment are nonzero. This covers the two pure projective coefficient
choices as well and proves the terminal-six case of multivariate
\(\operatorname{SFC}(2)\).

Together with the
[low-terminal theorem](MULTIVARIATE_TWO_MONOMIAL_NONDIVISIBILITY.md)
and the [terminal-five theorem](MULTIVARIATE_TERMINAL_FIVE_ANGULAR_SEPARATION.md),
this gives every terminal from two through six.
