# Binomial weighting can destroy interlaced-ratio coprimality

Consider the two positive rational arrays
\[
     a=(1,3/5,3/5,1),\qquad
     b=(1,7/12,121/225,7/12,1).
\]
The degree-three and degree-four polynomials with those exact
coefficients satisfy the strict interlaced ratio hypothesis of
the [endpoint descent theorem](MULTIVARIATE_BOTH_PRIVATE_EUCLIDEAN_DESCENT.md).
Both coefficient arrays are truncated moments of positive
measures supported on \((0,\infty)\). The forward and reversed
constant-normalized coefficient quotient arrays are also
truncated positive moments. Nevertheless the binomial-weighted
polynomials share an irreducible quadratic.

The successive ratios are
\[
  (3/5,1,5/3) \quad\text{and}\quad
                (7/12,484/525,525/484,12/7).
\]
They strictly alternate with an upper ratio first. Both arrays
are reciprocal, so their forward and reversed quotients are
the same list,
\[
                 (1,\ 36/35,\ 135/121,\ 12/7).
\]
For the two degree-three arrays \(a\) and \(a/b\), the leading
two-by-two Hankel matrix and its shift by one are strictly
positive definite. For \(b\), the three-by-three Hankel matrix
and the shifted two-by-two matrix are also strictly positive
definite. To see this for the larger matrix, write
\(u=7/12,\ v=121/225\).
The subspace of the vector \((1,0,-1)\) has eigenvalue \(1-v>0\).
On its orthogonal complement the matrix has diagonal entries
\(1+v,v\) and off-diagonal entry \(\sqrt2 u\). Its determinant
\(v(1+v)-2u^2\) is positive. The shifted matrix has eigenvalues
\(u\pm v>0\).

Here the elementary sufficient truncated moment criterion
can be seen directly. For an odd endpoint index \(2r+1\),
let \(H=(m_{i+j})_{0\le i,j\le r}\) and
\(K=(m_{i+j+1})_{0\le i,j\le r}\) be positive definite.
On real polynomials of degree at most \(r\) use the inner
product with Gram matrix \(H\). Define the self-adjoint positive
operator \(J\) by
\(\langle J f,g\rangle=\sum f_i g_j m_{i+j+1}\).
Then \(J^k1=x^k\) for \(0\le k\le r\), by testing multiplication
against all polynomials in the space. The spectral measure of
\(J\) at the vector \(1\) has moments \(m_0,\ldots,m_{2r+1}\).
All its support points are positive. For an even endpoint
index \(2r\), the shifted smaller Hankel matrix is the leading
block of \(K\). Choose the last unknown moment sufficiently
large to make \(K\) positive definite by the Schur complement.
This reduces to the odd case and proves the sufficiency claims.

However the binomial-weighted arrays give the two polynomials
\[
\begin{aligned}
 \sum_{k=0}^3\binom3k a_k z^k
    &=(1+z)(1+(4/5)z+z^2),\\
 \sum_{k=0}^4\binom4k b_k z^k
    &=(1+(23/15)z+z^2)(1+(4/5)z+z^2).
\end{aligned}
\]
Thus interlaced ratios and positive Mellin quotients alone
cannot justify reintroducing the binomial weights in the
both-private descent. The example concerns relaxed positive
moment arrays. Its arrays are not asserted to be an ordinary
factorial pair. The missing structural step for the general
binomial case must use further information from the actual
row products.
