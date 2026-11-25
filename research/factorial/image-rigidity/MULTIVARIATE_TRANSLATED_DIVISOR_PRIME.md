# A divisor-prime criterion for translated factorial binomials

Fix \(N\ge3\) and a nonempty list of exponent rows with smaller
entries \(A_i\ge0\), gaps \(d_i\ge1\), and orientations. At integer
translation \(T\ge0\), the entries of row \(i\) are \(T+A_i\) and
\(T+A_i+d_i\), in its prescribed order. A row with equal entries
contributes a nonzero scalar at every order and may be dropped.
For \(j=N-1,N\), put \(\ell_{j,i}(k)=k\) in an increasing row and
\(\ell_{j,i}(k)=j-k\) in a decreasing row. Dividing the ordinary
moment by \(\prod_i (j(T+A_i))!\) gives
\[
 F_j(T,z)=\sum_{k=0}^j\binom jk z^k
         \prod_i\prod_{h=1}^{d_i\ell_{j,i}(k)} (j(T+A_i)+h)
         \in\mathbb Z[T,z].
\]

Choose one of the adjacent orders \(J\), let \(J'\) denote the
other, and set
\[
 E_J(k)=\max_i\{JA_i+d_i\ell_{J,i}(k)\}.
\]
Choose an integer \(x>J\max_i A_i\) such that \(J\nmid x\) and the
set \(I_x=\{k\in\{0,\ldots,J\}:E_J(k)<x\}\) is nonempty with
diameter at most one. Define the nonzero integer
\[
 B=N!\prod_i\prod_{\substack{1\le h\le Jd_i\\ JA_i+h\ne x}}
       (JA_i+h-x)\,
    \prod_i\prod_{h=1}^{J'd_i}\bigl(J(J'A_i+h)-J'x\bigr).
 \tag{1}
\]

**Theorem.** At every \(T\ge0\) for which \(JT+x\) has a prime
divisor not dividing \(B\), the ordinary factorial moments of
orders \(N-1,N\) are coprime over every characteristic-zero field.
There are at most \(N(N-1)\sum_i d_i\) integer translations at
which the pair has a nonconstant common factor. Every such
translation has \(JT+x\) supported on the prime divisors of \(B\).

## Proof

The last factors in (1) are nonzero because \(\gcd(J,J')=1\) and
\(J\nmid x\). If \(p\mid JT+x\) and \(p\nmid B\), then \(p>N\),
and every factor \(J(T+A_i)+h\) is a unit modulo \(p\), except when
\(JA_i+h=x\). That exceptional value of \(h\) is positive by the
choice of \(x\), and it occurs in row \(i\) at index \(k\) precisely
when \(JA_i+d_i\ell_{J,i}(k)\ge x\). Therefore the reduction of
\(F_J\) is supported exactly on \(I_x\), since the binomial
coefficients are units. The reduction of \(F_{J'}\) has unit
endpoints. In fact its factors, multiplied by the unit \(J\), are
congruent modulo \(p\) to the last factors of (1).

A rational common divisor can be chosen integral and primitive at
\(p\). Its endpoint coefficients are units. To see this, choose a
primitive quotient in \(\mathbb Z_{(p)}[z]\) for its division into
\(F_{J'}\). Primitivity of the product follows by reduction in the
domain \(\mathbb F_p[z]\), so the residual scalar in the product
is a unit. The product has unit endpoints, which forces both
endpoints of the divisor to be units. Its reduction hence retains
its degree and divides the reduction of \(F_J\) after the initial
power of \(z\) is removed. That remaining polynomial has degree
at most one. One of the original ordinary moments has even order.
For real \(z\) it is the integral of an even power of a nonzero
real binomial against the positive product exponential measure, so
it has no real zero. Thus a nonconstant rational common divisor
would have even degree, which proves coprimality. Rational Bézout
extends it to every characteristic-zero coefficient field.

For any prime \(p\nmid B\), the congruence \(JT+x\equiv0\pmod p\)
has a nonnegative integer solution. At that translation the
resultant is nonzero. Consequently
\(\operatorname{Res}_z(F_{N-1},F_N)\in\mathbb Z[T]\) is nonzero.
Write \(D_+\) and \(D_-\) for the sums of the increasing and
decreasing gaps. The coefficient at \(z^k\) in \(F_j\) has
\(T\)-degree \(jD_-+k(D_+-D_-)\). The resultant is homogeneous of
degrees \(N\) and \(N-1\) in the two coefficient sets. It is also
weighted homogeneous of degree \(N(N-1)\) when a coefficient at
index \(k\) is given weight \(k\). The latter follows by replacing
\(z\) by \(sz\) in both polynomials in the product expression for
the resultant. Consequently each resultant monomial has
\(T\)-degree at most \(N(N-1)(D_++D_-)\). This proves the bound on
exceptional translations. The first part shows that only prime
divisors of \(B\) can occur in \(JT+x\) at those translations.

## When the narrow interval exists

Let \(m=\min_k E_J(k)\). The maximum over the increasing rows
is strictly increasing with steps at least one, and the maximum
over the decreasing rows is strictly decreasing with steps at
least one, when the respective rows are present. Their maximum
has at most two minimizing integer indices, which are adjacent.
Indeed each of the two branches can attain \(m\) at most once,
and the sublevel sets of the maximum are intervals. Also
\(m\ge J\max_i A_i\). Thus \(x=m+1\) meets all the conditions
except possibly \(J\nmid x\).

For every \(J\ge4\), a choice of \(x\) exists in the following
profile classes.

- **All gaps are one.** If increasing and decreasing rows both
  occur, let \(A_+,A_-\) be the largest bases of the respective
  orientations. Then \(E_J(k)=\max(JA_++k,JA_-+J-k)\).
  The branch with the larger base dominates throughout when the
  bases differ, giving \(m\equiv0\pmod J\). If the bases agree then
  \(m\equiv\lceil J/2\rceil\pmod J\). A missing orientation also
  gives \(m\equiv0\pmod J\). In each case \(x=m+1\) works.
- **Every gap is at least two.** Choose a minimizer \(k_0\) and an
  active row there. If it is increasing, all indices to the right
  have \(E_J\ge m+2\). For indices left of \(k_0\), each increasing
  row has value less than \(m\), so the decreasing branch supplies
  \(E_J\). At every index at most \(k_0-2\) that branch is at least
  \(m+2\), by comparison with \(k_0-1\). If \(k_0-1\) is missing
  there are no indices to its left to check. The argument with an
  active decreasing row is the reversal. Thus \(x=m+2\) also gives
  width at most one, and one of \(m+1,m+2\) is not divisible by \(J\).
- **The rows occur in mirror pairs.** In this case write
  \(E_J(k)=\max_i(JA_i+d_i\max(k,J-k))\) using one representative
  of each pair. For even \(J\), the value \(m\) at \(k=J/2\) is
  congruent to either \(0\) or \(J/2\) modulo \(J\), so \(m+1\)
  works. For odd \(J\), let \(q=(J+1)/2\). The two central indices
  give \(m=\max_i(JA_i+d_iq)\). If \(J\mid m+1\), any row pair
  attaining this maximum has \(d_i\equiv-2\pmod J\), hence \(d_i\ge3\).
  At every noncentral index this pair gives \(E_J\ge m+3\). Then
  \(x=m+2\) works.

The criterion permits either adjacent order and arbitrarily many
rows contributing to the support cut. The searched primes are
divisors of a linear function of the translation. They may be
smaller than both orders' maximal factorial arguments. This scope
complements the physical intervals in the
[single marked-row](MULTIVARIATE_MARKED_ROW_VARIABLE_PRIME.md) and
[opposite marked-row](MULTIVARIATE_OPPOSITE_MARKED_ROW_PRIME.md)
criteria.

## Examples

At terminal \(17\), take the two mirror pairs with \((A,d)=(2,3)\)
and \((1,6)\). The choice \(x=72\) gives \(I_x=\{8,9\}\) at
\(J=17\). At \(T=1717\), \(17T+72=29\cdot1009\). The normalized
order-\(17\) reduction at \(1009\) is \(260z^8+260z^9\), and both
order-\(16\) endpoint coefficients are \(336\) modulo \(1009\).
The maximal factorial argument is \(29308\), so this prime
lies well below both physical marked-row windows.

The narrow interval hypothesis can fail at both orders. For odd
\(N=2d+1\ge5\), take the increasing gap \(d\) and decreasing gap one
with the same base \(A\). At order \(N\) the smallest envelope value
is \(NA+N-1\), attained at indices \(1,2\). The envelope at index
zero equals \(NA+N\). At order \(N-1=2d\), the minimum is
\((N-1)A+N-2\) at index \(1\), and both neighboring values equal
\((N-1)A+N-1\). In either order the only nonempty width-one
cutoff is \(x=j(A+1)\), which is divisible by \(j\).
