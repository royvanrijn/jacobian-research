# All-degree coprimality for dense distributions of factorial gaps

Let \(a,b\) be distinct ordinary exponent columns, and omit equal rows.
Select a group of unequal rows containing both orientations. Let \(X\)
be the largest selected gap and \(M\) the largest selected upper exponent.
The selected rows have gaps \(1\le d_i\le X\) and upper exponents
\(M-e_i\). For \(\sigma\in\{+,-\}\), define \(A_\sigma(t)\) as the number
of selected rows of orientation \(\sigma\) with \(d_i\le t\).
Counts include multiplicities. Fix \(J\ge4\) at least the
[absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
Put \(H=2J^2\).

Suppose there are real densities \(\lambda_+,\lambda_->0\), and put
\(\lambda_*=\min(\lambda_+,\lambda_-)\). Assume
\[
 \max(\lambda_+,\lambda_-)\le2\lambda_*,\qquad
 0\le e_i\le\frac{X}{256J^4},\qquad
 |A_\sigma(t)-\lambda_\sigma t|
       \le\frac{\lambda_* X}{128J^3}
       \quad(0\le t\le X,\ \sigma\in\{+,-\}).                 \tag{1}
\]
In particular the selected gaps can occupy the whole range from \(1\)
to \(X\), with arbitrarily many distinct scales and multiplicities.

For each retained row, let \(M_i=\max(a_i,b_i)\) and \(d_i=|a_i-b_i|\).
Let \(\mathcal L\) be the sum of \(d_i\log(JM_i)\) over those retained
rows with \(JM_i>HX\), omitting any such rows whose two exponents
lie in \([M-X/(4J),M]\).

**Theorem.** If
\[
                 X\ge\max(144,8J),\qquad
        \log X\ge576JH,\qquad
                 \mathcal L\le\frac{X\log M}{256J^3},          \tag{2}
\]
then the original profile satisfies SFC(2) in every degree over
every characteristic-zero field. The result concerns coprimality
in the physical variable. The numerical value of \(J\) is unevaluated.

## Counting all copies of a factor

Write \(m=N\) and \(m'=N-1\) for a terminal \(4\le N<J\).
Strip the smaller-end factorial from each row in each degree.
For \(u\in[X/4,X/2]\) with \(u\equiv1\bmod m\), mark the
literal factor \(x_u=mM-u\). A selected increasing row contributes
one mark at index \(k\) exactly when
\[
                 me_i+(m-k)d_i\le u<m(e_i+d_i).
\]
For a decreasing row the condition is \(me_i+kd_i\le u<m(e_i+d_i)\).
These conditions follow by comparing \(x_u\) with the first and last
sampled factors in that row's stripped coefficient.

Write \(C_k(u)\) for the number of marks at index \(k\), from all
selected rows. Temporarily omit the strict upper inequalities in the
two conditions above, and call the resulting count \(Q_k(u)\).
Then
\[
 Q_k(u)-C_k(u)=\#\{i:m(e_i+d_i)\le u\}.
\]
Every row counted on the right satisfies its lower threshold at every
index. This difference therefore has no dependence on \(k\).
In particular, \(C_k\) and \(Q_k\) minimize at the same indices.

For \(1\le k\le m-1\), condition (1) gives
\[
 Q_k(u)=u\left(\frac{\lambda_+}{m-k}+
                    \frac{\lambda_-}{k}\right)+R_k,\qquad
 |R_k|\le\frac{\lambda_* X}{32J^3}.                           \tag{3}
\]
For example, the increasing count is between
\(A_+((u-m\max_i e_i)/(m-k))\) and \(A_+(u/(m-k))\).
Both arguments lie in \([0,X]\). The estimation error on that
orientation is at most
\(\lambda_* X/(128J^3)+\lambda_+m\max_i e_i\le\lambda_* X/(64J^3)\).
The same estimate with \(k\) in the denominator applies to the
other orientation.

The second difference of \(1/k\) is \(2/[k(k^2-1)]\).
At \(2\le k\le m-2\), the second difference of the main term in
(3) is at least \(4u\lambda_*/J^3\ge X\lambda_*/J^3\).
The second difference of \(R_k\) has absolute value at most
\(X\lambda_* /(8J^3)\). Thus \(Q_k\) is strictly convex on the
internal integer grid. For \(k=\lfloor m/2\rfloor\) its main term
is at most \(3X\lambda_*/4\), because both denominators are at
least two. At \(k=0\) the count includes all the decreasing
rows, and at \(k=m\) all the increasing rows, since \(m e_i\le u\).
Either total is at least \(\lambda_* X-\lambda_* X/(128J^3)\).
Consequently the full minimizing index set is contained in the
interior and consists of one index or two consecutive indices.
This uses only the mark counts and also works when every coefficient
has at least one mark.

## Prime factors and the two Newton polygons

If the total absolute curvature \(S\) is at least twelve, use the
[absolute curvature theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md).
Otherwise a selected row of gap \(X\) gives \(M>X^2/12\).

For a prime \(p>HX\) dividing \(x_u\), no selected factor different
from \(x_u\) at order \(m\) is divisible by \(p\). In fact its
difference from \(x_u\) is \(u-me_i-h\) for \(0\le h<md_i\),
of magnitude at most \(HX\). At order \(m'\) we have
\[
 m(m'(M-e_i)-h)=m' x_u+(m'u-mm'e_i-mh),\qquad 0\le h<m'd_i.
\]
The offset is nonzero modulo \(m\) and has magnitude at most \(HX\).
All selected factors at order \(m'\) are therefore units.
Also the omitted retained rows make no contributions at this prime.
Those with \(JM_i\le HX\) only contain smaller factors. Those in the
top strip have order-\(m\) factors strictly greater than
\(m(M-X/(4J))\ge x_u\). Their two offset bounds hold with
\(e_i+d_i\le X/(4J)\).

For \(j=m,m'\) let \(B_j\) be \(j!\) times the product of full
stripped factorial ramps of the retained rows counted in
\(\mathcal L\). Write \(E=v_p(x_u)\), \(U=v_p(B_m)\) and
\(L=v_p(B_{m'})\). The marked coefficient valuations are
\(EC_k(u)\) plus residual ordinates in \([0,U]\). Subtract
\(E\min_i C_i(u)\) from every ordinate. Outside the minimizing
index set the adjusted ordinate has at least \(E\) from the marks.
The other polygon has all slopes in \([-L,L]\). If \(E>U+mL\),
the [deep mark comparison](MULTIVARIATE_CENTRAL_DEEP_PRIME_DIVISOR_GATE.md#marks-in-arbitrary-blocks)
forces all possible common slopes into edges of total width at
most one. The even moment is strictly positive on the real axis,
so any common rational factor has even degree. This excludes a
common root for \(E>U+mL\).

It follows that at a possible common root the large-prime part of the
product of all \(x_u\) divides \((B_m B_{m'}^m)_{>HX}\), since a prime
\(p>HX\) can divide only one mark. The \(s\) marks form an arithmetic
progression of step \(m\), coprime to that step. As in the
[progression product argument](MULTIVARIATE_LARGE_FACTORIAL_MARK_PRIME_BUDGET.md#a-product-bound-at-a-possible-common-root),
delete one term of maximal valuation for each prime at most \(HX\)
that divides a mark. The remaining small-prime part divides
\((s-1)!\). Hence a common root would give
\[
       \bigl(\min_u x_u\bigr)^{s-\pi(HX)}
                     \le (s-1)!\,(B_mB_{m'}^m)_{>HX}.       \tag{4}
\]
Here \(s\ge X/(8J)\), \(s\le X\), and every mark is greater than \(M\).
By the elementary estimate \(\pi(y)\le9y/\log y\) used in the
[large-prime bound](MULTIVARIATE_LARGE_FACTORIAL_MARK_PRIME_BUDGET.md#uniform-estimates-for-the-certificate),
the log condition in (2) ensures \(\pi(HX)\le X/(64J)\le s/8\).
Since \(X\ge144\), we have \(\log s\le\log X\le(2/3)\log M\).
There are at most \(m^2d_i\) ramp factors per retained row in
\(B_mB_{m'}^m\); they are at most \(JM_i\). The factorials \(j!\)
are units at every prime above \(HX\), so
\[
 \log (B_mB_{m'}^m)_{>HX}\le J^2\,\mathcal L
                                   \le s\log M/32.
\]
The logarithms of the two sides of (4) now contradict
\(7/8>2/3+1/32\). The first two pairs use
[ordinary low-order separation](MULTIVARIATE_TWO_MONOMIAL_NONDIVISIBILITY.md),
and terminals \(N\ge J\) use the absolute cutoff. This finishes the proof
of the physical statement and its extension by the rational Bézout identity.

## Example and the underlying count criterion

For each \(d=1,\ldots,X\), take \((M-d,M),(M,M-d)\), and include
the fixed shared rows \((1,3),(3,2),(5,4)\).
Take \(M=\lfloor X^3/9\rfloor\), or any integer \(M\) between
\(X^3/12\) and \(X^3/6\). The total degree difference is zero.
The selected distributions have unit density and discrepancy at most
one, so (1) holds for sufficiently large \(X\). The fixed rows make
no contribution to \(\mathcal L\). Thus every indicated profile is
regular in every degree once \(X\) is sufficiently large.
The selected gaps have no empty interval \([D/12,D)\) for
integer \(2\le D\le X\).

The proof of (4) only requires marks in the indicated progression
with a minimizing count set of width at most one. More generally,
one can restrict to a subset of \(s\) marks with this property.
After removal of one term per small prime, the small-prime part
then divides \((s_{\rm full}-1)!\) using the span of the full
progression with \(s_{\rm full}\) terms. This also gives an exact
certificate for other integer gap distributions when the count
condition holds on sufficiently many marks.
