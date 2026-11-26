# Single-order rigidity of ordinary factorial row profiles

Let \(A,A'\) be ordinary two-monomial profiles, ignoring rows on
which the two exponents agree. Write \(g_+(A)\) for the sum of
gaps on increasing rows, and \(g_-(A)\) for the sum on decreasing
rows. Define
\[
 D=\max\{ g_+(A)+g_-(A'),\ g_-(A)+g_+(A')\}.
\]

**Theorem.** If the degree-\(m\) factorial polynomials of \(A,A'\)
agree up to a positive scalar and positive argument change, and
\(m>D\), then the multisets of their ordered unequal rows are
equal. In particular the normalized polynomial of a single
profile is reciprocal at order \(m>2\max(g_+,g_-)\) exactly when
the row profile is mirrored. All share and private
configurations are allowed. The bound depends on gap sums and is
independent of the common exponents in the rows.

## Exact signed quotient of adjacent coefficients

In order \(m\), expand each ordinary unequal row of smaller
exponent \(\ell\) and gap \(d\) into its Gauss slots with bases
\(b_h=(m\ell+h)/d\), \(1\le h\le d\). For a coefficient with the
binomial factor removed, the ratio at indices \(k+1,k\) is, up
to a positive constant independent of \(k\), the evaluation at
\(x=k\) of
\[
                 \mathcal R_A(x)=
   \frac{\prod_{+\text{ slots of } A}(b_h+x)}
        {\prod_{-\text{ slots of } A}(b_h+m-1-x)} .
\]
This follows directly by grouping the \(d\) successive factorial
factors on each row. It holds for \(0\le k<m\).
[Exact matched Gauss grid](MULTIVARIATE_MATCHED_TANH_INTEGER_MULTIPLICITY.md)

Positive equivalence of the two polynomials gives
\(\mathcal R_A(k)=K\mathcal R_{A'}(k)\), \(0\le k<m\), for a fixed
\(K>0\). Clearing denominators yields a polynomial of degree at
most \(D\) vanishing at these \(m\) distinct points. Thus it
vanishes identically. The zeros of increasing slot factors are
strictly to the left of zero. Those of decreasing slot factors
are strictly to the right of \(m-1\). Equating multiplicities
on those two sides shows equality of the two slot multisets in
each orientation separately.

To recover rows from their slot multiset in a fixed orientation,
consider the signed row difference. Choose the largest gap
\(d\) with a nonzero difference at some smaller exponent \(\ell\).
Among the \(d\) numerators \(m\ell+1,\ldots,m\ell+d\), choose
one coprime to \(d\). Its reduced denominator is exactly \(d\).
No smaller gap can contribute it. Rows of that same gap with
different smaller exponents have disjoint slot intervals because
\(m\ge d\). The selected atom therefore has nonzero signed weight,
contrary to equality of the slot multisets. This proves row
recovery. The hypothesis \(m>D\) implies \(m\ge d\) at every
row under consideration. Taking \(A'\) to be the reflection of
\(A\) gives the reciprocal statement.

## A sign-variation refinement for reciprocity

At a fixed order \(m\ge2\), write \(G>0\) for the positive
endpoint normalization scale of a single profile. Group the
bases \(b_h\) at this order without distinguishing rows. For
each distinct \(v=b_h(b_h+m-1)\), let \(c(v)\) be the number of
increasing slots minus the number of decreasing slots at that
value. Delete zero entries, order the remaining values \(v\)
in increasing order, and let \(s\) be the number of sign changes
of the resulting list. If the normalized polynomial is reciprocal,
then either all \(c(v)=0\) or
\[
                         s\ge\left\lfloor\frac{m+1}{2}\right\rfloor-1.
\]
Hence reciprocity forces exact slot pairing whenever the signed
list has too few alternations. If additionally \(m\ge d\) on all
unequal rows, it forces the rows to be mirrored.

In fact, reciprocity requires the two endpoint-normalized
successive coefficient ratios at \(k\) and \(m-1-k\) to multiply
to one. The increasing pair of slot factors multiplies to
\(v+x\), where \(x=k(m-1-k)\), and the decreasing slot gives
the inverse. Thus
\[
                         L(x)=\sum_v c(v)\log(x+v)
\]
takes a common value at \(\lfloor(m+1)/2\rfloor\) distinct
nonnegative \(x\)'s. Between these points are at least
\(\lfloor(m+1)/2\rfloor-1\) distinct zeros of \(L'\).
For a nonzero signed list, \(L'(x)=\sum_v c(v)/(x+v)\) has
at most \(s\) distinct nonnegative zeros. To see this, suppose
there are \(s+1\) such zeros \(x_j\). Choose a real polynomial
\(P(t)\) with a simple root in each of the \(s\) gaps between
successive \(v\)'s where the signs change. Choose its overall
sign to match all \(c(v)\). The rational function
\(P(t)/\prod_j(t+x_j)\) has a partial fraction expansion in the
functions \(1/(t+x_j)\). Its sum against \(c(v)\) would be
zero by the equations \(L'(x_j)=0\). But all terms of this
sum are strictly positive, because \(v>0\) and \(x_j\ge0\).
This is a contradiction. The argument also shows that if all
unequal rows have one orientation, a degree-\(m\) reciprocal
normalization is impossible for \(m\ge3\).

When all signed slot counts vanish, the factorization in that
degree gives simple roots on a single circle by the
[mirrored fixed-order product theorem](MULTIVARIATE_MIRRORED_ROOT_CIRCLE.md#mirrored-profiles-at-every-fixed-order),
even before recovering the rows. More generally, if an ordinary
degree-\(m\) array is reciprocal up to positive scaling then
\(G^2\in\mathbb Q_{>0}\), since the endpoint-neighbor coefficient
ratios \(C_1/C_0\) and \(C_m/C_{m-1}\) have product \(G^2\).
Here \(C_k\) may be the coefficients after removing binomial
factors, which are themselves reciprocal. The endpoint ratio
also gives \(G^m\in\mathbb Q_{>0}\), so \(G\) is rational for
odd \(m\). For two adjacent slot-paired arrays a physical
common root would force their normalization scales equal, as
their original roots lie on circles of those inverse radii.
Since one order is odd, the resulting common circle then has
rational radius, permitting rational reciprocal common-factor
tests after one exact rescaling.

By contrast, the [row replacement theorem](MULTIVARIATE_EXACT_FACTORIAL_ROW_REPLACEMENTS.md)
constructs distinct mirrored and nonmirrored profiles whose
raw factorial polynomials coincide through any preassigned
finite cutoff. A bound on row complexity is necessary for
unconditional recovery from one specified order.
