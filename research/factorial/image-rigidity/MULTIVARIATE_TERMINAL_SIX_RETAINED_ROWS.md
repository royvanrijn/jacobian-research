# Terminal six with large retained Gauss reciprocals

Let \(a,b\) be distinct ordinary exponent columns whose unequal rows
are all shared. In each row put \(\ell_i=\min(a_i,b_i)\),
\(d_i=|a_i-b_i|\), and \(s_i=d_i^2/(\ell_i+d_i)\). Define
\[
 S=\sum_i s_i,\qquad x_i=\frac{d_i}{5\ell_i+1}.
\]
Suppose \(S<4\) and at least one row has \(x_i>1/8\).
Retain all rows with \(x_i>1/8\), with their orientations.
If the retained list differs from a single row of unordered exponents
\((1,2)\), the order-five and order-six polynomials have disjoint root
arguments under independent positive rescalings. Arbitrarily many
other rows and unbounded exponents are allowed.

The [checker](verify_multivariate_terminal_six_retained_rows.py)
uses rational arithmetic for the enumeration, coefficient bounds,
and ray tests. It uses
[FLINT's certified polynomial root isolation](https://flintlib.org/doc/acb_poly/index.html)
through python-flint for the finite set of anchors.
The replay covers 822 oriented retained lists up to column reversal,
2,060 tail intervals, and 129,186 radial cells.

## Retain the same integer rows in both orders

If \(x_i>t\) and \(s_i\le B\), then
\[
 \ell_i<\frac{d_i}{5t}-\frac15,\qquad
 d_i < B\left(1+\frac1{5t}\right).
\]
The second bound follows by inserting the first into \(d_i^2\le
B(\ell_i+d_i)\). The checker enumerates the finite integer rectangle
given by these inequalities and then applies the exact predicates.
For any marked row of \(x_i>1/8\), reverse both columns if needed so
that this row is increasing. Enumerate the signed multiplicities
of all eligible other retained rows under the remaining curvature
budget. Each eligible row has positive cost, so there are finitely
many choices, including repetitions. Retain the sorted list up to
simultaneous column reversal, which sends every nonzero root of an
endpoint-normalized polynomial to its reciprocal and preserves the
disjointness condition. The relaxed budget permits equality on the
core after marking and thus includes all lists with \(S<4\).
After column reversal it generates 823 lists. The one exceptional
list is the unordered row \((1,2)\) by itself.

For order \(j\), remove the binomial coefficients and normalize the
two endpoints to one. If a retained row is increasing, the \(j\)-th
power of its normalized multiplier at \(k\) is
\[
 \frac{((j\ell+kd)!)^j}{((j\ell)!)^{j-k}
                         ((j(\ell+d))!)^k}.
\]
Reverse the array for a decreasing row. Multiply these rational
numbers down the list and take the positive \(j\)-th root. In this
way both orders use the same exact oriented core and its two
different arrays of factorials.

## Bound the unretained rows by their common outer mean

Write \(v\) for the unretained order-five outer mean, equal to the
square root of the product of the two outer normalized multipliers.
Use \(v=1\) if the remainder is empty. A row with \(x_i\le1/8\)
has \(8d_i\le5\ell_i+1\), hence \(d_i/\ell_i\le5/7\). If
\(r=d_i/\ell_i\), that row's contribution \(U_i\) to \(-\log v\)
satisfies
\[
 \frac{U_i}{s_i}\le
 \frac{1+r}{2r}\log(1+4r/5)\le\frac65\log(11/7).
\]
For the first bound use
\(U_i=\frac12\sum_{h=1}^{d_i}\log(1+4d_i/(5\ell_i+h))\).
For the second, the middle expression is
\(\frac25\int_0^1(1+r)/(1+4rt/5)\,dt\) and is increasing in \(r\).
If \(B=4-\sum_{\rm retained}s_i\), it follows that
\[
                  (7/11)^{6B/5}\le v\le1.
\]
Each factor in both Gauss arrays of the remainder is now controlled
by the [linked outer-mean envelopes](MULTIVARIATE_TERMINAL_SIX_JOINT_RAY_BANDS.md#coupling-the-outer-means)
at the reciprocal cap \(1/8\). The checker bounds every actual
interior multiplier on each interval of \(v\), and multiplies by the
enclosed core multipliers and the binomial coefficients. Both orders
use the same interval of \(v\), with the order-six exponents enclosing
the ratio of the two outer logarithms for all ordinary remainder
rows.

## Certify separation for every tail interval

The constant and leading coefficients remain one. All remainder
multipliers are at most one by log convexity. At the inner radius
\(1/8\) and outer radius \(8\), the checker uses the retained
multipliers alone to prove strict dominance of the constant and
leading terms, respectively. Thus every root lies between these
two radii.

For an interval of \(v\) with midpoint \(v_0\), take the anchor
obtained by multiplying each retained multiplier at \(j,k\) by
\(v_0^{5k(j-k)/(4j)}\). The checker verifies each exponent belongs
to the corresponding enclosure for the remainder. It encloses all
anchor roots with certified complex balls and orders the upper
roots by the cosine of their argument. The certified order in
increasing argument is \(6,5,6,5,6\). The quintic also has one
negative real root, by the count of both nonreal halves and
positivity of its coefficients.

In each gap between those upper anchor roots the checker finds
a rational cosine strictly between the two root enclosures.
For that ray it proves both polynomial arrays are nonzero at
every allowed radius, using the [rational Bernstein ray
test](MULTIVARIATE_TERMINAL_SIX_JOINT_RAY_BANDS.md#moving-rays-and-exact-replay)
on the whole coefficient box for each order. Sorting approximate
midpoints and proposing each cosine use floating arithmetic.
Acceptance requires strict certified comparisons for the anchors
and strict rational signs for each radial cell. If an interval
fails either test, the checker bisects it and tries again.
A depth cap makes an unresolved interval an error.

The straight segment from each anchor to its actual array stays
in those coefficient boxes, preserving both root counts in every
sector bounded by the four rays and their conjugates. The first
and last upper sectors are counted together with their lower
reflections across the real axis. The last joined sector contains
one quintic root, so it stays real by conjugation. The actual
sextic is strictly positive on the real line because it integrates
a nonzero real binomial raised to the sixth power against a
positive product exponential measure. This accounts for the two
sextic roots in the last joined sector, both nonreal, and proves
disjointness of the arguments of all roots. The original factorial
polynomials have rational coefficients, and their rational
Bézout identity gives coprimality over every characteristic-zero
field.
