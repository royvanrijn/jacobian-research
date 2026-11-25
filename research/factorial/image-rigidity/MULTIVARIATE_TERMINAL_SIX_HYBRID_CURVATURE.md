# Fully shared terminal-six separation for curvature at least four

Let \(a,b\) be distinct ordinary exponent columns whose unequal rows
are all shared. Their order-five and order-six factorial moment
polynomials have disjoint root arguments under independent positive
rescalings whenever
\[
 S:=\sum_i \frac{(b_i-a_i)^2}{\max(a_i,b_i)}\ \ge\ 4.
\]
The absolute exponent sizes and the dimension are unrestricted.

For \(S\ge6\) use the existing
[linked high-curvature proof](MULTIVARIATE_TERMINAL_SIX_PRIVATE_ANGLES.md#total-curvature-at-least-six).
We treat \(4\le S<6\) by intersecting total-curvature and adjacent attenuation caps with
the [coupled outer-mean envelopes](MULTIVARIATE_TERMINAL_SIX_JOINT_RAY_BANDS.md).
The intersection controls the endpoints and the middle coefficients
using the actual shared shape.

## Bounds forced by the integer rows

In a shared row put \(\ell=\min(a_i,b_i)>0\) and
\(d=|a_i-b_i|\ge1\). Its Gauss reciprocals are
\(x_j=d/(j\ell+h)\), \(1\le h\le d\).
For \(S<6\), the row curvature \(d^2/(\ell+d)\) gives
\(d\le6,7,8\) at \(\ell=1,2,3\) respectively. For
\(\ell\ge4\) it gives \(d\le(5/2)\ell\), since the
curvature at \(d=(5/2)\ell\) is \(25\ell/14>6\).

At order five, a row contributes the outer-mean square
\[
 u_{\ell,d}^2=\prod_{h=1}^d\frac{5\ell+h}{5\ell+4d+h}.
\]
At fixed \(\ell\) this strictly decreases with \(d\), since
every existing factor decreases and an additional factor is
strictly less than one. The whole outer mean \(u\) is the
product of the row means. The verifier checks the rational
inequalities
\[
 u_{1,6}<.020,\quad u_{1,5}<.045,\quad
 u_{1,4}<.101,\quad u_{2,6}<.101 .
\]
Thus, on the intervals of \(u\) in the table below, these
rows and their larger gaps are successively excluded.
Combining the exclusions with the curvature bounds gives
\[
\begin{array}{c|ccc}
u&\text{upper }x_5&\text{upper }x_6&
       \text{upper }d/\ell\\ \hline
[0,.020]&1&6/7&6\\
[.020,.045]&5/6&5/7&5\\
[.045,.101]&2/3&4/7&4\\
[.101,.205]&1/2&3/7&3 .
\end{array}
\]
The displayed endpoints are exact rationals. The endpoint
\(u=0\) is used only as a limiting anchor.

## Intersecting the two coefficient bounds

After endpoint normalization, the total-curvature estimate
bounds the order-\(j\) coefficient \(c_{j,k}\) by
\[
 c_{j,k}\le\binom jk
     \exp\!\left(-4\frac{k(j-k)}{2j}\right).
 \tag{1}
\]
In particular \(u\le\exp(-8/5)<.205\).
Independently, in each row of the parameter table, use
the coupled outer-mean envelopes at the separate caps
for \(x_5,x_6\). The upper log-ratio bound is
\((5/4)(5+R)/(6+R)\) with \(R\) the cap on \(d/\ell\);
the lower log-ratio bound is \(L(t)\) from the
[coupling formula](MULTIVARIATE_TERMINAL_SIX_JOINT_RAY_BANDS.md#coupling-the-outer-means)
with \(t\) the cap on \(x_5\). In each case the checker
verifies the skew and defect constants by rational
logarithmic bounds. Cap the coefficient upper bounds from
those envelopes by (1) while retaining their lower bounds.

For \(1\le k\le5\), apply the
[adjacent attenuation theorem](MULTIVARIATE_ORDINARY_ADJACENT_COEFFICIENT_ATTENUATION.md#an-all-order-endpoint-inequality)
from each end of the order-six coefficient array. If \(b_k\)
is the current upper bound on the normalized order-five
multiplier after removing \(\binom5k\), this gives the
additional upper bound
\[
 c_{6,k}\le\binom 6k\min\!\left\{
          b_k(65/66)^{4k^2},\
          b_{k-1}(65/66)^{4(6-k)^2}\right\}.
 \tag{2}
\]
Here \(S\ge4\) permits replacing \(S\) by four in the
attenuation exponent. The right endpoint term uses the
same distance \(6-k\) at both orders, so its index at
order five is \(k-1\). This comparison uses the two
orders on the same parameter interval.

## Exact ray checks

The
[checker](verify_multivariate_terminal_six_hybrid_curvature.py)
interpolates two ray cosines from its exact rational knots,
and keeps the second and fourth cosines at \(2/5\) and
\(-7/10\). The rays stay ordered. The constant and leading
terms dominate on radii \(1/8\) and \(8\). For each subinterval
of \(u\), the checker encloses the ray projections in rational
Bernstein polynomials on this radial interval, as in the
[joint-ray certificate](MULTIVARIATE_TERMINAL_SIX_JOINT_RAY_BANDS.md#moving-rays-and-exact-replay).
It bisects the radius or the outer parameter until one
projection is strictly signed in each terminal box. Each
exponential or rational-power inequality is enclosed by
the stated exact arithmetic. The accepted radial totals
are \(468,165,714,183\) at order five, and
\(588,173,435,181\) at order six.

At \(u=0\) the polynomials are \(1+z^5\) and \(1+z^6\),
whose nonreal upper arguments alternate around the four rays.
Follow the canonical coefficient exponents \(1,3/2,3/2,1\)
at order five and \((25/24)(1,8/5,9/5,8/5,1)\) at order six
from zero up to the actual \(u\). This path stays inside the
outer-mean envelopes. It also respects (1), because
\(u\le\exp(-8/5)\) along the path and each exponent equals
\(\frac{5}{4j}k(j-k)\). For (2), the differences between the
canonical exponent in order six at index \(k\) and the
canonical exponents in order five at indices \(k,k-1\) are
\(k^2/24\) and \((6-k)^2/24\), respectively. Hence its two
order ratios are bounded by
\(\exp(-k^2/15)\) and \(\exp(-(6-k)^2/15)\).
These are at most the attenuation factors in (2), since
\(\log(66/65)<1/65<1/60\). At the target \(u\), interpolate
straight to the actual coefficients inside the intersection.
Since no root crosses a ray along either path, the sector
counts stay fixed. The last full conjugate sector contains
one quintic root and it is real, while the actual sextic
factorial moment is positive on the real axis.
Thus no root arguments agree.
