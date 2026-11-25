# Terminal-seven separation for one private linear row

Let \(a,b\) be distinct ordinary exponent columns with a single
unequal private row \((0,1)\) or \((1,0)\). All other unequal rows
are shared, with both exponents positive. The sixth and seventh
factorial moment polynomials have disjoint root arguments after
independent positive rescalings, and are coprime over every
characteristic-zero field.

Together with the
[larger private-degree proof](MULTIVARIATE_TERMINAL_SEVEN_PRIVATE_ROWS.md),
this covers every ordinary profile with at least one private row.
The [exact checker](verify_multivariate_terminal_seven_private_linear.py)
uses rational ray and anchor inequalities. The shared rows may have
unbounded exponents and their number is unrestricted.

## Shared reciprocal caps

Exchange columns if necessary so that the private row is \((0,1)\).
Write \(S_0=\sum_i d_i^2/(\ell_i+d_i)\) over shared unequal rows,
where \(\ell_i\) is the smaller exponent and \(d_i\) is the gap.
After endpoint normalization, the binomial coefficient times the
private multiplier is
\[
 B_{j,k}=\binom jk\frac{k!}{(j!)^{k/j}}.
\]
In the range \(S_0\ge11/2\), multiply it by the shared upper bound
\(\exp(-S_0 k(j-k)/(2j))\) from the
[shared log-curvature inequality](MULTIVARIATE_TERMINAL_SEVEN_PRIVATE_ROWS.md#private-factors-and-shared-loss).
The boxes from zero to these interior caps pass the fixed rays of the
[high-curvature certificate](MULTIVARIATE_TERMINAL_SEVEN_HIGH_CURVATURE.md).
The checker makes this additional high-strength test at \(S_0=11/2\),
which bounds every larger value.

For \(S_0\le11/2\), let \(u,w\) be the order-six and order-seven
outer geometric means of the shared multipliers, respectively.
The common Gauss factors of a shared row at the two orders have
reciprocals
\[
 x=\frac{d}{6\ell+h},\qquad
 x'=\frac{d}{7\ell+h},\qquad 1\le h\le d.
\]
The row's squared contribution to \(u\) is
\[
 u_{\ell,d}^{\,2}=
       \prod_{h=1}^d\frac{6\ell+h}{6\ell+5d+h}.
\]
Hence \(u\le u_{\ell,d}\) for every shared row, regardless of its
orientation.

Set \(t=1/40\). A row with \(d/(6\ell+1)>t\) satisfies
\(\ell<(d/t-1)/6\), hence
\[
 d < \frac{11}{2}\left(1+\frac1{6t}\right)
\]
by its curvature bound. The checker enumerates that finite
integer region and retains exactly the rows satisfying both
inequalities \(d/(6\ell+1)>t\) and \(d^2/(\ell+d)\le11/2\).
For an outer-mean interval \([u_-,u_+]\) it uses only candidates
with \(u_{\ell,d}^{\,2}\ge u_-^2\).
The reciprocal cap \(c_j\) at each order \(j\) is the larger of
\(t\) and the largest \(d/(j\ell+1)\) among these candidates.
Rows absent because their sextic reciprocal is at most \(t\)
also have septic reciprocals at most \(t\).
In particular these caps apply on each interval without a bound
on the size or number of the other rows. The complete list has
2,262 candidate rows, and every checked reciprocal is at most one.

## The two normalized coefficient arrays

For one increasing Gauss factor at order \(j\), put \(m=j-1\),
\(g_h=\log(1+h x_j)\), \(L_k=\sum_{h=0}^{k-1}g_h\), and
\(s=\log(1+m x_j)/2\). Its normalized log multiplier is
\(L_k-(k/j)L_j\). Put \(r=\min(k,j-k)\).
The average of the two reflected logarithmic losses is
\[
 \frac12\sum_{h=0}^{r-1}(g_{m-h}-g_h).
\]
The difference \(g_{m-h}-g_h\) is
\(2\operatorname{artanh}((1-2h/m)\tanh s)\).
For \(0<a\le1\) the function
\(\operatorname{artanh}(a\tanh s)\) is concave on \(s\ge0\),
since its derivative
\(a(1-\tanh^2s)/(1-a^2\tanh^2s)\) decreases.
Consequently the ratio of the averaged loss to \(s\) is
between its value at the reciprocal cap \(c_j\) and its
limit \(k(j-k)/(j-1)\) at \(s=0\).

The local skew magnitude is
\[
 G_{j,r}(x_j)=\sum_{h=0}^{j-1}
       \frac{\mathbf 1_{\{r\le h<j-r\}}-(j-2r)/j}{2}\,g_h .
\]
Indeed the signed log skew in an increasing row is the negative
of \(G_{j,r}\), by subtracting the reflected normalized logarithms.
Write the rational weights in this formula as \(\omega_h\).
For \(\tau=\log(1+m x_j)\), the second derivative of \(G_{j,r}\)
in \(\tau\) has the sign of the polynomial
\[
 \sum_{h=1}^{m-1}\omega_h h(m-h)
       \prod_{\substack{1\le l<m\\l\ne h}}(1+l x_j)^2 .
\]
This follows because the second derivative of \(g_h\) is
\((1+m x_j)h(m-h)/(m^2(1+h x_j)^2)\).
The checker expands the numerator exactly. Where a coefficient
\(b\) of degree \(i\ge2\) is negative, it substitutes
\(b x_j^{i-2}\le b x_j^i\) on \([0,1]\); the resulting
lower polynomial has nonnegative coefficients in each case.
Also \(G_{j,r}(0)=0\) and its first derivative is zero there
because \(\sum_h h\omega_h=0\). Thus \(G_{j,r}\ge0\) and its
ratio to \(s\) increases. Reversing a row changes the skew sign.
The triangle inequality gives a bound on the total skew by its
capped ratio times the total outer logarithm.

For the outer-log ratio between orders, observe that
\[
 \frac{6x}{7}\le x'=\frac{6x}{7-xh/d}\le\frac{6x}{7-x},
 \qquad 1+\frac{36}{7}x\le1+6x'\le\frac{1+5x}{1-x/7}.
\]
The ratio \(\log(1+(36/7)x)/\log(1+5x)\) decreases, by concavity
of \(s\mapsto\log(1+(36/35)(e^s-1))\) and its value zero at zero.
The ratio \(1-\log(1-x/7)/\log(1+5x)\) increases, since
\(-\log(1-x/7)\) is convex and \(\log(1+5x)\) is concave,
with both vanishing at zero. Evaluating at the sextic reciprocal
cap therefore gives lower and upper bounds on \(\log w/\log u\).
The same bounds hold after summing the logarithms over factors.
When there are no shared unequal rows both means are one and every
shared coefficient multiplier is one.

The checker combines these comparisons to enclose every shared
multiplier as \([u^b,u^a]\) for rational exponent bounds
\(0<a\le b\). At the sextic use the reflected mean-loss and
skew bounds directly. At the septic multiply the lower length
bound by the lower outer-log ratio and the upper length bound by
its upper ratio. Multiplication by the private coefficient
\(B_{j,k}\), enclosed by outward rational \(j\)-th roots, gives
the final coefficient boxes. The checker verifies that the
reference exponents
\[
 \frac{k(6-k)}5\quad\text{and}\quad
 \frac{36}{35}\frac{k(7-k)}6
\]
belong to the respective intervals at each order and cap used.
The constant and leading coefficients are always one.
The logarithm bounds use the
[positive power series with geometric tail](MULTIVARIATE_TERMINAL_SIX_JOINT_RAY_BANDS.md#moving-rays-and-exact-replay).
The checker certifies all rational power enclosures by integer
power comparisons.

## Ray exclusions and the reference path

The five ray cosines vary continuously by rational affine
interpolation on the checker knots for \(0\le u\le1\).
The private-only coefficients bound the boxes from above. They
give strict dominance of the constant and leading terms on
the circles of radii \(1/8\) and \(8\), respectively. The
[Bernstein ray test](MULTIVARIATE_TERMINAL_SIX_JOINT_RAY_BANDS.md#moving-rays-and-exact-replay)
then excludes roots throughout the annulus at each ray and both
orders. Each accepted radial interval has a strict rational
projection sign. If needed, subdivide the outer-mean interval
and recompute its ordinary row caps.

At fixed \(u\) deform the actual arrays by a straight segment to
the arrays consisting of \(B_{j,k}\) times the displayed reference
powers of \(u\). Both endpoints and the segment lie in the
coefficient boxes. Along these reference arrays, decrease \(u\)
to zero using the continuous rays. All reference exponents are
positive at interior indices, so the anchor polynomials at zero
are \(1+t^6,1+t^7\). The checker verifies their sector counts by
the [rational binomial sign test](MULTIVARIATE_TERMINAL_SEVEN_HIGH_CURVATURE.md).
The high-strength branch has the same anchors and sector counts
on its fixed rays. In both branches the six successive upper
sectors have root orders \(7,6,7,6,7,6\).

When continuing the counts, join the first and last sectors with
their respective lower reflections. The last joined sector
contains one real septic root, by conjugation, and two sextic
roots. The actual sextic is strictly positive for every real
coefficient ratio, since it is the Gamma integral of an even
power of a nonzero binomial. Therefore the sextic roots in the
last sector are nonreal, and root arguments are disjoint.
Positive normalizations preserve arguments, and rational Bézout
transfers the coprimality to all characteristic-zero fields.

The replay verifies 52 high-strength radial cells and 5,551
cells in the outer-mean intervals.
