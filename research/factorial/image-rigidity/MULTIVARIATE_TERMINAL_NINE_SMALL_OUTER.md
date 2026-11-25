# Terminal nine for small outer mean and light shared rows

Let every unequal row be shared and write \(a_i\) for its
smaller exponent and \(d_i>0\) for its gap. Let \(u\) be the
order-eight outer geometric mean of the endpoint-normalized
factorial multipliers with the binomial coefficient removed.
If \(d_i/(8a_i+1)\le1/8\) in every row and \(0<u\le19/200\),
then the order-eight and order-nine moments have disjoint
root arguments under independent positive rescalings. In
particular they are coprime in characteristic zero.

## Joint bounds from the same outer mean

Let \(F_{j,k}\) be the normalized multiplier without the
binomial coefficient. The
[general mean and skew formulas](MULTIVARIATE_TERMINAL_SEVEN_SHARED_HIGH_CURVATURE.md#joint-coefficient-bounds)
bound \(-\log F_{j,k}\) as a multiple of \(-\log u\).
The paired reciprocals are \(x,x'=8x/(9-\gamma x)\) with
\(0<\gamma\le1\). The checker encloses the corresponding
outer-log ratio by the exact rational log bounds from that
proof. Its Gaussian anchor ratio is \(64/63\).
For the skew bound it verifies nonnegative second derivatives
at both orders in \(s=\tfrac12\log(1+(j-1)x)\), by checking
the sign of the rational derivative expressions on \(0\le x\le1\).
The reflected mean bounds follow from the concavity of the
single-factor artanh expression.

Integrality and the reciprocal cap give \(d_i\le a_i\).
For \(r_i=d_i/a_i\), the row's contribution \(U_i\) to
\(-\log u\) satisfies
\[
 \frac{U_i}{d_i^2/(a_i+d_i)}
       \le\frac{1+r_i}{2r_i}\log(1+7r_i/8)
       \le\log(15/8).
\]
The middle expression equals
\(\frac7{16}\int_0^1(1+r_i)/(1+7r_it/8)\,dt\),
whose integrand increases with \(r_i\).
Thus \(S\ge-\log u/\log(15/8)\). The exact comparison
\(\log(15/8)<16/25\) and the
[total curvature cap](MULTIVARIATE_FIXED_TERMINAL_CURVATURE.md#large-curvature-and-the-endpoint-radii)
give \(F_{j,k}\le u^{25k(j-k)/(32j)}\).
The [adjacent endpoint attenuation](MULTIVARIATE_ORDINARY_ADJACENT_COEFFICIENT_ATTENUATION.md#an-all-order-endpoint-inequality)
and \((153/152)^{96}>15/8\) also give
\[
 F_{9,k}\le u^{k^2/96} F_{8,k},\qquad
 F_{9,9-k}\le u^{k^2/96} F_{8,8-k}.
\]
The checker intersects all these bounds before testing each ray.

## Gaussian anchors and exact ray replay

On a tested \(u\)-interval choose its midpoint \(u_0\) and
symmetric reference multipliers
\(u_0^{k(8-k)/7}\) and \(u_0^{(64/63)k(9-k)/8}\).
These exponents are in the joint envelopes and meet the
total curvature cap. For each fixed distance \(k\) from
an endpoint, the exponent difference is \(k^2/63\), which
passes the adjacent attenuation bound. Thus each
reference array and actual array is in the same convex
coefficient box for its order.
The outer caps give strict endpoint dominance at radii
\(1/8\) and \(8\) for the boxes.

The [checker](verify_multivariate_terminal_nine_small_outer.py)
uses certified complex balls to isolate all midpoint reference
roots. Their upper arguments alternate \(9,8,9,8,9,8,9,8\).
Choose rational ray cosines inside the seven gaps. The
[rational Bernstein projection tests](MULTIVARIATE_TERMINAL_SIX_JOINT_RAY_BANDS.md#moving-rays-and-exact-replay)
give strict ray exclusion over the entire admissible annulus.
If those tests are inconclusive, halve the \(u\)-interval.
There are 28 accepted outer-mean cells with 3,164 radial
cells. Straight continuation inside the coefficient boxes
preserves the sector counts. The unmatched ninth-order
real root is negative. The eighth-order moment is strictly
positive on the real axis by its even positive-factorial
integral, which completes the argument separation at the
two real boundary sectors.
