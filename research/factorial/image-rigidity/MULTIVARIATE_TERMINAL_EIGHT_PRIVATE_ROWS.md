# Terminal-eight separation with total private degree at least two

Let \(a,b\) be distinct ordinary exponent columns. Their private
degrees are \(\sum_{b_i=0}a_i\) and \(\sum_{a_i=0}b_i\).
If the sum of these degrees is at least two, the seventh and eighth
factorial moment polynomials have disjoint root arguments under
independent positive rescalings, and are coprime over every
characteristic-zero field. The number and size of the shared
rows are unrestricted.

The [exact checker](verify_multivariate_terminal_eight_private_rows.py)
uses the [private-factor and shared-curvature method](MULTIVARIATE_TERMINAL_SEVEN_PRIVATE_ROWS.md)
at orders seven and eight. Its numerical acceptance tests are
rational, including the anchors.

## Coefficient bounds and coverage

Normalize the endpoints to one at each order. A private row
\((0,d)\) contributes \(f_{j,k}(d)=(dk)!/((dj)!)^{k/j}\)
at index \(k\), and the reverse row contributes \(f_{j,j-k}(d)\).
Log convexity of \((dk)!/(k!)^d\) in \(k\) gives
\(f_{j,k}(d)\le f_{j,k}(1)^d\). The multiplier decreases with \(d\)
by monotonicity of the digamma function in its Gamma interpolation.
Let \(S_0=\sum d_i^2/(\ell_i+d_i)\) over the shared unequal rows.
Each normalized shared multiplier is at most
\(\exp[-S_0 k(j-k)/(2j)]\).

The following high-strength boxes already prove separation for all
profiles with larger private degree at least four, with a private
exponent at least three, or with both private degrees positive and
at least one of them at least three. For these, retain respectively
the upper caps from four private unit rows, one gap of three, or
three unit gaps and an opposite unit gap. Both private degrees at
least two are also covered directly by retaining two unit gaps in
each orientation. These four boxes require no shared curvature.
The checker uses the
[terminal-eight binomial rays](MULTIVARIATE_TERMINAL_EIGHT_HIGH_CURVATURE.md)
except that the third cosine for the single gap-three cap is \(7/25\).
It rechecks the binomial anchor signs for that choice.

Up to column reversal the other allocations occur in the following
table. For each allocation the high-strength certificate applies
at \(S_0\ge S_*\), with the exact retained private factors and the
shared upper bound at \(S_*\). A box uses independent coefficient
intervals from zero to these upper bounds for the interior indices,
and has \(1+t^j\) as its anchor.

| Increasing private gaps | Decreasing private gaps | \(S_*\) |
| --- | --- | ---: |
| \((1,1,1)\) | \(()\) | \(5/4\) |
| \((1,2)\) | \(()\) | \(3/4\) |
| \((1,1)\) | \(()\) | \(13/4\) |
| \((2)\) | \(()\) | \(11/4\) |
| \((1)\) | \((1)\) | \(3\) |
| \((1,1)\) | \((1)\) | \(1\) |
| \((2)\) | \((1)\) | \(1/2\) |

For \(S_0\le S_*\) use both sides of a bound on the shared
normalized multiplier \(F_{j,k}\):
\[
 e^{-L S_0 k(j-k)/(2j)}\le F_{j,k}\le e^{-S_0 k(j-k)/(2j)}.
\]
The values of \(L\) are \(3/2\) for \(S_0\le1\), \(9/5\)
for \(1<S_0<9/4\), \(201/100\) for \(9/4\le S_0\le3\),
and \(221/100\) for \(3<S_0\le13/4\).
For a single shared row write \(r=d/\ell\) and \(s=d^2/(\ell+d)\).
The second-difference argument in the
[terminal-seven proof](MULTIVARIATE_TERMINAL_SEVEN_PRIVATE_ROWS.md#private-factors-and-shared-loss)
bounds the ratio of its normalized log loss to \(s k(j-k)/(2j)\)
by \(1+r\). This suffices if \(r\le L-1\).
Otherwise \(s>\ell(L-1)^2/L\). At each of the four stated
curvature bounds this forces \(\ell\le6\). Then \(d\ge7\)
would give \(s\ge49/13>13/4\), so \(d\le6\).
For the remaining rows the checker compares the reciprocal of
the \(j\)-th power of the normalized multiplier with the
positive Taylor sum through degree 64 for
\(\exp(Ls\,k(j-k)/2)\). In the second range it tests only
\(s<9/4\). The next range includes the endpoint \(9/4\).
Summing the row losses proves the two-sided bound.

## Ray projections and anchors

For a lower-strength allocation bound the retained coefficients by
outward rational \(j\)-th roots and combine them with the two-sided
shared bound on each curvature interval. The six ray cosines are
continuous piecewise affine functions on the quarter-unit knots in
the checker. In checking the closed radial rectangles the program
subdivides in \(S_0\) as necessary. At \(S_0=9/4\) use the larger
loss bound in the interval beginning at that point. The reference
path traverses the knot because its coefficients satisfy both
enclosures there.

The checker verifies dominance of the constant and leading terms
at radii \(1/8\) and \(8\), respectively. On the annulus the
[rational Bernstein projection test](MULTIVARIATE_TERMINAL_SIX_JOINT_RAY_BANDS.md#moving-rays-and-exact-replay)
gives a strict sign on every accepted radial cell.
In each lower-strength branch, the straight coefficient segment
to the private array times the shared factors
\(\exp[-S_0 k(j-k)/(2j)]\) stays in the certified box. Move
\(S_0\) to zero along these reference arrays and their continuous
rays. The anchor at zero is the exact private array.
Disjoint rational Rouché disks of radius \(1/1000\) isolate
its upper roots and their conjugates, and the extra negative real
root at order seven. The checker also bounds the distance from
each disk center to the six anchor rays.

In both branches, the upper sectors have root orders
\(8,7,8,7,8,7,8\). The remaining septic root is negative real.
Continuing the counts joins the first and last sectors with their
respective lower reflections. In the last joined sector the
single septic root stays real by conjugation. The two eighth-order
roots are nonreal because the eighth ordinary factorial moment
is a strictly positive Gamma integral on the real axis.
Hence the root arguments are disjoint, and rational Bézout transfers
coprimality to every characteristic-zero field.

The replay verifies 299 row-index comparisons, 762 high-strength
radial cells, and 14,982 lower-strength radial cells. Run it from the
repository root with

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_eight_private_rows.py
~~~
