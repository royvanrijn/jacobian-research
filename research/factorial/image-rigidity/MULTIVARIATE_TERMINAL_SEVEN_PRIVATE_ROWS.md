# Terminal-seven separation with total private degree at least two

Let \(a,b\) be distinct nonnegative integer exponent columns. Their
private degrees are \(\sum_{b_i=0}a_i\) and \(\sum_{a_i=0}b_i\).
If their sum is at least two, then the ordinary factorial moments at
orders six and seven are coprime over every characteristic-zero field.
In fact their root arguments are disjoint under independent positive
rescalings. Equal rows and any number of shared unequal rows are allowed.

The proof retains the private factorial contributions and bounds the
shared contribution in the same curvature interval at both orders.
The [exact checker](verify_multivariate_terminal_seven_private_rows.py)
uses only rational arithmetic during replay.

## Private factors and shared loss

Normalize both polynomial endpoints to one by dividing by the
constant coefficient and rescaling the argument positively. For a
private row \((0,d)\), \(d\ge1\), its normalized multiplier at order
\(j\) and index \(k\) is
\[
 f_{j,k}(d)=\frac{(dk)!}{((dj)!)^{k/j}}.
\]
For \((d,0)\) the index is \(j-k\). Each multiplier is at most one.
The sequence \((dk)!/(k!)^d\) is log convex in \(k\), since the
successive ratios are \(\prod_{h=1}^d(dk+h)/(k+1)\) and each factor
increases weakly with \(k\). Its endpoint chord gives
\[
 f_{j,k}(d)\le f_{j,k}(1)^d .
 \tag{1}
\]
Also \(f_{j,k}(d)\) decreases in \(d>0\) under Gamma interpolation.
Indeed the derivative of its logarithm is
\(k[\psi(dk+1)-\psi(dj+1)]<0\) for \(0<k<j\), by positivity
of the trigamma function.

For one shared unequal row put \(\ell=\min(a_i,b_i)\),
\(d=|b_i-a_i|\), \(s=d^2/(\ell+d)\). Its normalized multiplier
at index \(k\) has logarithmic loss between
\(s\,k(j-k)/(2j)\) and \((1+d/\ell)s\,k(j-k)/(2j)\).
To see this for the increasing orientation, the second difference
of the raw log factorial coefficients is
\[
 \sum_{h=1}^d
  \log\!\left(1+\frac{d}{j\ell+(k-1)d+h}\right).
\]
At every interior index it lies between \(s/j\) and \(d^2/(j\ell)\),
using \(x/(1+x)\le\log(1+x)\le x\). The quadratic that vanishes at
the endpoints and has the corresponding constant second difference
then bounds the chord-normalized logarithm. Reversing the row preserves
both bounds.

Write \(S_0\) for the sum of \(s\) over all shared unequal rows.
For orders six and seven and \(S_0\le3\), the checker improves the
upper loss bound uniformly. The shared factor \(F_{j,k}\) satisfies
\[
 \exp\!\left(-L S_0\frac{k(j-k)}{2j}\right)
 \le F_{j,k}\le \exp\!\left(-S_0\frac{k(j-k)}{2j}\right),
 \qquad
 \begin{array}{c|ccc}
   \text{upper bound on }S_0&1&2&3\\ \hline
   L&3/2&7/4&2
 \end{array}
 \tag{2}
\]
in the respective ranges. For a single row with \(d/\ell\le L-1\)
the preceding second-difference estimate suffices. Otherwise
\[
 s=\ell\,\frac{(d/\ell)^2}{1+d/\ell}>
       \ell\,\frac{(L-1)^2}{L}.
\]
For each of the three pairs of bounds in (2), \(\ell\ge7\)
would contradict \(s\le S_0\). When \(\ell\le6\), \(d\ge7\)
would give \(s\ge49/13>3\). For the remaining rows the exact
test compares the reciprocal of the \(j\)-th power of the
normalized factor with the lower Taylor sum for
\(\exp(Ls\,k(j-k)/2)\) through degree 64. This proves the
upper loss bound for each row. The bounds then add over rows.

## Complete coverage by the private allocations

Exchange the columns if necessary so that the larger private degree
is on the increasing side. The checker treats the following
private allocations at the indicated cutoff \(S_*\).
An empty allocation is written as \(()\).

| Increasing gaps | Decreasing gaps | \(S_*\) |
| --- | --- | ---: |
| \((1,1,1)\) | \(()\) | \(1\) |
| \((1,2)\) | \(()\) | \(1/2\) |
| \((1,1)\) | \(()\) | \(3\) |
| \((2)\) | \(()\) | \(3\) |
| \((1)\) | \((1)\) | \(3\) |
| \((1,1)\) | \((1)\) | \(2\) |
| \((2)\) | \((1)\) | \(2\) |

For \(S_0\ge S_*\) the checker keeps only the upper coefficient
caps, using the right side of (2), whose proof holds without an
upper bound on \(S_0\). The one-sided allocation \((2)\) uses
the \((1,1)\) cap by (1). The same proof covers every profile with
one private degree at least four by keeping the fourth power of
\(f_{j,k}(1)\) and discarding all other curvature. It also covers
every profile with a single private exponent at least three by
keeping \(f_{j,k}(3)\). These tests use the five rational rays and
the binomial anchors \(1+t^j\) of the
[terminal-seven high-curvature proof](MULTIVARIATE_TERMINAL_SEVEN_HIGH_CURVATURE.md).
The cap is checked by raising to the \(j\)-th power, applying a
lower exponential Taylor sum to the reciprocal exponential, and
rounding the resulting root bound outward. Every ray is certified
on the independent boxes from zero to the caps.

For \(S_0\le S_*\) retain the exact private allocation and enclose
its normalized coefficients, including the binomial factor, by
outward rational \(j\)-th roots. On each interval of shared
curvature use both sides of (2) and multiply by the retained
coefficient intervals. The rational bounds on the exponential use
the positive Taylor sum through degree 64 and its geometric tail.
Each of the seven allocations has its own five ray cosines,
continuous and piecewise affine in \(S_0\), given by the quarter-unit
knots in the checker. The two endpoint coefficients stay one.

The cases with both private degrees at least two are also covered
by the [all-order theorem](MULTIVARIATE_TWO_MONOMIAL_NONDIVISIBILITY.md#all-orders-when-both-private-exponent-sums-are-at-least-two).
When the private degrees are \(1,3\), the private degree-three
upper caps apply. If that side has no exponent at least three,
use its \((1,1,1)\) or \((1,2)\) allocation and count the
opposite private linear row in the curvature lower bound, since
that row has strength one and obeys the same lower loss estimate.
This exhausts the cases with total private degree at least two.

## Ray checks and sector counts

The checker uses the
[Bernstein projection test](MULTIVARIATE_TERMINAL_SIX_JOINT_RAY_BANDS.md#moving-rays-and-exact-replay)
for each coefficient box at each order. At radii \(1/8\) and \(8\)
the constant and leading terms, respectively, dominate strictly.
These inequalities use the private core alone in the lower-curvature
branches because the shared normalized multipliers are at most one.
On the annulus, the checker projects onto the directions obtained
from the sine-ratio recurrence. A radial cell is accepted only
when some projected polynomial has strictly signed Bernstein
coefficients. An unresolved interval of shared curvature is
bisected, with a depth-limit failure reported as an error.

In a lower-curvature branch, at the actual \(S_0\) first deform
each polynomial along a straight coefficient segment to the
private core times the multipliers \(\exp(-S_0 k(j-k)/(2j))\).
Both arrays and their segment satisfy the certified box constraints.
Then decrease \(S_0\) to zero along these reference arrays and
their continuous rays. The anchor at zero is the private core.
The checker isolates its upper roots and their conjugates by
disjoint disks of rational radius \(1/1000\), with one additional
negative real disk for the septic. For the polynomial \(P\) and
disk center \(c\), its exact interval checks give
\[
 |P'(c)|r>|P(c)|+\sum_{m=2}^{j}|P^{(m)}(c)|r^m/m!,
 \qquad r=1/1000 .
\]
Rouché's theorem gives one root in every disk. The certified
distance from each upper center to the anchor rays puts its disk
in the indicated sector. The successive upper orders are
\(7,6,7,6,7,6\), with exactly one root per sector. The higher-curvature
anchors have the same assignment by the binomial sign test.

Along either deformation the first and last sectors are counted
together with their respective lower reflections. At the actual
polynomials, the last joined sector has one septic root, which
stays real by conjugation, and two sextic roots. The actual
sextic has no real roots because it is the positive Gamma integral
of the sixth power of a nonzero real binomial. This proves
disjointness of the arguments, including in the last joined sector.
The original factorial polynomials have rational coefficients and
their Bézout identity proves coprimality over characteristic zero.

The replay checks 176 row-index inequalities, 404 radial cells
in the upper-cap branches, and 6,288 radial cells in the seven
lower-curvature branches.
