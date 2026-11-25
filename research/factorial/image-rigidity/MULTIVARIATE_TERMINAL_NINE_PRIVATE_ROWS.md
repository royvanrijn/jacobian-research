# Terminal-nine separation with total private degree at least two

Let \(a,b\) be distinct ordinary exponent columns. The private degree on
each side is the sum of its entries in rows where the opposite entry is zero.
When the total private degree is at least two, the eighth and ninth
factorial moment polynomials have disjoint root arguments under independent
positive rescalings. They are coprime over every characteristic-zero field.
There is no restriction on the shared rows.

## Private and shared coefficient bounds

Normalize the constant and leading coefficients to one in each order.
A private increasing row of gap \(d\) contributes
\[
 f_{j,k}(d)=\frac{(dk)!}{((dj)!)^{k/j}},\qquad j\in\{8,9\}.
\]
A decreasing row contributes \(f_{j,j-k}(d)\). The sequence
\((dk)!/(k!)^d\) is log convex in \(k\), since its consecutive ratios
are products of the increasing factors \((dk+h)/(k+1)\) for
\(1\le h\le d\). It follows that \(f_{j,k}(d)\le f_{j,k}(1)^d\).
Each private multiplier is at most one and decreases as \(d\) increases,
as follows by differentiating its Gamma interpolation and using the
increasing digamma function.

For the shared unequal rows with lower base \(\ell\) and gap \(d\), put
\(S_0=\sum d^2/(\ell+d)\). Their combined normalized coefficient
multiplier is at most \(\exp[-S_0 k(j-k)/(2j)]\) by the
[shared second-difference bound](MULTIVARIATE_TERMINAL_SEVEN_PRIVATE_ROWS.md#private-factors-and-shared-loss).
In the high-strength test each interior coefficient varies independently
between zero and its cap formed from this shared bound, the retained
private factors, and \(\binom jk\). The following three zero-curvature
caps cover all profiles with a larger private degree at least four,
both private degrees at least two, or a private degree at least three
with a positive opposite private degree. Use four increasing unit
factors, two units in each direction, or three increasing units and
one decreasing unit, respectively. All these caps allow discarding
additional private factors.

Up to column reversal, the remaining allocations and their exact
high-strength cutoffs are as follows. The table also records
separate certificates for the two mixed allocations of degrees three and one.

| Increasing gaps | Decreasing gaps | High-strength cutoff on \(S_0\) | Below the cutoff |
| --- | --- | ---: | --- |
| \(3\) | — | \(1/8\) | Shared curvature |
| \(1,1,1\) | — | \(1\) | Shared curvature |
| \(1,2\) | — | \(1/2\) | Shared curvature |
| \(1,1\) | — | \(25/8\) | Shared outer mean |
| \(2\) | — | \(21/8\) | Shared curvature |
| \(1\) | \(1\) | \(3\) | Shared outer mean |
| \(1,1\) | \(1\) | \(3/4\) | Shared curvature |
| \(2\) | \(1\) | \(1/4\) | Shared curvature |

## Two-sided shared-curvature boxes

In the six branches marked shared curvature, use
\[
 e^{-L S_0 k(j-k)/(2j)}\le F_{j,k}\le e^{-S_0 k(j-k)/(2j)}
\]
for the combined shared multiplier \(F_{j,k}\). The successive
upper strength endpoints are
\(1/8,1/4,3/8,7/8,5/4,17/8,25/8\), and the associated \(L\)'s
are \(11/10,6/5,5/4,29/20,3/2,9/5,21/10\).
For a row of strength \(s=d^2/(\ell+d)\), the corresponding
upper loss-to-strength ratio is at most \(1+d/\ell\) by the same
second-difference bound. If \(d/\ell>L-1\), then
\(s>\ell(L-1)^2/L\), so only bases \(\ell\le S_0L/(L-1)^2\)
need extra checks. The gap loop stops when \(d^2/(\ell_{\max}+d)\)
exceeds the applicable strength endpoint. This expression increases
with \(d\), so the loop exhausts the remaining gaps. For each
remaining row and index, the checker bounds the reciprocal
\(j\)-th power of the multiplier by the positive Taylor sum of
\(\exp[Ls\,k(j-k)/2]\) through degree 64. Summing the row losses
gives the displayed bound.

The reference array uses the exact private coefficients times
\(\exp[-S_0 k(j-k)/(2j)]\). The endpoint-normalized profile and
reference array are in the same coefficient box at the actual
\(S_0\), so the straight coefficient segment stays in the box.
Move \(S_0\) to zero along the reference arrays. Their limit is
the private array, with the roots isolated by exact rational
Rouché disks.

## Shared outer-mean boxes

For the two remaining degree-two branches, let \(u,w\) be the
shared outer reflected means at orders eight and nine.
On a row, the Gauss reciprocals are
\[
 x=\frac{d}{8\ell+h},\qquad x'=\frac{d}{9\ell+h},\qquad
 1\le h\le d .
\]
The squared contribution of the row to \(u\) is
\(\prod_h(8\ell+h)/(8\ell+7d+h)\). Every other row factor is at most
one, so the row can occur in an interval with lower endpoint \(u_-\)
only when its squared contribution is at least \(u_-^2\).
To bound the reciprocals exactly, use the floor \(t=1/35\).
A row with \(d/(8\ell+1)>t\) satisfies
\(\ell<(d/t-1)/8\), and the condition \(S_0\le S_*\) gives
\(d<S_*(1+1/(8t))\). Enumerating the resulting finite rows, and
taking the maximum of \(t\) and each retained row's corresponding
reciprocal, bounds \(x,x'\) on the whole interval.

For \(u<1\), put \(\rho=(-\log w)/(-\log u)\). The identity
\(x'=8x/(9-xh/d)\) gives
\[
 1+\frac{64}{9}x\le 1+8x'\le\frac{1+7x}{1-x/9}.
\]
With \(C\) an upper bound on \(x\), the corresponding bounds on
\(\rho\) are
\[
 \frac{\log(1+(64/9)C)}{\log(1+7C)}
 \le \rho\le 1-\frac{\log(1-C/9)}{\log(1+7C)}.
\]
For the lower bound the numerator logarithm is a concave function
of \(\log(1+7x)\) and vanishes at zero. For the upper bound
\(-\log(1-x/9)\) is convex and \(\log(1+7x)\) is concave, each
with value zero at zero. Quotient monotonicity therefore gives
the two bounds termwise and then after summation. The resulting
interval contains \(64/63\). At \(u=1\) all shared unequal rows
are absent.

To obtain the remaining coefficient exponents in \(u\), use
the [reflected mean and skew formulas](MULTIVARIATE_TERMINAL_EIGHT_PRIVATE_LINEAR.md#reflected-means-and-skew-links)
with \(j=8,9\), \(m=j-1\), and with the reciprocal at the
corresponding order. For one Gauss factor set \(s=\log(1+mx)/2\)
and \(g_h=\log(1+hx)\). At distance \(r=\min(k,j-k)\) from the
endpoints, the reflected mean loss divided by \(s\) is
\[
 \sum_{h=0}^{r-1}\frac{\operatorname{artanh}((1-2h/m)\tanh s)}{s}.
\]
Each summand is nonincreasing in \(s\), because its numerator is
concave in \(s\) and vanishes at zero. The ratio ranges between
its value at the reciprocal cap and \(e_{j,k}=k(j-k)/(j-1)\).
The local skew up to orientation is \(\sum_h\omega_hg_h\), with
\(\omega_h=(\mathbf 1_{\{r\le h<j-r\}}-(j-2r)/j)/2\).
As a function of \(\log(1+mx)\) its first derivative starts at
zero. After positive factors are removed, the numerator of its
second derivative is
\[
 \sum_{h=1}^{m-1}\omega_h h(m-h)
   \prod_{\substack{1\le l<m\\ l\ne h}}(1+lx)^2 .
\]
The checker proves it nonnegative on \([0,1]\) by replacing each
negative higher-degree monomial \(b x^i\) by its lower bound
\(b x^{i-2}\) and checking the remaining coefficients.
The midpoint skew at order eight is zero. Convexity bounds
the absolute skew divided by \(s\) by its value at the reciprocal
cap. Adding these bounds to the reflected mean bounds controls
the exponents in \(u\), with the extra factor \(\rho\) at order nine.

The checker confirms that the coefficient boxes contain the reference
array with exact private factors times \(u^{e_{8,k}}\) or
\(u^{(64/63)e_{9,k}}\), respectively. The actual array connects to
the reference by the straight coefficient segment at its \(u\).
Move \(u\) to zero on the reference arrays, whose endpoints at
zero are \(1+z^8\) and \(1+z^9\).

## Ray projections and root counts

The [exact replay](verify_multivariate_terminal_nine_private_rows.py)
uses fixed rays in each high-strength box and seven rays with
continuous piecewise affine cosines along each lower-strength path.
It verifies the roots are between \(1/8\) and \(8\) by constant and
leading-term domination. Within that annulus it uses the
[rational Bernstein projection test](MULTIVARIATE_TERMINAL_SIX_JOINT_RAY_BANDS.md#moving-rays-and-exact-replay)
on each closed radial cell. Subdivision of strength or outer mean
uses closed subintervals so the displayed coefficient bounds apply
at their common endpoints. The checker establishes the angular
counts at the binomial and private anchors by exact sign and disk
tests, respectively.

In increasing upper-half-plane argument the eight sectors contain
one root each, of orders \(9,8,9,8,9,8,9,8\). The remaining
ninth-order root is negative real at the anchors. During continuation
join the first and last sectors to their conjugates. In the last joined
sector the single ninth-order root stays real by conjugation and stays
negative since the constant coefficient is one. The two actual
eighth-order roots in that sector are nonreal by strict positivity of
the even factorial moment on the real axis, from its positive Gamma
integral. Thus the actual arguments cannot coincide. Rational Bézout
transfers coprimality to every characteristic-zero field.

Replay from the repository root with

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_nine_private_rows.py
~~~

The replay verifies 420 shared row-index comparisons, 7 skew
convexity inequalities, 848 high-band radial cells, 6,376
shared-curvature radial cells, and 8,077 shared outer-mean radial
cells. The two outer-mean enumerations retain 248 and 227 row
types.
