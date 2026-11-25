# Terminal-ten angular separation with private degree at least two

For two distinct ordinary exponent columns, delete the equal rows and
write \(p,q\) for the sums of the positive exponents on their two
respective sides in rows with a zero on the opposite side.
Write \(S_0=\sum d^2/(\ell+d)\) over the unequal shared rows, with
smaller exponent \(\ell\ge1\) and gap \(d>0\). Put \(S=p+q+S_0\).

**Theorem.** The ninth and tenth ordinary factorial moment polynomials
have disjoint root arguments, including after independent positive
rescalings, in each of these cases.

- \(p+q\ge2\).
- \(S\ge77/10\).
- \(p+q=1\) and \(S_0\ge57/10\).

Their coprimality holds over every characteristic-zero field.

## Fixed coefficient boxes

Normalize the two endpoint coefficients to one at each order.
A private row of increasing gap \(d\) contributes
\[
             f_{j,k}(d)=\frac{(dk)!}{((dj)!)^{k/j}}.
\]
A decreasing row replaces \(k\) by \(j-k\). All normalized multipliers
are at most one. The [private-factor bound](MULTIVARIATE_TERMINAL_NINE_PRIVATE_ROWS.md#private-and-shared-coefficient-bounds)
\(f_{j,k}(d)\le f_{j,k}(1)^d\) follows from log convexity of the ratios
\((dk)!/(k!)^d\). The [shared curvature bound](MULTIVARIATE_TERMINAL_SEVEN_PRIVATE_ROWS.md#private-factors-and-shared-loss)
bounds the product of the shared multipliers by
\(\exp[-S_0 k(j-k)/(2j)]\). Together these give the interior coefficient
caps used in the fixed ray checks. The total-curvature test uses
\(\binom jk \exp[-S k(j-k)/(2j)]\) without distinguishing private rows.

After column reversal take \(p\ge q\). The three zero-strength fixed
boxes retain four increasing unit factors, two unit factors on
each side, or three increasing unit factors and one decreasing unit
factor. They cover \(p\ge4\), \(q\ge2\), and \(p\ge3,\ q\ge1\),
respectively. The remaining private configurations with \(p+q\ge2\)
have the following cutoffs. The singleton private configuration has
the cutoff \(57/10\).

| Increasing gaps | Decreasing gaps | Fixed-ray cutoff on \(S_0\) | Lower range |
|---|---|---:|---|
| \(3\) | none | \(1/4\) | Shared curvature |
| \(1,1,1\) | none | \(9/8\) | Shared curvature |
| \(1,2\) | none | \(5/8\) | Shared curvature |
| \(1,1\) | none | \(27/8\) | Shared outer mean |
| \(2\) | none | \(11/4\) | Shared curvature |
| \(1\) | \(1\) | \(13/4\) | Shared outer mean |
| \(1,1\) | \(1\) | \(7/8\) | Shared curvature |
| \(2\) | \(1\) | \(3/8\) | Shared curvature |

For each fixed box the [exact checker](verify_multivariate_terminal_ten_private_degrees.py)
verifies leading and constant coefficient domination on \(|z|=8\)
and \(|z|=1/8\). It encloses the coefficients from zero to their
rational caps, with endpoints fixed. The binomial array \(1+z^j\)
belongs to the box, and its eight ray levels have exact alternating
cosine signs. The Bernstein projection tests exclude zeros along
these rays throughout the annulus and the straight segments from
the binomial array to the actual arrays.

## Two-sided shared-curvature boxes

For the six lower-strength branches marked shared curvature the
shared multipliers \(F_{j,k}\) satisfy
\[
 e^{-L S_0 k(j-k)/(2j)}\le F_{j,k}
                  \le e^{-S_0 k(j-k)/(2j)}.
\]
The strength endpoints are \(1/8,1/4,3/8,7/8,5/4,17/8,25/8\),
with factors \(11/10,6/5,313/250,29/20,3/2,9/5,21/10\).
For a row of strength \(s=d^2/(\ell+d)\), the logarithmic loss
comparison follows from the shared bound if \(1+d/\ell\le L\).
Otherwise \(s>\ell(L-1)^2/L\), so only
\(\ell\le S_*L/(L-1)^2\) need be retained at strength endpoint \(S_*\).
The checker enumerates gaps until \(d^2/(\ell_{\max}+d)>S_*\)
and proves the upper log loss at both orders by integer factorial
comparison against a positive exponential Taylor lower bound.
Summation proves the inequality over all shared rows.

At actual strength \(S_0\), both coefficient arrays and the reference
formed from their exact private factors times \(e^{-S_0k(j-k)/(2j)}\)
are included in the coefficient boxes. Straight coefficient segments
join them. Sliding the reference along \(S_0\) down to zero ends at the
private arrays, with roots isolated in exact rational Rouché disks.
All eight moving rays are strictly excluded throughout the annulus.

## Shared outer-mean boxes

In the two remaining degree-two configurations, let \(u,w\) be the
shared reflected outer means at orders nine and ten after removing
the binomial coefficients. The row with base \(\ell\) and gap \(d\)
contributes squared lower mean
\[
                  \prod_{h=1}^d\frac{9\ell+h}{9\ell+8d+h}
\]
and its lower reciprocal is \(x=d/(9\ell+h)\). The upper reciprocal
is \(x'=9x/(10-xh/d)\). On each outer interval, the checker enumerates
rows above reciprocal floor \(1/35\) subject to the strength cutoff,
and discards a row if its squared mean is below the squared interval
minimum. Rows at or below the floor have upper reciprocal at most
\(9(1/35)/(10-1/35)\). The enumeration is finite because
\(x>c\) implies \(\ell<d/(9c)-1/9\), and then
\(d<S_*(1+1/(9c))\).

Let \(C\) bound the lower reciprocals. For \(u<1\), the outer ratio
\(\rho=\log w/\log u\) is in
\[
 \left[\frac{\log(1+(81/10)C)}{\log(1+8C)},
       1-\frac{\log(1-C/10)}{\log(1+8C)}\right].
\]
The identity for \(x'\) gives \(1+(81/10)x\le1+9x'
\le(1+8x)/(1-x/10)\). Taking the two logarithmic ratios termwise
and then summing proves the interval. Its lower endpoint follows
from concavity in \(\log(1+8x)\), and its upper endpoint uses that
\(-\log(1-x/10)\) is convex and \(\log(1+8x)\) is concave.

The [reflected mean and skew formulas](MULTIVARIATE_TERMINAL_NINE_PRIVATE_ROWS.md#shared-outer-mean-boxes)
apply with \(j=9,10\). They bound the reflected loss divided by the
outer log loss between its value at the reciprocal cap and
\(k(j-k)/(j-1)\). The skew divided by the outer loss is bounded by
its absolute value at the reciprocal cap. The exact checker proves
the skew convexity inequalities for these orders on \(0\le x\le1\),
which includes every retained reciprocal. All positive logarithm
and power intervals are outward enclosures.

The resulting coefficient boxes contain the actual arrays and the
comparison arrays with the same private cores and shared factors
\(u^{k(9-k)/8}\), \(u^{(81/80)k(10-k)/9}\). The comparison arrays slide
continuously down to \(u=0\), where their binomial anchors have exact
cosine signs. Straight coefficient segments connect the pairs to
these comparison arrays at the actual outer mean. The checker
certifies the whole path by exact Bernstein projection on the eight
moving rays and domination at the two annulus boundaries.

## From sector counts to coprimality

The nine upper sectors contain one root each, in order
\(10,9,10,9,10,9,10,9,10\), and there is one further negative real
nonic root at the anchors. Join the first and last sectors to their
lower reflections across the real axes. In the last joined sector
the sole nonic root stays real by conjugacy. The even moment is
strictly positive for every real argument because it is the integral
of the tenth power of a nonzero real binomial against a positive
Gamma measure. Thus the two decic roots in that sector are nonreal,
so the arguments are disjoint there as well. Rational Bézout gives
the characteristic-zero claim.

The checker replays 1,160 fixed-ray cells, 476 row inequalities,
14,213 curvature cells, eight skew polynomials and 21,432 outer-mean
cells. Run it with SymPy available.

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_ten_private_degrees.py
~~~
