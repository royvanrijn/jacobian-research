# Terminal-eight separation for three marked singleton profiles

Let all unequal rows of two distinct ordinary exponent columns be shared.
Suppose one row has unordered exponents \((1,3)\), \((2,4)\), or \((3,6)\).
On every remaining unequal row assume \(d/(7\ell+1)\le1/8\),
where \(\ell\) is its lower exponent and \(d\) its gap.
If the sum \(S_R=\sum d^2/(\ell+d)\) over the remaining
rows is at most \(15/4\), the order-seven and order-eight
polynomials have disjoint root arguments after independent
positive rescalings. Thus they are coprime over every
characteristic-zero field.

This includes all three singleton exceptions to the retained-row
test with total curvature less than \(19/4\). The marked
row contributes \(4/3\), \(1\), or \(3/2\), respectively.
The [exact checker](verify_multivariate_terminal_eight_marked_singletons.py)
uses the paired coefficient model of the
[terminal-eight private-linear proof](MULTIVARIATE_TERMINAL_EIGHT_PRIVATE_LINEAR.md).

## The retained coefficients and residual bounds

Orient the marked row as \((a,a+d)\), with \((a,d)=(1,2),(2,2),(3,3)\), by
exchanging the columns if necessary. Keep its factorial contribution
to each coefficient. At order \(j\) and index \(k\), including the
binomial coefficient, the marked contribution is
\[
 B_{j,k}=\binom jk\frac{(ja+dk)!}
             {((ja)!)^{1-k/j}((j(a+d))!)^{k/j}} .
\]
The checker encloses it by integer comparisons of \(j\)-th powers.

Let \(v\) be the order-seven outer geometric mean of the residual
shared arrays. On one residual row put \(r=d/\ell\) and
\(s_i=d^2/(\ell+d)\). Its reciprocal cap gives
\(8d\le7\ell+1\le7\ell+d\), so \(r\le1\).
If \(U_i\) is its contribution to \(-\log v\), then
\[
 \frac{U_i}{s_i}
 \le \frac{1+r}{2r}\log(1+6r/7)\le\log(13/7).
\]
The first bound uses \(d/(7\ell+h)\le r/7\) in each
outer-mean factor. The middle function is
\(\frac37\int_0^1(1+r)/(1+6rt/7)\,dt\) and increases in \(r\)
because its integrand does. Therefore
\[
                    v\ge(7/13)^{15/4}>9/100 .
\]
The replay verifies the final comparison by raising both sides
to the fourth power.

For each \(v\) interval use the candidate row caps from the
[private-linear reciprocal bound](MULTIVARIATE_TERMINAL_EIGHT_PRIVATE_LINEAR.md#ordinary-reciprocal-caps).
That list assumes shared curvature at most \(21/4\), so it includes
every residual row that exceeds its reciprocal floor here.
Intersect its order-seven cap with \(1/8\).
If that cap is \(c\), the matching order-eight reciprocals are
at most \(7c/(8-c)\); intersect this with the candidate row cap
at order eight. These are the caps used in the linked skew,
mean defect and outer-log ratio intervals of the private-linear
proof.

## The two paired ray replays

In the paired projection model use \(B_{j,k}\) in place of the
private coefficient of that proof. It multiplies all shared
contributions by the same marked row throughout the deformation.
In particular the shared reference powers are
\(v^{e_{7,k}}\) and \(v^{(49/48)e_{8,k}}\), where
\(e_{j,k}=k(j-k)/(j-1)\). The proof keeps each shared scalar and
each quadratic product together in the projected polynomials.
Only the cubic exponential remainder is enclosed independently.

Each marked row has six continuously moving ray cosines, given by
rational piecewise affine functions on the residual mean interval.
The marked-only caps give strict dominance on \(|t|=1/10\) and
\(|t|=10\), since all actual and reference residual multipliers
are at most one. On this annulus each accepted radial cell has a
strict rational Bernstein projection sign. The checker bisects
the \(v\) interval when needed.

At fixed \(v\) take the straight coefficient segment to the marked
array times the symmetric reference residual. The coefficient
model contains both arrays and their segment by its independent
scalar, product and remainder intervals. Next increase \(v\) to
one along the reference arrays. The canonical outer-log ratio
\(49/48\) lies in every ratio interval, and \(v=1\) gives the
marked-only anchors.

For each anchor the checker verifies disjoint rational Rouché disks
at all roots and checks that the upper disks lie strictly in the
prescribed sectors. The upper sector counts are one each of
orders \(8,7,8,7,8,7,8\), and there is one negative real septic
root. These counts continue along the two paths, with the first
and last sectors joined to their respective lower reflections.
The sole septic root in the last joined sector stays real.
The actual eighth moment is strictly positive for real coefficient
ratios by the Gamma integral of the eighth power of a nonzero
binomial. Thus its two roots in that sector are nonreal.
Rational Bézout transfers the coprimality to characteristic zero.

The exact replay verifies 2,989, 9,883, and 3,150 radial cells,
respectively, for the three listed marked rows.

Replay from the repository root with

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_eight_marked_singletons.py
~~~
