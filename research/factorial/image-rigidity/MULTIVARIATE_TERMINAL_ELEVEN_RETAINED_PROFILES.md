# Terminal-eleven separation for retained shared rows

Let \(a,b\) be distinct columns of nonnegative integer exponents.
Delete equal rows and write the smaller exponent and the gap of each
unequal shared row as \(\ell,d\). Its strength is
\(s=d^2/(\ell+d)\), and \(S_0\) is the sum of these strengths. Its
largest order-ten Gauss reciprocal is \(d/(10\ell+1)\).

**Theorem.** The tenth and eleventh factorial moment polynomials
have disjoint root arguments under independent positive rescalings
for profiles with exactly one private row of positive exponent one,
and for fully shared profiles with a Gauss reciprocal exceeding
\(1/8\). In particular they are coprime over every
characteristic-zero field. The [fixed-ray theorem](MULTIVARIATE_TERMINAL_ELEVEN_HIGH_CURVATURE.md)
already separates the first class for \(S_0\ge6\) and the second
for \(S_0\ge8\). Below these cutoffs the proof keeps every shared
row above the reciprocal threshold exactly and uses a joint
enclosure for the residual factors.

## Finite retention and bounds for the residual outer mean

Retain every row with \(d/(10\ell+1)>1/8\).
More generally, if \(x=d/(10\ell+1)>c\), then
\(\ell<d/(10c)-1/10\). Therefore a row of strength at most \(B\)
satisfies \(d<B(1+1/(10c))\).
The [checker](verify_multivariate_terminal_eleven_retained_profiles.py)
enumerates all those integer bases and gaps and all their
multiplicities under the two strength cutoffs. It records the
orientations of the retained rows. For fully shared profiles a
simultaneous column reversal is allowed, and the retained list
is nonempty. The tuples are sorted in the same lexicographic order
before comparing a list to its reversal. For one private row its orientation is fixed by
column reversal, so the retained profiles include both shared
orientations and the empty list.

For a retained list write \(B\) for the remaining strength budget.
Let \(v=(F_{10,1}F_{10,9})^{1/2}\), where \(F_{j,k}\) denotes
the endpoint-normalized residual multiplier after removing the
binomial coefficient. The Gauss product gives
\[
 -\log v\le \sum_{\rm residual\ rows}\frac d2\log(1+9x)
 \le\frac9{10}\log(17/8)B .
\]
In fact \(d/s=(\ell+d)/d<1+1/(10x)\), and
\((1+1/(10x))\log(1+9x)\) is increasing. To check the latter
assertion put \(z=9x\). Its derivative has the sign of
\(z(1+(10/9)z)-(1+z)\log(1+z)\), whose derivative is
\((20/9)z-\log(1+z)>0\) for \(z>0\).
The maximum at \(x=1/8\) gives the stated constant.
Thus the residual runs over \([(8/17)^{9B/10},1]\), with the
bottom endpoint rounded downward by interval arithmetic.
Every residual reciprocal satisfies \(x_h\le(v_-^{-2}-1)/9\)
on a cell starting at \(v_-\), since \(v^{-2}\) is the product
of the factors \(1+9x_h\). The checker also uses the
[ordinary cap enumeration](verify_multivariate_terminal_eleven_shared_light.py)
at \(v_-\) to reduce this bound, rounding upward before that
enumeration.

## Joint coefficient comparison

Use endpoint normalization and let \(F_{j,k}\) be a residual
coefficient multiplier at order \(j\) after removing its binomial
coefficient. A retained core multiplies these factors by known
positive numbers. The
[endpoint attenuation inequality](MULTIVARIATE_ORDINARY_ADJACENT_COEFFICIENT_ATTENUATION.md#an-all-order-endpoint-inequality)
applies to the residual because each unequal shared row has maximal
exponent at least two. For \(1\le k\le10\) its factors bound
\[
 F_{11,k}\le F_{10,k}\,v^{k^2/160},\qquad
 F_{11,k}\le F_{10,k-1}\,v^{(11-k)^2/160}.
\]
Indeed at orders ten and eleven the strength attenuation is at
least \(\log(231/230)\ge1/231\). The bound on \(-\log v\) and
the exact comparison \(231(9/10)\log(17/8)<160\) give the
displayed exponents.

For the comparison Gaussians, put \(e_{j,k}=k(j-k)/(j-1)\).
The order-ten residual has multipliers \(v^{e_{10,k}}\), and
the order-eleven residual uses \(v^{\rho e_{11,k}}\), where
\(\rho\) is the successive-order outer log ratio. Matching the
reciprocals gives
\[
 \frac{\log(47/22)}{\log(17/8)}\le \rho\le
  1+\frac{\log(1/(1-c/11))}{\log(1+9c)}
  \quad(c=1/8).
\]
For a matched factor, \(x'=10x/(11-\gamma x)\) with
\(\gamma=h/d\in(0,1]\). It follows that
\[
 1+(100/11)x\le1+10x'\le\frac{1+9x}{1-x/11}.
\]
The resulting lower quotient is nonincreasing in \(x\), by
concavity of \(y\mapsto\log(1+(100/99)(e^y-1))\) and its zero
at \(y=0\). The upper quotient increases because
\(-\log(1-x/11)/x\) increases and \(\log(1+9x)/x\) decreases.
Evaluate the quotients at the cap and sum with the positive
weights \(\log(1+9x)\). The lower ratio is greater than
\(1+1/160\), as checked by rational logarithm enclosures.
For the attenuation of the Gaussians,
\[
 (1+1/160)e_{11,k}-e_{10,k}-k^2/160
   =k(k-1)\left(1/90-11/(10\cdot160)\right)\ge0 .
\]
Reflection gives the same bound at the other endpoint. The
reference ratio \(100/99\) lies in the allowed interval.

The checker intersects its independent coefficient boxes with
the two attenuation bounds. It also uses the grouped
[matched tanh coefficient comparison](MULTIVARIATE_TERMINAL_TEN_SHARED_LIGHT.md#grouped-moment-and-moving-ray-comparison)
with the order-ten and eleven constants. In particular the skew
transfer has midpoint coefficient \((100/99)^2(1+\eta^2)\)
and normalized error coefficient \(2(100/99)^2\eta\) with
\(\eta=c/(22+98c)\). The residual signed parameter lies in
\([-1,1]\), because the retained orientations are fixed.
Each grouped vector coordinate is multiplied by an
outward enclosure of the corresponding retained core coefficient.

At the actual residual mean, scale the matched log skew and mean
defects jointly to zero. Their signed moment bounds are preserved,
including the Jensen lower bounds by \(t^{2l}\le t\) for \(0\le t\le1\).
The logarithmic coefficient boxes and attenuation inequalities
are also preserved, since they hold at both endpoints of this
logarithmic interpolation. Then move the Gaussian outer log
ratio to \(100/99\) and move the residual mean to the cell
midpoint. The same bounds hold throughout these Gaussian
segments. At residual mean one the comparison is the retained
core itself.

## Ray counts for the comparison paths

The checker verifies constant and leading coefficient domination
on the circles of radii \(1/12\) and \(12\) from the retained
core times the binomial coefficients alone. All normalized
residual multipliers along the interpolation are at most one,
which gives the same annulus for every array on the path.

For each mean cell, outward complex balls isolate all roots of
the two midpoint Gaussians with their fixed core. The checker
requires alternating upper root arguments of orders
\(11,10,11,10,11,10,11,10,11,10\), and chooses nine rational
ray cosines strictly between the isolated arguments. It proves
every ray exclusion along the whole coefficient homotopy by
strict interval Bernstein projections. The order-eleven root
in the sector joined across the negative real axis stays real
by conjugation and the root count of that sector. The actual
even order-ten moment is strictly positive on the real axis. Its
Gamma integral has nonnegative integrand and strict positivity
on an open set because the two monomials are distinct. Its two
roots in that joined sector are therefore nonreal. Disjoint upper and lower arguments give
characteristic-zero coprimality by rational Bézout.

Replay the two enumerations from the repository root.

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_eleven_retained_profiles.py --native private
python3 research/factorial/image-rigidity/verify_multivariate_terminal_eleven_retained_profiles.py --native shared
~~~

The full replay checks 333 one-private lists with 458 outer cells
and 25,665 radial cells. It checks 881 fully shared lists with
1,213 outer cells and 65,884 radial cells. The checker also accepts
specific indices after either family name. The native path uses
the same outward Bernstein arithmetic from the terminal-ten
checker. Order eleven uses an additional coefficient slot and
a uniformly rescaled outward Bernstein transform.
