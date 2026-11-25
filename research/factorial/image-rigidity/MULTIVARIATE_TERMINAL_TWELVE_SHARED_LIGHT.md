# Terminal-twelve light shared profiles through outer mean .66

Let two distinct ordinary exponent columns have only shared unequal
rows. Let \(u\) be their endpoint-normalized order-eleven reflected
outer mean after removal of the binomial factors. In each band below,
assume that every reciprocal \(d/(11\ell+1)\) is at most the given cap,
where \(\ell\) is the smaller exponent and \(d\) is the gap. The
eleventh and twelfth factorial moment polynomials then have disjoint
root arguments under independent positive rescalings. Their
coprimality follows over every characteristic-zero field.

| Outer mean \(u\) | Reciprocal cap |
| --- | ---: |
| \([0,.36]\) | \(1/8\) |
| \([.36,.40]\) | \(3/50\) |
| \([.40,.47]\) | \(3/40\) |
| \([.47,.48]\) | \(13/200\) |
| \([.48,.50]\) | \(3/40\) |
| \([.50,.66]\) | \(1/20\) |

## Exact row cutoffs and outer ratios

At order eleven a row has squared outer mean
\[
 u_{\ell,d}^2=\prod_{h=1}^d\frac{11\ell+h}{11\ell+10d+h}.
\]
If \(d/(11\ell+1)>t\), this is strictly below
\((1+10t/(1+t))^{-d}\), since \(11\ell+d<d(1+1/t)\).
Each row mean is at least the product mean. Thus at a slab starting
at \(u_->0\) only finitely many integer rows with reciprocal above
\(t=1/40\) can occur. The checker intersects the declared cap with
the maximum of the floor and the reciprocals of these rows.
The extra knot \(u=.221\) is justified by this calculation. The
rows of maximal reciprocal \(2/17\) cannot occur there.

For successive Gauss factors \(x=d/(11\ell+h)\) and
\(x'=11x/(12-hx/d)\), a lower-order cap \(c\) gives upper-order
cap \(11c/(12-c)\). Put \(U=-\log u\), \(W=-\log w\), and
\(\rho=W/U\). The comparison assumes
\[
 R=\frac{121}{120}\le \rho
 \le 1+\frac{\log(1/(1-c/12))}{\log(1+10c)}.
\]
The lower bound is the [matched outer-ratio inequality](MULTIVARIATE_MATCHED_OUTER_RATIO_MONOTONICITY.md).
The upper bound follows termwise from
\(1+11x'\le(1+10x)/(1-x/12)\) and the monotone logarithm ratio,
as in the [terminal-eleven light proof](MULTIVARIATE_TERMINAL_ELEVEN_SHARED_LIGHT.md).
The coarser lower logarithm ratio from that proof is also sufficient
in the independent annulus boxes.

In the tighter models the checker bounds \(\rho\) above by the
maximum of the upper bound at the reciprocal floor and the exact
outward log-ratio bounds on the finitely retained individual rows.
This uses the fact that \(W/U\) is the weighted average of their
individual ratios, with weights given by the lower-order outer
losses. The same row enumeration used for the cap proves that
the upper ratio is valid throughout each slab.

## Grouped coefficient enclosure

The ten separating rays use rational cosines affine in \(u\)
within each slab. As in the
[preceding grouped proof](MULTIVARIATE_TERMINAL_ELEVEN_SHARED_LIGHT.md#ray-exclusions-and-root-counts),
the logarithmic coefficient model consists of a Gaussian direction,
signed even-power skew corrections, and unsigned odd-power reflected
mean defects. The exact directions at orders eleven and twelve
are obtained by expanding the row logarithms in the matched tanh
coordinates. The signed second moment is \(P=\sum_i\epsilon_i T_i^2/2\).
A common column reversal permits a nonnegative signed normalized
parameter in the fully shared enclosures.

For the coupled models the
[outer-loss and skew inequality](MULTIVARIATE_MATCHED_OUTER_RATIO_MONOTONICITY.md#coupling-to-the-lower-order-tanh-moments)
gives a lower bound on \(\rho-R\) from the minimum of \(|P|/U\) in
the signed subbox and the upper bound on \(U\).
Some slabs also use the
[excess-loss transport theorem](MULTIVARIATE_TANH_EXCESS_LOSS_TRANSPORT.md)
to bound \(P'-R^2P\) by the smaller of the separate matched envelope
and a rational multiple of \(U(\rho-R)\). These bounds apply to the
ordinary row factors at the actual profile.

For the tighter centered slabs the reference log ratio is the
midpoint of the admitted ratio interval. The reference powers
with these rational exponents are enclosed by outward logarithms
followed by a bounded exponential Taylor series and repeated
squaring. In this calculation the truncation tail is the one used
in the earlier grouped exponential comparison. The ordinary and
Gaussian paths still lie in the same coefficient model because
the Gaussian ratio \(R\) belongs to the zero-skew interval.

Quadratic or cubic exponential terms are evaluated as grouped
directions before closing the remainder interval. The angular
projection keeps the common first-order mean displacement in
the coefficient and ray directions. The exact native replay
uses dyadic interval arithmetic, with common outward positive
rescaling of Bernstein controls at order twelve so that signed
integer operations stay in range. A coarser enclosure may combine
small coefficient groups into the interval base, and the
checker tries the complete grouped model when that enclosure
does not prove the sign.

## Ray counts and continuation

At the first anchor the two arrays are \(1+z^{11}\) and \(1+z^{12}\).
Exact cosine polynomial signs certify the two alternating sector
counts there. Other anchors have all their reference roots
isolated by rational Rouché disks. The checker verifies
leading and constant coefficient domination or the strict
adjacent-ratio inequalities for the Eneström–Kakeya annulus,
using independent boxes valid throughout the band.

At fixed \(u,\rho\) scale all the logarithmic skew and reflected
mean defects to zero together. The Jensen lower bounds remain
valid because the lower sides shrink by powers of the scale
at least two, and the actual defects shrink linearly. The
outer-loss constraints and the centered transport errors are
also preserved. Move the Gaussian ratio to \(R\) and then the
mean to the lower anchor. The exact Bernstein signs keep each
ray free of roots throughout this comparison.

In increasing upper argument the eleven sectors have alternating
root orders \(12,11,12,11,12,11,12,11,12,11,12\). The negative
real order-eleven root remains in the last sector joined with
its reflected sector. Positivity of the real order-twelve
moment from its even Gamma integral separates it from the
order-twelve conjugate pair there. Rational Bézout gives the
conclusion over every characteristic-zero field.

The [exact checker](verify_multivariate_terminal_twelve_shared_light.py)
replays 520,781 radial cells across the rational ray and anchor tables.

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_twelve_shared_light.py --native
python3 research/factorial/image-rigidity/verify_multivariate_terminal_twelve_shared_light.py --native 0 .22 .47
~~~
