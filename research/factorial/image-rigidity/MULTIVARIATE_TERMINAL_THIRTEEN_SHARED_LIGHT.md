# Terminal-thirteen light shared profiles on certified outer bands

Let two distinct ordinary exponent columns have only shared unequal
rows. Let \(u\) be their endpoint-normalized order-twelve reflected
outer mean after removal of the binomial factors. On a shared row
write \(\ell\ge1\) for the smaller exponent and \(d\ge1\) for the gap.

**Theorem.** At every \(u\) in the following closed bands, if each
\(d/(12\ell+1)\) is at most the stated cap, then the order-twelve
and order-thirteen factorial moment polynomials have disjoint root
arguments. They are therefore coprime over every characteristic-zero
field.

| Outer mean \(u\) | Reciprocal cap |
| --- | ---: |
| \([0,.18]\) | \(1/8\) |
| \([.18,.34]\) | \(3/40\) |
| \([.34,.43]\) | \(3/50\) |
| \([.43,.47]\) | \(1/20\) |
| \([.47,.50]\) | \(1/22\) |
| \([.50,.52]\) | \(1/20\) |
| \([.52,.54]\) | \(1/22\) |
| \([.54,.55]\) | \(1/30\) |
| \([.55,.56]\) | \(1/20\) |
| \([.56,.57]\) | \(1/25\) |
| \([.57,.60]\) | \(1/40\) |
| \([.60,.61]\) | \(1/36\) |
| \([.61,.66]\) | \(1/40\) |
| \([.66,.70]\) | \(1/50\) |
| \([.71,.72]\) | \(1/50\) |

## Matched logarithmic enclosure

At order twelve the squared row mean is
\[
 u_{\ell,d}^{\,2}
   =\prod_{h=1}^d\frac{12\ell+h}{12\ell+11d+h}.
\]
If \(d/(12\ell+1)>t\), then every Gauss reciprocal in the same
row exceeds \(t/(1+t)\). Thus \(u_{\ell,d}^{\,2}\) is below
\((1+11t/(1+t))^{-d}\). Each row mean is at least the product mean.
On a slab with positive lower endpoint this restricts the gaps
of rows above the reciprocal floor \(1/40\). For each such gap the
condition on \(d/(12\ell+1)\) bounds \(\ell\). The checker uses this
finite list to tighten the reciprocal cap and, above the floor,
bounds the upper/lower log ratio by the individual factor products
of the same integer rows.

Write \(U=-\log u\) and \(W=-\log w\) for the two outer losses.
With a lower-order reciprocal cap \(c\), matched factors obey
\(x'=12x/(13-hx/d)\) and \(x'\le12c/(13-c)\).
The [outer-loss comparison](MULTIVARIATE_MATCHED_OUTER_RATIO_MONOTONICITY.md)
gives \(W/U\ge R=144/143\). The termwise upper ratio bound used in
the [terminal-twelve model](MULTIVARIATE_TERMINAL_TWELVE_SHARED_LIGHT.md#exact-row-cutoffs-and-outer-ratios)
becomes \(1+\log(1/(1-c/13))/\log(1+11c)\).
The exact integer list above the reciprocal floor sharpens this
upper ratio on each certified slab.

The logarithmic enclosure retains the Gaussian direction, signed
even-power and unsigned odd-power tanh corrections, then groups the
quadratic and cubic exponential terms before taking interval bounds.
The [excess-loss inequality](MULTIVARIATE_TANH_EXCESS_LOSS_TRANSPORT.md)
couples \(W-RU\) to the change in signed second moment. The
[linked lower bound on the excess](MULTIVARIATE_MATCHED_OUTER_RATIO_MONOTONICITY.md#coupling-to-the-lower-order-tanh-moments)
also uses the signed second moment within a subinterval of the model.
The [higher signed-moment bound](MULTIVARIATE_CENTERED_HIGHER_SIGNED_TANH.md)
absorbs exact multiples of that moment into the shared skew direction.
The checker verifies the rational centering constants before evaluation.

## Exact separating rays

The eleven ray cosines are rational and affine in \(u\) between
successive knots. The checker uses the small-mean coefficient boxes
on the initial band and grouped Bernstein projection thereafter.
Annulus bounds follow from the independent coefficient boxes, using
endpoint domination or adjacent coefficient ratio bounds. At the first
anchor the polynomials are \(1+z^{12}\) and \(1+z^{13}\). Exact cosine
signs certify their initial sector counts. The subsequent anchors
use disjoint rational Rouché disks on the same Gaussian references.

For any profile in a grouped band, scale its skew and reflected mean
defects to zero, keeping \(u\) and the loss ratio fixed. The Jensen
lower bounds remain valid because they shrink at least quadratically,
and the defects themselves shrink linearly. The signed transport
error and higher signed defects shrink linearly as well, so the same
enclosures apply on this path. Then move the Gaussian ratio to \(R\)
and the outer mean to the band's lower anchor. The certified rays
exclude roots throughout the path and preserve the angular sector
counts.

The twelve upper sectors alternate starting with order thirteen.
Join the first and last sectors to their conjugate sectors across
the two real rays. The last joined sector contains the order-twelve
conjugate pair and the unique negative-real order-thirteen root.
The even Gamma moment is strictly positive on the real axis, so the
order-twelve pair stays nonreal. This separates the arguments in
the last sector as well. Since both moment polynomials have integer
coefficients and nonzero endpoints, a rational Bézout identity
transfers coprimality to every characteristic-zero field.

The [exact checker](verify_multivariate_terminal_thirteen_shared_light.py)
replays 805,379 radial cells on the displayed bands. The native
Bernstein evaluator uses outward dyadic arithmetic with an explicit
smaller transform bound for the fourteen-coefficient case.

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_thirteen_shared_light.py --native
~~~
