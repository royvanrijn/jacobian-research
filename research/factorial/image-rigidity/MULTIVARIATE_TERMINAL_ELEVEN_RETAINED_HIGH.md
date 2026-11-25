# High and tightened retained bands at terminal eleven

Let two distinct ordinary exponent columns have only shared unequal
rows. Write \(u\) for their normalized order-ten reflected outer mean,
and \(x=d/(10\ell+1)\) for the largest Gauss reciprocal of a row
with smaller exponent \(\ell\) and gap \(d\). In a row of the table
assume \(u\) lies in the band, every reciprocal is at most the parent
cap, and at least one exceeds the light cap.

| Outer mean \(u\) | Parent cap | Light cap | Retained rows up to reversal |
| --- | ---: | ---: | --- |
| \([.76,.88]\) | \(1/8\) | \(1/55\) | 12 singletons, 20 pairs, and five triples |
| \([.88,.95]\) | \(1/8\) | \(1/70\) | One \((a,a+1)\), \(a=3,4,5,6\) |
| \([.42,.44]\) | \(1/8\) | \(93/1000\) | One \((2,4)\) |
| \([.57,.64]\) | \(9/100\) | \(3/50\) | One \((3,5)\) |
| \([.73,.76]\) | \(11/200\) | \(9/200\) | One \((2,3)\) |

The tenth and eleventh moment polynomials have disjoint root
arguments under independent positive rescalings in every stated
band. Rational Bézout gives coprimality over any
characteristic-zero field.

## Exact verification

The [middle-band enumeration](MULTIVARIATE_TERMINAL_ELEVEN_RETAINED_MID_BANDS.md#exhaustion-and-residual-ranges)
applies with the displayed parent and light caps. Retain each
row with reciprocal exceeding the light cap. Its squared outer
mean is at most \((1+9c/(1+c))^{-d}\) for light cap \(c\).
The product of retained means is at least the band minimum.
These inequalities give a finite enumeration of bases, gaps
and multiplicities. The checker computes the exact products,
restricts the retained reciprocals to the corresponding parent
cap, and identifies the two global orientations. This leaves
the indicated lists. It brackets the residual mean
\(v=u/u_C\) outward on the rational grid with denominator
10,000, using squared means and exact comparisons.

The [retained homotopy](MULTIVARIATE_TERMINAL_ELEVEN_RETAINED_PROFILES.md#joint-coefficient-comparison)
applies to each remaining residual profile. In the grouped
coefficient enclosure use the sharper
[ordinary outer-ratio inequality](MULTIVARIATE_MATCHED_OUTER_RATIO_MONOTONICITY.md)
\(\rho\ge100/99\). The homotopy holds the actual \(\rho\) fixed
while scaling log defects, then moves it to \(100/99\).
Thus both stages respect the stronger ratio bound, and the
Gaussian reference has the same ratio. The signed transfer
error, mean-defect and tail estimates are those in the
earlier matched comparison.

The retained multipliers have exact rational order-th powers.
The checker encloses their roots by integer-power comparisons
on a dyadic grid of denominator \(2^{55}\). It also computes
the Gaussian reference powers on a dyadic grid of denominator
\(2^{52}\) using integer-power comparisons. These higher
precision bounds control the small projections close to the
negative real axis. Both constructions are outward before
entering the Bernstein evaluation. The optional native replay
may combine small parameter-vector groups into a coefficient
interval in the base group. This only enlarges the enclosure,
and the uncompressed model is still available when the enlarged
one gives no strict sign.

The checker bounds every adjacent coefficient ratio, in both
directions, strictly below 12 using coefficient intervals for
the entire residual band. This gives an annulus by the
Eneström–Kakeya argument. In detail, for an upper ratio bound
\(R\), rescale the polynomial to \(Q(Rz)\), whose coefficients
\(a_k\) are then increasing. In \((1-z)Q(Rz)\), the term
\(-a_jz^{j+1}\) strictly dominates the sum of the other terms
when \(|z|>1\), as their nonnegative coefficients sum to
\(a_j\) and have degree at most \(j\). Apply the upper bound
to the reversed polynomial for the inner radius. The strict
inequalities enforce both circles throughout the homotopy.

At each accepted mean cell the midpoint core Gaussian has
alternating upper arguments certified by outward complex root
balls. Nine rational rays separate these arguments. Exact
Bernstein signs exclude every ray along the homotopy, so
the sector counts transfer to the actual moments. The
order-eleven negative root is the only root of its order
in the last conjugate-joined sector. The Gamma integral
makes the actual tenth moment strictly positive for real
arguments, which separates its last pair from the real root.

## Replay

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_eleven_retained_high.py --native
python3 research/factorial/image-rigidity/verify_multivariate_terminal_eleven_retained_high.py --native --upper-middle
python3 research/factorial/image-rigidity/verify_multivariate_terminal_eleven_retained_high.py --native --tightened
~~~

The first command checks the four high unit-row lists.
The second checks the 37 upper-middle lists. Its sorted list indices
0 through 36 may also be selected explicitly. The third checks the
three tightened bands.

The upper-middle replay checks 584 outer cells and 93,368 radial
cells. The high replay checks 57 outer cells and 11,764 radial cells.
The tightened replay checks 48 outer cells and 6,296 radial
cells, with per-band radial counts 877, 3,734, and 1,685.
