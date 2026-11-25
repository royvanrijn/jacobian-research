# Retained bands for shared terminal twelve

Let \(a,b\) be distinct ordinary exponent columns with every
unequal row shared. Let \(u\) be the order-eleven outer mean
after endpoint normalization and removal of binomial factors.
For any shared row of smaller exponent \(\ell\) and gap \(d>0\),
its largest order-eleven reciprocal is \(d/(11\ell+1)\).

**Theorem.** In each band below, if all largest reciprocals
are at most the listed upper cap and at least one exceeds
the listed lower cap, the eleventh and twelfth moment
polynomials have disjoint root arguments under independent
positive rescalings. They are therefore coprime in
characteristic zero.

| Range of \(u\) | Lower cap | Upper cap | Retained lists | Outer cells | Radial cells |
| --- | ---: | ---: | ---: | ---: | ---: |
| \([.35,.40]\) | \(17/200\) | \(1/8\) | 2 | 64 | 10,483 |
| \([.36,.40]\) | \(3/40\) | \(17/200\) | 5 | 112 | 16,495 |
| \([.36,.40]\) | \(3/50\) | \(3/40\) | 1 | 32 | 3,855 |
| \([.40,.56]\) | \(3/40\) | \(1/8\) | 6 | 260 | 38,542 |
| \([.46,.48]\) | \(13/200\) | \(3/40\) | 1 | 8 | 986 |
| \([.49,.56]\) | \(11/200\) | \(3/40\) | 1 | 32 | 4,143 |
| \([.50,.56]\) | \(1/20\) | \(11/200\) | 1 | 16 | 2,051 |
| \([.56,.70]\) | \(1/20\) | \(1/8\) | 2 | 114 | 23,434 |
| \([.66,.70]\) | \(1/25\) | \(1/20\) | 4 | 64 | 10,536 |
| \([.70,.76]\) | \(3/80\) | \(1/8\) | 2 | 48 | 9,257 |
| \([.70,.76]\) | \(3/100\) | \(3/80\) | 1 | 16 | 2,405 |
| \([.76,.88]\) | \(1/55\) | \(1/8\) | 14 | 291 | 61,394 |
| \([.855,.88]\) | \(1/60\) | \(1/55\) | 1 | 16 | 3,323 |
| \([.88,.95]\) | \(1/70\) | \(1/8\) | 3 | 80 | 18,341 |
| \([.92,.95]\) | \(1/90\) | \(1/8\) | 4 | 51 | 11,615 |

## Finite enumeration

Retain every row whose reciprocal exceeds the lower cap \(c\).
The squared row mean is
\[
 C_{\ell,d}^{\,2}
   =\prod_{h=1}^d\frac{11\ell+h}{11\ell+10d+h}.
\]
A retained row satisfies \(\ell<(d/c-1)/11\) and
\(11\ell+h<d(1+1/c)\). Hence
\[
 C_{\ell,d}^{\,2}<\left(1+\frac{10c}{1+c}\right)^{-d}.
\]
If the band begins at \(u_-\), each retained squared row mean
and their product must be at least \(u_-^2\). This bounds each
gap, then each smaller exponent, and finally each multiplicity.
The checker compares exact rational products while enumerating
the oriented tuples. It orders the tuple and its reverse-column
image lexicographically before reducing by whole-column
reversal.

For a retained list with product mean \(C\) the residual mean
is \(v=u/C\). The checker encloses
\([u_-/C,\min(1,u_+/C)]\) outward on a rational grid with
denominator 10,000 by squaring both sides of each endpoint
comparison. Every residual reciprocal is at most \(c\).
The [ordinary cap bound](verify_multivariate_terminal_twelve_shared_light.py)
further improves the residual cap on each subinterval.

## Exact angular separation

Apply the
[paired coefficient bounds and sector comparison](MULTIVARIATE_TERMINAL_TWELVE_RETAINED_PROFILES.md#coupled-coefficient-bounds)
with the retained list fixed. Its multipliers are bracketed
on a grid with denominator \(2^{55}\) by comparing integer
powers. The same endpoint attenuation, the independent
coefficient inequalities, and the grouped tanh inequalities
apply to the residual. The checker uses the refined power
enclosure and the
[monotone outer ratio bound](MULTIVARIATE_MATCHED_OUTER_RATIO_MONOTONICITY.md)
for its grouped projections.

The independent coefficient boxes enclose the entire
interpolation from each actual array to the Gaussian
midpoint array. Constant and leading coefficient domination,
or positive adjacent coefficient ratios, bound all roots
between radii \(1/12\) and \(12\). The strict
Eneström–Kakeya inequalities are checked when domination does not
certify the required annulus. Outward complex balls isolate
the midpoint roots, and strict Bernstein projections keep
both interpolations free of roots on the ten same rays.
The even-order Gamma positivity and conjugation argument
then separate the last joined sectors, as in the retained
profile theorem.

## Replay

The [checker](verify_multivariate_terminal_twelve_retained_bands.py)
enumerates these exact cases and verifies their rays. Select
a band to replay its full finite list. Optional numbers
following the band select indices in that list.

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_twelve_retained_bands.py --native low-tight
python3 research/factorial/image-rigidity/verify_multivariate_terminal_twelve_retained_bands.py --native low-narrow
python3 research/factorial/image-rigidity/verify_multivariate_terminal_twelve_retained_bands.py --native low-reduced
python3 research/factorial/image-rigidity/verify_multivariate_terminal_twelve_retained_bands.py --native middle-lower
python3 research/factorial/image-rigidity/verify_multivariate_terminal_twelve_retained_bands.py --native middle-lower-short
python3 research/factorial/image-rigidity/verify_multivariate_terminal_twelve_retained_bands.py --native middle-lower-tight
python3 research/factorial/image-rigidity/verify_multivariate_terminal_twelve_retained_bands.py --native middle-lower-narrow
python3 research/factorial/image-rigidity/verify_multivariate_terminal_twelve_retained_bands.py --native middle
python3 research/factorial/image-rigidity/verify_multivariate_terminal_twelve_retained_bands.py --native middle-narrow
python3 research/factorial/image-rigidity/verify_multivariate_terminal_twelve_retained_bands.py --native middle-upper
python3 research/factorial/image-rigidity/verify_multivariate_terminal_twelve_retained_bands.py --native middle-tight
python3 research/factorial/image-rigidity/verify_multivariate_terminal_twelve_retained_bands.py --native cayley-middle
python3 research/factorial/image-rigidity/verify_multivariate_terminal_twelve_retained_bands.py --native cayley-middle-tight
python3 research/factorial/image-rigidity/verify_multivariate_terminal_twelve_retained_bands.py --native cayley-upper
python3 research/factorial/image-rigidity/verify_multivariate_terminal_twelve_retained_bands.py --native cayley-upper-tight
~~~

The listed replays check 48 retained list occurrences in 1,204 outer
cells and 216,860 radial cells.
